import httpclient, json, algorithm, tables, sequtils, re, math
import os, times
import htmlparser
import q, xmltree # Parse XML
from streams import newStringStream
import docopt # Formating the command-line
import strutils # Split string
import hts
import lapper

# Local libs
import ./clinvcfpkg/gff

type
  ClinSig* = enum
    csBenign = "Benign",
    csBenignLikelyBenign = "Benign/Likely benign",
    csLikelyBenign = "Likely benign",
    csUncertainSignificance = "Uncertain significance",
    csLikelyPathogenic = "Likely pathogenic",
    csPathogenicLikelyPathogenic = "BPathogenic/Likely pathogenic",
    csPathogenic = "Pathogenic",
    csUnknown = "not provided",
    csDrugResponse = "drug response",
    csRiskFactor = "risk factor",
    csAffects = "Affects",
    csAssociation = "association",
    csProtective = "protective",
    csConflictingDataFromSubmitters = "conflicting data from submitters",
    csOther = "other"

  RevStat* = enum
    rsNoAssertion = "no assertion provided",
    rsNoAssertionCriteria = "no assertion criteria provided",
    rsNoAssertionVariant = "no assertion for the individual variant",
    rsSingleSubmitter = "criteria provided, single submitter",
    rsConflicting = "criteria provided, conflicting interpretations",
    rsMutlipleSubmitterNoConflict = "criteria provided, multiple submitters, no conflicts",
    rsExpertPanel = "reviewed by expert panel",
    rsPracticeGuideline = "practice guideline"

  Submission* = ref object
    clinical_significance: ClinSig
    review_status: RevStat
    submitter_id: int

  MolecularConsequence* = ref object
    description: string
    so_term : string

  ClinVariant* = ref object
    variant_id: int32
    allele_id: int32
    rsid: int32
    chrom: string # we could use an 8bit/16bit integer ...
    pos: int32
    ref_allele: string
    alt_allele: string
    molecular_consequences: seq[MolecularConsequence]
    submissions: seq[Submission]

var
  acmg_clinsig = @[csBenign, csLikelyBenign, csUncertainSignificance, csLikelyPathogenic, csPathogenic]
  non_acmg_clinsig = @[csDrugResponse, csRiskFactor, csAffects, csAssociation, csProtective, csConflictingDataFromSubmitters, csOther]

proc findNodes(n: XmlNode, tag: string): seq[XmlNode] =
  for xref_node in n:
    if xref_node.kind == xnElement:
      if xref_node.tag == tag:
        result.add(xref_node)

proc quantile*(xs: seq[float], q: float): float =
  # CODE TAKEN FROM R quantile function
  # index <- 1 + (n - 1) * probs 1 + (6 - 1) 
  # lo <- floor(index)
  # hi <- ceiling(index)
  # x <- sort(x, partial = unique(c(lo, hi)))
  # qs <- x[lo]
  # i <- which(index > lo)
  # h <- (index - lo)[i] # > 0	by construction
  # ##	    qs[i] <- qs[i] + .minus(x[hi[i]], x[lo[i]]) * (index[i] - lo[i])
  # ##	    qs[i] <- ifelse(h == 0, qs[i], (1 - h) * qs[i] + h * x[hi[i]])
  # qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
  var ys = xs
  sort(ys, system.cmp[float])
  let 
    index = 1.0 + float((len(xs) - 1)) * q
    lo = int(floor(index))
    hi = int(ceil(index))
    qs = ys[lo - 1]
    h = (index - float(lo))

  result = (1.0 - h) * qs + h * ys[hi - 1]

proc IQRoutlierBounds*(xs: seq[float]): tuple[min_val: float, max_val: float] =
  ## Use InterQuartil Range procedure to remove outlier
  ## (Same procedure used in boxplot to show outliers)
  var ys = xs
  sort(ys, system.cmp[float])

  # Compute quartiles
  var
    q1 = quantile(xs, 0.25)
    q3 = quantile(xs, 0.75)
    iqr      = q3 - q1
    min_val  = float(q1) - (float(iqr) * 1.5) # Min acceptable value
    max_val  = float(q3) + (float(iqr) * 1.5) # Max acceptable value

  result.min_val = min_val
  result.max_val = max_val

proc `$`*(mc: MolecularConsequence): string =
  return mc.so_term & "|" & mc.description

proc nbStars*(rs: RevStat): int =
  case rs:
    of rsNoAssertion:
      result = 0
    of rsNoAssertionCriteria:
      result = 0
    of rsNoAssertionVariant:
      result = 0
    of rsSingleSubmitter:
      result = 1
    of rsConflicting:
      result = 1
    of rsMutlipleSubmitterNoConflict:
      result = 2
    of rsExpertPanel:
      result = 3
    of rsPracticeGuideline:
      result = 4

proc clnsigToFloat*(cs: ClinSig): float =
  case cs:
    of csBenign:
      result = 1
    of csBenignLikelyBenign:
      result = 1.5
    of csLikelyBenign:
      result = 2
    of csUncertainSignificance:
      result = 3
    of csLikelyPathogenic:
      result = 4
    of csPathogenicLikelyPathogenic:
      result = 4.5
    of csPathogenic:
      result = 5
    else:
      result = -1

proc selectACMGsubmissions(subs: seq[Submission]): seq[Submission] =
  for sub in subs:
    if sub.clinical_significance in acmg_clinsig:
      result.add(sub)

proc chromToInt*(chrom: string): int =
  ## Return -1 if chrom is not an integer (eg: X, Y)
  try:
    result = chrom.parseInt()
  except:
    result = -1

proc cmpVariant*(x, y: ClinVariant): int =
  ## Cmp based on genomic order
  let cmp_chrom = cmp(x.chrom, y.chrom)
  if cmp_chrom == 0:
    let cmp_pos = cmp(x.pos,y.pos)
    if cmp_pos == 0:
      let cmp_ref = cmp(x.ref_allele,y.ref_allele)
      if cmp_ref == 0:
        return cmp(x.alt_allele,y.alt_allele)
      else:
        return cmp_ref
    else:
      return cmp_pos
  else:
    let 
      x_chrom_int = x.chrom.chromToInt()
      y_chrom_int = y.chrom.chromToInt()
    # Both chromosomes are integer, we return the smallest one
    if x_chrom_int != -1 and y_chrom_int != -1:
      return cmp(x_chrom_int,y_chrom_int)
    # Only X ins an integer, it is the smallest one
    elif x_chrom_int != -1:
      return -1
    # Only Y is an integer, it is the smalles one
    elif y_chrom_int != -1:
      return 1
    # Neither are integer, we return the lexicographic order
    else:
      return cmp_chrom

let ncbi_conversion_regex = re(r"^Converted during submission to (.*)\.$")
proc parseNCBIConversionComment*(comment: string): ClinSig =
  ## Parse 'Converted during submission to Likely pathogenic.' to ClinSig csLikelyPathogenic, return csUnknown if parsing failed
  var arr: array[1, string]
  if match(comment, ncbi_conversion_regex, arr, 0):
    result = parseEnum[ClinSig](arr[0], csUnknown)
  else:
    result = csUnknown

proc aggregateReviewStatus*(revstat_count: TableRef[RevStat, int], nb_submitters: int, has_conflict = false): RevStat =
  if revstat_count.hasKey(rsPracticeGuideline):
    result = rsPracticeGuideline
  elif revstat_count.hasKey(rsExpertPanel):
    result = rsExpertPanel
  elif revstat_count.hasKey(rsSingleSubmitter):
    if has_conflict: 
      result = rsConflicting
    elif nb_submitters > 1 and revstat_count.hasKey(rsSingleSubmitter):
      result = rsMutlipleSubmitterNoConflict
    else:
      result = rsSingleSubmitter
  elif revstat_count.hasKey(rsNoAssertionCriteria):
    result = rsNoAssertionCriteria
  elif revstat_count.hasKey(rsNoAssertionVariant):
    result = rsNoAssertionVariant
  else:
    result = rsNoAssertion

proc selectElligibleSubmissions*(submissions: seq[Submission]): seq[Submission] =
  # Count submissions with one star or more
  var 
    has_one_star_sub : bool = false
    has_three_star_sub : bool = false
  for sub in submissions:
    let nb_stars = sub.review_status.nbStars()
    if nb_stars >= 3:
      has_three_star_sub = true
      break
    elif nb_stars >= 1:
      has_one_star_sub = true
  
  # Select eligible submissions depending on the submission with the highest number of stars
  for sub in submissions:
    # When there is a submission from an Expert panel or from a group providing practice guidelines, 
    # only the interpretation from that group is reported in the aggregate record, 
    # even if other submissions provide different interpretations. 
    # If we have submissions(s) with one star or more, we only uses this submissions in the aggregation
    # Otherwise we use all submissions
    if (has_three_star_sub and (sub.review_status == rsExpertPanel or sub.review_status == rsPracticeGuideline)) or (not has_three_star_sub and sub.review_status.nbStars >= 1) or (not has_one_star_sub and not has_three_star_sub):
      result.add(sub)

proc isConflicting*(clinsig_count: TableRef[ClinSig, int]): bool =
  result = ((clinsig_count.hasKey(csPathogenic) or clinsig_count.hasKey(csLikelyPathogenic) or clinsig_count.hasKey(csBenign) or clinsig_count.hasKey(csLikelyBenign)) and clinsig_count.hasKey(csUncertainSignificance)) or ((clinsig_count.hasKey(csPathogenic) or clinsig_count.hasKey(csLikelyPathogenic)) and (clinsig_count.hasKey(csBenign) or clinsig_count.hasKey(csLikelyBenign)))

proc countSubmissions*(submissions: seq[Submission]): tuple[clinsig_count: TableRef[ClinSig, int], revstat_count: TableRef[RevStat, int], submitter_ids: seq[int]] =
  result.clinsig_count = newTable[ClinSig, int]()
  result.revstat_count = newTable[RevStat, int]()
  result.submitter_ids = newSeq[int]()
  for sub in submissions:
    if result.clinsig_count.hasKey(sub.clinical_significance):
      inc(result.clinsig_count[sub.clinical_significance])
    else:
      result.clinsig_count[sub.clinical_significance] = 1
    
    if result.revstat_count.hasKey(sub.review_status):
      inc(result.revstat_count[sub.review_status])
    else:
      result.revstat_count[sub.review_status] = 1

    if sub.submitter_id notin result.submitter_ids:
      result.submitter_ids.add(sub.submitter_id)

proc removeOutlyingSubmissions*(acmg_submissions: seq[Submission], submissions: seq[Submission]): seq[Submission] =
  let 
    cs_values = map(acmg_submissions, proc (x: Submission): float = x.clinical_significance.clnsigToFloat())
    (min_val, max_val) = cs_values.IQRoutlierBounds()
  if min_val != -1 and max_val != -1:
    #var corrected_retained_submissions = newSeq[Submission]()
    for sub in submissions:
      if sub.clinical_significance in acmg_clinsig:
        let clnsig_float =  sub.clinical_significance.clnsigToFloat() 
        if clnsig_float >= min_val and clnsig_float <= max_val:
          result.add(sub)
      else:
        result.add(sub)

proc aggregateSubmissions*(submissions: seq[Submission], autocorrect_conflicts = false): tuple[clinsig: string, revstat: string, old_clinsig: string] =
  var
    retained_submissions = submissions.selectElligibleSubmissions()
    (clinsig_count, revstat_count, submitter_ids) = retained_submissions.countSubmissions()
  
  # Correct conflicting submissions
  let is_conflicting = clinsig_count.isConflicting()
  #
  if autocorrect_conflicts and is_conflicting:
    let acmg_submissions = retained_submissions.selectACMGsubmissions()
    if acmg_submissions.len() >= 5:
        retained_submissions = removeOutlyingSubmissions(acmg_submissions, retained_submissions)
        # Update counts with outlier removed
        (clinsig_count, revstat_count, submitter_ids) = retained_submissions.countSubmissions()
  
  # Filter acmg_only values:
  var
    nb_acmg_tags = 0
    acmg_tag : ClinSig
  
  # Need refacto
  for tag in clinsig_count.keys:
    if tag in acmg_clinsig:
      inc(nb_acmg_tags)
      acmg_tag = tag

  # Case #1, agreement between all submissions
  if nb_acmg_tags == 1 and (not is_conflicting or (is_conflicting and (acmg_tag == csPathogenic or acmg_tag == csLikelyPathogenic))):
    result.clinsig = $acmg_tag
    result.revstat = $revstat_count.aggregateReviewStatus(submitter_ids.len(), false)
  # Case #2 Patho and Likely Patho (only)
  elif nb_acmg_tags == 2 and clinsig_count.hasKey(csPathogenic) and clinsig_count.hasKey(csLikelyPathogenic):
    result.clinsig = "Pathogenic/Likely pathogenic"
    result.revstat = $revstat_count.aggregateReviewStatus(submitter_ids.len(), false)
  # Case #3, Only patho entries
  elif nb_acmg_tags == 2 and clinsig_count.hasKey(csBenign) and clinsig_count.hasKey(csLikelyBenign) and not is_conflicting:
    result.clinsig = "Benign/Likely benign"
    result.revstat = $revstat_count.aggregateReviewStatus(submitter_ids.len(), false)
  # Case #4, Conflict !!!
  # TODO: Do some desambiguiations !!!
  elif clinsig_count.isConflicting() or is_conflicting:
      result.clinsig = "Conflicting interpretations of pathogenicity"
      result.revstat = $revstat_count.aggregateReviewStatus(submitter_ids.len(), true)
  
  # Add the info that the variant was in fact reclassified
  if is_conflicting and result.clinsig != "Conflicting interpretations of pathogenicity":
    result.old_clinsig = "Conflicting interpretations of pathogenicity"
  
  # Add non-ACMG values to the end of clinsig
  # If ClinVar aggregates submissions from groups that provided a standard term not recommend by ACMG/AMP ( e.g. drug response), then those values are reported after the ACMG/AMP-based interpretation (see the table below).
  var additional_cstags : seq[string]
  for cstag in non_acmg_clinsig:
    if clinsig_count.hasKey(cstag):
      additional_cstags.add($cstag)
  if additional_cstags.len() > 0:
    # Sort the additional cstags:
    additional_cstags.sort(system.cmp[string])
    if result.clinsig == "":
      result.clinsig = additional_cstags.join(", ")
    else:
      result.clinsig.add(", " & additional_cstags.join(", "))
  
  # Handle default values
  if result.clinsig == "": 
    result.clinsig = $csUnknown
  if result.revstat == "":
    result.revstat = $revstat_count.aggregateReviewStatus(submitter_ids.len(), false)

iterator nextClinvarSet*(file: var BGZ): string =
  var chunk: string
  for line in file:
    # This is the start of a new record
    if line == "</ClinVarSet>":
      chunk.add(line & "\n")
      yield chunk
      chunk = ""
    if line.startsWith("<ClinVarSet"):
      chunk = ""
      chunk.add(line & "\n")
    elif chunk != "":
      chunk.add(line & "\n")

proc loadVariants*(clinvar_xml_file: string, genome_assembly: string): tuple[variants: TableRef[int, ClinVariant], filedate: string] =
  result.variants = newTable[int, ClinVariant]()

  var 
    file : BGZ
    submitters_hash = initTable[string, int]()
    i = 0

  file.open(clinvar_xml_file, "r")
  
  # Parse headers
  # <?xml version="1.0" encoding="UTF-8" standalone="yes"?>
  # <ReleaseSet Dated="2019-12-31" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Type="full" xsi:noNamespaceSchemaLocation="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/clinvar_public_1.59.xsd">
  for line in file:
    if line.startsWith("<ClinVarSet"):
      break
    elif line.startsWith("<ReleaseSet"):
      let 
        # Add closing ReleaseSet tag as it is at the end of the file
        new_line = line & "</ReleaseSet>\n"
        root = parseHtml(newStringStream(new_line))
      for releaseset_node in root.findNodes("releaseset"):
        result.filedate = releaseset_node.attr("Dated")
  
  file.close()
  file.open(clinvar_xml_file, "r")

  for clinvarset_string in file.nextClinvarSet():
    # TODO: Add some kind of loader ever 10K parsed variants
    if clinvarset_string != "" and clinvarset_string.startsWith("<ClinVarSet"):
      let 
        doc = q(clinvarset_string)
        reference_clinvar_assertion_nodes = doc.select("referenceclinvarassertion")
      
      if reference_clinvar_assertion_nodes.len() > 0:
        # Skipping Het-compond variants
        # TODO: We should have a more pretty way to do it as the MeasureSet in this case is bellow the GenotypeSet:
        # <GenotypeSet Type="CompoundHeterozygote" ID="424779" Acc="VCV000424779" Version="1">
        #   <MeasureSet Type="Variant" ID="928" Acc="VCV000000928" Version="2" NumberOfChromosomes="1">
        if reference_clinvar_assertion_nodes[0].select("genotypeset").len() > 0:
          continue
          
        let
          measureset_nodes = reference_clinvar_assertion_nodes[0].select("measureset")
        
        if measureset_nodes.len() > 0:
          let 
            measureset_node = measureset_nodes[0]
            variant_id = measureset_node.attr("ID").parseInt()
            measure_nodes = measureset_nodes[0].select("measure")

          # Only parse measure node to extract variant position if we do not have seen this variant
          # Already
          if not result.variants.hasKey(variant_id) and measure_nodes.len() > 0:
            let
              measure_node = measure_nodes[0] 
              measure_relationship_nodes = measure_node.select("measurerelationship")
      
            for sequence_loc in measure_node.select("sequencelocation"):
              if sequence_loc.attr("Assembly") == genome_assembly:
                # <SequenceLocation Assembly="GRCh38" AssemblyAccessionVersion="GCF_000001405.38" AssemblyStatus="current" Chr="2" Accession="NC_000002.12" start="219469373" stop="219469408" display_start="219469373" display_stop="219469408" variantLength="36" positionVCF="219469370" referenceAlleleVCF="ATGACACAGTGTACGTGTCTGGGAAGTTCCCCGGGAG" alternateAlleleVCF="A"/>
                let
                  allele_id = measure_node.attr("ID").parseInt
                  chrom = sequence_loc.attr("Chr")
                  pos_string = sequence_loc.attr("positionVCF")
                  ref_allele = sequence_loc.attr("referenceAlleleVCF")
                  alt_allele = sequence_loc.attr("alternateAlleleVCF")
                
                var
                  pos : int = -1
                
                if pos_string != "":
                  pos = pos_string.parseInt()
                              
                # Parse dbSNP rsid
                # FIXME: Use this kind of loop to replace q calls and only explore first line childs in loops !!!
                # <XRef Type="rs" ID="846664" DB="dbSNP"/>
                var rsid: int = -1
                for xref_node in measure_node:
                  if xref_node.kind == xnElement:
                    if xref_node.tag == "xref":
                      if xref_node.attr("Type") == "rs" and xref_node.attr("DB") == "dbSNP":
                        rsid = xref_node.attr("ID").parseInt()
                
                var 
                  variant = ClinVariant(
                    chrom: chrom,
                    pos: cast[int32](pos),
                    variant_id: cast[int32](variant_id),
                    allele_id: cast[int32](allele_id),
                    rsid: cast[int32](rsid),
                    ref_allele: ref_allele,
                    alt_allele: alt_allele,
                  )
                result.variants[variant_id] = variant

                # Parse Molecular Consequence
                # <AttributeSet>
                #   <Attribute Type="MolecularConsequence">missense variant</Attribute>
                #   <XRef ID="SO:0001583" DB="Sequence Ontology"/>
                #   <XRef ID="NM_000059.3:c.241T&gt;A" DB="RefSeq"/>
                # </AttributeSet>
                for attribute_set_node in measure_node.select("attributeset"):
                  let 
                    attribute_nodes = attribute_set_node.select("attribute")
                  if attribute_nodes.len > 0:
                    let attribute_node = attribute_nodes[0]
                    var
                      description : string
                      so_term: string
                    if attribute_node.attr("Type") == "MolecularConsequence":
                      description = attribute_node.innerText
                      for xref_node in attribute_set_node.select("xref"):
                        if xref_node.attr("DB") == "Sequence Ontology":
                          so_term = xref_node.attr("ID")
                      var 
                        mc_new = MolecularConsequence(description: description, so_term: so_term)
                        found_mc = false
                      # Do not add if the MC is already there TODO: Need refactor (simplified code)
                      for mc in variant.molecular_consequences:
                        if mc.so_term == mc_new.so_term:
                          found_mc = true
                      if not found_mc:
                        variant.molecular_consequences.add(mc_new)

                break # We found our "sequenceLocation"

          # Not lets add the submissions
          if result.variants.hasKey(variant_id):
            for clinvar_assertion_node in doc.select("clinvarassertion"):
              let 
                clinsig_nodes = clinvar_assertion_node.select("clinicalsignificance")
                clinvar_submission_id_nodes = clinvar_assertion_node.select("clinvarsubmissionid")
                measure_relationship_nodes = clinvar_assertion_node.select("measurerelationship")

              var 
                submitter_id = -1
              
              # Extract Submitter ID
              if clinvar_submission_id_nodes.len() > 0:
                let submitter_name = clinvar_submission_id_nodes[0].attr("submitter")
                if submitters_hash.hasKey(submitter_name):
                  submitter_id = submitters_hash[submitter_name]
                else:
                  # Add the new submitter to the submutter hash
                  submitter_id = submitters_hash.len()
                  submitters_hash[submitter_name] = submitter_id
                  
              if clinsig_nodes.len() > 0: # FIXME: Should not be > to 1 ...
                var
                  clinical_significance : ClinSig = csUnknown
                  review_status : RevStat = rsNoAssertion
                
                if clinsig_nodes.len() > 0:
                  var
                    desc_nodes = clinsig_nodes[0].select("description")
                    revstat_nodes = clinsig_nodes[0].select("reviewstatus")
                    comment_nodes = clinsig_nodes[0].select("comment")
                  # <ClinicalSignificance>
                  #   <ReviewStatus>no assertion criteria provided</ReviewStatus>
                  #   <Description>likely pathogenic - adrenal pheochromocytoma</Description>
                  #   <Comment Type="ConvertedByNCBI">Converted during submission to Likely pathogenic.</Comment>
                  # </ClinicalSignificance>
                  # Here we handle the clinsig that were automatically converted by NCBI and has to be
                  # extracted with a regex from the comment node
                  for comment in comment_nodes:
                    let parse_clnsig = parseNCBIConversionComment(comment.innerText)
                    if parse_clnsig != csUnknown:
                      clinical_significance = parse_clnsig
                  if clinical_significance == csUnknown and desc_nodes.len() > 0:
                    clinical_significance = parseEnum[ClinSig](desc_nodes[0].innerText, csUnknown)
                  if revstat_nodes.len() > 0:
                   review_status = parseEnum[RevStat](revstat_nodes[0].innerText, rsNoAssertion)
                   
                  # Add the submission to the variant record
                  result.variants[variant_id].submissions.add(Submission(
                    clinical_significance: clinical_significance,
                    review_status: review_status,
                    submitter_id: submitter_id
                  ))

proc formatVCFString*(vcf_string: string): string =
  result = vcf_string.replace(' ', '_')

proc printVCF*(variants: seq[ClinVariant], genome_assembly: string, filedate: string, genes_index: TableRef[string, Lapper[GFFGene]], coding_priority : bool) =
  echo "##fileformat=VCFv4.1"
  if filedate != "":
    echo "##fileDate=" & filedate # TODO: Get date from XML headers
  echo "##source=ClinVar"
  echo "##reference=" & genome_assembly
  echo "##ID=<Description=\"ClinVar Variation ID\">"
  # ##INFO=<ID=AF_ESP,Number=1,Type=Float,Description="allele frequencies from GO-ESP">
  # ##INFO=<ID=AF_EXAC,Number=1,Type=Float,Description="allele frequencies from ExAC">
  # ##INFO=<ID=AF_TGP,Number=1,Type=Float,Description="allele frequencies from TGP">
  echo "##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description=\"the ClinVar Allele ID\">"
  # ##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
  # ##INFO=<ID=CLNDNINCL,Number=.,Type=String,Description="For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB">
  # ##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN">
  # ##INFO=<ID=CLNDISDBINCL,Number=.,Type=String,Description="For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN">
  # ##INFO=<ID=CLNHGVS,Number=.,Type=String,Description="Top-level (primary assembly, alt, or patch) HGVS expression.">
  echo "##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description=\"ClinVar review status for the Variation ID\">"
  echo "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical significance for this single variant\">"
  echo "##INFO=<ID=OLD_CLNSIG,Number=.,Type=String,Description=\"Clinical significance was deciphered and this value is the original one given by ClinVar aggregation method\">"
  # ##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Conflicting clinical significance for this single variant">
  # ##INFO=<ID=CLNSIGINCL,Number=.,Type=String,Description="Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance.">
  # ##INFO=<ID=CLNVC,Number=1,Type=String,Description="Variant type">
  # ##INFO=<ID=CLNVCSO,Number=1,Type=String,Description="Sequence Ontology id for variant type">
  # ##INFO=<ID=CLNVI,Number=.,Type=String,Description="the variant's clinical sources reported as tag-value pairs of database and variant identifier">
  # ##INFO=<ID=DBVARID,Number=.,Type=String,Description="nsv accessions from dbVar for the variant">
  echo "##INFO=<ID=GENEINFO,Number=1,Type=String,Description=\"Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)\">"
  echo "##INFO=<ID=MC,Number=.,Type=String,Description=\"comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence\">"
  # ##INFO=<ID=ORIGIN,Number=.,Type=String,Description="Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 3
  # 2 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other">
  echo "##INFO=<ID=RS,Number=.,Type=String,Description=\"dbSNP ID (i.e. rs number)\">"
  # ##INFO=<ID=SSR,Number=1,Type=Integer,Description="Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 
  # - 1kg_failed, 1024 - other">
  echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  
  var 
    nb_corrections = 0
    nb_variants = 0
  for v in variants:
    inc(nb_variants)
    let (clinsig, revstat, old_clinsig) = v.submissions.aggregateSubmissions(true) # Autocorrect conflicts
    var info_fields : seq[string] = @["ALLELEID=" & $v.allele_id]
    
    if not genes_index.isNil():
      let gene_info = genes_index.getInfoString(v.chrom, int(v.pos), int(v.pos) + v.ref_allele.len() - 1, coding_priority)
      if gene_info != "":
        info_fields.add("GENEINFO=" & gene_info)

    info_fields.add("CLNSIG=" & clinsig.formatVCFString())
    info_fields.add("CLNREVSTAT=" & revstat.formatVCFString())

    if old_clinsig != "":
      info_fields.add("OLD_CLNSIG=" & old_clinsig.formatVCFString())
      inc(nb_corrections)

    if v.molecular_consequences.len > 0:
      var formated_consequences = map(v.molecular_consequences, proc (x: MolecularConsequence): string = formatVCFString($x))
      info_fields.add("MC=" & join(formated_consequences, ","))
    
    if v.rsid != -1:
      info_fields.add("RS=" & $v.rsid)
  
    if v.chrom != "" and v.ref_allele != "" and v.alt_allele != "":
      echo [
        v.chrom,
        $v.pos,
        $v.variant_id,
        v.ref_allele,
        v.alt_allele,
        ".",
        ".",
        info_fields.join(";")
      ].join("\t")

  stderr.writeLine("[Log] " & $nb_variants & " variants have been extracted from the XML")
  stderr.writeLine("[Log] " & $nb_corrections & " variants had a conflicting interpretation deciphering")

proc main*(argv: seq[string]) =

  # TODO: Create a usage and expose api_keys as options
  let doc = format("""
Usage: clinvcf [options] <clinvar.xml.gz>

Options:
  --genome <version>              Genome assembly to use [default: GRCh37]
  --gff <file>                    NCBI GFF to annotate variations with genes
  --coding-first                  Give priority to coding gene in annotation (even if intronic and exonic for another gene)
  """)

  let 
    args = docopt(doc)
    genome_assembly = $args["--genome"]
    clinvar_xml_file = $args["<clinvar.xml.gz>"]
    coding_priority = args["--coding-first"]
  
  var 
    variants_hash: TableRef[int, ClinVariant]
    variants_seq: seq[ClinVariant]
    filedate: string
    genes_index: TableRef[string, Lapper[GFFGene]]

  # Load variants from XML
  stderr.writeLine("[Log] Parsing variants from " & clinvar_xml_file)
  (variants_hash, filedate) = loadVariants(clinvar_xml_file, genome_assembly)

  if args["--gff"]:
    let gff_file = $args["--gff"]
    stderr.writeLine("[Log] Load genes coordinates from " & gff_file)
    genes_index = loadGenesFromGFF(gff_file)
  
  # Sort variants by genomic order
  stderr.writeLine("[Log] Sorting variants")
  variants_seq = toSeq(variants_hash.values())
  variants_seq.sort(cmpVariant)
  
  # Print VCF of STDOUT
  stderr.writeLine("[Log] Printing variants")
  printVCF(variants_seq, genome_assembly, filedate, genes_index, coding_priority)

when isMainModule:
  main(commandLineParams())