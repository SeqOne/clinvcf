import httpclient, json, algorithm, tables, sequtils, re
import os, times
import q, xmltree # Parse XML
import docopt # Formating the command-line
import strutils # Split string
import hts

type
  ClinSig* = enum
    csBenign = "Benign",
    csLikelyBenign = "Likely benign",
    csUncertainSignificance = "Uncertain significance",
    csLikelyPathogenic = "Likely pathogenic",
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
    gene_id: int32
    gene_symbol: string
    molecular_consequences: seq[MolecularConsequence]
    submissions: seq[Submission]

var
  acmg_clinsig = @[csBenign, csLikelyBenign, csUncertainSignificance, csLikelyPathogenic, csPathogenic]
  non_acmg_clinsig = @[csDrugResponse, csRiskFactor, csAffects, csAssociation, csProtective, csConflictingDataFromSubmitters, csOther]

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

proc aggregateSubmissions*(submissions: seq[Submission]): tuple[clinsig: string, revstat: string] =
  var 
    clinsig_count = newTable[ClinSig, int]()
    revstat_count = newTable[RevStat, int]()
    submitter_ids = newSeq[int]()
    at_least_one_star_subs: seq[Submission]

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
      if clinsig_count.hasKey(sub.clinical_significance):
        inc(clinsig_count[sub.clinical_significance])
      else:
        clinsig_count[sub.clinical_significance] = 1
      
      if revstat_count.hasKey(sub.review_status):
        inc(revstat_count[sub.review_status])
      else:
        revstat_count[sub.review_status] = 1
      if sub.submitter_id notin submitter_ids:
        submitter_ids.add(sub.submitter_id)
  
  # Filter acmg_only values:
  var
    nb_acmg_tags = 0
    acmg_tag : ClinSig

  for tag in clinsig_count.keys:
    if tag in acmg_clinsig:
      inc(nb_acmg_tags)
      acmg_tag = tag

  # Case #1, agreement between all submissions
  if nb_acmg_tags == 1:
    result.clinsig = $acmg_tag
    result.revstat = $revstat_count.aggregateReviewStatus(submitter_ids.len(), false)
  # Case #2 Patho and Likely Patho (only)
  elif nb_acmg_tags == 2 and clinsig_count.hasKey(csPathogenic) and clinsig_count.hasKey(csLikelyPathogenic):
    result.clinsig = "Pathogenic/Likely pathogenic"
    result.revstat = $revstat_count.aggregateReviewStatus(submitter_ids.len(), false)
  # Case #3, Only patho entries
  elif nb_acmg_tags == 2 and clinsig_count.hasKey(csBenign) and clinsig_count.hasKey(csLikelyBenign):
    result.clinsig = "Benign/Likely benign"
    result.revstat = $revstat_count.aggregateReviewStatus(submitter_ids.len(), false)
  # Case #4, Conflict !!!
  # TODO: Do some desambiguiations !!!
  elif (clinsig_count.hasKey(csPathogenic) or clinsig_count.hasKey(csLikelyPathogenic) or clinsig_count.hasKey(csBenign) or clinsig_count.hasKey(csLikelyBenign)) and clinsig_count.hasKey(csUncertainSignificance):
    result.clinsig = "Conflicting interpretations of pathogenicity"
    result.revstat = $revstat_count.aggregateReviewStatus(submitter_ids.len(), true)
  elif (clinsig_count.hasKey(csPathogenic) or clinsig_count.hasKey(csLikelyPathogenic)) and (clinsig_count.hasKey(csBenign) or clinsig_count.hasKey(csLikelyBenign)):
    result.clinsig = "Conflicting interpretations of pathogenicity"
    result.revstat = $revstat_count.aggregateReviewStatus(submitter_ids.len(), true)
  
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
    if line == "":
      yield chunk
      chunk = ""
    else:
      chunk.add(line & "\n")
  yield chunk

proc loadVariants*(clinvar_xml_file: string, genome_assembly: string): tuple[variants: TableRef[int, ClinVariant], filedate: string] =
  result.variants = newTable[int, ClinVariant]()

  var 
    file : BGZ
    submitters_hash = initTable[string, int]()
    i = 0

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
                  gene_symbol: string
                  gene_id : int = -1
                  pos : int = -1
                
                if pos_string != "":
                  pos = pos_string.parseInt()
                
                if measure_relationship_nodes.len() > 0:
                  let gene_symbol_nodes = measure_relationship_nodes[0].select("symbol elementvalue")
                  if gene_symbol_nodes.len() > 0:
                    gene_symbol = gene_symbol_nodes[0].innerText
                  for xref_node in measure_relationship_nodes[0].select("xref"):
                    if xref_node.attr("DB") == "Gene":
                      gene_id = xref_node.attr("ID").parseInt()
                      # <XRef ID="672" DB="Gene" />
                      # <XRef Type="MIM" ID="113705" DB="OMIM" />
                              
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
                    gene_id: cast[int32](gene_id),
                    gene_symbol: gene_symbol
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

              var 
                submitter_id = -1

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
                    submitter_id: submitter_id,
                  ))
    # THese are the headers to be parsed to retrieve the date
    # and we want to avoid the closing </ReleaseSet> at the end of the file
    elif not clinvarset_string.startsWith("</"):
      let 
        doc = q(clinvarset_string)
        releaseset_nodes = doc.select("releaseset")
      if releaseset_nodes.len() > 0:
        result.filedate = releaseset_nodes[0].attr("Dated")

proc formatVCFString*(vcf_string: string): string =
  result = vcf_string.replace(' ', '_')

proc printVCF*(variants: seq[ClinVariant], genome_assembly: string, filedate: string) =
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

  for v in variants:
    let (clinsig, revstat) = v.submissions.aggregateSubmissions()
    var info_fields : seq[string] = @["ALLELEID=" & $v.allele_id]
    
    if v.gene_id > 0:
      info_fields.add("GENEINFO=" & v.gene_symbol & ":" & $v.gene_id)

    info_fields.add("CLNSIG=" & clinsig.formatVCFString())
    info_fields.add("CLNREVSTAT=" & revstat.formatVCFString())

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

proc main*(argv: seq[string]) =

  # TODO: Create a usage and expose api_keys as options
  let doc = format("""
Usage: clinvcf [options] <clinvar.xml.gz>

Options:
  --genome <version>              Genome assembly to use [default: GRCh37]
  """)

  let 
    args = docopt(doc)
    genome_assembly = $args["--genome"]
    clinvar_xml_file = $args["<clinvar.xml.gz>"]
  
  var 
    variants_hash: TableRef[int, ClinVariant]
    variants_seq: seq[ClinVariant]
    filedate: string

  # Load variants from XML
  stderr.writeLine("[Log] Parsing variants from " & clinvar_xml_file)
  (variants_hash, filedate) = loadVariants(clinvar_xml_file, genome_assembly)
  
  # Sort variants by genomic order
  stderr.writeLine("[Log] Sorting variants")
  variants_seq = toSeq(variants_hash.values())
  variants_seq.sort(cmpVariant)
  
  # Print VCF of STDOUT
  stderr.writeLine("[Log] Printing variants")
  printVCF(variants_seq, genome_assembly, filedate)

when isMainModule:
  main(commandLineParams())