import algorithm, tables, sequtils, re, math
import os
import htmlparser
import q, xmltree # Parse XML
from streams import newStringStream
import docopt # Formating the command-line
import strutils # Split string
import strformat
import hts
from regex import match, RegexMatch, groupFirstCapture
import logging

# Local libs
from ./clinvcfpkg/utils import logger
import ./clinvcfpkg/gff
import ./clinvcfpkg/lapper
import ./clinvcfpkg/hgnc

type
  ClinSig* = enum
    csNA = "",
    csBenign = "Benign",
    csBenignLikelyBenign = "Benign/Likely benign",
    csLikelyBenign = "Likely benign",
    csUncertainSignificance = "Uncertain significance",
    csLikelyPathogenic = "Likely pathogenic",
    csPathogenicLikelyPathogenic = "Pathogenic/Likely pathogenic",
    csPathogenic = "Pathogenic",
    csUnknown = "not provided",
    csDrugResponse = "drug response",
    csRiskFactor = "risk factor",
    csAffects = "Affects",
    csAssociation = "association",
    csProtective = "protective",
    csConflictingDataFromSubmitters = "conflicting data from submitters",
    csOther = "other",
    csConflictingInterpretation = "Conflicting interpretations of pathogenicity"

  RevStat* = enum
    rsNA = "",
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
    variant_in_gene: string

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
    pathologies: TableRef[string, seq[string]]

const ignoredPathoTag = @["not specified", "see cases", "not provided", "variant of unknown significance"]

var
  acmg_clinsig = @[csBenign, csLikelyBenign, csUncertainSignificance, csLikelyPathogenic, csPathogenic]
  non_acmg_clinsig = @[
    csDrugResponse, csRiskFactor, csAffects, csAssociation, csProtective, csConflictingDataFromSubmitters, csOther
  ]
  clinicalPathoType: seq[string] = @[]

method correctGeneAlias(cv: ClinVariant, hgnc: HgncIndex) =
  ## Loop on submission and switch alias gene to official gene name
  for i, submission in cv.submissions:
    if hgnc.alias.hasKey(submission.variant_in_gene):
      var newGeneName = hgnc.alias[submission.variant_in_gene]
      submission.variant_in_gene = newGeneName
      cv.submissions[i] = submission
    elif not hgnc.entrez.hasKey(submission.variant_in_gene):
      # logger.log(lvlInfo, fmt"{submission.variant_in_gene}")
      submission.variant_in_gene = ""
      cv.submissions[i] = submission

proc findNodes(n: XmlNode, tag: string): seq[XmlNode] =
  for xref_node in n:
    if xref_node.kind == xnElement:
      if xref_node.tag == tag:
        result.add(xref_node)

proc quantile*(xs: seq[float], q: float): float =
  ## CODE TAKEN FROM R quantile function
  ## This correspond to the default implemented method for quantile calculation
  ## index <- 1 + (n - 1) * probs 1 + (6 - 1)
  ## lo <- floor(index)
  ## hi <- ceiling(index)
  ## x <- sort(x, partial = unique(c(lo, hi)))
  ## qs <- x[lo]
  ## i <- which(index > lo)
  ## h <- (index - lo)[i] # > 0	by construction
  ##	    qs[i] <- qs[i] + .minus(x[hi[i]], x[lo[i]]) * (index[i] - lo[i])
  ##	    qs[i] <- ifelse(h == 0, qs[i], (1 - h) * qs[i] + h * x[hi[i]])
  ## qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
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
    of rsNA:
      result = 0
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
  ## Parse 'Converted during submission to Likely pathogenic.' to ClinSig csLikelyPathogenic,
  ## return csUnknown if parsing failed
  var arr: array[1, string]
  if match(comment, ncbi_conversion_regex, arr, 0):
    result = parseEnum[ClinSig](arr[0], csUnknown)
  else:
    result = csUnknown

proc parseClinicalPathologies*(pathoType: string, pathoList: seq[string]): string =
  ## Take a patho type (disease, finding etc) and return a formatted INFO field with all 
  ## pathologies associated to a variant for a specific type
  ## EXAMPLE FINDING=PATHO1|PATHO2|PATHO3 ...
  result = "CLN" & pathoType.toUpperAscii & "="
  for i, patho in pathoList:
    # pred() gives the len - 1 value
    if i == pathoList.len.pred:
      result = result & patho
    else:
      result = result & patho & '|'

proc aggregateReviewStatus*(revstat_count: TableRef[RevStat, int], nb_submitters: int,
  has_conflict = false): RevStat =
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
    var cond = (
      (
        has_three_star_sub and (
          sub.review_status == rsExpertPanel or sub.review_status == rsPracticeGuideline
        )
      ) or (
        not has_three_star_sub and sub.review_status.nbStars >= 1
      ) or (
        not has_one_star_sub and not has_three_star_sub
      )
    )
    if cond:
      result.add(sub)

proc isConflicting*(clinsig_count: TableRef[ClinSig, int]): bool =
  result = (
    (
      (
        clinsig_count.hasKey(csPathogenic) or clinsig_count.hasKey(csLikelyPathogenic) or 
        clinsig_count.hasKey(csBenign) or clinsig_count.hasKey(csLikelyBenign)
      ) and clinsig_count.hasKey(csUncertainSignificance)
    ) or (
      (
        clinsig_count.hasKey(csPathogenic) or clinsig_count.hasKey(csLikelyPathogenic)
      ) and (
        clinsig_count.hasKey(csBenign) or clinsig_count.hasKey(csLikelyBenign)
      )
    )
  )

proc countSubmissions*(submissions: seq[Submission]): tuple[clinsig_count: TableRef[ClinSig, int],
  revstat_count: TableRef[RevStat, int], submitter_ids: seq[int]] =
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

proc removeOutlyingSubmissions*(submissions: seq[Submission]): seq[Submission] =
  let
    acmg_submissions = submissions.selectACMGsubmissions()
    cs_values = map(acmg_submissions, proc (x: Submission): float = x.clinical_significance.clnsigToFloat())
    (min_val, max_val) = cs_values.IQRoutlierBounds()
  if min_val != -1 and max_val != -1:
    #var corrected_retained_submissions = newSeq[Submission]()
    for sub in submissions:
      if sub.clinical_significance in acmg_clinsig:
        let clnsig_float = sub.clinical_significance.clnsigToFloat()
        if clnsig_float >= min_val and clnsig_float <= max_val:
          result.add(sub)
      else:
        result.add(sub)

proc aggregatSubmissionsClinvar*(submissions: seq[Submission]): tuple[clinsig: ClinSig, revstat: RevStat] =
  result.clinsig = csNA
  result.revstat = rsNA
  # Update counts with outlier removed
  let
    (clinsig_count, revstat_count, submitter_ids) = submissions.countSubmissions()
    is_conflicting = clinsig_count.isConflicting()
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
  if nb_acmg_tags == 1:
    result.clinsig = acmg_tag
    result.revstat = revstat_count.aggregateReviewStatus(submitter_ids.len(), false)
  # Case #2 Patho and Likely Patho (only)
  elif nb_acmg_tags == 2 and clinsig_count.hasKey(csPathogenic) and clinsig_count.hasKey(csLikelyPathogenic):
    result.clinsig = csPathogenicLikelyPathogenic
    result.revstat = revstat_count.aggregateReviewStatus(submitter_ids.len(), false)
  # Case #3, Only patho entries
  elif nb_acmg_tags == 2 and clinsig_count.hasKey(csBenign) and clinsig_count.hasKey(csLikelyBenign):
    result.clinsig = csBenignLikelyBenign
    result.revstat = revstat_count.aggregateReviewStatus(submitter_ids.len(), false)
  # Case #4, Conflict !!!
  elif is_conflicting:
    result.clinsig = csConflictingInterpretation
    result.revstat = revstat_count.aggregateReviewStatus(submitter_ids.len(), true)

proc aggregateVariantInGene(submissions: seq[Submission], hgnc: HgncIndex): string =
  result = ""
  var
    genes: seq[string] = @[]
    genesToReturn: seq[string] = @[]
  for submission in submissions:
    if submission.variant_in_gene != "":
      if submission.variant_in_gene notin genes:
        genes.add(submission.variant_in_gene)
        try:
          genesToReturn.add(fmt"{submission.variant_in_gene}:{hgnc.entrez[submission.variant_in_gene]}")
        except KeyError:
          logger.log(lvlInfo, fmt"{submission.variant_in_gene} is not found in hgnc table.")

  if len(genesToReturn) > 0:
    result = genesToReturn.join("|")
  
proc aggregateSubmissions*(submissions: seq[Submission], hgncIndex: HgncIndex,
  autocorrect_conflicts = false): tuple[clinsig: string, revstat: string, old_clinsig: string,
  nb_reclassification_stars: int, geneInfo: string] =

  var
    retained_submissions = submissions.selectElligibleSubmissions()
    (clinsig, revstat) = retained_submissions.aggregatSubmissionsClinvar()
    (clinsig_count, revstat_count, submitter_ids) = retained_submissions.countSubmissions()

  result.geneInfo = retained_submissions.aggregateVariantInGene(hgncIndex)
  
  if clinsig != csNA:
    result.clinsig = $clinsig
  
  if revstat != rsNA:
    result.revstat = $revstat

  result.nb_reclassification_stars = -1

  # Correct conflicting submissions
  if autocorrect_conflicts and clinsig == csConflictingInterpretation:
    let acmg_submissions = retained_submissions.selectACMGsubmissions()
    if acmg_submissions.len() >= 4:
      let
        retained_submissions_without_outliers = removeOutlyingSubmissions(retained_submissions)
        (
          clinsig_without_outliers, revstat_without_outliers
        ) = retained_submissions_without_outliers.aggregatSubmissionsClinvar()
        # (
        #   clinsig_count_without_outliers, revstat_count_without_outliers, submitter_ids_without_outliers
        # ) = retained_submissions_without_outliers.countSubmissions()
        (clinsig_count_without_outliers, _, _) = retained_submissions_without_outliers.countSubmissions()
      # We were "conflicting" but re-assigned the clinsig to a pathogenic tag after outlier removal
      # Now we try to see if we have 1, 2 or 3 stars for reclassification
      # - 1 star : default
      # - 2 stars : reclassification remains even if we add a virtual VUS submission
      # - 3 stars : 2 stars requirements and at least 1 pathogenic classification
      if clinsig_without_outliers.clnsigToFloat() >= 4:
        var
          submissions_with_one_vus = retained_submissions
          nb_reclassification_stars = 1
        submissions_with_one_vus.add(
          Submission(clinical_significance: csUncertainSignificance, review_status: rsSingleSubmitter)
        )
        let
          #acmg_submissions_with_one_vus = submissions_with_one_vus.selectACMGsubmissions()
          submissions_with_one_vus_without_outlier = removeOutlyingSubmissions(submissions_with_one_vus)
          (
            clinsig_with_one_vus, revstat_with_one_vus
          ) = submissions_with_one_vus_without_outlier.aggregatSubmissionsClinvar()

        # Debug lines for this scary code section
        # stderr.writeLine("[Log] submissions: " & $map(submissions, proc (x: Submission): string = $x.clinical_significance.clnsigToFloat()).join(","))
        # stderr.writeLine("[Log] retained_submissions: " & $map(retained_submissions, proc (x: Submission): string = $x.clinical_significance.clnsigToFloat()).join(","))
        # stderr.writeLine("[Log] retained_submissions_without_outliers: " & $map(retained_submissions_without_outliers, proc (x: Submission): string = $x.clinical_significance.clnsigToFloat()).join(","))
        # stderr.writeLine("[Log] submissions_with_one_vus: " & $map(submissions_with_one_vus, proc (x: Submission): string = $x.clinical_significance.clnsigToFloat()).join(","))
        # stderr.writeLine("[Log] submissions_with_one_vus_without_outlier: " & $map(submissions_with_one_vus_without_outlier, proc (x: Submission): string = $x.clinical_significance.clnsigToFloat()).join(","))

        if clinsig_with_one_vus.clnsigToFloat() >= 4:
          if clinsig_count_without_outliers.hasKey(csPathogenic):
            nb_reclassification_stars = 3
          else:
            nb_reclassification_stars = 2
          # If we have 2/3 stars we use the classification with the virtual "VUS"
          result.clinsig = $clinsig_with_one_vus
          result.revstat = $revstat_with_one_vus
        else:
          # Otherwise we keep our original classification as the virtual "VUS" got us back to conflicting status
          result.clinsig = $clinsig_without_outliers
          result.revstat = $revstat_without_outliers
        result.old_clinsig = $csConflictingInterpretation
        result.nb_reclassification_stars = nb_reclassification_stars

  # Add non-ACMG values to the end of clinsig
  # If ClinVar aggregates submissions from groups that provided a standard term not recommend by ACMG/AMP ( e.g. drug response), 
  # then those values are reported after the ACMG/AMP-based interpretation (see the table below).
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
  # <ReleaseSet Dated="2019-12-31" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
  # Type="full" xsi:noNamespaceSchemaLocation="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/clinvar_public_1.59.xsd">
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

  var clinvarSetCount = 0
  for clinvarset_string in file.nextClinvarSet():
    inc(clinvarSetCount)
    var moduloClinvarSetCount = clinvarSetCount mod 100000
    if moduloClinvarSetCount == 0:
      logger.log(lvlInfo, fmt"[main] clinvarSet={clinvarSetCount}")
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

        let measureset_nodes = reference_clinvar_assertion_nodes[0].select("measureset")

        if measureset_nodes.len() > 0:
          let
            measureset_node = measureset_nodes[0]
            variant_id = measureset_node.attr("ID").parseInt()
            measure_nodes = measureset_nodes[0].select("measure")

          # Only parse "variant" and skip "Haplotype"
          if measureset_node.attr("Type") != "Variant":
            continue

          # Only parse measure node to extract variant position if we do not have seen this variant
          # Already
          if not result.variants.hasKey(variant_id) and measure_nodes.len() > 0:
            let measure_node = measure_nodes[0]
              # if measure_relationship_nodes.attr("Type") == "variant in gene":
                # element_value_node = measure_relationship_nodes.select("symbol").select("ElementValue")
            for sequence_loc in measure_node.select("sequencelocation"):
              if sequence_loc.attr("Assembly") == genome_assembly:
                # <SequenceLocation Assembly="GRCh38" AssemblyAccessionVersion="GCF_000001405.38" AssemblyStatus="current"
                # Chr="2" Accession="NC_000002.12" start="219469373" stop="219469408" display_start="219469373"
                # display_stop="219469408" variantLength="36" positionVCF="219469370" referenceAlleleVCF="ATGACACAGTGTACGTGTCTGGGAAGTTCCCCGGGAG"
                # alternateAlleleVCF="A"/>
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
                    pathologies: newTable[string, seq[string]]()
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

          # Now lets add the submissions and pathology
          if result.variants.hasKey(variant_id):
            for clinvar_assertion_node in doc.select("clinvarassertion"):
              let
                clinsig_nodes = clinvar_assertion_node.select("clinicalsignificance")
                clinvar_submission_id_nodes = clinvar_assertion_node.select("clinvarsubmissionid")
                measure_relationship_nodes = clinvar_assertion_node.select("measurerelationship")
                traitSetNodes = clinvar_assertion_node.select("traitset")

              # Extract gene
              var variant_in_gene = ""
              for measure_relationship in measure_relationship_nodes:
                if measure_relationship.attr("Type") == "variant in gene":
                  for symbol in measure_relationship.select("symbol"):
                    for elementvalue in symbol.select("elementvalue"):
                      if elementvalue.attr("Type") == "Preferred":
                        if len(elementvalue) > 0:
                          # Gene is present in the submission
                          variant_in_gene = elementvalue[0].innerText

              # Extract pathologies
              for trait in traitSetNodes[0].select("trait"):
                # Check if there is informations about pathology
                if trait.select("elementvalue").len() > 0:
                  # <traitset Type="Finding">
                  #   <trait Type="Finding">
                  #     <name>
                  #       <elementvalue Type="Preferred">nuclear cataracts</elementvalue>
                  #     </name>
                  #   </trait>
                  #   <trait Type="Finding">
                  #     <name>
                  #       <elementvalue Type="Preferred">microcornea</elementvalue>
                  #     </name>
                  #   </trait>
                  # </traitset>
                  # <traitset Type="Disease">
                  #   <trait Type="Disease">
                  #    <name>
                  #      <elementvalue Type="Preferred">Cataract 1</elementvalue>
                  #    </name>
                  #    <xref DB="OMIM" Type="MIM" ID="116200" />
                  #   </trait>
                  # </traitset>
                  let pathoType = trait.attr("Type")
                  # Patho to skip : all values in ingnoredPathoTag
                  if trait.select("elementvalue")[0].innerText.toLowerAscii in ignoredPathoTag:
                    continue
                  let pathology = trait.select("elementvalue")[0].innerText.toLowerAscii
                  # Add result inside pathology table:
                  #   key: pathology type : (Disease, Finding...)
                  #   value: a list contaning pathology's names
                  if result.variants[variant_id].pathologies.hasKey(pathoType):
                    if pathology in result.variants[variant_id].pathologies[pathoType]:
                      continue
                    result.variants[variant_id].pathologies[pathoType].add(pathology)
                  else:
                    result.variants[variant_id].pathologies[pathoType] = @[]
                    result.variants[variant_id].pathologies[pathoType].add(pathology)
                  if pathoType in clinicalPathoType:
                    continue
                  clinicalPathoType.add(pathoType)

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
                  var submission = Submission(
                    clinical_significance: clinical_significance,
                    review_status: review_status,
                    submitter_id: submitter_id,
                    variant_in_gene: variant_in_gene
                  )                 
                  result.variants[variant_id].submissions.add(submission)

proc formatVCFString*(vcf_string: string): string =
  result = vcf_string.replace(' ', '_')

proc formatPathoString*(pathoString: string): string =
  # First : remove all non-word char after "="
  # e.g CLNDISEASE= cancer|
  result = pathoString.replace(re"=\W+", "=")
  # Second : remove all non-word chars at the end of the string (or first)
  # and non-words after/before pipe ('|')
  result = result.replace(re"^\W+|\W+$|\W+?(?=\|)|(?<=\|)\W+")
  # Third replace all non-word characters by '_'
  result = result.replace(re"[^\w\||=]+", "_")

proc printVCF*(variants: seq[ClinVariant], genome_assembly: string, filedate: string,
  genes_index: TableRef[string, Lapper[GFFGene]], coding_priority : bool, hgncIndex: HgncIndex) =
  # Commented lines correspond to the NCBI Clinvar original header that are not currently supported by clinVCF
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
  echo "##INFO=<ID=CLNRECSTAT,Number=1,Type=Integer,Description=\"3-levels stars confidence for automatic reclassfication of conflicting variants\">"
  for pathoType in clinicalPathoType:
    echo "##INFO=<ID=CLN" & pathoType.toUpperAscii & ",Number=.,Type=String,Description=\"Clinical pathology(ies) ranked as " & pathoType & " referenced for a variant\">"
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
    let (clinsig, revstat, old_clinsig, nb_reclassification_stars, rawGeneInfo) = v.submissions.aggregateSubmissions(
      hgncIndex,
      true
    ) # Autocorrect conflicts
    var info_fields : seq[string] = @["ALLELEID=" & $v.allele_id]

    if not genes_index.isNil() and rawGeneInfo == "":
      let gene_info = genes_index.getInfoString(v.chrom, int(v.pos), int(v.pos) + v.ref_allele.len() - 1, coding_priority)
      if gene_info != "":
        info_fields.add("GENEINFO=" & gene_info)
    
    elif rawGeneInfo != "":
      info_fields.add("GENEINFO=" & rawGeneInfo)

    info_fields.add("CLNSIG=" & clinsig.formatVCFString())
    info_fields.add("CLNREVSTAT=" & revstat.formatVCFString())

    if old_clinsig != "":
      info_fields.add("OLD_CLNSIG=" & old_clinsig.formatVCFString())
      info_fields.add("CLNRECSTAT=" & $nb_reclassification_stars)
      inc(nb_corrections)

    # Add each pathologies per variants in info_field
    if v.pathologies.len > 0:
      for pathoType in v.pathologies.keys:
        # Edit pathology string to fit with header
        # Example : given Disease return CLNDISEASE
        let pathology = pathoType.parseClinicalPathologies(v.pathologies[pathoType])
        info_fields.add(pathology.formatPathoString)

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
        if v.ref_allele == v.alt_allele: "." else: v.alt_allele,
        ".",
        ".",
        info_fields.join(";")
      ].join("\t")

  logger.log(lvlInfo, fmt"[printVCF] {nb_variants} variants have been extracted from the XML")
  logger.log(lvlInfo, fmt"[printVCF] {nb_corrections} variants had a conflicting interpretation deciphering")

proc main*(argv: seq[string]) =

  # TODO: Create a usage and expose api_keys as options
  let doc = format("""
Usage: clinvcf [options] --hgnc <table> --genome <version> <clinvar.xml.gz>

Arguments:
  --genome <version>              Genome assembly to use
  --hgnc <table>                  HGNC table used for gene name alias corrections

Options:
  --filename-date                 Use xml filename date instead of inner date which may differ

Gene annotation:
  --gff <file>                    NCBI GFF to annotate variations with genes
  --coding-first                  Give priority to coding gene in annotation (even if intronic and exonic for another gene)
  --gene-padding <int>            Padding to annotation upstream/downstream genes (not applied for MT) [default: 5000]
  """)

  let
    args = docopt(doc)
    genome_assembly = $args["--genome"]
    clinvar_xml_file = $args["<clinvar.xml.gz>"]
    coding_priority = args["--coding-first"]
    gene_padding = parseInt($args["--gene-padding"])
    filename_date = args["--filename-date"]

  var
    variants_hash: TableRef[int, ClinVariant]
    variants_seq: seq[ClinVariant]
    filedate: string
    genes_index: TableRef[string, Lapper[GFFGene]]

  stderr.writeLine("GENOME = " & genome_assembly)


  # Load variants from XML
  logger.log(lvlInfo, "[main] Parsing variants from " & clinvar_xml_file)
  (variants_hash, filedate) = loadVariants(clinvar_xml_file, genome_assembly)
  if filename_date:

    var m: RegexMatch
    var r = regex.re".*ClinVarFullRelease_(?P<date>[0-9]{4}-[0-9]{2}).xml.gz"

    if clinvar_xml_file.match(r, m):
      filedate = m.groupFirstCapture("date", clinvar_xml_file) & "-01"

  # Gene name submission correction
  var hgncIndex: HgncIndex
  if args["--hgnc"]:
    let hgnc_file = $args["--hgnc"]
    hgncIndex = initHgncDbfromFile(hgnc_file)
    for id, clinvariant in variants_hash:
      correctGeneAlias(clinvariant, hgncIndex)

  if args["--gff"]:
    let gff_file = $args["--gff"]
    logger.log(lvlInfo, "[main] Load genes coordinates from " & gff_file)
    genes_index = loadGenesFromGFF(gff_file, gene_padding)

  # Sort variants by genomic order
  logger.log(lvlInfo, "[main] Sorting variants")
  variants_seq = toSeq(variants_hash.values())
  variants_seq.sort(cmpVariant)

  # Print VCF of STDOUT
  logger.log(lvlInfo, "[main] Printing variants")
  printVCF(variants_seq, genome_assembly, filedate, genes_index, coding_priority, hgncIndex)

when isMainModule:
  main(commandLineParams())
