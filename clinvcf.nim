import httpclient, json, tables, sequtils
import os, times
import q, xmltree # Parse XML
import docopt # Formating the command-line
import strutils # Split string
import hts

type
  ClinSig = enum
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

  RevStat = enum
    rsNoAssertion = "no assertion provided",
    rsNoAssertionCriteria = "no assertion criteria provided",
    rsNoAssertionVariant = "no assertion for the individual variant",
    rsSingleSubmitter = "criteria provided, single submitter",
    rsMultipleSubmitterConflicting = "criteria provided, conflicting interpretations",
    rsMutlipleSubmitterNoConflict = "criteria provided, multiple submitters, no conflicts",
    rsExpertPanel = "reviewed by expert panel",
    rsPracticeGuideline = "practice guideline"

  Submission = ref object
    clinical_significance: ClinSig
    review_status: RevStat

  MolecularConsequence = ref object
    description: string
    so_term : string

  ClinVariant = ref object
    variant_id: int
    allele_id: int
    chrom: string
    pos: int
    ref_allele: string
    alt_allele: string
    gene_id: int
    gene_symbol: string
    molecular_consequences: seq[MolecularConsequence]
    submissions: seq[Submission]

var
  acmg_clinsig = @[csBenign, csLikelyBenign, csUncertainSignificance, csLikelyPathogenic, csPathogenic]
  non_acmg_clinsig = @[csDrugResponse, csRiskFactor, csAffects, csAssociation, csProtective, csConflictingDataFromSubmitters, csOther]

proc `$`*(mc: MolecularConsequence): string =
  return mc.so_term & "|" & mc.description

proc aggregateReviewStatus*(revstat_count: TableRef[RevStat, int], total: int, has_conflict = false): RevStat =
  if total > 1 and revstat_count.hasKey(rsSingleSubmitter):
    if has_conflict:
      result = rsMultipleSubmitterConflicting
    else:
      result = rsMutlipleSubmitterNoConflict
  elif revstat_count.hasKey(rsSingleSubmitter):
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
    total = 0
    no_aggregation_needed = false

  for sub in submissions:
    # When there is a submission from an Expert panel or from a group providing practice guidelines, only the interpretation from that group is reported in the aggregate record, even if other submissions provide different interpretations. 
    if sub.review_status == rsExpertPanel or sub.review_status == rsPracticeGuideline:
      result.clinsig = $sub.clinical_significance
      result.revstat = $sub.review_status
      no_aggregation_needed = true
      break

    if clinsig_count.hasKey(sub.clinical_significance):
      inc(clinsig_count[sub.clinical_significance])
    else:
      clinsig_count[sub.clinical_significance] = 1
    
    if revstat_count.hasKey(sub.review_status):
      inc(revstat_count[sub.review_status])
    else:
      revstat_count[sub.review_status] = 1
    inc(total)
  
  # If we have found and rsExpertPanel or rsPracticeGuideline we do not perform aggreagtion of submissions
  if not no_aggregation_needed:
    
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
      result.revstat = $revstat_count.aggregateReviewStatus(total, false)
    # Case #2 Patho and Likely Patho (only)
    elif nb_acmg_tags == 2 and clinsig_count.hasKey(csPathogenic) and clinsig_count.hasKey(csLikelyPathogenic):
      result.clinsig = "Pathogenic/Likely pathogenic"
      result.revstat = $revstat_count.aggregateReviewStatus(total, false)
    # Case #3, Only patho entries
    elif nb_acmg_tags == 2 and clinsig_count.hasKey(csBenign) and clinsig_count.hasKey(csLikelyBenign):
      result.clinsig = "Benign/Likely benign"
      result.revstat = $revstat_count.aggregateReviewStatus(total, false)
    # Case #4, Conflict !!!
    # TODO: Do some desambiguiations !!!
    elif (clinsig_count.hasKey(csPathogenic) or clinsig_count.hasKey(csLikelyPathogenic) or clinsig_count.hasKey(csBenign) or clinsig_count.hasKey(csLikelyBenign)) and clinsig_count.hasKey(csUncertainSignificance):
      result.clinsig = "Conflicting interpretations of pathogenicity"
      result.revstat = $revstat_count.aggregateReviewStatus(total, true)
    elif (clinsig_count.hasKey(csPathogenic) or clinsig_count.hasKey(csLikelyPathogenic)) and (clinsig_count.hasKey(csBenign) or clinsig_count.hasKey(csLikelyBenign)):
      result.clinsig = "Conflicting interpretations of pathogenicity"
      result.revstat = $revstat_count.aggregateReviewStatus(total, true)
    
    # Add non-ACMG values to the end of clinsig
    # If ClinVar aggregates submissions from groups that provided a standard term not recommend by ACMG/AMP ( e.g. drug response), then those values are reported after the ACMG/AMP-based interpretation (see the table below).
    var additional_cstags : seq[string]
    for cstag in non_acmg_clinsig:
      if clinsig_count.hasKey(cstag):
        additional_cstags.add($cstag)
    if additional_cstags.len() > 0:
      if result.clinsig == "":
        result.clinsig = additional_cstags.join(", ")
      else:
        result.clinsig.add(", " & additional_cstags.join(", "))
    
    # Handle default values
    if result.clinsig == "": 
      result.clinsig = $csUnknown
    if result.revstat == "":
      result.revstat = $revstat_count.aggregateReviewStatus(total, false)

proc nextClinvarSet*(file: BGZ): string =
  for line in file:
    if line == "":
      break
    else:
      result.add(line & "\n")

proc loadVariants*(clinvar_xml_file: string, genome_assembly: string): TableRef[int, ClinVariant] =
  result = newTable[int, ClinVariant]()

  stderr.writeLine("[Log] Parsing variants from " & clinvar_xml_file)

  var 
    file : BGZ
    i = 0

  file.open(clinvar_xml_file, "r")
  
  var found_clinvarset = true
  while found_clinvarset:
    # TODO: Add some kind of loader ever 10K parsed variants
    let clinvarset_string = file.nextClinvarSet()
    if clinvarset_string != "" and clinvarset_string.startsWith("<ClinVarSet"):
      #echo clinvarset_string
      let 
        doc = q(clinvarset_string)
        reference_clinvar_assertion_nodes = doc.select("referenceclinvarassertion")
      
      if reference_clinvar_assertion_nodes.len() > 0:
        let 
          #clinsig_nodes = reference_clinvar_assertion_nodes[0].select("clinicalsignificance")
          measureset_nodes = reference_clinvar_assertion_nodes[0].select("measureset")
        
        if measureset_nodes.len() > 0:
          let 
            measureset_node = measureset_nodes[0]
            variant_id = measureset_node.attr("ID").parseInt()
            measure_nodes = measureset_nodes[0].select("measure")

          # Only parse measure node to extract variant position if we do not have seen this variant
          # Already
          if not result.hasKey(variant_id) and measure_nodes.len() > 0:
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
                
                var 
                  variant = ClinVariant(
                    chrom: chrom,
                    pos: pos,
                    variant_id: variant_id,
                    allele_id: allele_id,
                    ref_allele: ref_allele,
                    alt_allele: alt_allele,
                    gene_id: gene_id,
                    gene_symbol: gene_symbol
                  )
                result[variant_id] = variant

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
                      variant.molecular_consequences.add(MolecularConsequence(description: description, so_term: so_term))

                break # We found our "sequenceLocation"

          # Not lets add the submissions
          if result.hasKey(variant_id):
            for clinvar_assertion_node in doc.select("clinvarassertion"):
              let clinsig_nodes = clinvar_assertion_node.select("clinicalsignificance")
              if clinsig_nodes.len() > 0: # FIXME: Should not be > to 1 ...
                var
                  clinical_significance : ClinSig = csUnknown
                  review_status : RevStat = rsNoAssertion
                
                if clinsig_nodes.len() > 0:
                  var
                    desc_nodes = clinsig_nodes[0].select("description")
                    revstat_nodes = clinsig_nodes[0].select("reviewstatus")
                  if desc_nodes.len() > 0:
                    clinical_significance = parseEnum[ClinSig](desc_nodes[0].innerText, csUnknown)
                  if revstat_nodes.len() > 0:
                   review_status = parseEnum[RevStat](revstat_nodes[0].innerText, rsNoAssertion)
                   
                  # Add the submission to the variant record
                  result[variant_id].submissions.add(Submission(clinical_significance: clinical_significance, review_status: review_status))
    else:
      # We did not found a clinvarset entry, it can either mean that we are at the end of the file
      # Or that we are parsing headers of the xml file
      if clinvarset_string == "":
        found_clinvarset = false
      # TODO: THese are the headers to be parsed to retrieve the date
      #else:
      

proc formatVCFString*(vcf_string: string): string =
  result = vcf_string.replace(' ', '_')

proc printVCF*(variants: TableRef[int, ClinVariant], genome_assembly: string) =
  echo "##fileformat=VCFv4.1"
  ##fileDate=2019-12-23 # TODO: Get date from XML headers
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
  # ##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
  # ##INFO=<ID=SSR,Number=1,Type=Integer,Description="Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 
  # - 1kg_failed, 1024 - other">
  echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

  for v in variants.values():
    let (clinsig, revstat) = v.submissions.aggregateSubmissions()
    var info_fields : seq[string] = @["ALLELEID=" & $v.allele_id]
    
    if v.gene_id > 0:
      info_fields.add("GENEINFO=" & v.gene_symbol & ":" & $v.gene_id)

    info_fields.add("CLNSIG=" & clinsig.formatVCFString())
    info_fields.add("CLNREVSTAT=" & revstat.formatVCFString())

    if v.molecular_consequences.len > 0:
      var formated_consequences = map(v.molecular_consequences, proc (x: MolecularConsequence): string = $x)
      info_fields.add("MC=" & join(formated_consequences, ","))
  
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
  
  var variants: TableRef[int, ClinVariant]
    #variation_allele_file = $args["<variation_allele.txt.gz>"]
    #allele_variant_table = loadAlleleVariantTable(variation_allele_file)

  variants = loadVariants(clinvar_xml_file, genome_assembly)
  printVCF(variants, genome_assembly)

when isMainModule:
  main(commandLineParams())