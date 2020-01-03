import httpclient, json, tables
import os, times
import q, xmltree # Parse XML
import docopt # Formating the command-line
import strutils # Split string
import hts

proc nextClinvarSet*(file: BGZ): string =
  for line in file:
    if line == "":
      break
    else:
      result.add(line & "\n")

proc formatVCFString*(vcf_string: string): string =
  result = vcf_string.replace(' ', '_')

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
    #variation_allele_file = $args["<variation_allele.txt.gz>"]
    #allele_variant_table = loadAlleleVariantTable(variation_allele_file)

  # TODO: Print VCF headers
  stderr.writeLine("[Log] Parsing variants from " & clinvar_xml_file)

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
  # ##INFO=<ID=MC,Number=.,Type=String,Description="comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence">
  # ##INFO=<ID=ORIGIN,Number=.,Type=String,Description="Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 3
  # 2 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other">
  # ##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID (i.e. rs number)">
  # ##INFO=<ID=SSR,Number=1,Type=Integer,Description="Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 
  # - 1kg_failed, 1024 - other">
  echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

  var 
    file : BGZ
    parsed_variants = initTable[string, int]()
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
          clinsig_nodes = reference_clinvar_assertion_nodes[0].select("clinicalsignificance")
          measureset_nodes = reference_clinvar_assertion_nodes[0].select("measureset")
        
        if measureset_nodes.len() > 0:
          let 
            measureset_node = measureset_nodes[0]      
            measure_nodes = measureset_nodes[0].select("measure")

          if measure_nodes.len() > 0:
            let
              measure_node = measure_nodes[0] 
              measure_relationship_nodes = measure_node.select("measurerelationship")
      
            for sequence_loc in measure_node.select("sequencelocation"):
              if sequence_loc.attr("Assembly") == genome_assembly:
                # <SequenceLocation Assembly="GRCh38" AssemblyAccessionVersion="GCF_000001405.38" AssemblyStatus="current" Chr="2" Accession="NC_000002.12" start="219469373" stop="219469408" display_start="219469373" display_stop="219469408" variantLength="36" positionVCF="219469370" referenceAlleleVCF="ATGACACAGTGTACGTGTCTGGGAAGTTCCCCGGGAG" alternateAlleleVCF="A"/>
                let
                  variant_id = measureset_node.attr("ID")#.parseInt()
                  allele_id = measure_node.attr("ID")
                  chrom = sequence_loc.attr("Chr")
                  pos = sequence_loc.attr("positionVCF")#.parseInt()
                  ref_allele = sequence_loc.attr("referenceAlleleVCF")
                  alt_allele = sequence_loc.attr("alternateAlleleVCF")

                if parsed_variants.hasKey(variant_id):
                  break
                
                var
                  clinical_significance : string
                  review_status : string
                  info_fields : seq[string] = @["ALLELEID=" & $allele_id]
                  gene_symbol: string
                  gene_id : int = -1
                
                if clinsig_nodes.len() > 0:
                  clinical_significance = clinsig_nodes[0].select("description")[0].innerText
                  review_status = clinsig_nodes[0].select("reviewstatus")[0].innerText
                  info_fields.add("CLNSIG=" & clinical_significance.formatVCFString())
                  info_fields.add("CLNREVSTAT=" & review_status.formatVCFString())
                
                if measure_relationship_nodes.len() > 0:
                  let gene_symbol_nodes = measure_relationship_nodes[0].select("symbol elementvalue")
                  if gene_symbol_nodes.len() > 0:
                    gene_symbol = gene_symbol_nodes[0].innerText
                  for xref_node in measure_relationship_nodes[0].select("xref"):
                    if xref_node.attr("DB") == "Gene":
                      gene_id = xref_node.attr("ID").parseInt()
                      # <XRef ID="672" DB="Gene" />
                      # <XRef Type="MIM" ID="113705" DB="OMIM" />

                #var info_fields = @["ALLELEID=" & $allele_id, "CLNSIG=" & clinical_significance.replace(' ', '_'), "CLNREVSTAT=" & review_status.replace(' ', '_')]
              
                if gene_id > 0:
                  # TODO: handle multi-gene variants
                  let gene_info : string = "GENEINFO=" & gene_symbol & ":" & $gene_id
                  info_fields.add(gene_info)
                
                if chrom != "" and pos != "" and ref_allele != "" and alt_allele != "":
                  parsed_variants[variant_id] = 1
                  echo [
                    chrom,
                    $pos,
                    $variant_id,
                    ref_allele,
                    alt_allele,
                    ".",
                    ".",
                    info_fields.join(";")
                  ].join("\t")
                break
      else:
        found_clinvarset = false

when isMainModule:
  main(commandLineParams())