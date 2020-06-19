import httpclient, json, tables
import os, times
import xmltree # Parse XML
import htmlparser
import docopt # Formating the command-line
import strutils # Split string
from streams import newStringStream
import hts

iterator nextClinvarSet*(file: var BGZ): string =
  var chunk: string
  for line in file:
    if line == "":
      yield chunk
      chunk = ""
    else:
      chunk.add(line & "\n")
  yield chunk

proc formatVCFString*(vcf_string: string): string =
  result = vcf_string.replace(' ', '_')

proc findNodes(n: XmlNode, tag: string): seq[XmlNode] =
  for xref_node in n:
    if xref_node.kind == xnElement:
      if xref_node.tag == tag:
        result.add(xref_node)

proc main*(argv: seq[string]) =

  # TODO: Create a usage and expose api_keys as options
  let doc = format("""
Usage: extractClinvarSet <clinvar.xml.gz> <variant_id>

  """)

  let
    args = docopt(doc)
    searched_id = $args["<variant_id>"]
    clinvar_xml_file = $args["<clinvar.xml.gz>"]
    #variation_allele_file = $args["<variation_allele.txt.gz>"]
    #allele_variant_table = loadAlleleVariantTable(variation_allele_file)

  # TODO: Print VCF headers
  stderr.writeLine("[Log] Parsing variants from " & clinvar_xml_file)

  var
    file : BGZ
    parsed_variants = initTable[string, int]()
    i = 0

  file.open(clinvar_xml_file, "r")

  for clinvarset_string in file.nextClinvarSet():
    if clinvarset_string != "" and clinvarset_string.startsWith("<ClinVarSet"):
      let
        root = parseHtml(newStringStream(clinvarset_string))
      for clinvarset_node in root.findNodes("clinvarset"):
        for reference_clinvar_assertion_nodes in clinvarset_node.findNodes("referenceclinvarassertion"):
          for measureset_node in reference_clinvar_assertion_nodes.findNodes("measureset"):
            let variant_id = measureset_node.attr("ID")
            if variant_id == searched_id:
              echo clinvarset_string

when isMainModule:
  main(commandLineParams())
