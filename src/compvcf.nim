import httpclient, json, tables, os, times
import docopt
import strutils # Split string
import hts

type

  ClinVariant* = ref object
    variant_id: int
    clinsig: string
    revstat: string

proc loadVariantsFromVCF(filename: string): TableRef[int, ClinVariant] =
  var 
    file : BGZ
    nb_reclassif = 0
  
  result = newTable[int, ClinVariant]()
  file.open(filename, "r")
  for line in file:
    if line.len() == 0 or line[0] == '#':
      continue
    else:
      let 
        v = line.split('\t')
        f = v[7].split(';')
      var 
        h = initTable[string, string]()
      for field in f:
        let a = field.split('=')
        h[a[0]] = a[1]
      # We ignore re-classified variants because they have a different clinsig
      if h.hasKey("OLD_CLNSIG"):
        inc(nb_reclassif)
        var variant = ClinVariant(variant_id: v[2].parseInt(), clinsig: "reclassified", revstat: "reclassified")
        result[variant.variant_id] = variant
        #continue
      elif h.hasKey("CLNSIG") and h.hasKey("CLNREVSTAT"):
        var variant = ClinVariant(variant_id: v[2].parseInt(), clinsig: h["CLNSIG"], revstat: h["CLNREVSTAT"])
        result[variant.variant_id] = variant
  if nb_reclassif > 0:
    stderr.writeLine("[Log] Found " & $nb_reclassif & " reclassified variants in " & filename)


proc main*(argv: seq[string]) =

  # TODO: Create a usage and expose api_keys as options
  let doc = format("""
Usage: compVCF <1.vcf> <2.vcf>

  """)

  let 
    args = docopt(doc)
    vcf1 = $args["<1.vcf>"]
    vcf2 = $args["<2.vcf>"]
  
  stderr.writeLine("[Log] Loading variants from " & vcf1)
  var variants1 = vcf1.loadVariantsFromVCF()
  stderr.writeLine("[Log] " & $variants1.len() & " variant loaded")
  stderr.writeLine("[Log] Loading variants from " & vcf2)
  var variants2 = vcf2.loadVariantsFromVCF()
  stderr.writeLine("[Log] " & $variants2.len() & " variant loaded")

  var
    nb_wrong_clinsig = 0
    nb_wrong_revstat = 0
    nb_missing_variant_v1 = 0
    nb_missing_variant_v2 = 0

  for vid, v1 in variants1:
    if variants2.hasKey(vid):
      if variants2[vid].clinsig == "reclassified" or v1.clinsig == "reclassified":
        continue
      if variants2[vid].clinsig != v1.clinsig:
        inc(nb_wrong_clinsig)
        echo "DIFF OF CLNSIG for variant " & $vid & " : " & vcf1 & " = " & v1.clinsig & " <-> " & vcf2 & " = " & variants2[vid].clinsig
      if variants2[vid].revstat != v1.revstat:
        inc(nb_wrong_revstat)
        echo "DIFF OF REVSTAT for variant " & $vid & " : " & vcf1 & " = " & v1.revstat & " <-> " & vcf2 & " = " & variants2[vid].revstat
    else: 
      echo "MISSING variant " & $vid & " in " & vcf2
      inc(nb_missing_variant_v2)
  
  for vid, v2 in variants2:
    if not variants1.hasKey(vid):
      inc(nb_missing_variant_v1)
      echo "MISSING variant " & $vid & " in " & vcf1
  
  stderr.writeLine("[Stats] NB_WRONG_CLINSIG " & $nb_wrong_clinsig)
  stderr.writeLine("[Stats] NB_WRONG_REVSTAT " & $nb_wrong_revstat)
  stderr.writeLine("[Stats] NB_MISSING_VARIANT_VCF1 " & $nb_missing_variant_v1)
  stderr.writeLine("[Stats] NB_MISSING_VARIANT_VCF2 " & $nb_missing_variant_v2)

when isMainModule:
  main(commandLineParams())