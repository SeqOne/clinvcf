import logging, strutils, re

var logger* = newConsoleLogger(fmtStr="[$datetime] - $appname - $levelname : ", useStderr=true)

proc contains*[T](s: seq[T], value: T): bool =
    for elem in s:
        if elem == value:
            return true
    return false

proc formatVCFString*(vcf_string: string): string =
  result = vcf_string.replace(' ', '_')

let
  rePathoTrim = re"^\W+|\W+$"
  rePathoReplace = re"\W+"

proc formatPathoString*(pathoString: string): string =
  # First : remove all non-word chars at the end of the string (or first)
  result = pathoString.replace(rePathoTrim)
  # Second replace all non-word characters by '_'
  result = result.replace(rePathoReplace, "_")