# Package

version       = "0.0.1"
author        = "Jérôme Audoux"
description   = "Generate a clean Clinvar VCF"
license       = "SEQONE"


# Dependencies

requires "hts >= 0.2.20", "q", "docopt"
srcDir = "src"
installExt = @["nim"]

bin = @["clinvcf"]

skipDirs = @["tests"]

import ospaths,strutils

task test, "run the tests":
  exec "bash tests/functional-tests.sh"
  exec "nim c  -d:useSysAssert -d:useGcAssert --lineDir:on --debuginfo --lineDir:on --debuginfo -r --threads:on tests/all"