# Package

version       = "0.0.1"
author        = "Jérôme Audoux, Sacha Beaumeunier"
description   = "Generate a clean Clinvar VCF"
license       = "SEQONE"


# Dependencies

requires "hts >= 0.2.20 & <= 0.2.23", "q", "docopt"#, "lapper"
requires "https://github.com/GULPF/tiny_sqlite#head"
requires "regex >= 0.13"
srcDir = "src"
installExt = @["nim"]

bin = @["clinvcf", "extractClinvarSet"]

skipDirs = @["tests"]

import ospaths,strutils

task test, "run the tests":
  exec "nim c  -d:useSysAssert -d:useGcAssert --lineDir:on --debuginfo --lineDir:on --debuginfo -r --threads:on tests/all"
  exec "bash tests/functional-tests.sh"
