import unittest
import tables
import clinvcfpkg/hgnc

suite "test HGNC functions":

  test "test initHgncDbfromFile":
    var f = "tests/files/hgnc_entrez.tsv"
    var hgncIndex = initHgncDbfromFile(f)
    # assert alias give the same entrez ID
    check hgncIndex.alias["PPP1R53"] == "BRCA1"
    check hgncIndex.alias["KIAA0644"] == "TRIL"
    check hgncIndex.alias["DBA"] == "RPS19"
    check hgncIndex.alias["LFS1"] == "TP53"
