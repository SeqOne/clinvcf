import unittest
import tables
import clinvcfpkg/hgnc

suite "test HGNC functions":

  test "test initHgncDbfromFile":
    var f = "tests/files/hgnc_toy.tsv"
    var hgncIndex = initHgncDbfromFile(f)
    # assert alias give the same entrez ID
    check hgncIndex.alias["LFS1"] == "TP53"
