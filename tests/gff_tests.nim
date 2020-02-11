import unittest
import clinvcfpkg/gff

suite "test GFF functions":

  test "test interval":
    var 
      a = Region(chrom: "A", start: 1, stop: 8)
      b = Region(chrom: "A", start: 5, stop: 12)
      c = Region(chrom: "A", start: 1, stop: 4)
      d = Region(chrom: "A", start: 1, stop: 4)

    check a.isOverlapping(b) == true
    check a.isOverlapping(c) == true
    check b.isOverlapping(c) == false

    check (d == c) == true
    a.merge(b)
    check a.start == 1
    check a.stop == 12
    check b.start == 5
    check b.stop == 12
