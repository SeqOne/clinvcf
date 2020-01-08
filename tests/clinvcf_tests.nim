import unittest, tables, hts
import clinvcf

suite "test utils functions":

  test "test nbStars":
    check rsNoAssertion.nbStars() == 0
    check rsNoAssertionCriteria.nbStars() == 0
    check rsNoAssertionVariant.nbStars() == 0
    check rsSingleSubmitter.nbStars() == 1
    check rsConflicting.nbStars() == 1
    check rsMutlipleSubmitterNoConflict.nbStars() == 2
    check rsExpertPanel.nbStars() == 3
    check rsPracticeGuideline.nbStars() == 4
  
  test "test parseNCBIConversionComment":
    check parseNCBIConversionComment("Converted during submission to Likely pathogenic.") == csLikelyPathogenic
    check parseNCBIConversionComment("Converted during submission to Benign.") == csBenign


  test "test IQRoutlierBounds":
    let 
      d = @[3.0,3.0,4.0,4.0,4.0,4.0]
      (min_val, max_val) = d.IQRoutlierBounds()

    check quantile(d, 0.25) == 3.25
    check quantile(d, 0.75) == 4

    var filtered_d : seq[float32]
    for v in d:
      if v >= min_val and v <= max_val:
        filtered_d.add(v)
    check filtered_d.len() == 6
