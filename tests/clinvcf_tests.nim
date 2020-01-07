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