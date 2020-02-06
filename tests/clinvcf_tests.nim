import unittest, tables, hts, strutils
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

  test "test clinsig conversion":
    check parseEnum[ClinSig]("Benign") == csBenign
    check parseEnum[ClinSig]("Benign/Likely benign") == csBenignLikelyBenign
    check parseEnum[ClinSig]("Likely benign") == csLikelyBenign
    check parseEnum[ClinSig]("Uncertain significance") == csUncertainSignificance
    check parseEnum[ClinSig]("Likely pathogenic") == csLikelyPathogenic
    check parseEnum[ClinSig]("Pathogenic") == csPathogenic
    check parseEnum[ClinSig]("Pathogenic/Likely pathogenic") == csPathogenicLikelyPathogenic
    check parseEnum[ClinSig]("drug response") == csDrugResponse
    check parseEnum[ClinSig]("association") == csAssociation
    check parseEnum[ClinSig]("risk factor") == csRiskFactor
    check parseEnum[ClinSig]("protective") == csProtective
    check parseEnum[ClinSig]("Affects") == csAffects
    check parseEnum[ClinSig]("conflicting data from submitters") == csConflictingDataFromSubmitters
    check parseEnum[ClinSig]("other") == csOther
    check parseEnum[ClinSig]("not provided") == csUnknown
  
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
