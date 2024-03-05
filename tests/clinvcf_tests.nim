import unittest, tables, hts, strutils
import strformat
import clinvcf
import compvcf

suite "test compare vcf files versions":
  let
    vcf1 = "tests/files/ClinVarFullRelease_5MB_Rand_00-latest.vcf.gz"
    vcf2 = "tests/files/ClinVarRCVRelease_5MB_Rand_00-latest.vcf.gz"
  var vcf1_tab = vcf1.loadClinvarVariantsFromVCF()
  var vcf2_tab = vcf2.loadClinvarVariantsFromVCF()

  var
    missing_rows_from_vcf1 = 0
    missing_rows_from_vcf2 = 0
  for id1, v1 in vcf1_tab:
    if vcf2_tab.hasKey(id1):
      var
         v2 = vcf2_tab[id1]
      # Target file is supposed to contain more INFO fields than the base file. In other words,
      # base file INFO fields must be a subset of the target file ones.
      check v1.info_fields in v2.info_fields
      # The 7 main fields(CHROM,POS,ID,REF,ALT,QUAL,FILTER) from both files must be equal.
      check v1.main_fields == v2.main_fields
    else:
      inc(missing_rows_from_vcf2)
  check missing_rows_from_vcf1 == 0
  for id2, v2 in vcf2_tab:
    if not vcf1_tab.hasKey(id2):
      inc(missing_rows_from_vcf1)
  check missing_rows_from_vcf2 == 0


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


  test "test pathology string format":
    check formatPathoString(" Factor X Deficiency ") == "Factor_X_Deficiency"
    check formatPathoString("Factor (X) Deficiency, pathology") == "Factor_X_Deficiency_pathology"
    check formatPathoString("Factor, (X) Deficiency, pathology,|cancer/") == "Factor_X_Deficiency_pathology|cancer"
    check formatPathoString("Factor, X, Deficiency      ,pathology/| cancer") == "Factor_X_Deficiency_pathology|cancer"
    check formatPathoString("  , Factor X,,Deficiency/pathology| , ,cancer/") == "Factor_X_Deficiency_pathology|cancer"
    check formatPathoString("CLNDISEASE=  , Factor X,,Deficiency/pathology| , ,cancer/") == "CLNDISEASE=Factor_X_Deficiency_pathology|cancer"


  test "test clinical pathology parsing":
    check parseClinicalPathologies("DISEASE", @["coagulation_x_deficiency", "factor_x_deficiency"]) == "CLNDISEASE=coagulation_x_deficiency|factor_x_deficiency"
