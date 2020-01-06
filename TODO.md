[x] Packaging with Nimble
[x] Sort VCF before
[ ] Gene stats module
[ ] Comp Clinvar VCF and ours
[ ] Unit testing
[x] Functional testing
[x] Gitlab CI
[ ] Implement new aggregation function + output stats of reclassification
[ ] Add RS tag (dbSNP rsid)
[x] Add date header
[ ] Add gnomad annot (using API calls and cache)
[ ] Add progressiong bar when loading variants
[ ] Create README file
[ ] Add NB_STARS tag 
[ ] Handle GenotypeSet (These variants should not be included in the VCF (see case 424779))
[x] Handle automatic clinsig conversions from NCBI
[x] Multiple submitter is not multiple submissions !!! (see case 307134)

<ClinicalSignificance>
      <ReviewStatus>no assertion criteria provided</ReviewStatus>
      <Description>likely pathogenic - adrenal pheochromocytoma</Description>
      <Comment Type="ConvertedByNCBI">Converted during submission to Likely pathogenic.</Comment>
    </ClinicalSignificance>

Correct discrepancies :

CLINSIG ERRORS 

A submission has a non harmonized clinsig values (ex: likely pathogenic - adrenal pheochromocytoma) that
is said to be converted to Likely pathogenic, but information is not really there in XML (or is it ?)
EX: DIFF OF CLNSIG for variant 109 : clinvar_2020-01.vcf = risk_factor <-> clinvar_20191223.vcf.gz = Likely_pathogenic,_risk_factor
EX: DIFF OF CLNSIG for variant 1365 : clinvar_2020-01.vcf = Pathogenic <-> clinvar_20191223.vcf.gz = Pathogenic/Likely_pathogenic
EX: DIFF OF CLNSIG for variant 1762 : clinvar_2020-01.vcf = not_provided <-> clinvar_20191223.vcf.gz = Benign

https://www.ncbi.nlm.nih.gov/clinvar/variation/1762/

It looks like we do not filter-out the 0-star "Pathogenic" submission from this one.
DIFF OF CLNSIG for variant 928 : clinvar_2020-01.vcf = Pathogenic/Likely_pathogenic <-> clinvar_20191223.vcf.gz = Likely_pathogenic
EX: DIFF OF CLNSIG for variant 1274 : clinvar_2020-01.vcf = Pathogenic/Likely_pathogenic <-> clinvar_20191223.vcf.gz = Pathogenic

REVSTAT ERRORS

MISSING VARIANTS
