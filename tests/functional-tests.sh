#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

nim c -d:debug  -d:useSysAssert -d:useGcAssert --lineDir:on --debuginfo --boundChecks:on -x:on src/clinvcf
exe=./src/clinvcf

run simple_parsing $exe --gff tests/files/BRCA2.gff tests/files/37785.xml
assert_exit_code 0
assert_in_stdout "##fileDate=2019-12-31"
assert_in_stdout "13	32893387	37785	T	A"
assert_in_stdout "CLNSIG=Conflicting_interpretations_of_pathogenicity"
assert_in_stdout "ALLELEID=46341"
assert_in_stdout "GENEINFO=BRCA2:675"
assert_in_stdout "CLNREVSTAT=criteria_provided,_conflicting_interpretations"
assert_in_stdout "MC=SO:0001583|missense_variant"
assert_in_stdout "RS=80358507"

# Check integration of NCBI clinsig conversion
run ncbi_clnsig_conversion $exe tests/files/109.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Likely_pathogenic,_risk_factor"

# Multiple submission from same submitter
run mutli_subs_from_same_submitter $exe tests/files/307134.xml
assert_exit_code 0
assert_in_stdout "CLNREVSTAT=criteria_provided,_single_submitter"

# Multiple submission from same submitter
run skip_het_compound $exe tests/files/928.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Likely_pathogenic"

# Conflicting variants should always has a ReviewStatus conflicting
# Even if all submission are from the same submitter
run same_submitter_conflict $exe tests/files/1166.xml
assert_exit_code 0
assert_in_stdout "CLNREVSTAT=criteria_provided,_conflicting_interpretations"

# Multiple 3-4 stars subs, take them all !!! (see case 7108)
run run_multiple_3_4_star_subs $exe tests/files/7108.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Pathogenic,_drug_response"
assert_in_stdout "CLNREVSTAT=reviewed_by_expert_panel"

# Sort non-ACMG clnsig lexicographically
run sort_non_acmg_cnlsig_tags $exe tests/files/5333.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Affects,_risk_factor"

# Sort non-ACMG clnsig lexicographically
run expert_panel $exe tests/files/582.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Pathogenic;"

# Correction of conflicting interpretation
run conflict_deciphering $exe tests/files/9.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Pathogenic"
assert_in_stdout "OLD_CLNSIG=Conflicting_interpretations_of_pathogenicity"

# Handle multiple gene and select the prefered one from HGVS
run multi_gene_selection $exe --gff tests/files/TREX1.gff tests/files/225499.xml
assert_exit_code 0
assert_in_stdout "GENEINFO=TREX1:11277"

# Handle mutliple gene and select the prefered based on submissions (HGVS has no gene)
run multi_gene_selection $exe --gff tests/files/CFTR.gff tests/files/618897_2019-05.xml
assert_exit_code 0
assert_in_stdout "GENEINFO=CFTR:1080|"
