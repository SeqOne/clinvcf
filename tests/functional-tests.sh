#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

nim c -d:debug  -d:useSysAssert -d:useGcAssert --lineDir:on --debuginfo --boundChecks:on -x:on src/clinvcf
exe=./src/clinvcf

run simple_parsing $exe tests/files/37785.xml
assert_exit_code 0
assert_in_stdout "##fileDate=2019-12-31"
assert_in_stdout "13	32893387	37785	T	A"
assert_in_stdout "CLNSIG=Conflicting_interpretations_of_pathogenicity"
assert_in_stdout "ALLELEID=46341"
assert_in_stdout "GENEINFO=BRCA2:675"
assert_in_stdout "CLNREVSTAT=criteria_provided,_conflicting_interpretations"
assert_in_stdout "MC=SO:0001583|missense variant"

# Check integration of NCBI clinsig conversion
run ncbi_clnsig_conversion $exe tests/files/109.xml
assert_in_stdout "CLNSIG=Likely_pathogenic,_risk_factor"

# Multiple submission from same submitter
run mutli_subs_from_same_submitter $exe tests/files/307134.xml
assert_in_stdout "CLNREVSTAT=criteria_provided,_single_submitter"