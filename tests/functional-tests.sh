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