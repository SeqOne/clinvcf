#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

nim c -d:debug  -d:useSysAssert -d:useGcAssert --lineDir:on --debuginfo --boundChecks:on -x:on src/clinvcf
grch37_version="--genome GRCh37"
exe=./src/clinvcf

run simple_parsing $exe --gff tests/files/BRCA2.gff $grch37_version tests/files/37785.xml
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
run ncbi_clnsig_conversion $exe $grch37_version tests/files/109.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Likely_pathogenic,_risk_factor"

# Multiple submission from same submitter
run mutli_subs_from_same_submitter $exe $grch37_version tests/files/307134.xml
assert_exit_code 0
assert_in_stdout "CLNREVSTAT=criteria_provided,_single_submitter"

# Multiple submission from same submitter
run skip_het_compound $exe $grch37_version tests/files/928.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Likely_pathogenic"

# Conflicting variants should always has a ReviewStatus conflicting
# Even if all submission are from the same submitter
run same_submitter_conflict $exe $grch37_version tests/files/1166.xml
assert_exit_code 0
assert_in_stdout "CLNREVSTAT=criteria_provided,_conflicting_interpretations"

# Multiple 3-4 stars subs, take them all !!! (see case 7108)
run run_multiple_3_4_star_subs $exe $grch37_version tests/files/7108.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Pathogenic,_drug_response"
assert_in_stdout "CLNREVSTAT=reviewed_by_expert_panel"

# Sort non-ACMG clnsig lexicographically
run sort_non_acmg_cnlsig_tags $exe $grch37_version tests/files/5333.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Affects,_risk_factor"

# Sort non-ACMG clnsig lexicographically
run expert_panel $exe $grch37_version tests/files/582.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Pathogenic;"

# Correction of conflicting interpretation
run conflict_deciphering $exe $grch37_version tests/files/9.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Pathogenic"
assert_in_stdout "OLD_CLNSIG=Conflicting_interpretations_of_pathogenicity"

# Handle multiple gene and select the prefered one from HGVS
run multi_gene_selection $exe --gff tests/files/TREX1.gff $grch37_version tests/files/225499.xml
assert_exit_code 0
assert_in_stdout "GENEINFO=TREX1:11277"

# Handle mutliple gene and select the prefered based on submissions (HGVS has no gene)
run multi_gene_selection $exe --gff tests/files/CFTR.gff $grch37_version tests/files/618897_2019-05.xml
assert_exit_code 0
assert_in_stdout "GENEINFO=CFTR:1080|"

# Mitochondrial annotations
run mito_anno $exe --gff tests/files/MT.gff $grch37_version tests/files/9618.xml
assert_exit_code 0
assert_in_stdout "GENEINFO=TRNE:4556"

# Coding first option for gene anno to force using protein coding annotation
run coding_first_control $exe --gff tests/files/ADORA2A.gff $grch37_version tests/files/225974.xml
assert_exit_code 0
assert_in_stdout "GENEINFO=ADORA2A-AS1:646023"
run coding_first_option $exe --gff tests/files/ADORA2A.gff --coding-first $grch37_version tests/files/225974.xml
assert_exit_code 0
assert_in_stdout "GENEINFO=ADORA2A:135"

# Consider close exonic regions (20bp padding) as exonic in gene priorization module
# In this case variant is in FTCD-AS1 (protein-coding) exon but 9bp away from FTCD.
# We discriminate these two gene using the gene_id of FTCD that is smaller that FTCD-AS1
run close_exonic_region $exe --gff tests/files/FTCD.gff $grch37_version tests/files/340430.xml
assert_exit_code 0
assert_in_stdout "GENEINFO=FTCD:10841|FTCD-AS1:100861507"

# For antivariant (same as the reference)
# We use the "." for the alternate allele representation
run antivariant $exe $grch37_version tests/files/242771.xml
assert_exit_code 0
assert_in_stdout "22	42523943	242771	A	."

# Haplotypes should not be exported in the VCF
run haplotype $exe $grch37_version tests/files/16895.xml
assert_equal "$(grep -v '^#' $STDOUT_FILE)" ""

# 3-stars reclassification system
run three_star_reclassification $exe $grch37_version tests/files/184976.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Pathogenic/Likely_pathogenic"
assert_in_stdout "CLNRECSTAT=3"

run two_star_reclassification $exe $grch37_version tests/files/140866.xml
assert_exit_code 0
assert_in_stdout "CLNSIG=Likely_pathogenic"
assert_in_stdout "CLNRECSTAT=2"

run one_star_reclassification $exe $grch37_version tests/files/182965.xml
assert_exit_code 0
assert_in_stdout "CLNRECSTAT=1"

# Pathology parsing
run pathology_field_parsing $exe $grch37_version tests/files/109.xml
assert_exit_code 0
assert_in_stdout "CLNDISEASE=pheochromocytoma_susceptibility_to|pheochromocytoma"
