#!/bin/bash

################### THIS SCRIPT INVOLVES A VCF FILE AND PHENOTYPE FILE FROM PLINK TO RUN RVTESR


############ DIRECTORIES
work_dir=/home/data/mirna_snp_ms
out_dir=${work_dir}/rare_epilepsy/output
rvtest_dir=${work_dir}/mirna_dge/rvtest/executable
plink_test_dir=${work_dir}/epilepsy/output/test_pipeline/plink

############ RUN THE TEST SCRIPT
${rvtest_dir}/rvtest --inVcf ${plink_test_dir}/test.vcf.gz \
	--pheno ${plink_test_dir}/pt2_pheno_test.pheno \
	--pheno-name phenotype \
	--out ${out_dir}/wald_rvtest_out \
	--single wald
#	--single exact
	--out ${out_dir}/rvtest_out \
