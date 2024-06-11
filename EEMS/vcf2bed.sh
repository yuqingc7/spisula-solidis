
vcftools --keep ../structure_20231004/sol_A_n110.txt --vcf ../vcf_20231004/SNP_20231004_n540_LDprune.recode.vcf --recode --out SNP_sol_A_n110_LDprune
vcftools --keep ../structure_20231004/sol_B_n415.txt --vcf ../vcf_20231004/SNP_20231004_n540_LDprune.recode.vcf --recode --out SNP_sol_B_n415_LDprune

/programs/plink2_linux_avx2_20230721/plink2 \
--vcf SNP_sol_A_n110_LDprune.recode.vcf --allow-extra-chr --make-bed \
--out SNP_sol_A_n110_LDprune --threads 24 --vcf-half-call m --max-alleles 2
# 'missing'/'m': Treat half-calls as missing.

/programs/plink2_linux_avx2_20230721/plink2 \
--vcf SNP_sol_B_n415_LDprune.recode.vcf --allow-extra-chr --make-bed \
--out SNP_sol_B_n415_LDprune --threads 24 --vcf-half-call m --max-alleles 2
# 'missing'/'m': Treat half-calls as missing.

bcftools view -Ov -S sol_A_n109.txt SNP_sol_A_n110_LDprune.recode.vcf > SNP_sol_A_n109_LDprune.recode.vcf

/programs/plink2_linux_avx2_20230721/plink2 \
--vcf SNP_sol_A_n109_LDprune.recode.vcf --allow-extra-chr --make-bed \
--out SNP_sol_A_n109_LDprune --threads 24 --vcf-half-call m --max-alleles 2
