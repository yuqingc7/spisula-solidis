# Files we need: 
# vcf file (/workdir/yc2644/Spisula/vcf_20231004/SNP_20231004_n540_basic.recode.vcf)
# bam files for each individual (soft linked from /workdir/yc2644/Spisula/dDocent_sol_20231004) - need to also link index bai files
# Reference genome (soft linked from /workdir/yc2644/Spisula/dDocent_sol_20231004)

rad_haplotyper.pl -v SNP_20231004_n540_basic.recode.vcf -r reference.fasta -p popmap_n540_OTU -x 16 -o haps_20231004_n540_OTU.vcf  --genepop haps_20231004_n540_OTU.gp
rad_haplotyper.pl -v SNP_20231004_n540_basic.recode.vcf -r reference.fasta -p popmap_n540_grp -x 16 -o haps_20231004_n540_grp.vcf  --genepop haps_20231004_n540_grp.gp
