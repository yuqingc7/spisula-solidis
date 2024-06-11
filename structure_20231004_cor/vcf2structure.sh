#!/bin/bash -l

## start date
start=`date +%s`

VCF=SNP_20231004_n540_LDprune.recode.vcf
PREFIX=SNP_20231004_n540

/programs/plink-1.9-x86_64-beta5/plink --vcf $VCF \
--allow-extra-chr --make-bed \
--out $PREFIX

/programs/plink-1.9-x86_64-beta5/plink --bfile $PREFIX \
--allow-extra-chr \
--recode structure --out "threader_input/"$PREFIX

# end date
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
