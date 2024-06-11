## subsample to N=23
# cat popmap_n540_grp | shuf | head -23 > total.keep
# awk '$2 == "GBE"' popmap_n540_grp | shuf | head -23 > GBE.keep
# awk '$2 == "CCB"' popmap_n540_grp | shuf | head -23 > CCB.keep
# awk '$2 == "SCC"' popmap_n540_grp | shuf | head -23 > SCC.keep
# awk '$2 == "SLI"' popmap_n540_grp | shuf | head -23 > SLI.keep
# awk '$2 == "NAT"' popmap_n540_grp | shuf | head -23 > NAT.keep
# awk '$2 == "NJ"' popmap_n540_grp | shuf | head -23 > NJ.keep
# awk '$2 == "DEL"' popmap_n540_grp | shuf | head -23 > DEL.keep
# awk '$2 == "A"' popmap_n540_OTU | shuf | head -23 > A.keep
# awk '$2 == "B"' popmap_n540_OTU | shuf | head -23 > B.keep

for f in total GBE CCB SCC SLI NAT NJ DEL A B
do
	vcftools --vcf /workdir/yc2644/Spisula/vcf_20231004/SNP_20231004_n540_basic.recode.vcf --keep /workdir/yc2644/Spisula/subsample_n23/$f.keep --mac 1 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out /workdir/yc2644/Spisula/subsample_n23/SpC_sol_$f
	#Filtering my min/max alleles is only about the maximum number of alleles, not how many are actually present in this subset of samples.
	#So minor allele count must be at least 1
done

# Clear the results file before writing new data
> results.txt

for f in total GBE CCB SCC SLI NAT NJ DEL A B
do
    grep -v "^#" /workdir/yc2644/Spisula/subsample_n23/SpC_sol_$f.recode.vcf | cut -f 1 | sort | uniq -c | awk '{print $1}' > SpC_sol_$f.txt
    awk '{ total += $1; count++ } END { print total/count }' SpC_sol_$f.txt > tmpmean
    awk '{x+=$1;y+=$1^2}END{print sqrt(y/NR-(x/NR)^2)}' SpC_sol_$f.txt > tmpstd
    echo "sol $f $(cat tmpmean) $(cat tmpstd)" >> results.txt
done

# Cleanup temporary files
rm tmpmean tmpstd
