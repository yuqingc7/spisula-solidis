vcftools --gzvcf TotalRawSNPs.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out raw.g5mac3
# After filtering, kept 582 out of 582 Individuals
# Outputting VCF file...
# After filtering, kept 696245 out of a possible 2898136 Sites

vcftools --vcf raw.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.g5mac3dp3 
# After filtering, kept 582 out of 582 Individuals
# Outputting VCF file...
# After filtering, kept 696245 out of a possible 696245 Sites

curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/filter_missing_ind.sh
chmod +x filter_missing_ind.sh
./filter_missing_ind.sh raw.g5mac3dp3.recode.vcf raw.g5mac3dplm
#                                           Histogram of % missing data per individual
#        70 +---------------------------------------------------------------------------------------------------------+
#           |         +  **    +         +         +        +         +        +         +         +        +         |
#           |            **                                     'totalmissing' using (bin($1,binwidth)):(1.0) ******* |
#        60 |-+          **                                                                                         +-|
#           |            **                                                                                           |
#           |            **                                                                                           |
#           |            **                                                                                           |
#        50 |-+          **                                                                                         +-|
#           |            **                                                                                           |
#           |            **                                                                                           |
#        40 |-+          **                                                                                         +-|
#           |            **            **                                                                             |
#           |            **            **                                                                             |
#        30 |-+          **           ***                                                                           +-|
#           |            **           *****                                                                           |
#           |            **           ******                                                                          |
#        20 |-+         ***           ******                                                                        +-|
#           |           ****        *********                                                                         |
#           |           ****  **    ***********                                                                       |
#           |           ****  **    ****************                                                                  |
#        10 |-+        ***********  *******************                                                             +-|
#           |          ********************************                                                               |
#           |         +*********************************    +**       +      *******     +         + ***** **         |
#         0 +---------------------------------------------------------------------------------------------------------+
#           0        0.1      0.2       0.3       0.4      0.5       0.6      0.7       0.8       0.9       1        1.1
#                                                        % of missing data

# The 85% cutoff would be 0.381249
# Would you like to set a different cutoff, yes or no
# yes
# Please enter new cutoff
# 0.5
# After filtering, kept 555 out of 582 Individuals
# Outputting VCF file...
# After filtering, kept 696245 out of a possible 696245 Sites

vcftools --vcf raw.g5mac3dplm.recode.vcf --max-missing 0.85 --maf 0.05 --recode --recode-INFO-all --out DP3g95maf05 --min-meanDP 20
# Note that hannah didn’t do --maf 0.05 here! (see alternative in filter_20231004_nomaf)
# After filtering, kept 555 out of 555 Individuals
# Outputting VCF file...
# After filtering, kept 38776 out of a possible 696245 Sites

# make a popmap by region+year
# 413 (NJ 1999)		413
# 536 (GBE 1999)	536
# CCB			BAR, BRW, CGC, PLY, PT, PVT, WV (### WV should be SCC)
# GBE (GBE 2019)	GBE
# SLI			BLP, MCX, SHN, SsLI			
# SCC			CTM, MW
# NAT			NAN
# FNJ (NJ 2019)	FNJ
# MAB (NJ 2022)	MAB22
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/pop_missing_filter.sh
chmod +x pop_missing_filter.sh
./pop_missing_filter.sh DP3g95maf05.recode.vcf popmap 0.1 9 DP3g95p5maf05
# Usage is pop_missing_filter vcffile popmap percent_missing_per_pop number_of_pops_for_cutoff name_for_output
# After filtering, kept 555 out of 555 Individuals
# Outputting VCF file...
# After filtering, kept 38760 out of a possible 38776 Sites

./dDocent_filters_hh DP3g95p5maf05.recode.vcf DP3g95p5maf05
# Usage is dDocent_filters_hh VCF_file Output_prefix
# Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
#  4931 of 38760
# Are reads expected to overlap?  In other words, is fragment size less than 2X the read length?  Enter yes or no.
# no
# Number of additional sites filtered based on overlapping forward and reverse reads
#  1912 of 33829
# Is this from a mixture of SE and PE libraries? Enter yes or no.
# no
# Number of additional sites filtered based on properly paired status
#  830 of 31917
# Number of sites filtered based on high depth and lower than 2*DEPTH quality score
#  1285 of 31087
    #                                           Histogram of mean depth per site
    #   350 +---------------------------------------------------------------------------------------------------------+
    #       |    +    +     +    +    +    +     +    +    +    +    +     +    +    +    +     +    +    +    +    + |
    #       |                                               'meandepthpersite' using (bin($1,binwidth)):(1.0) ******* |
    #   300 |-+                                  *  *                                                               +-|
    #       |                                ** ** ****                                                               |
    #       |                           ***  ***** ****                                                               |
    #       |                   ** *    *** ************     *                                                        |
    #   250 |-+              ** ** *   ********************  *                                                      +-|
    #       |             ** ******* ********************** ***** **  **                                              |
    #       |             *************************************** *** **                                              |
    #   200 |-+        ** ******************************************* **                                            +-|
    #       |          ** *************************************************                                           |
    #       |       ********************************************************                                          |
    #   150 |-+     ***********************************************************                                     +-|
    #       |       ************************************************************  *                                   |
    #       |     ***************************************************************** ** **                             |
    #   100 |-+   **************************************************************************   *                    +-|
    #       |    ****************************************************************************  *                      |
    #       |    ********************************************************************************** *                 |
    #       |    **************************************************************************************               |
    #    50 |-+  *******************************************************************************************  *     +-|
    #       |    *****************************************************************************************************|
    #       |    *****************************************************************************************************|
    #     0 +---------------------------------------------------------------------------------------------------------+
    #       10   20   30    40   50   60   70    80   90  100  110  120   130  140  150  160   170  180  190  200  210
    #                                                       Mean Depth
# If distrubtion looks normal, a 1.645 sigma cutoff (~90% of the data) would be 128230.785
# The 95% cutoff would be 194
# Would you like to use a different maximum mean depth cutoff than 194, yes or no
# yes
# Please enter new cutoff
# 200
# Number of sites filtered based on maximum mean depth
#  1550 of 31087
# Number of sites filtered based on within locus depth mismatch
#  25 of 29537
# Total number of sites filtered
#  9248 of 38760
# Remaining sites
#  29512
# Filtered VCF file is called DP3g95p5maf05.FIL.recode.vcf
# Filter stats stored in DP3g95p5maf05.filterstats
# Filtered VCF file is called DP3g95p5maf05.FIL.recode.vcf
# Filter stats stored in DP3g95p5maf05.filterstats

vcfallelicprimitives DP3g95p5maf05.FIL.recode.vcf --keep-info --keep-geno > DP3g95p5maf05.prim.vcf
vcftools --vcf DP3g95p5maf05.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.DP3g95p5maf05
# After filtering, kept 555 out of 555 Individuals
# Outputting VCF file...
# After filtering, kept 30552 out of a possible 31565 Sites

./filter_hwe_by_pop.pl -v SNP.DP3g95p5maf05.recode.vcf -p popmap_hwe -o SNP.DP3g95p5maf05.HWE -h 0.001
# Processing population: 413 (54 inds)
# Processing population: 536 (26 inds)
# Processing population: BAR (10 inds)
# Processing population: BLP (15 inds)
# Processing population: BRW (5 inds)
# Processing population: CGC (12 inds)
# Processing population: CTM (1 inds)
# Processing population: FNJ (49 inds)
# Processing population: GBE (46 inds)
# Processing population: MAB22 (156 inds)
# Processing population: MCX (34 inds)
# Processing population: MW (3 inds)
# Processing population: NAN (48 inds)
# Processing population: PLY (10 inds)
# Processing population: PT (23 inds)
# Processing population: PVT (10 inds)
# Processing population: SHN (5 inds)
# Processing population: SsLI (53 inds)
# Processing population: WV (22 inds)
# Outputting results of HWE test for filtered loci to 'filtered.hwe'
# Kept 30137 of a possible 30552 loci (filtered 415 loci)

# SNP.DP3g95p5maf05.HWE.recode.vcf will be renamed to basic vcf
cp SNP.DP3g95p5maf05.HWE.recode.vcf SNP_20231004_basic.recode.vcf

# Filter to only one per radtag by minor allele frequency
vcftools --vcf SNP.DP3g95p5maf05.HWE.recode.vcf --freq --out freqlist
#Open freqlist.fq in excel
		#Data by tabs and :
		#Before sorting, measure the minimum between the two (because it isn't minor allele, it's reference vs alternate)
            # =MIN(VALUE(MID(F2,3,LEN(F2)-2)),VALUE(MID(G2,3,LEN(G2)-2)))
		#Sort by max minor allele frequency then sort by unique contig
			#By Data -> remove duplicates > expand selection > only select column A
            # 20398 duplicate values found and removed; 9738 unique values remain. Note that counts may include empty cells, spaces, etc.
		#Export by keeping column A and B without headers
		freq_keep.txt #9739  freq_keep.txt
			#1fprad = 1 (by frequency) per radtag
vcftools --vcf SNP.DP3g95p5maf05.HWE.recode.vcf --positions freq_keep.txt --recode --recode-INFO-all --out SNP_20231004_1fprad 
# After filtering, kept 555 out of 555 Individuals
# Outputting VCF file...
# After filtering, kept 9739 out of a possible 30137 Sites
# Run Time = 44.00 seconds

p=20231004
bcftools view -H SNP_"$p"_1fprad.recode.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > "$p"_1fprad.chrom-map.txt
vcftools --vcf SNP_"$p"_1fprad.recode.vcf --chrom-map "$p"_1fprad.chrom-map.txt --out SNP_"$p"_1fprad --plink
awk '{$1=""}1' SNP_"$p"_1fprad.ped > SNP_"$p"_1fprad2.ped
mv SNP_"$p"_1fprad2.ped SNP_"$p"_1fprad.ped
plink --file SNP_"$p"_1fprad --r2 --ld-window-r2 0.80 --noweb --no-sex --allow-no-sex --out "$p"_1fp_pLD
awk '{ print $6 }' "$p"_1fp_pLD.ld | awk 'BEGIN { FS = ":" } ; { print $1,$2 }' | tail -n+2 > "$p"_1fp_pLD_chr.txt
vcftools --vcf SNP_"$p"_1fprad.recode.vcf --exclude-positions "$p"_1fp_pLD_chr.txt --recode --out SNP_"$p"_LDprune
# After filtering, kept 555 out of 555 Individuals
# Outputting VCF file...
# After filtering, kept 9634 out of a possible 9739 Sites

# three outputs
# SNP_20231004_basic.recode.vcf
# SNP_20231004_1fprad.recode.vcf
# SNP_20231004_LDprune.recode.vcf
