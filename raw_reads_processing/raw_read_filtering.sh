#!/bin/bash -l

#SBATCH --partition=regular
#SBATCH --job-name=trim_seq_10t_40G
#SBATCH --output=s_solidistrim_sup.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=yc2644@cornell.edu

# cd solidis_sup_seq_raw
# ls *fastq.gz > /workdir/yc2644/Spisula/GBS/sol_sup_sample_temp.txt
# cd ..
# sed 's/_R[12]_001.fastq.gz//' sol_sup_sample_temp.txt | sort -u > sol_sup_sample.txt
# rm sol_sup_sample_temp.txt
# 224 sol_sup_sample.txt
# vim to remove 24 samples with reads below 1M 
# 200 sol_sup_sample.txt

export PATH=/programs/seqtk:$PATH
export PATH=/programs/seqkit-0.15.0:$PATH

export PYTHONPATH=/programs/cutadapt-4.1/lib/python3.9/site-packages:/programs/cutadapt-4.1/lib64/python3.9/site-packages
export PATH=/programs/cutadapt-4.1/bin:$PATH

INDIR="solidis_sup_seq_raw"
OUTDIR="solidis_sup_clean"

while read p; do
echo "$p"

#remove padding on R1 #possible that I should do this for R2 and just add back in the no cut sites but seems unnessecary
perl gbstrim.pl --enzyme1 psti --enzyme2 mspi --fastqfile $INDIR"/"$p"_R1_001.fastq.gz" --read R1 --outputfile $OUTDIR"/padtrim/"$p"_R1.trim.fastq" --verbose --threads 10 --minlength 50 >& $OUTDIR"/gbstrim_stats/"$p"_R1.log"

#remove adapter and padding on R2 and quality check
java -jar /programs/trimmomatic/trimmomatic-0.39.jar SE -phred33 -threads 8 $INDIR"/"$p"_R2_001.fastq.gz" $OUTDIR"/ommatic/"$p"_R2.fastq" ILLUMINACLIP:/programs/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 >& $OUTDIR"/ommatic_logs/"$p"_R2_trimmo.log"

seqtk trimfq -b 8 -e 8 $OUTDIR"/ommatic/"$p"_R2.fastq" > $OUTDIR"/ommatic_phase2/"$p"_R2.fastq"

#quality check R1
java -jar /programs/trimmomatic/trimmomatic-0.39.jar SE -phred33 -threads 8 $OUTDIR"/padtrim/"$p"_R1.trim.fastq" $OUTDIR"/trim_R1s_phase2/"$p"_R1.trim.fastq" SLIDINGWINDOW:4:15 MINLEN:50 >& $OUTDIR"/ommatic_logs/"$p"_R1_trimmo.log"

#resync
perl resync.pl $OUTDIR"/trim_R1s_phase2/"$p"_R1.trim.fastq" $OUTDIR"/ommatic_phase2/"$p"_R2.fastq" $OUTDIR"/sync_trim_6/"$p"_R1.trim.sync.fastq" $OUTDIR"/sync_trim_6/"$p"_R2.trim.sync.fastq" >& $OUTDIR"/sync_trim_6/"$p"_resync6.log"

done <sol_sup_sample_n200.txt 
