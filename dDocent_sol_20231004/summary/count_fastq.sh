#!/bin/bash

## start date
start=`date +%s`

## This script is used to count number of bases in raw, adapter clipped, and quality filtered fastq files. The result of this script will be stored in a nohup file.
BASEDIR=/workdir/yc2644/Spisula/dDocent_sol_20231004 # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories. An example for the Greenland cod data is: /workdir/cod/greenland-cod/
SAMPLELIST=$BASEDIR/sample_list_n540.txt # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table. An example of such a sample list is /workdir/cod/greenland-cod/sample_lists/sample_list_pe_1.tsv
FASTQDIR=/workdir/yc2644/Spisula/dDocent_sol_20231004 # Path to raw fastq files. An example for the Greenland cod data is: /workdir/backup/cod/greenland_cod/fastq/

# Create headers for the output
printf 'sample_seq_id\traw_reads\n'

# Loop over each sample in the sample table
for SAMPLEFILE in `cat $SAMPLELIST`; do
	FASTQFILES=$FASTQDIR'/'$SAMPLEFILE'*F.fq.gz'  # The input path and file prefix
	
	# Count the number of reads in raw fastq files. We only need to count the forward reads, since the reverse will contain exactly the same number of reads. 
    # fastq files contain 4 lines per read, so the number of total reads will be half of this line number. 
	RAWREADS=`zcat $FASTQFILES | wc -l`
	
	# Count the number of bases in raw fastq files. We only need to count the forward reads, since the reverse will contain exactly the same number of bases. The total number of reads will be twice this count. 
	#RAWBASES=`zcat $FASTQFILES | awk 'NR%4==2' | tr -d "\n" | wc -m` 

	# Write the counts in appropriate order.
	printf "%s\t%s\n" $SAMPLEFILE $((RAWREADS/4)) #$RAWBASES
	
done

# end date
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"

