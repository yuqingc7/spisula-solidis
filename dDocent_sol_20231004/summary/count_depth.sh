#!/bin/bash

## start date
start=`date +%s`

BASEDIR=/workdir/yc2644/Spisula/dDocent_sol_20231004
BAMLIST=$BASEDIR/sample_list_n540.txt
# Path to a list of bam files.

cat $BAMLIST | parallel -j 32 "samtools depth -aa $BASEDIR/{}-RG.bam | cut -f 3 | gzip > {}-RG.bam.depth.gz"

# end date
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
