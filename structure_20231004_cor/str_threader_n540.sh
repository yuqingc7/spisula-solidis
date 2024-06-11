#!/bin/bash -l

## start date
start=`date +%s`

PREFIX=SNP_20231004_n540

## set environment
export PYTHONPATH=/programs/structure_threader-1.3.10/lib/python3.9/site-packages
export PATH=/programs/structure_threader-1.3.10/bin:$PATH

mkdir $PREFIX"_out"

# run software
structure_threader run -i "/workdir/yc2644/Spisula/structure_20231004_cor/threader_input/"$PREFIX".recode.strct_in" \
-o "/workdir/yc2644/Spisula/structure_20231004_cor/"$PREFIX"_out" -st /programs/structure_2_3_4/bin/structure -K 8 -R 20 -t 8 --log true

# end date
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
