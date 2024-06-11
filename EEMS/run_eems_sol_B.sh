DIR=sol_B
DATE=20231025
IND=415 #415
SITES=9597

/programs/eems/runeems_snps/src/runeems_snps --datapath $DIR \
--mcmcpath $DIR"_"$DATE \
--nDemes 700 --nIndiv $IND --nSites $SITES \
--numMCMCIter 10000000  --numBurnIter 1000000  --numThinIter 999 > $DIR"_"$DATE/run_eems.log
