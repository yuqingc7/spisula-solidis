DIR=sol_A_109
DATE=20231025
IND=109 #415
SITES=9597

/programs/eems/runeems_snps/src/runeems_snps --datapath $DIR \
--mcmcpath $DIR"_"$DATE \
--nDemes 700 --nIndiv $IND --nSites $SITES \
--numMCMCIter 10000000  --numBurnIter 1000000  --numThinIter 999 > $DIR"_"$DATE/run_eems.log
