setwd("/local/workdir/yc2644/Spisula")

#install.packages("/workdir/yc2644/Spisula/PopGenKit_1.0.tar.gz", repos = NULL, type = "source")

library(PopGenKit)
#mv haps_20240516_*/*.gen .
batchconvert(ndigit = 3)
# https://rdrr.io/cran/PopGenKit/man/batchconvert.html 

convert("haps_20231004_total_n23.gen",ndigit=3)

# convert("./haps_20240516_grp/haps_20231004_n540_grp.gen", ndigit = 3)
# convert("./haps_20240516_OTU/haps_20231004_n540_OTU.gen", ndigit = 3)

