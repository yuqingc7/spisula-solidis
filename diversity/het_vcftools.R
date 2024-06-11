setwd("/local/workdir/yc2644/Spisula/vcf_20231004")

library(tidyverse)

het <- read_tsv("SNP_sol_A_n110_basic.het") %>% 
  mutate(`E(HET)`=N_SITES - `E(HOM)`) %>% 
  mutate(He=`E(HET)`/N_SITES)

mean(het$He)
sd(het$He)

