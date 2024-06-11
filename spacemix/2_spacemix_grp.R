require(devtools)
#install_github("gbradburd/SpaceMix",build_vignettes=TRUE)
require(spam)
require(SpaceMix)
require(tidyverse)

#### Input File Creation ####
# YC - 2023.10.31
setwd("/local/workdir/yc2644/Spisula/spacemix/grp")

prefix<-"SNP_sol_B_n415_LDprune" #ped file prefix
dat <- prefix
# 
# # argv <- commandArgs(TRUE) 
# # prefix <- argv[1] #"SNP_sol_A_n109_LDprune"
# 
# sol_A_n109_grp_lon_lat <- left_join(read_tsv("sol_A_n109_grpmap",col_names=F) %>% rename(Sample=X1, Region=X2),
#            read_tsv("grp_lon_lat"))
# write_tsv(sol_A_n109_grp_lon_lat, paste0(prefix,".grp_lon_lat"))

sol_B_n415_grp_lon_lat <- left_join(read_tsv("sol_B_n415_grpmap",col_names=F) %>% rename(Sample=X1, Region=X2),
                                    read_tsv("grp_lon_lat"))
write_tsv(sol_B_n415_grp_lon_lat, paste0(prefix,".grp_lon_lat"))

grp_lon_lat <- read_tsv(paste0(prefix,".grp_lon_lat"))

# by individual
samp.siz.full <- read_delim(paste0("../samp.siz.full.data.",prefix), col_names = F)
allel.count.full <- read_delim(paste0('../allel.count.full.data.',prefix)) # table of allelic count

samp.siz.grp <- left_join(grp_lon_lat,samp.siz.full,join_by(Sample == X2)) %>% 
  dplyr::select(-X1) %>% 
  group_by(Region) %>% 
  summarise(across(X3:X9599, sum, .names = "sum_{.col}"))

write_tsv(samp.siz.grp, paste0("samp.siz.grp.",prefix),col_names = F)

allel.count.grp <- left_join(grp_lon_lat,allel.count.full,join_by(Sample == V2)) %>% 
  dplyr::select(-V1) %>% 
  group_by(Region) %>% 
  summarise(across(V3:V9599, sum, .names = "sum_{.col}"))

write_tsv(allel.count.grp, paste0("allel.count.grp.",prefix),col_names=F)

grp_lon_lat <- grp_lon_lat %>% 
  group_by(Region) %>% 
  summarise(long=first(long),
            lat=first(lat))

write_tsv(grp_lon_lat[-1], paste0("coord.grp.", prefix), col_names=F)

#### Run Spacemix ####
temp <- read.table(paste0("allel.count.grp.",prefix),header=F)
allele.counts <- as.matrix(temp[,-c(1)])
nrow(allele.counts)
ncol(allele.counts)

temp <- read.table(paste0("samp.siz.grp.",prefix),header=F)
sample.sizes <- as.matrix(temp[,-c(1)])
nrow(sample.sizes)
ncol(sample.sizes)

temp <- read.table(paste0("coord.grp.",prefix),header=F)
population.coordinates <- as.matrix(temp)

# Data option: allele counts and sample sizes
# Fast Model option: estimating geogenetic locations
# Long Model option: estimating geogenetic locations and 
#					 admixture source locations
# Spatial priors: default variance,
#					observed geographic sampling locations
run.spacemix.analysis(n.fast.reps = 10,
                      fast.MCMC.ngen = 1e5,
                      fast.model.option = "source_and_target",
                      long.model.option = "source_and_target",
                      data.type = "counts",
                      sample.frequencies=NULL,
                      mean.sample.sizes=NULL,
                      counts = allele.counts,
                      sample.sizes = sample.sizes,
                      sample.covariance=NULL,
                      target.spatial.prior.scale=NULL,
                      source.spatial.prior.scale=NULL,
                      spatial.prior.X.coordinates = population.coordinates[,1],
                      spatial.prior.Y.coordinates = population.coordinates[,2],
                      round.earth = FALSE,
                      long.run.initial.parameters=NULL,
                      k = nrow(allele.counts),
                      loci = ncol(sample.sizes),
                      ngen = 1e6,
                      printfreq = 1e2,
                      samplefreq = 1e3,
                      mixing.diagn.freq = 50,
                      savefreq = 1e5,
                      directory=NULL,
                      prefix = prefix)



