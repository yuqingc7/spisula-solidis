require(devtools)
#install_github("gbradburd/SpaceMix",build_vignettes=TRUE)
require(spam)
require(SpaceMix)

#### Input File Creation ####
# YC - 2023.10.31
setwd("/local/workdir/yc2644/Spisula/spacemix")

# https://github.com/QuentinRougemont/spacemix_workflow/tree/master/00.scripts

#Script by QR - 15-09-16
#Miniscript to reshape a classical ped file into a genotypic matrix (usefull for several programms)
#Input files: 1) ped file 
#Output: genotypic matrix (AC, TG, AA, CC, etc. with inds in row, mk in cols)
# prefix<-"SNP_sol_A_n109_LDprune" #ped file prefix
# dat <- prefix

argv <- commandArgs(TRUE) 
prefix <- argv[1] #"SNP_sol_A_n109_LDprune"
dat <- prefix

# note that R is reading some "T"s as TRUE, to solve this issue I set all column classes as "character"
dat2<-as.matrix(read.table(paste0(dat, ".ped"),h=F,colClasses = c("character"))) 
dat3<-dat2[,-c(1:6)] 
dat3[dat3 == "0"] <- NA

start <- seq(1, by = 2, length = ncol(dat3) / 2)
sdf <- sapply(start,function(i, dat3) paste(as.character(dat3[,i]),as.character(dat3[,i+1]), sep="") ,dat3 = dat3) 

write.table(cbind(dat2[,c(1:2)],as.data.frame(sdf)),paste0("genotypic.matrix.",prefix),quote=F,col.names=F,row.names=F)

system("sed -i 's/NANA/NA/g' genotypic.matrix*", wait=TRUE)

#!/usr/bin/env Rscript
#QR - 12-09-16
#Input: #A genotypic matrix containing individuals in rows and loci in colomns (AG,GG,GC,TC,...), NA allowed
# Output:
#1. allele frequency normalized covariance table (see Bradburg et al. 2016)
#2. mean sample sizes
#3. allele count for each inds and loci 
#4. sample size  for each inds and loci

dat<-paste0("genotypic.matrix.",prefix)
dat0=read.table(dat,h=F)

dat1<-dat0[,-c(1:2)]
genot<-t(as.matrix(dat1)) #set class to matrix
#print(genot)
m <- data.frame( N=rowSums(!is.na(genot)) ) #numbers of individuals successfully genottyped at each m

alleles <- apply(cbind(substr(genot,1,1),substr(genot,2,2)),1,unique) #grep the different alleles at each loci
if( is.matrix(alleles) ) { alleles <- lapply(apply(alleles,1,as.list),as.character) } #transform the matrice as a class list and as characters
alleles <- lapply(alleles,sort) #sort alleles, alphabetical order; ignore NA
m$numAlleles = sapply(alleles,length) 
#print(alleles)
#print(m$numAlleles)
if( any(m$numAlleles>2) ) { stop("genot contains more than two alleles.\n") } #check that no m have more than 2 different alleles

m$A1 <- NA
inds <- which(m$numAlleles>0) #ceux sans NA
m$A1[inds] <- sapply(alleles[inds],'[[',1)
m$A2 <- NA
inds <- which(m$numAlleles>1)
m$A2[inds] <- sapply(alleles[inds],'[[',2)

ref <- m$A1
alt <- NA
inds <- which(ref==m$A1); alt[inds] <- m$A2[inds]
inds <- which(ref==m$A2); alt[inds] <- m$A1[inds]

if( any(ref!=m$A1 & ref!=m$A2) ) { warning("ref allele not present in genot for some m. Conversions for these m cannot be performed and will be coerced to NA.\n") }

m$G2 = paste(ref,ref,sep="")	#2 copies of ref
m$G1.1 = paste(ref,alt,sep="")	#1 copy of ref, ref allele coded first
m$G1.2 = paste(alt,ref,sep="")	#1 copy of ref, reversed coding
m$G0 = paste(alt,alt,sep="")	#0 copy of ref
m$G2[is.na(ref)] <- NA #imputing Missing
m$G1.1[is.na(alt)] <- NA
m$G1.2[is.na(alt)] <- NA
m$G0[is.na(alt)] <- NA

genot.mat <- matrix( NA, ncol=ncol(genot), nrow=nrow(genot), dimnames=dimnames(genot) ) #create final genotypic matrix
genot.mat[genot==m$G2] <- 2 #impute values
genot.mat[genot==m$G1.1 | genot==m$G1.2] <- 1
genot.mat[genot==m$G0] <- 0

#Get Sample Size
SampleSize <- matrix( NA, ncol=ncol(genot), nrow=nrow(genot), dimnames=dimnames(genot) )
SampleSize[genot.mat==2] <- 2
SampleSize[genot.mat==1] <- 2
SampleSize[genot.mat==0] <- 2
SampleSize[is.na(SampleSize)] <-0
SampSizful=t(SampleSize)
write.table(cbind((dat0[,1:2]),SampSizful), paste0("samp.siz.full.data.",prefix),quote=F,row.names=F,col.names=F) #write table of sample size

samp.siz.1<-cbind(dat0[,1:2],SampSizful)

genot.mat[is.na(genot.mat)] <-0

write.table(cbind((dat0[,c(1:2)]),t(genot.mat)),paste0('allel.count.full.data.',prefix),quote=F,row.names=F) #write table of allelic count
