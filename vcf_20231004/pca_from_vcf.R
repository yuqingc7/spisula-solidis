setwd("/local/workdir/yc2644/Spisula/vcf_20231004")

library(gdsfmt)
library(SNPRelate) # if there is something wrong with gfortran see link here https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/
library(tidyverse)
library(ggplot2)
library(SeqArray)
library(ggrepel)
library(RColorBrewer)

# PCA n=540 ---------------------------------------------------------------

VCF="SNP_20231004_n540_LDprune.recode.vcf"

# meta <- left_join(read_tsv("sample_list_n540.txt", col_names=F), 
#                   read_tsv("popmap_n543", col_names = F))
#write_tsv(meta, "popmap_n540", col_names = F)
#meta <- read_tsv("popmap_n540", col_names = F)
meta <- left_join(read_tsv("popmap_n540", col_names = F),read_tsv("metadata_n540", col_names=T),
                  by = join_by(X1 == sample)) %>% 
  select(X1, `analytical group`, `collection site`)

# Region			Site Codes
# Cape Cod Bay		BAR, BRW, CGC, PLY, PT, PVT
# Offshore	George’s Bank		GBE, 536
# Southern Long Island	BLP, MCX, SHN, SsLI			
# Southern Cape Cod		CTM, MW, WV
# Offshore	Nantucket Shoals		NAN
# Offshore	Delmarva			FNJ, 413
# MAB (NJ 2022)		MAB

names(meta)[1] <- "sample.id"
names(meta)[2] <- "grp"
names(meta)[3] <- "grp2"
unique(meta$grp) 

# vcf
vcf.fn <- paste0("./",VCF)
# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "pruned.noinvers.gds", method="biallelic.only")
# showfile.gds(closeall=TRUE)
# summary
snpgdsSummary("pruned.noinvers.gds")
# Open the GDS file
genofile <- snpgdsOpen("pruned.noinvers.gds")

pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = meta$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
dim(tab)
dim(meta)

tab_meta <- right_join(meta,tab, by="sample.id")

# regionsol<-c("#3BA4F4", # Cape Cod Bay (CCB)
#              "#D0679C", # GBE
#              "#FF2E17", # MAB
#              "#3D2C7D", # NAT
#              "#FFB357", # NJ/Delmarva
#              "#5BECAB", # Southern Cape Cod (SCC)
#              "#5ABB64","#5ABB64") # Southern Long Island (SLI)
#c("#D0679C","#3D2C7D","#3BA4F4","#5BECAB","#5ABB64","#FF2E17","#FFB357","#FFE017","#FF8A17")

regionsol<-c("#3BA4F4", # Cape Cod Bay (CCB)
             "#FFB357", # NJ/Delmarva
             "#D0679C", # GBE
             "#3D2C7D", # NAT
             "#FF2E17", # MAB
             "#5BECAB", # Southern Cape Cod (SCC)
             "#5ABB64","#5ABB64") # Southern Long Island (SLI)

ggplot(tab_meta, aes(x = EV1, y = EV2, color=grp)) + 
  geom_point(size=2) +
  scale_color_manual(values = regionsol, "Region") +
  scale_fill_manual(values = regionsol, "Region") +
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  #geom_label_repel(aes(label=label), max.overlaps = Inf, min.segment.length = 0)+
  theme_bw()+
  theme(legend.margin = margin(0.25, 0.25, 0.25, 0.25)) +
  # stat_ellipse(geom = "polygon",aes(fill = grp),
  #              type = "norm",level = 0.95,alpha=0.05) +
  coord_cartesian(xlim=c(-0.05, 0.105), ylim = c(-0.125,0.10))+
  guides(size = "none")

ggplot(tab_meta, aes(x = EV1, y = EV3, color=grp)) + 
  geom_point(size=2) +
  scale_color_manual(values = regionsol, "Region") +
  scale_fill_manual(values = regionsol, "Region") +
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC3 (",round(pc.percent[3],3),"%",")",sep=""))+ 
  #geom_label_repel(aes(label=label), max.overlaps = Inf, min.segment.length = 0)+
  theme_bw()+
  theme(legend.margin = margin(0.25, 0.25, 0.25, 0.25)) +
  # stat_ellipse(geom = "polygon",aes(fill = grp),
  #              type = "norm",level = 0.95,alpha=0.05) +
  coord_cartesian(xlim=c(-0.05, 0.105), ylim = c(-0.125,0.175))+
  guides(size = "none")

#OTU B
tab_meta %>% filter(grp=="MAB") %>% 
  filter(EV1<0.0) %>% 
  count()

#OTU A
tab_meta %>% filter(grp=="MAB") %>% 
  filter(EV1>0.0) %>% 
  count()
