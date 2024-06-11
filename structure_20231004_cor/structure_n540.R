setwd("/local/workdir/yc2644/Spisula/structure_20231004_cor")

library(RColorBrewer)
library(tidyverse)
library(pophelper)
library(viridis)
library(gridExtra)

### Read in data------------------------------------------------------------
# Q tables
n540_k2 <- readQ("./SNP_20231004_n540_cor_out/str_K2_rep8_f", filetype="auto")

# metadata
meta <- left_join(read_tsv("./popmap_n540", col_names = F),
                  read_tsv("../vcf_20231004/metadata_n540", col_names=T),
                  by = join_by(X1 == sample)) %>% 
  dplyr::select(X1, `analytical group`, `collection site`)

names(meta)[1] <- "sample.id"
names(meta)[2] <- "grp"
names(meta)[3] <- "grp2"

grps <- data.frame(meta[,2])
grps[,1] <- as.character(grps[,1])
colnames(grps) <- "Group"

# add sample IDs as rownames to qlist
rownames(n540_k2[[1]]) <- meta$sample.id


### Summarize A and B lists-----------------------------------------------------
str_n540_k2 <- cbind(n540_k2[[1]],grps)
str_n540_k2$Sample <- rownames(str_n540_k2)
rownames(str_n540_k2) <- 1:540
colnames(str_n540_k2) <- c("Cluster1","Cluster2","Group","Sample")
str_n540_k2 <- str_n540_k2[, c("Sample", "Group", "Cluster1", "Cluster2")]
site_codes <- read_tsv("../filter_20231004/popmap_hwe",col_names = F) %>% 
  dplyr::rename(Sample=X1, Site=X2)

str_n540_k2 <- left_join(str_n540_k2,site_codes)

#write_tsv(str_n540_k2,"str_n540_k2.txt", col_names = T)

list_A <- str_n540_k2$Sample[str_n540_k2$Cluster2 < 0.1]
length(list_A)

list_B <- str_n540_k2$Sample[str_n540_k2$Cluster1 < 0.1]
length(list_B)

#writeLines(list_A, "sol_A_n110.txt")
#writeLines(list_B, "sol_B_n415.txt")

# for report --------------------------------------------------------------
n540_k2
grps

# grps_order <- grps %>% mutate(Group = ifelse(Group == "NJ", "E", Group)) %>% #Delmarva
#   mutate(Group = ifelse(Group == "GBE", "B", Group)) %>%
#   mutate(Group = ifelse(Group == "CCB", "A", Group)) %>%
#   mutate(Group = ifelse(Group == "MAB", "D", Group)) %>% #New Jersey
#   mutate(Group = ifelse(Group == "NAT", "C", Group)) %>%
#   mutate(Group = ifelse(Group == "SCC", "F", Group)) %>%
#   mutate(Group = ifelse(Group == "SLI", "G", Group))

grps_order <- grps %>% mutate(Group = ifelse(Group == "Delmarva", "E", Group)) %>% #Delmarva
  mutate(Group = ifelse(Group == "George’s Bank", "B", Group)) %>%
  mutate(Group = ifelse(Group == "Cape Cod Bay", "A", Group)) %>%
  mutate(Group = ifelse(Group == "New Jersey", "D", Group)) %>% #New Jersey
  mutate(Group = ifelse(Group == "Nantucket Shoals", "C", Group)) %>%
  mutate(Group = ifelse(Group == "Southern Coast of Cape Cod", "F", Group)) %>%
  mutate(Group = ifelse(Group == "Southern Long Island - West", "G", Group)) %>% 
  mutate(Group = ifelse(Group == "Southern Long Island - East", "G", Group))

p1 <- plotQ(n540_k2, 
            clustercol=c("#F8766D","#00BFC4"),sortind="Cluster1",
            barsize=1,
            barbordercolour="white",barbordersize=0,
            splab=c("K=2"), splabsize = 12, splabangle = 90,
            showindlab=F, useindlab = F, 
            showgrplab=T, grplab=grps_order,
            grplabsize=3.5,pointsize=6,linesize=7,linealpha=0.2,
            pointcol="white",grplabpos=0.5,linepos=0.5,grplabheight=0.75,
            grplabspacer = 0, 
            ordergrp=T,
            grplabangle = 0,grplabjust=0.5,
            showyaxis=T,showticks=T, panelspacer = 0.4,
            basesize=15, exportplot=F, returnplot=T,
            divgrp="Group", divcol="black", divtype=2, divsize=0.25,
            showlegend=T, legendlab=c("A","B"),
            legendkeysize=10,legendtextsize=10,legendmargin=c(2,2,2,0),
            panelratio=c(9,1))

grid.arrange(p1$plot[[1]], nrow=1, widths=c(10,1))

# higher Ks --------------------------------------------------------------
# Q tables
n540_k3 <- readQ("./SNP_20231004_n540_cor_out/str_K3_rep1_f", filetype="auto")
n540_k4 <- readQ("./SNP_20231004_n540_cor_out/str_K4_rep1_f", filetype="auto")
n540_k5 <- readQ("./SNP_20231004_n540_cor_out/str_K5_rep1_f", filetype="auto")
n540_k6 <- readQ("./SNP_20231004_n540_cor_out/str_K6_rep1_f", filetype="auto")

# add sample IDs as rownames to qlist
rownames(n540_k3[[1]]) <- meta$sample.id
rownames(n540_k4[[1]]) <- meta$sample.id
rownames(n540_k5[[1]]) <- meta$sample.id
rownames(n540_k6[[1]]) <- meta$sample.id

# plot
p2 <- plotQ(c(n540_k3,n540_k4,n540_k5,n540_k6), imgoutput="join",
            #clustercol=c("#F8766D","#00BFC4"),
            sortind="Cluster1",
            barsize=1,
            barbordercolour="white",barbordersize=0,
            splab=c("K=3","K=4","K=5","K=6"), splabsize = 12, splabangle = 90,
            showindlab=F, useindlab = F, sharedindlab=F,
            showgrplab=T, grplab=grps_order,
            grplabsize=3.5,pointsize=6,linesize=7,linealpha=0.2,
            pointcol="white",grplabpos=0.5,linepos=0.5,grplabheight=0.75,
            grplabspacer = 0, 
            ordergrp=T,
            grplabangle = 0,grplabjust=0.5,
            showyaxis=T,showticks=T, panelspacer = 0.4,
            basesize=15, exportplot=F, returnplot=T,
            divgrp="Group", divcol="black", divtype=2, divsize=0.25,
            showlegend=F, #legendlab=c("A","B"),
            #legendkeysize=10,legendtextsize=10,legendmargin=c(2,2,2,0),
            panelratio=c(9,1))

grid.arrange(p2$plot[[1]], nrow=1, widths=c(10,1))
