setwd("/local/workdir/yc2644/Spisula")
library(ggforce)
library(scatterpie)
library(tidyverse)

# get percent admixture for each site
str_n540_k2 <- read_tsv("structure_20231004_cor/str_n540_k2.txt")

str_n540_k2 <- str_n540_k2 %>%
  mutate(Class = case_when(
    Cluster1 > 0.9 ~ "Class1",
    Cluster1 >= 0.65 & Cluster1 <= 0.9 ~ "Class2",
    Cluster1 >= 0.35 & Cluster1 < 0.65 ~ "Class3",
    Cluster1 >= 0.1 & Cluster1 < 0.35 ~ "Class4",
    Cluster1 < 0.1 ~ "Class5",
    TRUE ~ "Other"
  ))

unique(str_n540_k2$Class)

all_lat_lon <- read_tsv("all_lat_lon.txt")

df <- left_join(str_n540_k2,all_lat_lon)
missing_rows <- df[rowSums(is.na(df)) > 0, ]
df # lat and long for each exact site
#write_tsv(df[,6:7],"disperseNN/Input/SNP_sol_n540_basic.recode.locs",col_names = F)

df <- df %>%
  mutate(lat = round(lat, 2),
         long = round(long, 2))

df %>%
  distinct(lat, long) %>%
  n_distinct()

df <- df %>%
  dplyr::mutate(SiteNumber = dplyr::dense_rank(paste(lat, long))) %>% 
  dplyr::mutate(SiteNumber = case_when(
    Group == "SLI" ~ 100,
    Group == "SCC" ~ 101,
    Group == "CCB" ~ 102,
    #Group %in% c("SCC", "CCB") ~ 101,
    TRUE ~ SiteNumber
  ))

df_results <- df %>%
  group_by(SiteNumber, Class) %>%
  dplyr::summarize(tempCount = n()) 

df_results_wider <- df_results %>%
  group_by(SiteNumber) %>%
  dplyr::summarize(Count = sum(tempCount),
                   B = sum(tempCount[Class == "Class5"]),
                   Mix = sum(tempCount[Class == "Class3"]),
                   A = sum(tempCount[Class == "Class1"]),
                   `65-90 % A` = sum(tempCount[Class == "Class2"]),
                   `65-90 % B` = sum(tempCount[Class == "Class4"]))  %>%
  ungroup()
df_results_wider[, 3:7] <- df_results_wider[, 3:7]/df_results_wider$Count

long_lat <- df %>% dplyr::select(lat, long, SiteNumber) %>% 
  unique()
long_lat %>%
  filter(SiteNumber == 100) %>%
  summarize(mean_lat = mean(lat), mean_long = mean(long))
long_lat %>%
  filter(SiteNumber == 101) %>%
  summarize(mean_lat = mean(lat), mean_long = mean(long))
long_lat %>%
  filter(SiteNumber == 102) %>%
  summarize(mean_lat = mean(lat), mean_long = mean(long))

long_lat <- long_lat %>% 
  dplyr::mutate(lat = case_when(
    SiteNumber == 102 ~ 41.88167,
    SiteNumber == 101 ~ 41.5525,
    SiteNumber == 100 ~ 40.70286,
    TRUE ~ lat
)) %>% 
  dplyr::mutate(long = case_when(
    SiteNumber == 102 ~ -70.265,
    SiteNumber == 101 ~ -70.665,
    SiteNumber == 100 ~ -73.03762,
    TRUE ~ long
  )) %>% 
  unique()
df_results_wider <- left_join(df_results_wider,long_lat)

# # for report ------------------------------------------------------------

p <- ggplot(usa, aes(x = long, y = lat, group = group)) +
  coord_cartesian(xlim=c(-76, -65), ylim = c(36.5, 42.5))+
  geom_polygon(color = "black",fill = NA) +
  labs(x = "Longitude",y = "Latitude") +
  theme_minimal()+
  # Modifying text size for axis titles and ticks
  theme(axis.title = element_text(size = 16),  # Adjust the font size for axis titles
        axis.text = element_text(size = 14)    # Adjust the font size for axis ticks
  )

p + geom_scatterpie(data=df_results_wider,
                    aes(x=long, y=lat,r=0.1*log(Count)),
                    cols=c("A","65-90 % A","Mix","65-90 % B","B"),
                    color=NA)+
  scale_fill_manual(values=col_pal)+
  guides(fill = guide_legend(title = "Percent Admixture Per Sample between OTUs"))+
  theme(legend.position = "top")+
  geom_scatterpie_legend(0.1*log(df_results_wider$Count), n=3,
                         x=-67, y=38,labeller=function(x) 100*x)

df_results_major <- df_results %>%
  group_by(SiteNumber) %>%
  dplyr::summarize(Count = sum(tempCount),
                   `Majority OTU B` = sum(tempCount[Class == "Class5"])+sum(tempCount[Class == "Class4"]),
                   `Majority OTU A` = sum(tempCount[Class == "Class1"])+sum(tempCount[Class == "Class2"]))  %>%
  ungroup()
df_results_major[, 3:4] <- df_results_major[, 3:4]/df_results_major$Count
df_results_major <- left_join(df_results_major,long_lat)

p + geom_scatterpie(data=df_results_major,
                    aes(x=long, y=lat,r=log(Count)/10),
                    cols=c("Majority OTU A","Majority OTU B"),
                    color=NA, alpha=0.3)+
  scale_fill_manual(values=c("#ff6c6e","#00c6c9"))+
  guides(fill = guide_legend(title = "Geographic Distribution"))+
  theme(legend.position = c(0.55,0.28),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.background = element_rect(fill = "white", color = "black"))+
  geom_scatterpie_legend(log(df_results_wider$Count)/8, n=3,
                         x=-67, y=38,labeller=function(x) round(exp(8 * x), digits = 0))

grp_lon_lat <- read_tsv("/local/workdir/yc2644/Spisula/spacemix/grp/grp_lon_lat")

color=c(rep("#4ac0ff",3), #CCB1-3
        rep("#f9731a",5), #DEL1-5
        rep("#e20ffb",2), #GBE1-2
        "#8505ff", #NAT
        rep("#f4a261",2), #NJ
        rep("#4ac0ff",3), #PLY,PT1-2
        rep("#366bfd",2), #SCC1-2
        rep("#50ff35",4), #SLI1-4
        rep("#a5e12f",2) #SLI5-6
)

p + geom_scatterpie(data=df_results_major,
                    aes(x=long, y=lat,r=log(Count)/8),
                    cols=c("Majority A","Majority B"),
                    color=NA, alpha=0.3)+
  scale_fill_manual(values=c("#ff6c6e","#00c6c9"))+
  guides(fill = guide_legend(title = "Geographic Distribution"))+
  theme(legend.position = c(0.55,0.28),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.background = element_rect(fill = "white", color = "black"))+
  geom_scatterpie_legend(log(df_results_wider$Count)/8, n=3,
                         x=-67, y=38,labeller=function(x) round(exp(8 * x), digits = 0)) + 
  geom_label_repel(data=grp_lon_lat, aes(x=long,y=lat,label=Region,
                                           group=1,color=Region, fontface = "bold"),
                     min.segment.length = 0, point.padding = 0.3) +
  scale_color_manual(values=color)+
  guides(color="none")


# unsequenced MAB
MAB_unseq <- read_sheet("https://docs.google.com/spreadsheets/d/10B0XTQFUzH_Jyodvml-XC_yBuGemeQ3oTqGyKO2IfeA/edit?usp=sharing",sheet=2)

MAB_unseq <- MAB_unseq %>% select(`Sample ID`,`Yuqing's GBS result`,Latitude,Longitude) %>% 
  filter(`Yuqing's GBS result`=="failed sequencing")

MAB_unseq <- MAB_unseq %>%
  dplyr::mutate(SiteNumber = dplyr::dense_rank(paste(Latitude, Longitude)))

MAB_unseq_sum <- MAB_unseq %>%
  group_by(SiteNumber) %>%
  dplyr::summarize(Count = n()) 

long_lat_MAB_unseq <- MAB_unseq %>% dplyr::select(Latitude, Longitude, SiteNumber) %>% 
  unique()
MAB_unseq_sum <- left_join(MAB_unseq_sum ,long_lat_MAB_unseq) %>% 
  mutate(Size=0.1)

p + geom_scatterpie(data=df_results_wider,
                    aes(x=long, y=lat,r=log(Count)/8),
                    cols=c("A","65-90 % A","Mix","65-90 % B","B"),
                    color=NA, alpha=0.3)+
  scale_fill_manual(values=rep("#9b70f8",5))+
  guides(fill = "none")+
  geom_scatterpie_legend(log(df_results_wider$Count)/8, n=3,
                         x=-67, y=38,labeller=function(x) round(exp(8 * x), digits = 0))+
  geom_point(data=MAB_unseq_sum,aes(x=Longitude,y=Latitude),
             group=NA,shape=4,color="red",alpha=0.3,size=1,stroke=0.5)
  
p + geom_point(data=MAB_unseq_sum,aes(x=Longitude,y=Latitude,size=Count),
             group=NA,shape=4,color="red",alpha=0.3,stroke=0.5)

# 
# p + geom_point(data=df_results_wider,
#                aes(x=long, y=lat,group=NA,size=log(Count)),
#                fill="#9b70f8",shape=21,alpha=0.3)+
#   guides(size = guide_legend(title = "Sample size"))+
#   theme(legend.position = c(0.8, 0.3)) + 
#   geom_point(data=MAB_unseq,aes(x=Longitude,y=Latitude),
#              group=NA,shape=1,color="red",alpha=0.3)

