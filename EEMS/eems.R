## Part 0: Set up EEMS runs
setwd("/local/workdir/yc2644/Spisula/EEMS")

# all_lat_lon <- read_tsv("../all_lat_lon.txt")
# 
# sol_A <- read_tsv("sol_A/sol_A.order", col_names = F)
# sol_A <- lapply(sol_A, function(x) gsub("^0*\\s*", "", x))
# sol_A <- tibble(Sample = unlist(sol_A))
# sol_A_coord <- left_join(sol_A,all_lat_lon) %>% 
#   select(long, lat)
# write_tsv(sol_A_coord, "sol_A/sol_A.coord", col_names = F)
# 
# sol_B <- read_tsv("sol_B/sol_B.order", col_names = F)
# sol_B <- lapply(sol_B, function(x) gsub("^0*\\s*", "", x))
# sol_B <- tibble(Sample = unlist(sol_B))
# sol_B_coord <- left_join(sol_B,all_lat_lon) %>% 
#   select(long, lat)
# write_tsv(sol_B_coord, "sol_B/sol_B.coord", col_names = F)

## Part 1: Install rEEMSplots
## Check that the current directory contains the rEEMSplots source directory
library("devtools")
#install_github("dipetkov/reemsplots2")

## Possibly change the working directory with setwd()

## Part 2: Generate graphics
library("reemsplots2")
library("rgdal")
library("rworldmap")
library("rworldxtra")
library("broom")

mcmcpath <- "sol_B_20231025"
plots <- make_eems_plots(mcmcpath, longlat = TRUE,
                         add_demes=T, add_outline =T,
                         col_outline="black")

# check if converged
plots$pilogl01

# plot the effective migration surface
plots$mrates01

# plot the effective diversity surface
plots$qrates01

# "Tidy" the map so that each polygon is a "group"
map <- rworldmap::getMap(resolution = "high")
map <- broom::tidy(map)

plots$mrates01+ 
  geom_path(data=map, aes(x = long, y = lat, group = group), 
               color = "black", size=0.25) +
  coord_cartesian(xlim=c(-76, -65), ylim = c(36.5, 43))+
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "bottom",axis.title = element_text(size = 16),  # Adjust the font size for axis titles
        axis.text = element_text(size = 14))+
  scale_x_continuous(breaks = seq(-76, -65, by = 1), labels = seq(-76, -65, by = 1)) +
  scale_y_continuous(breaks = seq(36.5, 43, by = 1), labels = seq(36.5, 43, by = 1))

plots$mrates02+ 
  geom_path(data=map, aes(x = long, y = lat, group = group), 
            color = "black", size=0.25) +
  coord_cartesian(xlim=c(-76, -65), ylim = c(36.5, 43))+
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "bottom",axis.title = element_text(size = 16),  # Adjust the font size for axis titles
        axis.text = element_text(size = 14))+
  scale_x_continuous(breaks = seq(-76, -65, by = 1), labels = seq(-76, -65, by = 1)) +
  scale_y_continuous(breaks = seq(36.5, 43, by = 1), labels = seq(36.5, 43, by = 1))


plots$qrates01+ 
  geom_path(data=map, aes(x = long, y = lat, group = group), 
            color = "black", size=0.25) +
  coord_cartesian(xlim=c(-76, -65), ylim = c(36.5, 43))+
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "bottom",axis.title = element_text(size = 16),  # Adjust the font size for axis titles
        axis.text = element_text(size = 14))+
  scale_x_continuous(breaks = seq(-76, -65, by = 1), labels = seq(-76, -65, by = 1)) +
  scale_y_continuous(breaks = seq(36.5, 43, by = 1), labels = seq(36.5, 43, by = 1))
# blue indicates higher than average diversity and orange – lower than average diversity.

# ---------------

mcmcpath <- "sol_A_109_20231025"
plots_A <- make_eems_plots(mcmcpath, longlat = TRUE,
                         add_demes=T, add_outline =T,
                         col_outline="black")

# check if converged
plots_A$pilogl01

# plot the effective migration surface
plots_A$mrates01

# # "Tidy" the map so that each polygon is a "group"
# map <- rworldmap::getMap(resolution = "high")
# map <- broom::tidy(map)
# 
plots_A$mrates01+ 
  geom_path(data=map, aes(x = long, y = lat, group = group), 
            color = "black", size=0.25) +
  coord_cartesian(xlim=c(-75, -68.5), ylim = c(39.5, 42.5))+
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "bottom",axis.title = element_text(size = 16),  # Adjust the font size for axis titles
        axis.text = element_text(size = 14))+
  scale_x_continuous(breaks = seq(-75, -68.5, by = 1), labels = seq(-75, -68.5, by = 1)) +
  scale_y_continuous(breaks = seq(39.5, 42.5, by = 1), labels = seq(39.5, 42.5, by = 1))

plots_A$qrates01+ 
  geom_path(data=map, aes(x = long, y = lat, group = group), 
            color = "black", size=0.25) +
  coord_cartesian(xlim=c(-75, -68.5), ylim = c(39.5, 42.5))+
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "bottom",axis.title = element_text(size = 16),  # Adjust the font size for axis titles
        axis.text = element_text(size = 14))+
  scale_x_continuous(breaks = seq(-75, -68.5, by = 1), labels = seq(-75, -68.5, by = 1)) +
  scale_y_continuous(breaks = seq(39.5, 42.5, by = 1), labels = seq(39.5, 42.5, by = 1))


# ---------------
mcmcpath <- "sol_A_20231025"
plots_Am <- make_eems_plots(mcmcpath, longlat = TRUE,
                           add_demes=T, add_outline =T,
                           col_outline="black")

# check if converged
plots_Am$pilogl01

# plot the effective migration surface
plots_Am$mrates01

# # "Tidy" the map so that each polygon is a "group"
# map <- rworldmap::getMap(resolution = "high")
# map <- broom::tidy(map)
# 
plots_Am$mrates01+ 
  geom_path(data=map, aes(x = long, y = lat, group = group), 
            color = "black", size=0.25) +
  coord_cartesian(xlim=c(-76, -65), ylim = c(36.5, 43))+
  theme_minimal()+
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "bottom")+
  scale_x_continuous(breaks = seq(-76, -65, by = 1), labels = seq(-76, -65, by = 1)) +
  scale_y_continuous(breaks = seq(36.5, 43, by = 1), labels = seq(36.5, 43, by = 1))


