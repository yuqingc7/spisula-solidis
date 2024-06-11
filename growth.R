#Plot growth function
#from: https://danstich.github.io/we-r-nycafs/fishStats.html
library(FSAdata) # for data
library(dplyr)   # for filter(), mutate()
library(ggplot2)
library(tidyverse)
library(FSA)
library(nlstools)
library(plotrix)
library("readxl")

my_data <- read_excel("ShellLengthSolidis_by_age_taxon_HBWedit_yc.xlsx", sheet=1) %>% 
  filter(!is.na(`OTU (A, B, H)`)) %>% 
  filter(!is.na(AgeUpdate)) %>% 
  filter(AgeUpdate!="NA") %>% 
  filter(Length_mmUpdate!="NA") %>% 
  filter(!is.na(Length_mmUpdate)) %>% 
  mutate(across(c(Length_mmUpdate, AgeUpdate), as.numeric)) 

table(my_data$`OTU (A, B, H)`)

sol_data <- my_data %>% 
  filter(`OTU (A, B, H)` %in% c("B", "A")) %>% 
  select("Samples", "OTU (A, B, H)","Length_mmUpdate","AgeUpdate")

#Plot length by age
ggplot(sol_data) +
  geom_point(aes(x=AgeUpdate, y=`Length_mmUpdate`, color=`OTU (A, B, H)`)) +
  labs(title="Solidissima", x="Age (years)", y="Shell Length (mm)") + scale_shape_manual(values=c(16, 3)) +
  theme_minimal()

# subset B
lengths_B <- sol_data %>%
  filter(`OTU (A, B, H)`=="B")
vbmod <- `Length_mmUpdate` ~ Linf * (1 - exp(-K * (AgeUpdate - t0)))
starts_B <- vbStarts(formula = `Length_mmUpdate` ~ AgeUpdate, data = lengths_B)
mymod_B <- nls(vbmod, data = lengths_B, start = starts_B)
summary(mymod_B)
# Formula: Length_mmUpdate ~ Linf * (1 - exp(-K * (AgeUpdate - t0)))
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# Linf 141.8485     5.3209  26.659  < 2e-16 ***
#   K      0.1260     0.0412   3.059  0.00244 ** 
#   t0    -6.2396     2.8359  -2.200  0.02862 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 20.39 on 274 degrees of freedom
# 
# Number of iterations to convergence: 13 
# Achieved convergence tolerance: 9.721e-06

# subset A
lengths_A <- sol_data %>%
  filter(`OTU (A, B, H)`=="A")
vbmod <- `Length_mmUpdate` ~ Linf * (1 - exp(-K * (AgeUpdate - t0)))
starts_A <- vbStarts(formula = `Length_mmUpdate` ~ AgeUpdate, data = lengths_A)
mymod_A <- nls(vbmod, data = lengths_A, start = starts_A)
summary(mymod_A)
# Formula: Length_mmUpdate ~ Linf * (1 - exp(-K * (AgeUpdate - t0)))
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# Linf 147.83302    3.57900  41.306  < 2e-16 ***
#   K      0.23985    0.04031   5.950  2.8e-07 ***
#   t0    -0.41658    0.52173  -0.798    0.428    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 12.95 on 49 degrees of freedom
# 
# Number of iterations to convergence: 3 
# Achieved convergence tolerance: 1.795e-07

# functions for fitting (parameters estimated by nls)
f_A <- function(x) 147.83302*(1-exp(-0.23985*(x+0.41658)))
f_B <- function(x) 141.8485*(1-exp(-0.1260*(x+6.2396)))

ggplot(sol_data) +
  geom_point(aes(x=AgeUpdate, y=Length_mmUpdate, color=`OTU (A, B, H)`)) +
  labs(title="Solidissima Growth by OTU Genotype", x="Age (years)", y="Shell Length (mm)") + scale_shape_manual(values=c(16, 3)) +
  theme_bw() + geom_function(fun = f_A, aes(color = "A"))  +
  geom_function(fun = f_B, aes(color = "B")) +
  # Modify legend titles
  labs(color = "OTU") +
  scale_color_manual(labels = c("A", "B"), values = c("#F8766D", "#00BFC4"))
  #geom_label(aes(x=AgeUpdate, y=Length_mmUpdate,label=Samples))



