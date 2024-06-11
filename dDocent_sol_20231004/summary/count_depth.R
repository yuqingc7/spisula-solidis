#remove.packages(rlang)
#install.packages("rlang")
#install.packages("dplyr") #‘dplyr’

library(tidyverse)
library(data.table)
setwd("/workdir/yc2644/Spisula/dDocent_sol_20231004")
# ls *F.bam > LM.list
BAMLIST <- "sample_list_n540.txt"
basedir <- "/workdir/yc2644/Spisula/dDocent_sol_20231004" # Make sure to edit this to match your $BASEDIR
bam_list <- read_lines(paste0(basedir, "/",BAMLIST))

for (i in 1:length(bam_list)){
  
  bamfile = paste0(bam_list[i],'-RG.bam.depth.gz')
  bamhead = bam_list[i]
  # Compute depth stats
  depth <- read_tsv(paste0(basedir, '/', bamfile), col_names = F)$X1
  mean_depth <- mean(depth)
  sd_depth <- sd(depth)
  mean_depth_nonzero <- mean(depth[depth > 0])
  mean_depth_within2sd <- mean(depth[depth < mean_depth + 2 * sd_depth])
  median <- median(depth)
  presence <- as.logical(depth)
  proportion_of_reference_covered <- mean(presence)
  
  # Bind stats into dataframe and store sample-specific per base depth and presence data
  if (i==1){
    output <- data.frame(bamhead, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered)
    #total_depth <- depth
    #total_presence <- presence
  } else {
    output <- rbind(output, cbind(bamhead, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered))
    #total_depth <- total_depth + depth
    #total_presence <- total_presence + presence
  }
}
print(output)

output %>%
  mutate(across(where(is.numeric), round, 3))

write_tsv(output, paste0(BAMLIST, "_depth_per_position_per_sample_summary.tsv"))
#write_lines(total_depth, paste0(bam_list_prefix, "_depth_per_position_all_samples.txt"))
#write_lines(total_presence, paste0(bam_list_prefix, "_presence_per_position_all_samples.txt"))

# Plot the depth distribution
# tibble(total_depth = total_depth, position = 1:length(total_depth))  %>%
#   ggplot(aes(x = position, y = total_depth)) +
#   geom_point(size = 0.1)

# # Total depth per site across all individuals 
# total_depth_summary <- count(tibble(total_depth = total_depth), total_depth)
# total_presence_summary <- count(tibble(total_presence = total_presence), total_presence)
# total_depth_summary %>%
#   ggplot(aes(x = log(total_depth), y = n)) +
#   geom_point()
# total_presence_summary %>%
#   ggplot(aes(x = total_presence, y = n)) +
#   geom_col()

