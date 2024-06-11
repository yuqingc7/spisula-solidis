subset_genepop_individual <- function(GenePopData,indiv=NULL,path){
  
  # read in genepop file
  genepop <- read_tsv(GenePopData, col_names = F)

  #Remove first label
  header <- genepop[1,]
  colnames(header) <- "data"
  genepop <- genepop[-1,]
  colnames(genepop) <- "data"
  
  #ID the rows which flag the Populations
  Pops  <-  which(genepop$data == "POP")
  npops  <-  1:length(Pops)
  
  ## separate the data into the column headers and the rest
  ColumnData <- genepop$data[1:(Pops[1]-1)]
  ColumnData <- gsub("\r","",ColumnData)#remove any hidden carriage returns
  snpData <- genepop[Pops[1]:NROW(genepop),]
  
  ## now subset the individuals
  filtered_snpData <- snpData %>%
    filter(grepl(paste(paste(indiv, collapse = "|"),"|POP"), data))
  
  # Now recompile the GenePop format
  # Convert each item in the list to a data frame
  ColumnData_as_df <- lapply(ColumnData, function(x) data.frame(data = x))
  # Combine all data frames in the list into one data frame row-wise
  combined_df <- do.call(rbind, ColumnData_as_df)
  
  filtered <- as_tibble(rbind(header,rbind(combined_df,filtered_snpData)))
  filtered
  write_tsv(filtered, path, col_names = F)
  
} #End function


# retain the individuals listed
subset_genepop_individual(GenePopData="/workdir/yc2644/Spisula/haps_20231004_n540_grp.gen", 
                          indiv = readLines("/workdir/yc2644/Spisula/subsample_n23/grp_n161.keep"),
                          path = "/workdir/yc2644/Spisula/haps_20231004_grp_n161.gen")


subset_genepop_individual(GenePopData="/workdir/yc2644/Spisula/haps_20231004_n540_OTU.gen", 
                          indiv = readLines("/workdir/yc2644/Spisula/subsample_n23/OTU_A_B_n46.keep"),
                          path = "/workdir/yc2644/Spisula/haps_20231004_OTU_n46.gen")
# manually delete the last POP

subset_genepop_individual(GenePopData="/workdir/yc2644/Spisula/haps_20231004_n540_OTU.gen", 
                          indiv = readLines("/workdir/yc2644/Spisula/subsample_n23/total_n23.keep"),
                          path = "/workdir/yc2644/Spisula/haps_20231004_total_n23.gen")
# manually delete two POPs


