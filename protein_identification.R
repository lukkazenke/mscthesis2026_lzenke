#### MSc thesis ANK identification script 21 July 2025####

rm(list=ls(all=TRUE))
#IMPORT DATA (adjust strain name)
pdb <- read.table(file = 'prodigal.proteins.susy_bin_2.fa.DIAMOND.pdb.dmnd', sep = '\t', header = FALSE)
trembl <- read.table(file = 'prodigal.proteins.susy_bin_2.fa.DIAMOND.trembl.dmnd', sep = '\t', header = FALSE)
swissprot <- read.table(file = 'prodigal.proteins.susy_bin_2.fa.DIAMOND.swissprot.dmnd', sep = '\t', header = FALSE)
  #get proper column headers
header_names <- "Query_ID,Query_length,Subject_ID,Subject_length,Percentage_of_identical_matches, 
                     Alignment_length,Number_of_gap_openings,Start_of_alignment_in_query,End_of_alignment_in_query, 
                     Start_of_alignment_in_subject,End_of_alignment_in_subject,Query_coverage_per_HSP, 
                     Subject_coverage_per_HSP,Expected_value,Bit_score,Subject_title"
header_vector <- strsplit(header_names, ",")[[1]]
colnames(pdb) <- header_vector
colnames(trembl) <- header_vector
colnames(swissprot) <- header_vector

#QUALITY CONTROL: filtering for e-value
#filtering for e-value
df_list<- list(pdb, trembl, swissprot)
# Filter each data frame and store in a new list
filtered_list <- lapply(df_list, function(df) {
  # Ensure the e-value column is numeric
  df$Expected_value <- as.numeric(df$Expected_value)
  
  # Filter based on threshold 
  #(1.0e-05 is lower cutoff for homologous proteins, so I am using something a bit higher)
  df[df$Expected_value < 1.0e-5, ]
})

# Combine all filtered results into one big data frame
combined_quality_ok <- do.call(rbind, filtered_list)

#save protein file
write.csv(combined_quality_ok, file="SuSy_bin2_proteins_DIA")

#ANKYRIN SEARCH
library(dplyr)
ankyrin_df <- combined_quality_ok %>%
  filter(grepl("ankyrin", Subject_title, ignore.case = TRUE))
write.csv (ankyrin_df, file= "ANKs_pseudovibrioS9")
#ELP_df <- combined_quality_ok %>%
  #filter(grepl("leucin", Subject_title, ignore.case = TRUE))  
