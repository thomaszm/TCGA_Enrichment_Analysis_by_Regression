setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

#Directories
input_dir <- paste(getwd(), "LUAD_Data", sep = "/")
output_dir <- paste(input_dir, "regression_outputs", sep = "/")

mutation_file = "LUAD_mutation_reference_cbio-Xena-firebrowse.txt"
output_file = "LUAD_mutation_reference_transformed.txt"

library("tidyr")    #data
library("dplyr")    #data
library("stringr")  #Replace characters

#MUTATIONS
  
#Upload Sample Mutations
mutation_data <- read.table(paste(input_dir,
                                  mutation_file, 
                                  sep = "/"),
                            header = FALSE,
                            sep = "\t",
                            fill = TRUE,
                            comment.char = "")

#Replace Duplicate Column with Binary
mutation_data[, 3] <- TRUE

#Replace "-" with "_"
mutation_data[,2] <- gsub('-','_', mutation_data[,2])

#Shorten IDs
mutation_data[, 1] <- substr(mutation_data[, 1], 1, 12)

#Remove Duplicates
mutation_data <- mutation_data[!duplicated(mutation_data), ]

#Transform Data
mutation_data <- pivot_wider(mutation_data,
                             id_cols = "V1",
                             names_from = "V2",
                             values_from = "V3",
                             values_fill = FALSE)

#Convert to 0/1
mutation_data <- mutate_if(mutation_data,
                           is.logical,
                           function(x) as.numeric(x))

#Output
write.table(mutation_data, 
            file = paste(input_dir, output_file, sep = "/"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE, 
            sep = "\t")

#paste(input_dir, "MutationEnrichment_Results.txt", sep = "/")
#paste(input_dir, "Enrichment_Results_clin.var.txt", sep = "/")