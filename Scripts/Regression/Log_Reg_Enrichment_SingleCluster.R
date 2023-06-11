setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

library("tidyr")    #data
library("dplyr")    #data
library("stats")    #Regression
library("MASS")     #Regression
library("pscl")     #pR2
library("car")      #VIF
library("ggplot2")  #plotting
library("logistf")  #firth
library("corrr")    #Correlation
library("stringr")  #Replace characters

#SETTINGS
cancer <- "LUAD"            #Set Cancer for graphs
pvalue_threshold <- 0.05    #Sig P Value threshold for enrichment filtering
min_cluster_size <- 15      #Set minimum cluster size for regression
direction <- "both"         #Step direction: "both", "forward", "backward"
step_trace_value <- FALSE   #Perform step trace?
coef_threshold <- 10        #Coefficient Thresholds to perform Firth's Method
err_threshold <- 100        #Std. Error Thresholds to perform Firth's Method
vif_threshold <- 10         #VIF Threshold for severe correlation
pearson_threshold <- 0.9    #Pearson Threshold for severe correlation
min_mut_size <- 3           #Mutation n cutoff
select_cluster <- "R1.V8"   #Cluster to regress

#FILES
#Cluster assignment file
cluster_file <- "LUAD_MergedResult.txt"
#TCGA Mutations File
mutation_file <- "LUAD_mutation_reference_transformed.txt"
#TCGA Clinical Metadata file
clinical_file <- "TCGA-CDR-SupplementalTableS1.txt"
#Mutation Enrichments File
mut_enr_file <- "MutationEnrichment_Results.txt"
#Clinical Metadata Enrichments File
clin_enr_file <- "Enrichment_Results_clin.var.txt"

#Directories
input_dir <- paste(getwd(), paste(cancer, "Data", sep = "_"), sep = "/")
output_dir <- paste(input_dir, "regression_outputs", sep = "/")
coef_out <- paste(output_dir, "coefficients", sep = "/")
vif_out <- paste(output_dir, "vifs", sep = "/")
corr_out <- paste(output_dir, "correlation_matrix", sep = "/")
formulas_out <- paste(output_dir, "formulas", sep = "/")
step_pseudo_r2_out <- paste(output_dir, "step_pseudo_r2", sep = "/")
aic_out <- paste(output_dir, "aic", sep = "/")

#Create Directories if they don't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
if (!dir.exists(coef_out)) {
  dir.create(coef_out)
}
if (!dir.exists(corr_out)) {
  dir.create(corr_out)
}
if (!dir.exists(vif_out)) {
  dir.create(vif_out)
}
if (!dir.exists(formulas_out)) {
  dir.create(formulas_out)
}
if (!dir.exists(step_pseudo_r2_out)) {
  dir.create(step_pseudo_r2_out)
}
if (!dir.exists(aic_out)) {
  dir.create(aic_out)
}

#Function to upload cluster data
cluster_upload <- function(cluster_file) {
  #Upload File
  cluster_assignments <- read.table(cluster_file,
                                    header = TRUE,
                                    sep = "\t",
                                    fill = TRUE,
                                    comment.char = "")
  
  #Shorten IDs
  cluster_assignments[, 1] <- substr(cluster_assignments[, 1], 1, 12)
  
  #COMBINE SEVERAL SAMPLES FROM 1 PATIENT
  cluster_assignments <- cluster_assignments %>%
    #Group samples from same patient together
    group_by(across(1)) %>%
    #Obtain the max value for each cluster by group
    #Patient is assigned to a cluster if at least 1 sample is in the cluster
    summarise(across(all_of(names(.)[2:length(cluster_assignments)]), max))
  
  #REMOVE CLUSTERS WITH LOW N
  #Isolate numeric (cluster membership) and ID columns
  cluster_subset <- unlist(lapply(cluster_assignments, is.numeric))
  cluster_num <- cluster_assignments[, cluster_subset]
  cluster_ids <- cluster_assignments[, !cluster_subset]
  #Select only clusters larger than indicated minimum cluster size
  cluster_num <- cluster_num[, colSums(cluster_num) >= min_cluster_size]
  #Recombine IDs and cluster memberships
  cluster_assignments <- cbind(cluster_ids, cluster_num)
  
  #cluster_count <- colSums(cluster_num)
  
  #Output
  return(cluster_assignments)
  
}

#Function to upload clinical metadata
clinical_upload <- function(clinical_file) {
  #Upload Clinical Metadata
  metadata <- read.table(clinical_file,
                         header = TRUE,
                         sep = "\t",
                         fill = TRUE,
                         na.strings = c("", NA),
                         comment.char = "")
  
  #Shorten IDs
  metadata[ ,1] <- substr(metadata[ ,1], 1, 12)
  
  #Remove Patients not Assigned to Clusters
  metadata <- filter(metadata, metadata[, 1] %in% cluster_assignments[, 1])
  
  #Remove duplicates
  metadata <- metadata[!duplicated(metadata[ ,1]), ]
  
  #Replace unknown with NA
  metadata[metadata == "UNK"] <- NA
  metadata[metadata == "#N/A"] <- NA
  
  #Filter columns with only 1 value other than na
  metadata <- Filter(function(x) (length(unique(x[!is.na(x)])) > 1), metadata)
  
  #Remove columns where more than half of rows are NA
  w <- c()
  for (i in 1:ncol(metadata)) {
    w[i] <- sum(is.na(metadata[i])) > 0.5 * nrow(metadata)
  }
  metadata <- metadata[,!w]
  
  #Replace "-" with "_"
  for (a in 2:ncol(metadata)) {
    metadata[,a] <- gsub('-','_', metadata[,a])
  }
  
  #Output
  return(metadata)
}

#Function to upload mutation data
mutation_upload <- function(mutations_file) {
  #Upload Sample Mutations
  mutation_data <- read.table(mutations_file,
                              header = TRUE,
                              sep = "\t",
                              fill = TRUE,
                              comment.char = "")
  
  #Remove Duplicates
  mutation_data <- mutation_data[!duplicated(mutation_data), ]
  
  #Remove Patients Not Assigned to Clusters
  mutation_data <- filter(mutation_data,
                          mutation_data[, 1]
                          %in%
                            cluster_assignments[, 1])
  
  #REMOVE MUTATIONS WITH LOW N
  #Isolate numeric (cluster membership) and ID columns
  mutation_subset <- unlist(lapply(mutation_data, is.numeric))
  mutation_num <- mutation_data[, mutation_subset]
  mutation_ids <- mutation_data[, !mutation_subset]
  #Select only clusters larger than indicated minimum cluster size
  mutation_num <- mutation_num[, colSums(mutation_num) > min_mut_size]
  #Recombine IDs and cluster memberships
  mutation_data <- cbind(mutation_ids, mutation_num)
  
  #Output
  return(mutation_data)
}

#function to upload mutation enrichments
mut_enr_upload <- function(mutation_enrichment_file) {
  #Upload Data
  mut_enr_data <- read.table(mutation_enrichment_file,
                             header = TRUE,
                             sep = "\t",
                             fill = TRUE,
                             comment.char = "")

  #Keep Only Variables with Significant PValues
  mut_enr_data <- filter(mut_enr_data,
                         mut_enr_data$adjp.value < pvalue_threshold)
  
  #Replace "-" with "_"
  mut_enr_data <- mut_enr_data %>% 
    mutate_all(list(~str_replace_all(., "-", "_")))
  
  #Correct Cluster Names
  mut_enr_data$Cluster <- gsub("_", ".", mut_enr_data$Cluster)
  
  #Filter only clusters kept in Cluster Function
  mut_enr_data <- filter(mut_enr_data,
                         mut_enr_data$Cluster
                         %in%
                           colnames(cluster_assignments))
  
  #Remove duplicate cluster-category pairings
  mut_enr_data <- mut_enr_data[!duplicated(mut_enr_data[, 1:2]), ]
  
  #Add binary column
  mut_enr_data["Binary"] <- TRUE
  
  #Reshape enrichment data such that each row is a cluster
  mut_enr_data <- pivot_wider(mut_enr_data,
                              id_cols = "Cluster",
                              names_from = "Mutations",
                              values_from = "Binary",
                              values_fill = FALSE)
  
  #Output
  return(mut_enr_data)
}

#Function to upload clinical enrichmentsz
clin_enr_upload <- function(clinical_enrichment_file) {
  #Upload Data
  clin_enr_data <- read.table(clinical_enrichment_file,
                              header = TRUE,
                              sep = "\t",
                              fill = TRUE,
                              comment.char = "")
  
  #Keep Only Variables with Significant PValues
  clin_enr_data <- filter(clin_enr_data,
                          clin_enr_data$adjp.value < pvalue_threshold)
  
  #Replace "-" with "_"
  #clin_enr_data[,1] <- gsub("-", "_", clin_enr_data$Cluster)
  #clin_enr_data <- clin_enr_data[,1:2] %>% 
    #mutate_all(list(~str_replace_all(., "-", "_")))
  
  #Correct Cluster Names
  clin_enr_data$Cluster <- gsub("-", ".", clin_enr_data$Cluster)
  
  #Filter only clusters kept in Cluster Function
  clin_enr_data <- filter(clin_enr_data,
                          clin_enr_data$Cluster
                          %in%
                            colnames(cluster_assignments))
  
  #Change Values to Match Data
  clin_enr_data$Mutations <- toupper(clin_enr_data$Mutations)
  
  #REPLACE VARIABLE NAMES WITH CATEGORIES
  #List of variables
  sig_clin_vars <- clin_enr_data[, 1]
  #Initialize empty list of categories
  sig_clin_categories <- c()
  #Loop through variables
  for (i in sig_clin_vars) {
    #Obtain name of category containing the variable
    category <- grep(i, clin_data)
    
    #If category exists, append
    if (length(category) > 0) {
      #Append to list of categories
      sig_clin_categories <- append(sig_clin_categories,
                                    colnames(clin_data[category]))
    } else {
      #If variable was not in dataset (no cateogry), insert NA
      sig_clin_categories <- append(sig_clin_categories, NA)
    }
  }
  #Replace variables with categories
  clin_enr_data$Mutations <- sig_clin_categories
  #Remove duplicate cluster-category pairings
  clin_enr_data <- clin_enr_data[!duplicated(clin_enr_data[, 1:2]), ]
  #Remove NA Mutations
  clin_enr_data <- clin_enr_data[!is.na(clin_enr_data[,1]),]
  
  #Add binary column
  clin_enr_data["Binary"] <- TRUE
  
  #Reshape enrichment data such that each row is a cluster
  clin_enr_data <- pivot_wider(clin_enr_data,
                               id_cols = "Cluster",
                               names_from = "Mutations",
                               values_from = "Binary",
                               values_fill = FALSE)
  
  #Output
  return(clin_enr_data)
}

#Function to Merge Data
data_merge <- function(data1, data2) {
  data <- merge(data1,
                data2,
                by.x = 1,
                by.y = 1,
                all.x = FALSE,
                all.y = TRUE)
  #Output
  return(data)
}

#Function to detect correlated variables
correlation_matrix <- function(cdata) {
  
  #Create correlation matrix with 1 side & diagonal set to NAs
  cor_matrix <- shave(correlate(cdata[enr_index],
                                method = "pearson",
                                diagonal = NA))
  
  #Obtain ordered list of variables
  factors <- factor(cor_matrix$term,
                    levels <- c(cor_matrix$term),
                    ordered = TRUE)
  
  #Convert to long format for ggplot
  cor_matrix <- stretch(cor_matrix)
  
  #Set order of variables for ggplot
  cor_matrix$x <- factor(cor_matrix$x,
                         ordered = TRUE,
                         levels = factors)
  
  #Set order of variables for ggplot, reversed to graph lower triangle
  cor_matrix$y <- factor(cor_matrix$y,
                         ordered = TRUE,
                         levels = rev(factors))
  
  #Empty list to store removed variables
  vars_removed <- c()
  if (length(enr_index) > 1) {
  
    #Loop through all matrix entries
    for (j in 1:nrow(cor_matrix)) {
      #If value is real and above correlation threshold
      if (is.na(cor_matrix$r[j]) == FALSE && cor_matrix[j,3] > pearson_threshold) {
        #Save first variable in correlated pair for later removal
        vars_removed <- append(vars_removed, as.character(cor_matrix$x[j]))
        #Flag
        print(paste(cor_matrix$x[j],
                    "&",
                    cor_matrix$y[j],
                    ": severe correlation (Pearson Coefficient >",
                    pearson_threshold,
                    ")"))
        print(paste(cor_matrix$x[j], "will be removed from dataset"))
      }
    } 
  } else {
    vars_removed <- NA
  }
  
  #Create Heatmap
  heatmap <- ggplot(data = cor_matrix,
                    aes(x = x, y = y, fill = r)) +
    geom_tile(color = "white") +
    ggtitle(paste(cancer, "Cluster: ", select_cluster, "Mutation Correlations")) +
    #geom_text(aes(label = round(r, 3)),
    #          color = "black",
    #          size = 4) +
    scale_fill_gradient2(low = "blue",
                         high = "red",
                         mid = "white",
                         midpoint = 0,
                         limit = c(-1,1),
                         space = "Lab",
                         name="Pearson\nCorrelation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     size = 12,
                                     hjust = 1)) +
    coord_fixed()
  
  if (length(enr_index) > 1) {
    #Save & Remove Correlated Data
    removed_data <- cbind(data[,1], data[vars_removed])
    
    #Update index of enrichments to remove correlated variables
    enr_index <- enr_index[!enr_index %in% vars_removed]
  } else {
    removed_data <- "NA"
    enr_index <- enr_index
  }
  
  #Output
  return(list(matrix = cor_matrix, 
              heatmap = heatmap, 
              removed_variables = vars_removed,
              removed_data = removed_data,
              enr_index = enr_index))
}

#Function to calculate VIFs
calc_vifs <- function(log_model) {
  
  if (length(enr_index[[i]]) > 1) {
    
    #Caculate VIF
    vifs <- as.data.frame(vif(log_model))
    
    #Reformat for ggplot
    vifs <- cbind(rownames(vifs), data.frame(vifs, row.names=NULL))
    colnames(vifs) <- c('Variables', 'VIF')
    
    #Graph VIFs
    vifs_graph <- ggplot(data = vifs, aes(x = reorder(Variables, VIF), y = VIF)) +
      ggtitle(paste(cancer,
                    "Cluster: ",
                    select_cluster, "Variance Inflation Factors")) +
      geom_bar(stat = "identity") +
      labs(x = "Variables") +
      geom_hline(yintercept = 5,
                 linetype = "dashed",
                 color = "green") +
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 1,
                                       size = 12,
                                       hjust = 1)) +
      geom_hline(yintercept = 10,
                 linetype = "dashed",
                 color = "red")
    
    #Flag variables with severe correlation
    vifs_remove <- c()
    for (k in 1:nrow(vifs)) {
      if (vifs[k,2] > vif_threshold) {
        
        #FLag severe correlation
        print(paste(vifs[k,1], ": severe correlation (VIF >", vif_threshold, ")"))
        
        #Add variable to list for removal from data
        vifs_remove <- append(vifs_remove, as.character(vifs[i,1]))
      }
    }
    
    #Output
    return(list(vifs = vifs, 
                vifs_graph = vifs_graph, 
                removed_variables = vifs_remove))
  } else {
    vifs = NA
    vifs_graph = NA
    vifs_remove = NA
    
    #Output
    return(list(vifs = vifs, 
                vifs_graph = vifs_graph, 
                removed_variables = vifs_remove))
  }
}

#Upload Data if file entered above
if (is.na(cluster_file) == FALSE) {
  cluster_assignments <- cluster_upload(paste(input_dir, 
                                              cluster_file, 
                                              sep = "/"))
}
if (is.na(mutation_file) == FALSE) {
  mut_data <- mutation_upload(paste(input_dir, mutation_file, sep = "/"))
}
if (is.na(clinical_file) == FALSE) {
  clin_data <- clinical_upload(paste(input_dir, clinical_file, sep = "/"))
}
if (is.na(mut_enr_file) == FALSE) {
  mut_enr_data <- mut_enr_upload(paste(input_dir, mut_enr_file, sep = "/"))
}
if (is.na(clin_enr_file) == FALSE) {
  clin_enr_data <- clin_enr_upload(paste(input_dir, clin_enr_file, sep = "/"))
}

#Merge Raw Data
if (exists("mut_data") == TRUE) {
  #If mutation data exists, merge clusters & mutation data
  data <- data_merge(cluster_assignments, mut_data)
  
  if (exists("clin_data") == TRUE) {
    #If clinical data also exists, append to data
    data <- data_merge(data, clin_data)
    
  }
  
} else if (exists("clin_data") == TRUE) {
  #If only clinical data exists, merge clusters & clinical data
  data <- data_merge(cluster_assignments, clin_data)
  
}

if (exists("clin_enr_data") == TRUE & exists("mut_enr_data") == TRUE) {
  #Merge clinical & mutation enrichments
  enrichment_data <- data_merge(clin_enr_data, mut_enr_data)
  
} else if (exists("clin_enr_data") == FALSE & exists("mut_enr_data") == TRUE) {
  #Enrichment data is only mutation enrichments
  enrichment_data <- mut_enr_data
  
} else if (exists("clin_enr_data") == TRUE & exists("mut_enr_data") == FALSE) {
  #Enrichment data is only clinical enrichments
  enrichment_data <- clin_enr_data
  
}

#If enrichments file exists, use only enriched variables for that cluster. 
#Otherwise, use all variables
if (exists("enrichment_data") == TRUE) {
  #Find Single Cluster in Enrichment Data
  row_index <- which(enrichment_data[, 1] == select_cluster)
  
  #Obtain List of Enrichments for Cluster
  enr_index <- colnames(enrichment_data[which(enrichment_data[row_index, 
                                                              1:ncol(enrichment_data)] == TRUE)])
  
  #Remove indexed variables that were filtered out of data
  enr_index <- enr_index[(enr_index %in% colnames(data))]
  
} else {
  #Enrichments are all non-ID & cluster variables
  enr_index <- colnames(data[,!colnames(data) 
                             %in% 
                               colnames(cluster_assignments)])
}

#Correlation matrix, correlation heatmap, and correlated variables for removal
data_correlations <- correlation_matrix(data)

#Extract new data
enr_index <- data_correlations[[5]]

#Generate Formula for Cluster
full_formula <- paste(select_cluster,
                      " ~ ",
                      paste(colnames(data[enr_index]),
                            collapse = " + "))

#Full Logistic Regression Model
full_model <- glm(full_formula,
                  data = data,
                  family = "binomial",
                  na.action = na.omit)

#Stepwise Regression
step_model <- stepAIC(glm(full_formula,
                          data = data,
                          family = "binomial",
                          na.action = na.omit),
                      direction = direction,
                      trace = step_trace_value)

#Calculate pseudo_r2
full_pseudo_r2 <- pR2(full_model)
step_pseudo_r2 <- pR2(step_model)

#Output step formula
step_formula <- step_model$formula
step_aic <- step_model$aic


#Calculate Variance Inflation Factor
vifs <- calc_vifs(full_model)

#If any Std. Error or Coefficient is high, perform firth's method
large_coef <- any(abs(summary(step_model)$coefficients[,1]) > coef_threshold)
large_err <- any(abs(summary(step_model)$coefficients[,2]) > err_threshold)
if (large_coef == TRUE | large_err == TRUE) {
  
  #Perform firth's method on step model
  firth_reduced <- logistf(step_formula, data = data)
  
  #Save AIC for firth method
  firth_reduced_aic <- extractAIC(firth_reduced)
  
  #Save interesting firth results
  #Equations for Std. Error come from R/summary.logistf.R doc
  results <- tibble(Variables = firth_reduced$terms,
                    coefficients_log_odds = firth_reduced$coefficients,
                    standard_error = diag(firth_reduced$var)^0.5,
                    p_values = firth_reduced$prob)
  
} else {
  #Generate table of interesting stepAIC results
  results <- tibble(Variables = names(step_model$coefficients),
                    coefficients_log_odds = summary(step_model)$coefficients[,1],
                    standard_error = summary(step_model)$coefficients[,2],
                    p_values = summary(step_model)$coefficients[,4])
  
}

#Save interesting results to txt file
write.table(results, 
            file = paste(coef_out, 
                         paste(select_cluster, 
                               "coefficients.txt", 
                               sep = "_"), 
                         sep = "/"), 
            quote = FALSE, 
            row.names = FALSE, 
            sep = "\t")

#Save Correlation Matrix Graph
ggsave(file = paste(corr_out, 
                    paste(select_cluster, 
                          "pearson_matrix.png", 
                          sep = "_"), 
                    sep = "/"),
       device = png,
       plot = data_correlations[2])

#Save Correlation Matrix Data
write.table(data_correlations[1], 
            file = paste(corr_out, 
                         paste(select_cluster, 
                               "pearson_correlations.txt", 
                               sep = "_"), 
                         sep = "/"), 
            quote = FALSE, 
            row.names = FALSE, 
            sep = "\t")

#Save VIF graph
ggsave(file = paste(vif_out, 
                    paste(select_cluster, 
                          "vif_plot.png", 
                          sep = "_"), 
                    sep = "/"),
       device = png,
       plot = vifs[[2]])

#Save pseudo r2 tables
write.table(step_pseudo_r2, 
            file = paste(step_pseudo_r2_out, 
                         paste(select_cluster, 
                               "step_pseudo_r2.txt", 
                               sep = "_"), 
                         sep = "/"), 
            quote = FALSE, 
            row.names = TRUE, 
            sep = "\t")

#Save step formulas
writeLines(as.character(step_formula), con = paste(formulas_out, 
                                                    "step_formulas.txt",
                                                    sep = "/"))

#Save full formulas
writeLines(as.character(full_formula), con = paste(formulas_out, 
                                                    "full_formulas.txt",
                                                    sep = "/"))

#Save StepAIC table
write.table(step_aic, 
            file = paste(aic_out, 
                         "step_aic.txt", 
                         sep = "/"), 
            quote = FALSE, 
            row.names = FALSE, 
            sep = "\t")

#Save firthAIC table
write.table(firth_reduced_aic, 
            file = paste(aic_out, 
                         "firth_aic.txt", 
                         sep = "/"), 
            quote = FALSE, 
            row.names = FALSE, 
            sep = "\t")

rm(list = ls())
cat("\014")

