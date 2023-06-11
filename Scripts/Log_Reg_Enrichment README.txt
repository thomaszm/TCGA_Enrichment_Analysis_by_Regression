README

File: Log_Reg_Enrichment_V8.0.R
Version: V8.0
Created: 10 Jan 2023
Creator: Zachary Thomas
Last Edited: 23 Apr 2023
Last Editor: Zachary Thomas

Files:
	1. Log_Reg_Enrichment_V8.0_README.txt
	2. Log_Reg_Enrichment_V8.0.R
	3. Mutation_Transform_V1.0.R 
	4. Data files

Requirements/Dependencies:
	1. R 4.2.2+
	2. R Package "MASS"
	3. R Package "tidyverse 1.3.2+"
		3a. R Package "tidyr 1.2.1+"
		3b. R Package "ggplot2 3.4.0+" 
		3c. R Package "dplyr 1.0.10+" 
		3d. R Package "stringr 1.4.1+"
	4. R Package "logistf 1.24.1+"
	5. R Package "corrr 0.4.4+"
	6. R Package "stringr" 
	7. R Package "pscl 1.5.5+" 
	8. R Package "car 3.1-2+"
	9. RStudio (recommended)

Description: 
An R script to perform stepwise logistic regression on TCGA enrichment data. Aside from initial data and file entires, the whole script is automated. The accompanying script Mutation_Transform_V1.0.R can be used to transform the mutation data ahead of time if necessary. Correlated variables are removed by a given Pearson Correlation Coefficient threshold. All clusters for a given cancer are analyzed automatically in sequence; regression formulas are generated cluster-specific. A full model containing all candidate variables is generated first to obtain a baseline and calculate Variance Inflation Factors to monitor correlation in the model. Stepwise regression by StepAIC iteratively adds and removes variables to find the best model by minimum AIC. If correlation coefficients or error is large for any variable, firth's bias reduction method is performed to normalize this. If firth's method is performed those results are exported, elsewise, stepwise regression results are exported.

Steps to Run:
	1. Set main settings
		-Enter the name of cancer
		-p-value threshold for significance from enrichment results (fishers exact test)
		-minimum cluster size to model
		-stepwise regression direction (either forward, backward, or both)
		-stepwise trace (whether to view every step in the stepwise model)
		-coefficient threshold (coefficients in the stepwise model larger than this will trigger firth's regression)
		-standard error threshold (errors in the stepwise model larger than this will trigger firth's regression)
		-vif threshold (VIFs larger than this will be flagged for severe correlation)
		-pearson threshold (correlation coefficients larger than this will trigger the first variable in the pair to be removed from the model)
		-minimum mutation size (mutation #n less than this will be removed from the dataset)
	2. Enter files & directory
		-Enter file names in quotations and with extensions. If not using that filetype, replace the file name with NA
		-Input directory is currently assumed to be in the format of "cancer_Data". Change if desired
		-Output directory will be in the input directory
	3. Execute script

Regression Steps
	1. Packages, settings, files, directories, functions
	2. If specific files were provided, upload that data; if not, skip
	3. Create some empty lists for looping
	4. Merge raw data tables & enrichment data tables (if they exist)... glm can work with multiple data tables, but it's easier with 1.
	5. Main Regression
		1. Generate enrichment index for each cluster. This is a list of variables that are enriched, unique for that cluster
		2. Perform Pearson's Correlations on enrichments per cluster. If correlated variables are found, first in pair is removed 
		3. Enrichment index updated
		4. Generate formula for glm, unique for each cluster, using enrichment indices (all candidate variables)
		5. Perform regression on all candidate variables
		6. Calculate pseudo r2
		7. Calculate VIFs
		8. Perform stepwise regression
		9. Calculate pseudo r2
		10. Store stepwise formulas (list of variables that were retained in the model), and step AICs
		11. Determine if the absolute value of any variables coefficient or standard error is above threshold. If so, perform Firth's Regression
			1. Perform firth's regression, with step model as input
			2. Store firth AIC
			3. Generate titles
			4. Output results
			ELSE
			1. Generate titles
			2. Output results
		*If firth was performed, only firth results are output. If not, stepwise results are output

Outputs:
	1. A folder "regression outputs" in the input data folder
		1a. A folder "aic"
			1a.1 A text file "firth_aic" containing the final aic value for each firth model
			1a.2 A text file "step_aic" containing the final aic value for each stepwise model
		1b. A folder "coefficients"
			1b.1 A text file for each cluster containing a tab-separated table of key results from either the stepwise or firth's regression model; namely, the variable name, coefficient, standard error, and p-value
		1c. A folder "correlation_matrix"
			1c.1. A text file for each cluster containing a tab-separated table of Pearson's Correlation Coefficients for each variable pairing
			1c.2 A correlation matrix (png) heatmap for each variable
		1d. A folder "formulas"
			1d.1 A text file "full_formulas" containing the formulas for each cluster before stepwise regression
			1d.2 A text file "step_formulas" containing the stepwise formulas for each cluster (the variables that were retained in the modeled)
		1e. A folder "step_pseudo_r2"
			1e.1. A text file for each cluster containing a tab-separated table of pseudo r2 values for the model
		1f. A folder "vifs"
			1f.1. A text file for each cluster containing a tab-separated table of Variance Inflation Factors for each variable
			1f.2 A graph (png) of Variance Inflation Factors for each cluster
			
Known Errors:
-Glm: "glm.fit: fitted proabilities numerically 0 or 1 occurred" happens when predicted probabilities are 1 or 0. Can either 1) ignore it, 2) increase sample size, 3) remove outliers
-VIFs: "There are aliased coefficients in the model" happens when correlation coefficient of 1 is present

