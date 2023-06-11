README

File: Mutation_Transform_V1.0.R
Version: V1.0
Created: 11 Apr 2023
Creator: Zachary Thomas
Last Edited: 11 Apr 2023
Last Editor: Zachary Thomas

Files:
	1. Mutatoin_Transform_V1.0.R
	2. Data files

Requirements/Dependencies:
	1. R 4.2.2+
	2. R Package “tidyr 1.2.1+”
	3. R Package “dplyr 1.0.10+”
	4. R Package “stringr”

Accompaniment script to Log_Reg_Enrichment
This script was used to transform the LUAD dataset to look like the AML dataset I was provided. Note additional data processing is done within Log_Reg_Enrichment.R. Transformation Steps:
	1. Upload
	2. Replace Column 3 with all True Values (for transformation)
	3. Replace “.” In gene names with “_” (regression mis interprets period as “everything else”
	4. ID Correction
	5. Duplicate removal
	6. Transformation
	7. Convert T/F -> 1/0
	8. Output

Example:
	TCGA-55-6975-01	RFWD2	RFWD2
	TCGA-55-6975-01	PAPPA2	PAPPA2 

	Becomes

				RFWD2	PAPP2
	TCGA-55-6975	1		1
	

Outputs:
Tab separated text file containing transformed mutation data as 
