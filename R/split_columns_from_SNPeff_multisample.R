# Code snippets: Split columns from SNPeff (Alvaro Benítez, September 2022)

library(dplyr)
library(stringr)
source("C:/Users/varo_/Desktop/Álvaro/GitHub/Scripts-UAL/R/master_script.R")

# This script considers that the target vcf is modified to contain "SPLIT" or any other string in the place of "|" in the ANN field

# Loading the file

vcf <- readVCF(vcf="C:/Users/varo_/Desktop/Álvaro/TMP/VCF_cold_format.vcf")
colnames(vcf) <- paste("V", 1:ncol(vcf), sep = "")
#vcf <- read.delim("C:/path/to/file.vcf", header=FALSE)

# First, we create an object containing the number of columns to create after splitting the target
ncols <- max(stringr::str_count(na.omit(vcf$V8), fixed("ANN=")))+1

# Now it will generate the colnames for the new columns, they will be called "col" followed by a sequential number
colmn <- paste("col", 1:ncols, sep = "")

# Using separate, we finally split the column. Argument "sep" defines the delimiter
tmp <-
  tidyr::separate(
    data = vcf,
    col = V8,
    sep = "ANN=",
    into = colmn,
    remove = FALSE
  )

# "col2" is the new column storing all the annotation data
# Original delimiter "|" will cause trouble with R, so we change it before this step

ncols2 <- max(stringr::str_count(na.omit(tmp$col2), fixed("SPLIT")))+1
colmn2 <- paste("ann", 1:ncols2, sep = "")

tmp2 <-
  tidyr::separate(
    data = tmp,
    col = col2,
    sep = "SPLIT",
    into = colmn2,
    remove = FALSE
  )

# Many columns are usually generated. However, it is safe to consider that ann4 is the column storing the IDs of interest
# It is time to load the gene_description file to the environment

cpepo_gene_descr_ncbi <- read.delim("C:/Users/varo_/Desktop/Álvaro/UAL/resources/cpepo_ncbi_gene_descr_LOC.txt", header=FALSE)

# Join both dataframes based on the id

tmp3 <- left_join(tmp2, cpepo_gene_descr_ncbi, by = c("ann4"="V1"))

# Optional step: Removing all the undesired columns
final_vcf <- tmp3[,-c(26:(ncols2+10))]
final_vcf <- final_vcf[,-c(9,10)]
# Save to file
write.table(final_vcf, file = "VCF_cold_tolerance_columns_clean.csv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Selecting the QTL Region
final_vcf2 <- subset(final_vcf, final_vcf$V1 == "NC_036647.1")
write.table(final_vcf2, file = "VCF_cold_tolerance_columns_clean_CHR10.csv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
