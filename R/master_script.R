readVCF <- function(vcf) {
  #' Load VCF file
  #' 
  #' @description This function load a VCF file into R environment
  #'
  #' There is no need to remove the header lines as the function
  #' will handle them
  #' 
  #' @param vcf path. Path to VCF file to load
  #' @usage readVCF(vcf)
  #' @return The input VCF as a dataframe object.
  #' @note Header line should start with "#CHROM" for it to
  #' work properly
  #' @examples
  #' readVCF("path/to/file.vcf")
  tmp_vcf<-readLines(vcf)
  tmp_vcf_data<-read.table(vcf, stringsAsFactors = FALSE)
  tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
  vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
  names(tmp_vcf_data)<-vcf_names
  return(tmp_vcf_data)
}

readGeneDat <- function(file_path) {
  # Read the file skipping the first line and the first character of the second line
  tmp_data <- read.table(file_path, sep = "\t", header = TRUE, skip = 1, fill = TRUE,
                         colClasses = "character", comment.char = "")
  
  # If there are rows in the data frame and the second row is a character vector
  if (nrow(tmp_data) > 0 && is.character(tmp_data[2, ])) {
    # Remove the first character of the second row
    tmp_data[2, ] <- substr(tmp_data[2, ], 2, nchar(tmp_data[2, ]))
  }
  
  return(tmp_data)
}
