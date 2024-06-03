# Alvaro Benitez Mateo, 2024
# This script was initially created to gather into a single file the different CSV files resulting from the differential expression analysis of transcriptomic data.

library(dplyr)

# Define the directory containing the CSV files
csv_dir <- "path/to/directory"

# Get a list of all CSV files in the directory
csv_files <- list.files(path = csv_dir, pattern = "*.csv", full.names = TRUE)

# Read all CSV files into a list of dataframes
dfs <- lapply(csv_files, read.csv2)

# Name each dataframe based on the corresponding file name (without the path and extension)
names(dfs) <- basename(csv_files)
names(dfs) <- sub("\\.csv$", "", names(dfs))

# Alternatively, it is possible to make a list out of the dataframes available in the environment
dfs <- Filter(function(x) is.data.frame(x), mget(ls()))

# Select columns 1, 2, and 6 from each dataframe
dfs_clean <- lapply(dfs, function(df) df[, c(1,2,6)])

# Define a function to rename the columns of a dataframe with a suffix
rename_columns <- function(df, suffix) {
  colnames(df)[-1] <- paste0(colnames(df)[-1], "_", suffix)
  return(df)
}

# Apply the renaming function to each dataframe in dfs_clean
# Use the names of the dataframes as suffixes
dfs_clean <- mapply(rename_columns, dfs_clean, names(dfs_clean), SIMPLIFY = FALSE)

# Perform a left join on the dataframes based on the first column ("X")
result <- Reduce(function(x, y) merge(x, y, by = "X", all.x = TRUE), dfs_clean)
