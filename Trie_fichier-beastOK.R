#!/usr/bin/Rscript

# Clean up the environment
rm(list = ls())

# Load required libraries
library(ape)
library(ade4)
library(seqinr)
library(lubridate)

# Source custom functions
source("Temporal_signal_functions/mantelCounding.r")
source("Temporal_signal_functions/randRegression.r")
source("Temporal_signal_functions/tempSignalFunctions.r")

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]

# Extract subname from tree_file
subname1 <- gsub("(\\.*)\\.\\w+$", "\\1", tree_file)  # Corrected regex for file extension
subname <- gsub("Re_root_(.*)", "\\1", subname1)  # Simplified regex to capture the name
tree <- read.tree(file = tree_file)

# Extract and process tip labels for dates
labs <- tree$tip.label
labs <- gsub("/\\d+", "", labs)  # Remove trailing numbers after '/' in labels
tipDates <- numeric(length(labs))  # Initialize numeric vector for dates

for (i in seq_along(labs)) {
  date <- unlist(strsplit(labs[i], "\\|"))[2]  # Extract date portion after '|'
  
  if (!is.na(date)) {
    if (length(unlist(strsplit(date, "-"))) == 3) {
      tipDates[i] <- decimal_date(ymd(date))  # Full YYYY-MM-DD
    } else if (length(unlist(strsplit(date, "-"))) == 2) {
      tipDates[i] <- decimal_date(ymd(paste0(date, "-15")))  # Assume day = 15 for YYYY-MM
    } else if (length(unlist(strsplit(date, "-"))) == 1) {
      tipDates[i] <- as.numeric(date)  # Use as-is for single-year date
    }
  } else {
    stop(paste("Error: Invalid date format in tip label:", labs[i]))
  }
}

# Perform pathogen permutation test
test <- pathogen.permutation.test(
  phy = tree,
  dates = tipDates,
  use.clusters = FALSE,
  auto.round.dates = FALSE,
  nreps = 500
)

# Extract p-values and results
pValues <- test$p_value

# Write results to files
file_no <- paste0("PvaLUE_", subname, ".txt")
writeLines(
  text = paste(tree_file, pValues, test$r, sep = " "),
  con = file_no
)

if (pValues < 0.05) {
  # Significant p-value case
  file_ok <- paste0("Ok_", subname, "_tree.txt")
  writeLines(text = tree_file, con = file_ok)
} else {
  # Non-significant p-value case
  file_no <- paste0("No-Ok_", subname, "_tree.txt")
  writeLines(text = paste(tree_file, pValues, sep = " "), con = file_no)
}
