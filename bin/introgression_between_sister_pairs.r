#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if filename is provided
if (length(args) == 0) {
  stop("You must provide the name of the dataset file as a positional argument.")
}

# First argument is the filename
input_file <- args[1]

# Get the current working directory
current_dir <- getwd()

# Construct the full path to the input file
input_file <- file.path(current_dir, input_file)

# Check if the file exists
if (!file.exists(input_file)) {
  stop(paste("The file", input_file, "does not exist."))
}

# Read the dataset
data <- read.csv(input_file, sep = "\t", header = TRUE)

# Filter rows for Topology == "AB"
ab_data <- subset(data, Topology == "AB")

# Check if there's any data for Topology AB
if (nrow(ab_data) == 0) {
  stop("No rows found for Topology == 'AB'. Check your dataset.")
}

# Calculate skewness
library(e1071) # For skewness calculation
skewness_value <- skewness(ab_data$Distance)

# Print skewness results
cat("Skewness of 'Distance' for Topology AC: ")
cat(sprintf("Skewness: %.3f", skewness_value))

# Interpret skewness
if (skewness_value < 0) {
  cat(" The distribution is left-skewed (long tail on the left).\n")
} else if (skewness_value > 0) {
  cat(" The distribution is right-skewed (long tail on the right).\n")
} else {
  cat(" The distribution is symmetric.\n")
}

