#!/usr/bin/env Rscript

library(dplyr)
#install.packages("rstatix")
library(rstatix)
library(argparse)

# Argument parser
parser <- ArgumentParser(description = "Perform Chi-squared Test on AC vs BC counts.")
parser$add_argument("input_file", type = "character", help = "Path to the input TSV file")
args <- parser$parse_args()

# Load the input file
input_file <- args$input_file
data <- tryCatch({
  read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}, error = function(e) {
  cat("Error reading file:", e$message, "\n")
  quit(status = 1)
})



# Ensure 'Topology' column exists
if (!"Topology" %in% colnames(data)) {
  cat("Error: 'Topology' column not found in the file.\n")
  quit(status = 1)
}

# Count occurrences
# Count each topology
ac_count <- sum(data$Topology == "AC")
bc_count <- sum(data$Topology == "BC")
ab_count <- sum(data$Topology == "AB")
total_count <- ac_count + bc_count + ab_count

# Check if there's sufficient total data to perform the test
if (total_count > 0) {
  # Perform Two-Proportion Z-Test
  prop_test <- tryCatch(
    prop.test(c(ac_count, bc_count), c(total_count, total_count), correct = FALSE), 
    error = function(e) e  # Capture any errors
  )
  
  # Check if the test ran successfully
  if (inherits(prop_test, "error")) {
    result <- "Error: Unable to perform the two-proportion z-test due to an unexpected issue."
  } else {
    # Display results if the test is successful
    if (prop_test$p.value < 0.01) {
      result <- sprintf("The proportions of AC and BC are significantly different. Proportion Test: Chi-squared: %.4f, P-value: %.4e",
                        prop_test$statistic, prop_test$p.value)
    } else {
      result <- sprintf("The proportions of AC and BC are not significantly different. Proportion Test: Chi-squared: %.4f, P-value: %.4e",
                        prop_test$statistic, prop_test$p.value)
    }
  }
} else {
  result <- "Error: Insufficient data to perform the test. The total count of topologies is 0."
}

# Print the result
print(result)


sim1 <- data

# Perform the KS test for Out

ks_test <- tryCatch(
  ks.test(sim1$Out[sim1$Topology == "AC"], 
          sim1$Out[sim1$Topology == "AB"]), 
  error = function(e) NULL
)

if (is.null(ks_test)) {
  # If the KS test fails, assign the message to Out
  Out <- "Insufficient data to compare branch lengths."
} else {
  # Calculate the medians for each topology
  medians <- aggregate(Out ~ Topology, data = sim1, median)
  medians_filtered <- medians[medians$Topology %in% c("AC", "AB"), ]
  topology_with_highest_median <- medians_filtered[which.max(medians_filtered$Out), "Topology"]

  # Format the output with the KS test results
  Out <- sprintf(
    "The topology with the highest median is:%s; D-statistic:%.4f; P-value:%.4e",
    topology_with_highest_median, ks_test$statistic, ks_test$p.value
  )
}

cat(result, Out, sep = "~")
