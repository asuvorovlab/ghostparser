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
      result <- sprintf("The proportions of AC and BC are significantly different. ab_count: %.0f, ac_count: %.0f, bc_count: %.0f. Proportion Test: Chi-squared: %.4f, P-value: %.4e~",
                        ab_count, ac_count, bc_count,
                        prop_test$statistic, prop_test$p.value)
    } else {
      result <- sprintf("The proportions of AC and BC are not significantly different. ab_count: %.0f, ac_count: %.0f, bc_count: %.0f. Proportion Test: Chi-squared: %.4f, P-value: %.4e~",
                        ab_count, ac_count, bc_count,
                        prop_test$statistic, prop_test$p.value)
      
    }
  }
} else {
  result <- "Error: Insufficient data to perform the test. The total count of topologies is 0."
}

# Print the result
#print(result)


sim1 <- data


# Perform the KS test for tree_height



ks_test2 <- tryCatch(
  ks.test(sim1$tree_height[sim1$Topology == "AB"],
          sim1$tree_height[sim1$Topology == "BC"]),
  error = function(e) NULL
)

# … assume sim1, ks_test1 and ks_test2 are already defined …

# Compute medians
medians <- aggregate(tree_height ~ Topology, data = sim1, median)

# Prepare the output string
if (is.null(ks_test2)) {
  tree_height <- "Insufficient data to compare branch lengths."
} else {
  # Which median is higher?
  medians_filtered <- subset(medians, Topology %in% c("AB","BC"))

  # Extract medians for AB and BC (NA if a group is missing)
  ab_med <- medians$tree_height[match("AB", medians$Topology)]
  bc_med <- medians$tree_height[match("BC", medians$Topology)]


  # Now include BOTH D‐statistics and p‐values
    tree_height <- paste0(result, sprintf(
    " AB vs BC: D=%.4f (p=%.4e)",
    ks_test2$statistic, ks_test2$p.value
  ),"~",sprintf("%.6f\t%.6f", ab_med, bc_med))
  

}

# Print
cat(tree_height)

