#load libraries
library(dplyr)
library(e1071)
library(tidyr)


# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)



# Check if a path was provided
if (length(args) < 1) {
  stop("No path provided. Please provide a path as the first argument.")
}

# Assign the first command line argument to 'file_path'
file_path <- args[1]



calculate_summary_stats_by_topology <- function(df) {
  # Calculate the summary statistics grouped by Topology
  df_summary <- df %>%
    group_by(Topology) %>%
    summarise(
      # Means
      Out_mean = mean(Relative_Out, na.rm = TRUE),
      Internal_mean = mean(Relative_Internal, na.rm = TRUE),
      Distance_mean = mean(Relative_Distance, na.rm = TRUE),
      
      # Medians
      Out_median = median(Relative_Out, na.rm = TRUE),
      Internal_median = median(Relative_Internal, na.rm = TRUE),
      Distance_median = median(Relative_Distance, na.rm = TRUE),
      
      # Skewness and Kurtosis
      Out_skewness = skewness(Relative_Out, na.rm = TRUE),
      Internal_skewness = skewness(Relative_Internal, na.rm = TRUE),
      Distance_skewness = skewness(Relative_Distance, na.rm = TRUE),
      
      Out_kurtosis = kurtosis(Relative_Out, na.rm = TRUE),
      Internal_kurtosis = kurtosis(Relative_Internal, na.rm = TRUE),
      Distance_kurtosis = kurtosis(Relative_Distance, na.rm = TRUE),
      
      # Variance
      Out_variance = var(Relative_Out, na.rm = TRUE),
      Internal_variance = var(Relative_Internal, na.rm = TRUE),
      Distance_variance = var(Relative_Distance, na.rm = TRUE),
      
      # Range
      Out_range = max(Relative_Out, na.rm = TRUE) - min(Relative_Out, na.rm = TRUE),
      Internal_range = max(Relative_Internal, na.rm = TRUE) - min(Relative_Internal, na.rm = TRUE),
      Distance_range = max(Relative_Distance, na.rm = TRUE) - min(Relative_Distance, na.rm = TRUE),
      
      # Interquartile Range (IQR)
      Out_IQR = IQR(Relative_Out, na.rm = TRUE),
      Internal_IQR = IQR(Relative_Internal, na.rm = TRUE),
      Distance_IQR = IQR(Relative_Distance, na.rm = TRUE),
      
      # Mean Absolute Deviation (MAD)
      Out_MAD = mad(Relative_Out, na.rm = TRUE),
      Internal_MAD = mad(Relative_Internal, na.rm = TRUE),
      Distance_MAD = mad(Relative_Distance, na.rm = TRUE),
      
      # Percentiles (e.g., 10th and 90th percentiles)
      Out_10th_percentile = quantile(Relative_Out, 0.10, na.rm = TRUE),
      Out_90th_percentile = quantile(Relative_Out, 0.90, na.rm = TRUE),
      Internal_10th_percentile = quantile(Relative_Internal, 0.10, na.rm = TRUE),
      Internal_90th_percentile = quantile(Relative_Internal, 0.90, na.rm = TRUE),
      Distance_10th_percentile = quantile(Relative_Distance, 0.10, na.rm = TRUE),
      Distance_90th_percentile = quantile(Relative_Distance, 0.90, na.rm = TRUE)
    )
  
  # Spread the summary stats by topology (AB, BC, AC) into wide format
  df_wide <- df_summary %>%
    pivot_wider(names_from = Topology, 
                values_from = c(Out_mean, Internal_mean, Distance_mean, 
                                Out_median, Internal_median, Distance_median, 
                                Out_skewness, Internal_skewness, Distance_skewness, 
                                Out_kurtosis, Internal_kurtosis, Distance_kurtosis,
                                Out_variance, Internal_variance, Distance_variance,
                                Out_range, Internal_range, Distance_range,
                                Out_IQR, Internal_IQR, Distance_IQR,
                                Out_MAD, Internal_MAD, Distance_MAD,
                                Out_10th_percentile, Out_90th_percentile,
                                Internal_10th_percentile, Internal_90th_percentile,
                                Distance_10th_percentile, Distance_90th_percentile),
                names_sep = "_")
  
  return(df_wide)
}





file_path <- args[1]

data <- read.csv(file_path, header = T, sep = "\t")

combined_df <- calculate_summary_stats_by_topology(data)


# Interaction terms for AC-related features
combined_df$Internal_mean_AC_Distance_mean_AC <- combined_df$Internal_mean_AC * combined_df$Distance_mean_AC
combined_df$Internal_median_AC_Distance_median_AC <- combined_df$Internal_median_AC * combined_df$Distance_median_AC
combined_df$Distance_range_AC_Distance_skewness_AC <- combined_df$Distance_range_AC * combined_df$Distance_skewness_AC

#Interaction between different branches
combined_df$Distance_mean_AC_Distance_mean_AB <- combined_df$Distance_mean_AC * combined_df$Distance_mean_AB
combined_df$Distance_mean_AC_Distance_mean_BC <- combined_df$Distance_mean_AC * combined_df$Distance_mean_BC
combined_df$Distance_mean_AB_Distance_mean_BC <- combined_df$Distance_mean_AB * combined_df$Distance_mean_BC


combined_df$Out_mean_AC_Out_mean_AB <- combined_df$Out_mean_AC * combined_df$Out_mean_AB
combined_df$Out_mean_AC_Out_mean_BC <- combined_df$Out_mean_AC * combined_df$Out_mean_BC
combined_df$Out_mean_AB_Out_mean_BC <- combined_df$Out_mean_AB * combined_df$Out_mean_BC


combined_df$Internal_mean_AC_Internal_mean_AB <- combined_df$Internal_mean_AC * combined_df$Internal_mean_AB
combined_df$Internal_mean_AC_Internal_mean_BC <- combined_df$Internal_mean_AC * combined_df$Internal_mean_BC
combined_df$Internal_mean_AB_Internal_mean_BC <- combined_df$Internal_mean_AB * combined_df$Internal_mean_BC


# Function to calculate Shannon entropy
shannon_entropy <- function(x) {
  prob <- table(x) / length(x)
  -sum(prob * log(prob))
}

# Apply Shannon entropy to some key features (like AC)
combined_df$Entropy_Out_AC <- apply(combined_df[, c("Out_mean_AC", "Out_median_AC")], 1, shannon_entropy)
combined_df$Entropy_Distance_AC <- apply(combined_df[, c("Distance_mean_AC", "Distance_median_AC")], 1, shannon_entropy)
combined_df$Entropy_Out_BC <- apply(combined_df[, c("Out_mean_BC", "Out_median_BC")], 1, shannon_entropy)
combined_df$Entropy_Distance_BC <- apply(combined_df[, c("Distance_mean_BC", "Distance_median_BC")], 1, shannon_entropy)
combined_df$Entropy_Out_AB <- apply(combined_df[, c("Out_mean_AB", "Out_median_AB")], 1, shannon_entropy)
combined_df$Entropy_Distance_AB <- apply(combined_df[, c("Distance_mean_AB", "Distance_median_AB")], 1, shannon_entropy)


# Create distance ratio features
combined_df$Distance_ratio_mean_AC_AB <- combined_df$Distance_mean_AC / combined_df$Distance_mean_AB
combined_df$Distance_ratio_median_AC_AB <- combined_df$Distance_median_AC / combined_df$Distance_median_AB
combined_df$Distance_ratio_mean_BC_AB <- combined_df$Distance_mean_BC / combined_df$Distance_mean_AB


combined_df$Out_ratio_mean_AC_AB <- combined_df$Out_mean_AC / combined_df$Out_mean_AB
combined_df$Out_ratio_median_AC_AB <- combined_df$Out_median_AC / combined_df$Out_median_AB
combined_df$Out_ratio_mean_BC_AB <- combined_df$Out_mean_BC / combined_df$Out_mean_AB


combined_df$Internal_ratio_mean_AC_AB <- combined_df$Internal_mean_AC / combined_df$Internal_mean_AB
combined_df$Internal_ratio_median_AC_AB <- combined_df$Internal_median_AC / combined_df$Internal_median_AB
combined_df$Internal_ratio_mean_BC_AB <- combined_df$Internal_mean_BC / combined_df$Internal_mean_AB

combined_df$internal_distance_AC <- combined_df$Internal_mean_AC / combined_df$Distance_mean_AC
combined_df$out_distance_AC <- combined_df$Out_mean_AC / combined_df$Distance_mean_AC
combined_df$out_internal_AC <- combined_df$Out_mean_AC / combined_df$Internal_mean_AC


combined_df$internal_distance_AB <- combined_df$Internal_mean_AB / combined_df$Distance_mean_AB
combined_df$out_distance_AB <- combined_df$Out_mean_AB / combined_df$Distance_mean_AB
combined_df$out_internal_AB <- combined_df$Out_mean_AB / combined_df$Internal_mean_AB

combined_df$internal_distance_BC <- combined_df$Internal_mean_BC / combined_df$Distance_mean_BC
combined_df$out_distance_BC <- combined_df$Out_mean_BC / combined_df$Distance_mean_BC
combined_df$out_internal_BC <- combined_df$Out_mean_BC / combined_df$Internal_mean_BC



# Load the random forest model from file
# Load the random forest model from file
model_file <- args[2] 
load(model_file)

# Now the rf_model object is available again and you can make predictions
predictions <- predict(rf_tuned, newdata = combined_df, type = "prob")

output_dir <- args[3]

# Define the output file path
output_file <- file.path(output_dir, "ghostbuster_output.txt")

# Save the predictions to the file
write.table(predictions, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


