#Random forest predictions


# Load necessary packages

# Function to install and load libraries
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# List of required packages
required_packages <- c(
  "dplyr", "e1071", "tidyr", # Libraries you mentioned
  "ggplot2", "caret", "randomForest","parallel","gplots" # Add more libraries as needed
)

# Install and load the packages
install_and_load(required_packages)


#read in training dataset

args <- commandArgs(trailingOnly = TRUE)



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


####reading in "both"

# Set the path to the "both" directory
directory_path <- args[1]
#directory_path <- "~/Desktop/Ghostbuster_testing/ghostbuster_test/test_run_four/ghost_and_sampled_msprime_scripts/training_set"


# Get a list of all CSV files in the directory
file_list <- list.files(directory_path, pattern = "*.csv", full.names = TRUE)

file_list

# Read each file into a dataframe and store them in a list
df_list <- lapply(file_list, function(file) {
  read.csv(file, header = TRUE, sep = "\t")
})
df_list

# Apply the function to each dataframe in the list
summary_list <- lapply(df_list, calculate_summary_stats_by_topology)

# Combine the summaries into a single dataframe
ghost_and_sampled_large <- bind_rows(summary_list)

# Print the combined summary
ghost_and_sampled_large$classification <- as.factor("ghost_and_sampled")



#reading in sampled

# Set the path to the "both" directory
directory_path <- args[2]
#directory_path <- "~/Desktop/Ghostbuster_testing/ghostbuster_test/test_run_four/sampled_only_msprime_scripts/training_set"

# Get a list of all CSV files in the directory
file_list <- list.files(directory_path, pattern = "*.csv", full.names = TRUE)


# Read each file into a dataframe and store them in a list
df_list <- lapply(file_list, function(file) {
  read.csv(file, header = TRUE, sep = "\t")
})

# Apply the function to each dataframe in the list
summary_list <- lapply(df_list, calculate_summary_stats_by_topology)

# Combine the summaries into a single dataframe
sampled_only_large <- bind_rows(summary_list)

# Print the combined summary
sampled_only_large$classification <- as.factor("sampled_only")



#reading in ILS

# Set the path to the "both" directory
directory_path <- args[3]
#directory_path <- "~/Desktop/Ghostbuster_testing/ghostbuster_test/test_run_four/ils_only_msprime_scripts/training_set"

# Get a list of all CSV files in the directory
file_list <- list.files(directory_path, pattern = "*.csv", full.names = TRUE)


# Read each file into a dataframe and store them in a list
df_list <- lapply(file_list, function(file) {
  read.csv(file, header = TRUE, sep = "\t")
})

# Apply the function to each dataframe in the list
summary_list <- lapply(df_list, calculate_summary_stats_by_topology)

# Combine the summaries into a single dataframe
ILS_only_large <- bind_rows(summary_list)

# Print the combined summary
ILS_only_large$classification <- as.factor("ILS_only")

#reading in Ghost only


directory_path <- args[4]

#directory_path <- "~/Desktop/Ghostbuster_testing/ghostbuster_test/test_run_four/ghost_only_msprime_scripts/training_set"
# Get a list of all CSV files in the directory
file_list <- list.files(directory_path, pattern = "*.csv", full.names = TRUE)


# Read each file into a dataframe and store them in a list
df_list <- lapply(file_list, function(file) {
  read.csv(file, header = TRUE, sep = "\t")
})

# Apply the function to each dataframe in the list
summary_list <- lapply(df_list, calculate_summary_stats_by_topology)

# Combine the summaries into a single dataframe
ghost_only_large <- bind_rows(summary_list)

# Print the combined summary
ghost_only_large$classification <- as.factor("ghost_only")



# Combine all the dataframes into one
combined_df <- bind_rows(ghost_only_large, ILS_only_large, sampled_only_large, ghost_and_sampled_large)
#combined_df <- bind_rows(ghost_only, ILS_only, sampled_only)


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



combined_df <- na.omit(combined_df)
combined_df$classification <- as.factor(combined_df$classification)
control <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
results <- rfe(combined_df[, -1], combined_df$classification, sizes = c(10, 20, 30, 40), rfeControl = control)


# Check the selected features
selected_features <- results$optVariables

print(selected_features)

# Create a new data frame with only selected features and the target variable
final_data <- combined_df[, c(selected_features, "classification")]




set.seed(42)
trainIndex <- createDataPartition(final_data$classification, p = 0.8, list = FALSE, times = 1)
df_train <-  final_data[trainIndex, ]
df_test <- final_data[-trainIndex, ]



# Impute missing values using roughfix method from randomForest
df_train_imputed <- na.roughfix(df_train)

# Train the random forest model
set.seed(42)
rf_model <- randomForest(classification ~ ., 
                         data = df_train_imputed, 
                         ntree = 1500)




# Predict on the test set
predictions <- predict(rf_model, newdata = df_test)


# View model importance (optional)
#importance(rf_model)
#varImpPlot(rf_model)


# Evaluate the model
# Calculate the range of mtry based on the number of features
total_features <- ncol(df_train_imputed[, selected_features])  # Number of features in the dataset

# Define the range of mtry
tuneGrid <- expand.grid(mtry = seq(1, total_features, by = 1))  # Test all values from 1 to total_features

# Fine-tune the random forest model
rf_tuned <- train(classification ~ ., 
                  data = df_train_imputed[, selected_features],
                  method = "rf",
                  trControl = trainControl(method = "cv", number = 10),
                  tuneGrid = tuneGrid)


#importance(rf_tuned)
#varImpPlot(rf_tuned)

#predictions <- predict(rf_tuned, newdata = df_test)
#confusionMatrix(predictions, df_test$classification)

#confusionMatrix

pdf_file_path <- file.path(directory_path, "model_accuracy.pdf")

# Open the PDF device with larger dimensions
pdf(file = pdf_file_path, width = 25, height = 10)  # Adjust width and height

# Ensure a single plot per page
par(mfrow = c(1, 1)) 

# Save the importance table and varImpPlot
importance_values <- importance(rf_model)
print(importance_values)  # Ensures the importance table is printed in the console
varImpPlot(rf_model)  # Generate variable importance plot

# Generate predictions and save confusion matrix
predictions <- predict(rf_tuned, newdata = df_test)
conf_matrix <- confusionMatrix(predictions, df_test$classification)

# Print the confusion matrix output
# Use textOutput to capture and display textual information
library(gplots)  # For textplot

# Capture importance and confusion matrix in text
importance_text <- capture.output(print(importance_values))
confusion_text <- capture.output(print(conf_matrix))

# Display captured text in the PDF with proper alignment and font size
textplot(paste(importance_text, collapse = "\n"), halign = "left", valign = "top", cex = 0.8)  # Reduce font size with cex
textplot(paste(confusion_text, collapse = "\n"), halign = "left", valign = "top", cex = 0.8)  # Reduce font size with cex

# Close the PDF device
dev.off()

# Save the random forest models
# Save the tuned model
save(rf_tuned, file = file.path(directory_path, "tuned_model.RData"))

# Save the untuned model
save(rf_model, file = file.path(directory_path, "untuned_model.RData"))
