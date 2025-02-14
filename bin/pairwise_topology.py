import pandas as pd
import scipy.stats as stats
from scikit_posthocs import posthoc_dunn
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Compute means and perform pairwise comparisons for Topology data.")
parser.add_argument("input_file", type=str, help="Path to the input TSV file.")
parser.add_argument("output_file", type=str, help="Path to the output CSV file for mean values.")
args = parser.parse_args()

# Load the data
data = pd.read_csv(args.input_file, sep="\t")

# Calculate the mean for each value of Topology
means = data.groupby('Topology').mean()
print("Mean values for each topology:")
print(means)

# Perform Kruskal-Wallis H-test for each column
columns_to_test = ['Distance', 'Out', 'Internal']
results = {}

for col in columns_to_test:
    groups = [group[col].values for name, group in data.groupby('Topology')]
    h_stat, p_value = stats.kruskal(*groups)
    results[col] = (h_stat, p_value)
    print(f"\nKruskal-Wallis Test for {col}:")
    print(f"H-statistic: {h_stat}, P-value: {p_value}")

    posthoc = posthoc_dunn(data, val_col=col, group_col='Topology', p_adjust='bonferroni')
    print(f"\nDunn's Test Pairwise Comparisons for {col}:")
    for pair in [("AB", "BC"), ("AC", "BC"), ("AB", "AC")]:
        p_val = posthoc.loc[pair[0], pair[1]]
        larger_mean = pair[0] if means.loc[pair[0], col] > means.loc[pair[1], col] else pair[1]
        print(f"Comparison {pair[0]} vs {pair[1]}: P-value = {p_val}, Larger mean = {larger_mean}")

# Save results to the specified output file
means.to_csv(args.output_file, index=True)
print(f"\nMean values saved to '{args.output_file}'.")

