import argparse
import os
import random
import sys
from pathlib import Path
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate msprime simulation scripts with specified parameters and random values within given ranges."
    )

    # Define command-line arguments
    parser.add_argument("--Ne_min", type=int, help="Minimum effective population size for simulation.")
    parser.add_argument("--Ne_max", type=int, help="Maximum effective population size for simulation.")
    parser.add_argument("--outsplit_min", type=float, help="Minimum time for the split between the outgroup and the triplet in simulations.")
    parser.add_argument("--outsplit_max", type=float, help="Maximum time for the split between the outgroup and the triplet in simulations.")
    parser.add_argument("--Csplit_min", type=float, help="Minimum time for the split between C and (A,B) in simulations.")
    parser.add_argument("--Csplit_max", type=float, help="Maximum time for the split between C and (A,B) in simulations.")
    parser.add_argument("--AB_split_min", type=float, help="Minimum time for the split between A and B in simulations.")
    parser.add_argument("--AB_split_max", type=float, help="Maximum time for the split between A and B in simulations.")
    parser.add_argument("--num_sims", type=int, help="Total number of simulation scripts to generate.")
    parser.add_argument("--migration_rate", type=float, help="Migration rate between populations C and A.")
    parser.add_argument("--generation_time_min", type=float, help="Minimum generation time in years.")
    parser.add_argument("--generation_time_max", type=float, help="Maximum generation time in years.")
    parser.add_argument("--num_loci", type=float, help="Number of loci to generate per simulation replicate.")
    parser.add_argument("--introgression_length", type=float, help="Length of introgression event in generations.")
    parser.add_argument("--introgression_timing", type=float, help="Minimum time period for introgression to have occurred.")
    parser.add_argument("--working_directory", type=str, help="Path for working directory for ghostbuster to create.")
    parser.add_argument("--mode", type=str, help="Mode to run: either nucl_aligned, nucl_unaligned, aa_aligned, aa_unaligned, nucl_trees, or aa_trees.")
    parser.add_argument("--input_trees_file", type=str, help="Path to input trees file, containing one gene tree per line in newick format")
    parser.add_argument("--input_directory", type=str, help="Path to unaligned or aligned nucleotide or amino acid sequences. EACH SEQUENCE MUST BE IN ITS OWN FILE WITH A .FA EXTENSION")
    parser.add_argument("--A_taxon", type=str, help="Name of A taxon as found in tree or fasta files.")
    parser.add_argument("--B_taxon", type=str, help="Name of B taxon as found in tree or fasta files.")
    parser.add_argument("--C_taxon", type=str, help="Name of C taxon as found in tree or fasta files.")
    parser.add_argument("--Out_taxon", type=str, help="Name of Outgroup taxon as found in tree or fasta files.")

    args = parser.parse_args()

    # Check for missing arguments
    required_args = {
        "Ne_min": args.Ne_min,
        "Ne_max": args.Ne_max,
        "outsplit_min": args.outsplit_min,
        "outsplit_max": args.outsplit_max,
        "Csplit_min": args.Csplit_min,
        "Csplit_max": args.Csplit_max,
        "AB_split_min": args.AB_split_min,
        "AB_split_max": args.AB_split_max,
        "num_sims": args.num_sims,
        "migration_rate": args.migration_rate,
        "generation_time_min": args.generation_time_min,
        "generation_time_max": args.generation_time_max,
        "num_loci": args.num_loci,
        "introgression_length": args.introgression_length,
        "introgression_timing": args.introgression_timing,
        "working_directory": args.working_directory,
        "mode": args.mode,
        "A_taxon": args.A_taxon,
        "B_taxon": args.B_taxon,
        "C_taxon": args.C_taxon,
        "Out_taxon": args.Out_taxon
    }

    missing_args = [arg for arg, value in required_args.items() if value is None]
    if missing_args:
        print("Error: The following required arguments are missing:", ', '.join(missing_args))
        parser.print_help()
        exit(1)

#check if mode is valid
    valid_modes = ["nucl_aligned", "nucl_unaligned", "aa_aligned", "aa_unaligned", "nucl_trees", "aa_trees"]

    if args.mode not in valid_modes:
        print(f"Error: Invalid mode '{args.mode}'. Must be one of {', '.join(valid_modes)}.")
        parser.print_help()
        exit(1)

    # Ensure the correct arguments are provided based on the mode
    if args.mode in ["nucl_trees", "aa_trees"] and not args.input_trees_file:
        print("Error: --input_trees_file is required when mode is 'nucl_trees' or 'aa_trees'.")
        parser.print_help()
        exit(1)

    if args.mode not in ["nucl_trees", "aa_trees"] and not args.input_directory:
        print("Error: --input_directory is required when mode is not 'nucl_trees' or 'aa_trees'.")
        parser.print_help()
        exit(1)


    return args
####################Ghost_introgression_only#############################################
def generate_script_content_ghost_only(sim_id, Ne_min, Ne_max, outsplit_min, outsplit_max, Csplit_min, Csplit_max, AB_split_min, AB_split_max, migration_rate, generation_time_min, generation_time_max, introgression_timing, introgression_length):
    # Randomize parameters within specified ranges
    initial_size = random.randint(Ne_min, Ne_max)
    time1 = random.uniform(outsplit_min, outsplit_max)
    time3 = random.uniform(Csplit_min, Csplit_max)
    time4 = random.uniform(AB_split_min, AB_split_max)
    time2 = random.uniform(time1, time3)
    time5 = random.uniform(time4, introgression_timing)
    generation_time = random.uniform(generation_time_min, generation_time_max)
    time1 = time1 / generation_time
    time2 = time2 / generation_time
    time3 = time3 / generation_time
    time4 = time4 / generation_time
    time5 = time5 / generation_time
    time6 = time5 + introgression_length

    # Generate the "ghost only msprime scripts" content
    content = f"""
import msprime

# Define the demography
demography = msprime.Demography()
demography.add_population(name="Ghost", initial_size={initial_size})
demography.add_population(name="A", initial_size={initial_size})
demography.add_population(name="B", initial_size={initial_size})
demography.add_population(name="C", initial_size={initial_size})
demography.add_population(name="AB", initial_size={initial_size})
demography.add_population(name="ABC", initial_size={initial_size})
demography.add_population(name="All", initial_size={initial_size})
demography.add_population(name="Out", initial_size={initial_size})
demography.add_population(name="Ancestral", initial_size={initial_size})

# Define population splits with randomized times
demography.add_population_split(time={time1}, derived=["All", "Out"], ancestral="Ancestral")
demography.add_population_split(time={time2}, derived=["ABC", "Ghost"], ancestral="All")
demography.add_population_split(time={time3}, derived=["AB", "C"], ancestral="ABC")
demography.add_population_split(time={time4}, derived=["A", "B"], ancestral="AB")

# Migration rate changes
demography.add_migration_rate_change(time={time5}, source="Ghost", dest="A", rate={migration_rate})
demography.add_migration_rate_change(time={time6}, source="Ghost", dest="A", rate=0.0)
demography.add_migration_rate_change(time={time5}, source="A", dest="Ghost", rate={migration_rate})
demography.add_migration_rate_change(time={time6}, source="A", dest="Ghost", rate=0.0)

# Sort events in the demography
demography.sort_events()

# Simulate ancestry
ts = msprime.sim_ancestry(
    samples={{"A": 1, "B": 1, "C": 1, "Ghost": 1, "Out": 1}},
    demography=demography,
    recombination_rate=0,
    ploidy=1,
    sequence_length=1000
)

# Output the simulation result in Nexus format
print(ts.as_nexus(precision=3, include_alignments=False))
"""
    return content


###############################sampled introgression only#############################################

def generate_script_content_sampled_only(sim_id, Ne_min, Ne_max, outsplit_min, outsplit_max, Csplit_min, Csplit_max, AB_split_min, AB_split_max, migration_rate, generation_time_min, generation_time_max, introgression_timing, introgression_length):
    # Randomize parameters within specified ranges
    initial_size = random.randint(Ne_min, Ne_max)
    time1 = random.uniform(outsplit_min, outsplit_max)
    time3 = random.uniform(Csplit_min, Csplit_max)
    time4 = random.uniform(AB_split_min, AB_split_max)
    time2 = random.uniform(time1, time3)
    time5 = random.uniform(time4, introgression_timing)
    generation_time = random.uniform(generation_time_min, generation_time_max)
    time1 = time1 / generation_time
    time2 = time2 / generation_time
    time3 = time3 / generation_time
    time4 = time4 / generation_time
    time5 = time5 / generation_time
    time6 = time5 + introgression_length

    content = f"""
import msprime

# Define the demography
demography = msprime.Demography()
demography.add_population(name="Ghost", initial_size={initial_size})
demography.add_population(name="A", initial_size={initial_size})
demography.add_population(name="B", initial_size={initial_size})
demography.add_population(name="C", initial_size={initial_size})
demography.add_population(name="AB", initial_size={initial_size})
demography.add_population(name="ABC", initial_size={initial_size})
demography.add_population(name="All", initial_size={initial_size})
demography.add_population(name="Out", initial_size={initial_size})
demography.add_population(name="Ancestral", initial_size={initial_size})

# Define population splits with randomized times
demography.add_population_split(time={time1}, derived=["All", "Out"], ancestral="Ancestral")
demography.add_population_split(time={time2}, derived=["ABC", "Ghost"], ancestral="All")
demography.add_population_split(time={time3}, derived=["AB", "C"], ancestral="ABC")
demography.add_population_split(time={time4}, derived=["A", "B"], ancestral="AB")

# Migration rate changes
demography.add_migration_rate_change(time={time5}, source="C", dest="A", rate={migration_rate})
demography.add_migration_rate_change(time={time6}, source="C", dest="A", rate=0.0)
demography.add_migration_rate_change(time={time5}, source="A", dest="C", rate={migration_rate})
demography.add_migration_rate_change(time={time6}, source="A", dest="C", rate=0.0)

# Sort events in the demography
demography.sort_events()

# Simulate ancestry
ts = msprime.sim_ancestry(
    samples={{"A": 1, "B": 1, "C": 1, "Ghost": 1, "Out": 1}},
    demography=demography,
    recombination_rate=0,
    ploidy=1,
    sequence_length=1000
)

# Output the simulation result in Nexus format
print(ts.as_nexus(precision=3, include_alignments=False))
"""
    return content

################################################ILS_only########################################################

def generate_script_content_ils_only(sim_id, Ne_min, Ne_max, outsplit_min, outsplit_max, Csplit_min, Csplit_max, AB_split_min, AB_split_max, migration_rate, generation_time_min, generation_time_max, introgression_timing, introgression_length):
    # Randomize parameters within specified ranges
    initial_size = random.randint(Ne_min, Ne_max)
    time1 = random.uniform(outsplit_min, outsplit_max)
    time3 = random.uniform(Csplit_min, Csplit_max)
    time4 = random.uniform(AB_split_min, AB_split_max)
    time2 = random.uniform(time1, time3)
    time5 = random.uniform(time4, introgression_timing)
    generation_time = random.uniform(generation_time_min, generation_time_max)
    time1 = time1 / generation_time
    time2 = time2 / generation_time
    time3 = time3 / generation_time
    time4 = time4 / generation_time
    time5 = time5 / generation_time
    time6 = time5 + introgression_length

    content = f"""
import msprime

# Define the demography
demography = msprime.Demography()
demography.add_population(name="Ghost", initial_size={initial_size})
demography.add_population(name="A", initial_size={initial_size})
demography.add_population(name="B", initial_size={initial_size})
demography.add_population(name="C", initial_size={initial_size})
demography.add_population(name="AB", initial_size={initial_size})
demography.add_population(name="ABC", initial_size={initial_size})
demography.add_population(name="All", initial_size={initial_size})
demography.add_population(name="Out", initial_size={initial_size})
demography.add_population(name="Ancestral", initial_size={initial_size})

# Define population splits with randomized times
demography.add_population_split(time={time1}, derived=["All", "Out"], ancestral="Ancestral")
demography.add_population_split(time={time2}, derived=["ABC", "Ghost"], ancestral="All")
demography.add_population_split(time={time3}, derived=["AB", "C"], ancestral="ABC")
demography.add_population_split(time={time4}, derived=["A", "B"], ancestral="AB")

# Sort events in the demography
demography.sort_events()

# Simulate ancestry
ts = msprime.sim_ancestry(
    samples={{"A": 1, "B": 1, "C": 1, "Ghost": 1, "Out": 1}},
    demography=demography,
    recombination_rate=0,
    ploidy=1,
    sequence_length=1000
)

# Output the simulation result in Nexus format
print(ts.as_nexus(precision=3, include_alignments=False))
"""
    return content


##################Ghost and Sampled Introgression#######################
def generate_script_content_sampled_and_ghost(sim_id, Ne_min, Ne_max, outsplit_min, outsplit_max, Csplit_min, Csplit_max, AB_split_min, AB_split_max, migration_rate, generation_time_min, generation_time_max, introgression_timing, introgression_length):
    # Randomize parameters within specified ranges
    initial_size = random.randint(Ne_min, Ne_max)
    time1 = random.uniform(outsplit_min, outsplit_max)
    time3 = random.uniform(Csplit_min, Csplit_max)
    time4 = random.uniform(AB_split_min, AB_split_max)
    time2 = random.uniform(time1, time3)
    time5 = random.uniform(time4, introgression_timing)
    time7 = random.uniform(time4, introgression_timing)
    generation_time = random.uniform(generation_time_min, generation_time_max)
    time1 = time1 / generation_time
    time2 = time2 / generation_time
    time3 = time3 / generation_time
    time4 = time4 / generation_time
    time5 = time5 / generation_time
    time6 = time5 + introgression_length
    time7 = time7 / generation_time
    time8 = time7 + introgression_length

    content = f"""
import msprime

# Define the demography
demography = msprime.Demography()
demography.add_population(name="Ghost", initial_size={initial_size})
demography.add_population(name="A", initial_size={initial_size})
demography.add_population(name="B", initial_size={initial_size})
demography.add_population(name="C", initial_size={initial_size})
demography.add_population(name="AB", initial_size={initial_size})
demography.add_population(name="ABC", initial_size={initial_size})
demography.add_population(name="All", initial_size={initial_size})
demography.add_population(name="Out", initial_size={initial_size})
demography.add_population(name="Ancestral", initial_size={initial_size})

# Define population splits with randomized times
demography.add_population_split(time={time1}, derived=["All", "Out"], ancestral="Ancestral")
demography.add_population_split(time={time2}, derived=["ABC", "Ghost"], ancestral="All")
demography.add_population_split(time={time3}, derived=["AB", "C"], ancestral="ABC")
demography.add_population_split(time={time4}, derived=["A", "B"], ancestral="AB")

# Migration rate changes
demography.add_migration_rate_change(time={time5}, source="C", dest="A", rate={migration_rate})
demography.add_migration_rate_change(time={time6}, source="C", dest="A", rate=0.0)
demography.add_migration_rate_change(time={time5}, source="A", dest="C", rate={migration_rate})
demography.add_migration_rate_change(time={time6}, source="A", dest="C", rate=0.0)
demography.add_migration_rate_change(time={time5}, source="Ghost", dest="A", rate={migration_rate})
demography.add_migration_rate_change(time={time6}, source="Ghost", dest="A", rate=0.0)
demography.add_migration_rate_change(time={time5}, source="A", dest="Ghost", rate={migration_rate})
demography.add_migration_rate_change(time={time6}, source="A", dest="Ghost", rate=0.0)



# Sort events in the demography
demography.sort_events()

# Simulate ancestry
ts = msprime.sim_ancestry(
    samples={{"A": 1, "B": 1, "C": 1, "Ghost": 1, "Out": 1}},
    demography=demography,
    recombination_rate=0,
    ploidy=1,
    sequence_length=1000
)



# Output the simulation result in Nexus format
print(ts.as_nexus(precision=3, include_alignments=False))
"""
    return content


##############Function to rename treefiles######################

def replace_taxa_names(input_file, output_file, a_taxon, b_taxon, c_taxon, out_taxon):
    """
    Replace taxon names in Newick trees with A, B, C, and Out and save the result.

    Parameters:
    input_file (str): Path to the input file containing Newick trees.
    output_file (str): Path to save the modified Newick trees.
    a_taxon (str): Taxon name to replace with 'A'.
    b_taxon (str): Taxon name to replace with 'B'.
    c_taxon (str): Taxon name to replace with 'C'.
    out_taxon (str): Taxon name to replace with 'Out'.
    """
    try:
        # Read the input file
        with open(input_file, 'r') as infile:
            trees = infile.readlines()

        # Replace taxon names
        modified_trees = []
        for tree in trees:
            tree = tree.strip()
            tree = tree.replace(a_taxon, 'A')
            tree = tree.replace(b_taxon, 'B')
            tree = tree.replace(c_taxon, 'C')
            tree = tree.replace(out_taxon, 'Out')
            modified_trees.append(tree)

        # Write the modified trees to the output file
        with open(output_file, 'w') as outfile:
            outfile.write('\n'.join(modified_trees))

        print(f"Modified trees saved to {output_file}")

    except Exception as e:
        print(f"Error: {e}")




##################Running Simulations############################
def main():
    args = parse_args()

    # Create the working directory specified by --working_directory
    working_dir = Path(args.working_directory).resolve()
    working_dir.mkdir(parents=True, exist_ok=True)

    # Ghost Only
    # Create the output directory if it doesn't exist
    output_dir = Path(working_dir) / "ghost_only_msprime_scripts"
    if not output_dir.exists():
        os.makedirs(output_dir, exist_ok=True)

        # Generate specified number of simulation scripts
        for i in range(1, args.num_sims + 1):
            script_content = generate_script_content_ghost_only(
                sim_id=i,
                    Ne_min=args.Ne_min,
                    Ne_max=args.Ne_max,
                    outsplit_min=args.outsplit_min,
                    outsplit_max=args.outsplit_max,
                    Csplit_min=args.Csplit_min,
                    Csplit_max=args.Csplit_max,
                    AB_split_min=args.AB_split_min,
                    AB_split_max=args.AB_split_max,
                    migration_rate=args.migration_rate,
                    generation_time_min=args.generation_time_min,
                    generation_time_max=args.generation_time_max,
                    introgression_timing=args.introgression_timing,
                    introgression_length=args.introgression_length
            )
            script_path = os.path.join(output_dir, f"simulation_{i}.py")
            with open(script_path, "w") as script_file:
                script_file.write(script_content)

            print(f"Generated {args.num_sims} msprime simulation scripts in '{output_dir}' directory.")
    else:
        print(f"Output directory '{output_dir}' already exists. Skipping script generation.")

    # Sampled Introgression
    # Create the output directory if it doesn't exist
    output_dir = Path(working_dir) / "sampled_only_msprime_scripts"
    if not output_dir.exists():
        os.makedirs(output_dir, exist_ok=True)

        # Generate specified number of simulation scripts
        for i in range(1, args.num_sims + 1):
            script_content = generate_script_content_sampled_only(
                sim_id=i,
                    Ne_min=args.Ne_min,
                    Ne_max=args.Ne_max,
                    outsplit_min=args.outsplit_min,
                    outsplit_max=args.outsplit_max,
                    Csplit_min=args.Csplit_min,
                    Csplit_max=args.Csplit_max,
                    AB_split_min=args.AB_split_min,
                    AB_split_max=args.AB_split_max,
                    migration_rate=args.migration_rate,
                    generation_time_min=args.generation_time_min,
                    generation_time_max=args.generation_time_max,
                    introgression_timing=args.introgression_timing,
                    introgression_length=args.introgression_length
            )
            script_path = os.path.join(output_dir, f"simulation_{i}.py")
            with open(script_path, "w") as script_file:
                script_file.write(script_content)

            print(f"Generated {args.num_sims} msprime simulation scripts in '{output_dir}' directory.")
    else:
        print(f"Output directory '{output_dir}' already exists. Skipping script generation.")

    # ILS Only
    # Create the output directory if it doesn't exist
    output_dir = Path(working_dir) / "ils_only_msprime_scripts"
    if not output_dir.exists():
        os.makedirs(output_dir, exist_ok=True)

        # Generate specified number of simulation scripts
        for i in range(1, args.num_sims + 1):
            script_content = generate_script_content_ils_only(
                sim_id=i,
                    Ne_min=args.Ne_min,
                    Ne_max=args.Ne_max,
                    outsplit_min=args.outsplit_min,
                    outsplit_max=args.outsplit_max,
                    Csplit_min=args.Csplit_min,
                    Csplit_max=args.Csplit_max,
                    AB_split_min=args.AB_split_min,
                    AB_split_max=args.AB_split_max,
                    migration_rate=args.migration_rate,
                    generation_time_min=args.generation_time_min,
                    generation_time_max=args.generation_time_max,
                    introgression_timing=args.introgression_timing,
                    introgression_length=args.introgression_length
            )
            script_path = os.path.join(output_dir, f"simulation_{i}.py")
            with open(script_path, "w") as script_file:
                script_file.write(script_content)

            print(f"Generated {args.num_sims} msprime simulation scripts in '{output_dir}' directory.")
    else:
        print(f"Output directory '{output_dir}' already exists. Skipping script generation.")

    # Ghost and Sampled Introgression
    # Create the output directory if it doesn't exist
    output_dir = Path(working_dir) / "ghost_and_sampled_msprime_scripts"

    if not output_dir.exists():
        os.makedirs(output_dir, exist_ok=True)

        # Generate specified number of simulation scripts
        for i in range(1, args.num_sims + 1):
            script_content = generate_script_content_sampled_and_ghost(
                sim_id=i,
                    Ne_min=args.Ne_min,
                    Ne_max=args.Ne_max,
                    outsplit_min=args.outsplit_min,
                    outsplit_max=args.outsplit_max,
                    Csplit_min=args.Csplit_min,
                    Csplit_max=args.Csplit_max,
                    AB_split_min=args.AB_split_min,
                    AB_split_max=args.AB_split_max,
                    migration_rate=args.migration_rate,
                    generation_time_min=args.generation_time_min,
                    generation_time_max=args.generation_time_max,
                    introgression_timing=args.introgression_timing,
                    introgression_length=args.introgression_length
            )
            script_path = os.path.join(output_dir, f"simulation_{i}.py")
            with open(script_path, "w") as script_file:
                script_file.write(script_content)

            print(f"Generated {args.num_sims} msprime simulation scripts in '{output_dir}' directory.")
    else:
        print(f"Output directory '{output_dir}' already exists. Skipping script generation.")
    script_dir = os.path.dirname(os.path.abspath(__file__))


#different options based on mode
    if args.mode in ['nucl_aligned', 'nucl_unaligned', 'nucl_trees']:
    #tweak this to add indels DONE!
    # Assume 'num_loci' was provided to your Python script as an argument, e.g., --num_loci 100
        num_loci = sys.argv[sys.argv.index("--num_loci") + 1]

        # Define the first argument for the bash script
        #first_argument = Path(args.working_directory).resolve() /  "ghost_and_sampled_msprime_scripts"

        # Call the bash script using subprocess
        simulate_script = os.path.join(script_dir, "simulate_training_nucl_dataset.sh")
        
#        subprocess.run(["bash", simulate_script, first_argument, num_loci], check=True)


        # Define the first argument for the bash script
#        first_argument = Path(args.working_directory).resolve() / "sampled_only_msprime_scripts"

        # Call the bash script using subprocess
#        subprocess.run(["bash", simulate_script, first_argument, num_loci], check=True)


        # Define the first argument for the bash script
#        first_argument = Path(args.working_directory).resolve() / "ils_only_msprime_scripts"

        # Call the bash script using subprocess
#        subprocess.run(["bash", simulate_script, first_argument, num_loci], check=True)


        # Define the first argument for the bash script
#        first_argument = Path(args.working_directory).resolve() / "ghost_only_msprime_scripts"

        # Call the bash script using subprocess
#        subprocess.run(["bash", simulate_script, first_argument, num_loci], check=True)

#if args.mode in ['nucl_aligned', 'nucl_unaligned', 'nucl_trees']:
#    num_loci = sys.argv[sys.argv.index("--num_loci") + 1]

    # List of subdirectories to process
        directories = [
    	       "ghost_and_sampled_msprime_scripts",
     	       "sampled_only_msprime_scripts",
    	       "ils_only_msprime_scripts",
    	        "ghost_only_msprime_scripts"
        	    ]
        	
        for subdir in directories:
            first_argument = Path(args.working_directory).resolve() / subdir
            done_file = first_argument / "done.txt"  # Path to the checkpoint file

            if not done_file.exists():
                print(f"Starting simulation for {subdir}")
                try:
                    # Call the bash script
                    subprocess.run(["bash", simulate_script, first_argument, num_loci], check=True)
                    # If successful, create the 'done.txt' checkpoint file
                    with open(done_file, 'w') as f:
                        f.write("Simulation complete.\n")
                    print(f"Simulation for {subdir} complete. Marked as done.")
                except subprocess.CalledProcessError as e:
                    print(f"Error running simulation for {subdir}: {e}")
                    sys.exit(1)
            else:
                print(f"Skipping {subdir} - 'done.txt' exists, simulation already complete.")



        # Assume 'num_loci' was provided to your Python script as an argument, e.g., --num_loci 100
        num_sims = sys.argv[sys.argv.index("--num_sims") + 1]

        # Define the first argument for the bash script
        first_argument = Path(args.working_directory).resolve() /  "ghost_and_sampled_msprime_scripts"

        # Call the bash script using subprocess
        process_trees_script = os.path.join(script_dir, "process_trees.sh")
        subprocess.run(["bash", process_trees_script, first_argument, num_sims], check=True)


        # Define the first argument for the bash script
        first_argument = Path(args.working_directory).resolve() / "sampled_only_msprime_scripts"

        # Call the bash script using subprocess
        subprocess.run(["bash", process_trees_script, first_argument, num_sims], check=True)


        # Define the first argument for the bash script
        first_argument = Path(args.working_directory).resolve() / "ils_only_msprime_scripts"

        # Call the bash script using subprocess
        subprocess.run(["bash", process_trees_script, first_argument, num_sims], check=True)


        # Define the first argument for the bash script
        first_argument = Path(args.working_directory).resolve() / "ghost_only_msprime_scripts"

        # Call the bash script using subprocess
        subprocess.run(["bash", process_trees_script, first_argument, num_sims], check=True)

        train_ghostbuster_script = os.path.join(script_dir, "train_ghostbuster.sh")
        subprocess.run(["bash", train_ghostbuster_script, working_dir, num_sims], check=True)

        if args.mode in ['nucl_trees']:
            recoded_tree_file = os.path.join(args.working_directory, "recoded_trees.txt")
            replace_taxa_names(
            args.input_trees_file,
            recoded_tree_file,
            args.A_taxon,
            args.B_taxon,
            args.C_taxon,
            args.Out_taxon,
        )
            output_file=f"{args.working_directory}/input_tree_stats.csv"
            process_trees_python = os.path.join(script_dir, "relative_lengths.py")
            subprocess.run(["python", process_trees_python, recoded_tree_file, output_file])
            tuned_model_file = f"{working_dir}/tuned_model.RData"
            ghost_buster_apply = os.path.join(script_dir, "apply_ghost_buster_model.R")
            subprocess.run(["Rscript", ghost_buster_apply, output_file, tuned_model_file, working_dir], check=True)

        elif args.mode in ['nucl_aligned']:
            generate_trees_only = os.path.join(script_dir, "generate_alignment_and_trees_only.sh")
            subprocess.run(["bash", generate_trees_only, args.input_directory, working_dir], check=True)
            recoded_tree_file = os.path.join(args.working_directory, "recoded_trees.txt")
            input_trees_path = os.path.join(working_dir, "calculated_trees.txt")
            replace_taxa_names(
            input_trees_path,
            recoded_tree_file,
            args.A_taxon,
            args.B_taxon,
            args.C_taxon,
            args.Out_taxon,
            )
            output_file=f"{args.working_directory}/input_tree_stats.csv"
            process_trees_python = os.path.join(script_dir, "relative_lengths.py")
            subprocess.run(["python", process_trees_python, recoded_tree_file, output_file])
            tuned_model_file = f"{working_dir}/tuned_model.RData"
            ghost_buster_apply = os.path.join(script_dir, "apply_ghost_buster_model.R")
            subprocess.run(["Rscript", ghost_buster_apply, output_file, tuned_model_file, working_dir], check=True)

        elif args.mode in ['nucl_unaligned']:
            generate_alignment_trees = os.path.join(script_dir, "generate_alignment_and_trees_only.sh")
            subprocess.run(["bash", generate_alignment_trees, args.input_directory, str(working_dir)], check=True)
            recoded_tree_file = os.path.join(args.working_directory, "recoded_trees.txt")
            input_trees_path = os.path.join(working_dir, "calculated_trees.txt")
            replace_taxa_names(
            input_trees_path,
            recoded_tree_file,
            args.A_taxon,
            args.B_taxon,
            args.C_taxon,
            args.Out_taxon,
            )
            output_file=f"{args.working_directory}/input_tree_stats.csv"
            process_trees_python = os.path.join(script_dir, "relative_lengths.py")
            subprocess.run(["python", process_trees_python, recoded_tree_file, output_file])
            tuned_model_file = f"{working_dir}/tuned_model.RData"
            ghost_buster_apply = os.path.join(script_dir, "apply_ghost_buster_model.R")
            subprocess.run(["Rscript", ghost_buster_apply, output_file, tuned_model_file, working_dir], check=True)

    else:
        num_loci = sys.argv[sys.argv.index("--num_loci") + 1]

        # Define the first argument for the bash script
        first_argument = Path(working_dir) /  "ghost_and_sampled_msprime_scripts"

        # Call the bash script using subprocess
        simulate_script = os.path.join(script_dir, "simulate_training_aa_dataset.sh")

        subprocess.run(["bash", simulate_script, first_argument, num_loci], check=True)


        # Define the first argument for the bash script
        first_argument = Path(working_dir) / "sampled_only_msprime_scripts"

        # Call the bash script using subprocess
        subprocess.run(["bash", simulate_script, first_argument, num_loci], check=True)


        # Define the first argument for the bash script
        first_argument = Path(working_dir) / "ils_only_msprime_scripts"

        # Call the bash script using subprocess
        subprocess.run(["bash", simulate_script, first_argument, num_loci], check=True)


        # Define the first argument for the bash script
        first_argument = Path(working_dir) / "ghost_only_msprime_scripts"

        # Call the bash script using subprocess
        subprocess.run(["bash", simulate_script, first_argument, num_loci], check=True)


        # Assume 'num_loci' was provided to your Python script as an argument, e.g., --num_loci 100
        num_sims = sys.argv[sys.argv.index("--num_sims") + 1]

        # Define the first argument for the bash script
        first_argument = Path(working_dir) /  "ghost_and_sampled_msprime_scripts"

        # Call the bash script using subprocess
        process_trees_script = os.path.join(script_dir, "process_trees.sh")
        subprocess.run(["bash", process_trees_script, first_argument, num_sims], check=True)


        # Define the first argument for the bash script
        first_argument = Path(working_dir) / "sampled_only_msprime_scripts"

        # Call the bash script using subprocess
        subprocess.run(["bash", process_trees_script, first_argument, num_sims], check=True)


        # Define the first argument for the bash script
        first_argument = Path(working_dir) / "ils_only_msprime_scripts"

        # Call the bash script using subprocess
        subprocess.run(["bash", process_trees_script, first_argument, num_sims], check=True)


        # Define the first argument for the bash script
        first_argument = Path(working_dir) / "ghost_only_msprime_scripts"

        # Call the bash script using subprocess
        subprocess.run(["bash", process_trees_script, first_argument, num_sims], check=True)

        train_ghostbuster_script = os.path.join(script_dir, "train_ghostbuster.sh")
        subprocess.run(["bash", train_ghostbuster_script, working_dir, num_sims], check=True)

        if args.mode in ['aa_trees']:
            recoded_tree_file = os.path.join(args.working_directory, "recoded_trees.txt")
            replace_taxa_names(
            args.input_trees_file,
            recoded_tree_file,
            args.A_taxon,
            args.B_taxon,
            args.C_taxon,
            args.Out_taxon,
            )
            output_file=f"{args.working_directory}/input_tree_stats.csv"
            process_trees_python = os.path.join(script_dir, "relative_lengths.py")
            subprocess.run(["python", process_trees_python, recoded_tree_file, output_file])
            tuned_model_file = f"{working_dir}/tuned_model.RData"
            ghost_buster_apply = os.path.join(script_dir, "apply_ghost_buster_model.R")
            subprocess.run(["Rscript", ghost_buster_apply, output_file, tuned_model_file, working_dir], check=True)

        elif args.mode in ['aa_aligned']:
            generate_trees_only = os.path.join(script_dir, "generate_alignment_and_trees_only.sh")
            subprocess.run(["bash", generate_trees_only, args.input_directory, working_dir], check=True)
            recoded_tree_file = os.path.join(args.working_directory, "recoded_trees.txt")
            input_trees_path = os.path.join(working_dir, "calculated_trees.txt")
            replace_taxa_names(
            input_trees_path,
            recoded_tree_file,
            args.A_taxon,
            args.B_taxon,
            args.C_taxon,
            args.Out_taxon,
            )
            output_file=f"{args.working_directory}/input_tree_stats.csv"
            process_trees_python = os.path.join(script_dir, "relative_lengths.py")
            subprocess.run(["python", process_trees_python, recoded_tree_file, output_file])
            tuned_model_file = f"{working_dir}/tuned_model.RData"
            ghost_buster_apply = os.path.join(script_dir, "apply_ghost_buster_model.R")
            subprocess.run(["Rscript", ghost_buster_apply, output_file, tuned_model_file, working_dir], check=True)

        elif args.mode in ['aa_unaligned']:
            generate_alignment_trees = os.path.join(script_dir, "generate_alignment_and_trees_only.sh")
            subprocess.run(["bash", generate_alignment_trees, args.input_directory, working_dir], check=True)
            recoded_tree_file = os.path.join(args.working_directory, "recoded_trees.txt")
            input_trees_path = os.path.join(working_dir, "calculated_trees.txt")
            replace_taxa_names(
            input_trees_path,
            recoded_tree_file,
            args.A_taxon,
            args.B_taxon,
            args.C_taxon,
            args.Out_taxon,
            )
            output_file=f"{args.working_directory}/input_tree_stats.csv"
            process_trees_python = os.path.join(script_dir, "relative_lengths.py")
            subprocess.run(["python", process_trees_python, recoded_tree_file, output_file])
            tuned_model_file = f"{working_dir}/tuned_model.RData"
            ghost_buster_apply = os.path.join(script_dir, "apply_ghost_buster_model.R")
            subprocess.run(["Rscript", ghost_buster_apply, output_file, tuned_model_file, working_dir], check=True)

#run the simulations based off of the input
if __name__ == "__main__":
    main()

