#!/usr/bin/env python3
import argparse
import os
import re
import sys
import tempfile
import subprocess
import concurrent.futures
import threading
import itertools
from ete3 import Tree

# Global lock to serialize access to the shared tree_stats file.
lock = threading.Lock()

def usage_error():
    sys.stderr.write(
        "Usage: {} --out_taxa <file> --input_trees <file> --species_tree <file> "
        "--output_file <file> [--triplets <file>] [--threads <n>]\n".format(sys.argv[0])
    )
    sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--out_taxa",     required=True,
                        help="File listing each Outgroup taxon, one per line")
    parser.add_argument("--input_trees",  required=True,
                        help="File containing gene trees, one Newick per line")
    parser.add_argument("--species_tree", required=True,
                        help="Species-tree in rooted Newick format")
    parser.add_argument("--triplets",     help="Optional: CSV triplets file, one 't1,t2,t3' per line")
    parser.add_argument("--output_file",  required=True,
                        help="Path to write results")
    parser.add_argument("--threads", type=int, default=os.cpu_count(),
                        help="Parallel threads (default: all cores)")
    args, unknown = parser.parse_known_args()
    if unknown:
        sys.stderr.write("Error: Unknown option(s): {}\n".format(" ".join(unknown)))
        usage_error()
    return args

def read_file_lines(path):
    with open(path) as f:
        return [l.strip() for l in f if l.strip()]

def run_command(cmd, capture_output=False):
    try:
        if capture_output:
            out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            return out.decode("utf-8")
        else:
            subprocess.check_call(cmd)
            return ""
    except FileNotFoundError:
        sys.stderr.write("Error: command not found: {}\n".format(cmd[0]))
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Error running: {}\nStdout/Stderr:\n{}\n".format(" ".join(cmd), e.output.decode("utf-8") if hasattr(e, "output") and e.output else ""))
        # Do not hard exit here; let caller decide when appropriate.
        raise

def check_taxa_in_trees(taxa_file, trees_file):
    all_leaves = set()
    with open(trees_file) as gtf:
        for line in gtf:
            nw = line.strip()
            if not nw:
                continue
            gt = Tree(nw, format=1)
            all_leaves.update(gt.get_leaf_names())
    missing = [t for t in read_file_lines(taxa_file) if t not in all_leaves]
    if missing:
        sys.stderr.write("Error: The following taxa were not found in any gene tree:\n")
        for t in missing:
            sys.stderr.write(f"  - {t}\n")
        sys.exit(1)

def check_species_tree_vs_genetrees(species_tree_file, trees_file):
    st = Tree(open(species_tree_file).read().strip(), format=1)
    species_taxa = set(st.get_leaf_names())
    gene_taxa = set()
    with open(trees_file) as gtf:
        for line in gtf:
            nw = line.strip()
            if not nw:
                continue
            gt = Tree(nw, format=1)
            gene_taxa.update(gt.get_leaf_names())
    extra = gene_taxa - species_taxa
    if extra:
        sys.stderr.write("Error: these gene-tree taxa are missing from species tree:\n")
        for t in sorted(extra):
            sys.stderr.write(f"  - {t}\n")
        sys.exit(1)

def get_triplet_roles(triplet, species_tree_file, trees_file):
    st = Tree(open(species_tree_file).read().strip(), format=1)
    for t in triplet:
        if not st & t:
            sys.stderr.write(f"Error: '{t}' not in species tree\n")
            sys.exit(1)

    mrca = st.get_common_ancestor(*triplet)
    C = None
    for child in mrca.get_children():
        inter = set(child.get_leaf_names()).intersection(triplet)
        if len(inter) == 1:
            C = next(iter(inter))
            break
    if C is None:
        sys.stderr.write(f"Error: cannot determine outgroup for triplet {triplet}\n")
        sys.exit(1)

    t1, t2 = [x for x in triplet if x != C]

    # Decide B vs A using pruned 3-taxon gene trees; fall back gracefully
    prune_cmd = ["nw_prune", "-v", trees_file, C, t2, t1]  # keep exactly these 3; no duplicate C
    try:
        pruned = run_command(prune_cmd, capture_output=True).splitlines()
    except subprocess.CalledProcessError:
        pruned = []

    bc_count = ac_count = 0
    for tree_str in pruned:
        s = tree_str.strip()
        if not s:
            continue
        pt = Tree(s, format=1)
        sides = pt.get_children()
        if len(sides) != 2:
            continue
        for side in sides:
            leaves = set(side.get_leaf_names())
            if leaves == {t1, C}:
                bc_count += 1
            elif leaves == {t2, C}:
                ac_count += 1

    if bc_count == 0 and ac_count == 0:
        # Fallback: keep species-tree order (arbitrary but deterministic)
        B, A = t1, t2
    elif bc_count >= ac_count:
        B, A = t1, t2
    else:
        B, A = t2, t1

    return C, B, A

def exact_replace(text, mapping):
    # Replace whole-token taxon names only
    def repl(m):
        return mapping[m.group(0)]
    # Build a regex like r'\b(?:name1|name2|...)\b'
    names = sorted(mapping.keys(), key=len, reverse=True)
    pat = r'\b(?:' + "|".join(map(re.escape, names)) + r')\b'
    return re.sub(pat, repl, text)

def generate_combinations(args):
    out_taxa = read_file_lines(args.out_taxa)
    if args.triplets:
        raw = []
        for line in read_file_lines(args.triplets):
            parts = [p.strip() for p in line.split(",")]
            if len(parts) != 3:
                sys.stderr.write(f"Error: invalid triplet line: '{line}'\n")
                sys.exit(1)
            raw.append(tuple(parts))
    else:
        st = Tree(open(args.species_tree).read().strip(), format=1)
        taxa = sorted(st.get_leaf_names())
        raw = list(itertools.combinations(taxa, 3))

    combos = []
    seen = set()
    for trip in raw:
        C, B, A = get_triplet_roles(trip, args.species_tree, args.input_trees)
        key = tuple(sorted([A, B, C]))
        if key in seen:
            continue
        seen.add(key)
        for out in out_taxa:
            if out in (A, B, C):
                continue
            combos.append((C, B, A, out))
    return combos

def process_combination(line_tuple, args, script_dir):
    c_taxon, b_taxon, a_taxon, out_taxon = line_tuple
    combination_line = f"{c_taxon}\t{b_taxon}\t{a_taxon}\t{out_taxon}"

    # Keep only these 4 taxa
    nw_cmd = ["nw_prune", "-v", args.input_trees, out_taxon, a_taxon, b_taxon, c_taxon]
    try:
        nw_out = run_command(nw_cmd, capture_output=True)
    except subprocess.CalledProcessError:
        sys.stderr.write(f"[skip] nw_prune failed for {combination_line}\n")
        return None

    if not nw_out.strip():
        sys.stderr.write(f"[skip] No trees contained all of {out_taxon},{a_taxon},{b_taxon},{c_taxon}\n")
        return None

    triplet = exact_replace(
        nw_out,
        {
            a_taxon: "A",
            b_taxon: "B",
            c_taxon: "C",
            out_taxon: "Out",
        },
    )

    # Write temp triplet file
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False, prefix="triplet_", suffix=".txt")
    try:
        tmp.write(triplet)
    finally:
        tmp.close()
    triplet_fn = tmp.name

    stats_file = args.input_trees + ".tree_stats.txt"
    try:
        with lock:
            # Build TSV
            rc = 0
            try:
                run_command(["python", os.path.join(script_dir, "process_trees.py"), triplet_fn, stats_file])
            except subprocess.CalledProcessError:
                rc = 1

            # Verify TSV actually has data (header + at least one row)
            n_lines = 0
            if os.path.exists(stats_file):
                with open(stats_file) as f:
                    for n_lines, _ in enumerate(f, start=1):
                        pass

            if rc != 0 or n_lines <= 1:
                # No usable data — skip R, return NA row
                return f"{combination_line}\tNA\tNA\t0\t0\t0\tNA\tNA\tNA\tNo evidence of introgression (insufficient data)"

            # Run R only when we have data
            bl_stats = run_command(
                ["Rscript", os.path.join(script_dir, "branch_stats2.0.r"), stats_file],
                capture_output=True
            ).strip().replace("\n", " ")

        # Parse R output
        fields = bl_stats.split("~")
        chisq = fields[0].strip() if len(fields) >= 1 else ""
        chisqstat = ((m := re.search(r'Chi[\s\-–—]?squared\s*[:=]\s*([^,]+)', chisq)) and m.group(1).strip()) or chisq
        chisqp = ((m := re.search(r'(?i)\bP[\s\-–—]?value\s*[:=]\s*([^,;)\s]+)', chisq)) and m.group(1).strip()) or ""

        ABcount = int(((m := re.search(r'\bab_count\s*:\s*(\d+)', chisq)) and m.group(1)) or 0)
        ACcount = int(((m := re.search(r'\bac_count\s*:\s*(\d+)', chisq)) and m.group(1)) or 0)
        BCcount = int(((m := re.search(r'\bbc_count\s*:\s*(\d+)', chisq)) and m.group(1)) or 0)

        ks1_f = fields[1] if len(fields) >= 2 else ""
        medians = fields[2] if len(fields) >= 3 else ""

        p_match = re.search(r"\(p=([0-9.eE+-]+)\)", ks1_f)
        pval_str = p_match.group(1) if p_match else ""
        try:
            pval = float(pval_str) if pval_str else float("nan")
        except Exception:
            pval = float("nan")

        THTstat = float(((m:=re.search(r'\bD\s*=\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', ks1_f)) and m.group(1)) or 'nan')

        signif = "are significantly different" in chisq

        out_h = ""
        ab = bc = float("nan")
        if medians:
            parts = medians.split("\t")
            if len(parts) >= 2:
                try:
                    ab = float(parts[0]); bc = float(parts[1])
                except Exception:
                    pass
        if (ab == ab) and (bc == bc):  # both not NaN
            out_h = "AB" if ab > bc else "BC"

        if signif:
            if (pval_str and pval < 0.05 and out_h == "BC"):
                concl = "Evidence of unsampled introgression"
            elif (pval_str and pval < 0.05 and out_h == "AB"):
                concl = "Evidence of sampled introgression, possibly involving inflow."
            else:
                concl = "Evidence of sampled introgression, likely not involving inflow."
        else:
            concl = "No evidence of introgression"

        return f"{combination_line}\t{chisqstat}\t{chisqp}\t{ABcount}\t{BCcount}\t{ACcount}\t{THTstat}\t{pval_str}\t{medians}\t{concl}"

    finally:
        try:
            os.remove(triplet_fn)
        except OSError:
            pass
        # Clean up stats file after R runs or skip — optional:
        try:
            if os.path.exists(stats_file):
                os.remove(stats_file)
        except OSError:
            pass

def main():
    args = parse_args()

    check_taxa_in_trees(args.out_taxa, args.input_trees)
    check_species_tree_vs_genetrees(args.species_tree, args.input_trees)

    if args.triplets:
        for trip in [line.split(",") for line in read_file_lines(args.triplets)]:
            for t in trip:
                if not re.search(r'\b{}\b'.format(re.escape(t)), open(args.species_tree).read()):
                    sys.stderr.write(f"Error: '{t}' not in species tree\n")
                    sys.exit(1)
                if not re.search(r'\b{}\b'.format(re.escape(t)), open(args.input_trees).read()):
                    sys.stderr.write(f"Error: '{t}' not in gene trees\n")
                    sys.exit(1)

    with open(args.output_file, "w") as out_f:
        out_f.write("Taxon_C\tTaxon_B\tTaxon_A\tTaxon_out\tDCT_statistic\tDCT_p_value\tAB_count\tBC_count\tAC_count\tTHT_statistic\tTHT_p_value\tAB_median\tBC_median\tTest_conclusion\n")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    combinations = generate_combinations(args)

    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
        for res in pool.map(lambda c: process_combination(c, args, script_dir), combinations):
            if res:
                rows.append(res)

    with open(args.output_file, "a") as out_f:
        for line in rows:
            out_f.write(line + "\n")

if __name__ == "__main__":
    main()

