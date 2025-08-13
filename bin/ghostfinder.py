#!/usr/bin/env python3
"""
ghostfinder (step 5 + triplet CSV filter + finalize & run ghostparser)

Pipeline:
1) Build (A,B,C) from species tree using relational sister logic (A,B)|C.
2) Parallel-count gene-tree topologies rooted at the single Taxon_out.
3) Write ghostfinder_intermediate1.txt (TSV with AB/BC/AC counts).
4) Create ghostfinder_triplets.csv = rows with AC_count > BC_count (no header).
5) Remove ghostfinder_intermediate1.txt.
6) Write ghostfinder_OUT.txt containing the single Taxon_out.
7) Print the ghostparser command (python ghostparser.py ... --triplets ghostfinder_triplets.csv).
8) Run ghostparser.py from the same directory as this script.

Notes:
- --output_file is passed through to ghostparser.py (its output).
- --threads is used for parallel counting here; it is passed to ghostparser.py
  ONLY if provided on the command line.
"""

import argparse
import os
import sys
import re
import shlex
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Set, Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from ete3 import Tree
except ImportError:
    sys.stderr.write("[error] ete3 is required. Install with: pip install ete3\n")
    raise


# -------------------- CLI --------------------

def positive_int(x: str) -> int:
    try:
        v = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid integer: {x}")
    if v <= 0:
        raise argparse.ArgumentTypeError("threads must be a positive integer")
    return v


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Identify (A,B,C) triplets, count gene-tree topologies, and run ghostparser."
    )
    p.add_argument("--ghostparser_output", required=True,
                   help="Path to ghostparser output (TSV/whitespace).")
    p.add_argument("--species_tree", required=True,
                   help="Path to species tree (Newick).")
    p.add_argument("--gene_trees", required=True,
                   help="Path to gene trees (one Newick per line).")
    p.add_argument("--output_file", required=True,
                   help="Output file for ghostparser.py (passed through).")
    # Important: default=None so we can detect if user explicitly provided it.
    p.add_argument("--threads", type=positive_int, default=None,
                   help="Number of worker processes (and forwarded to ghostparser if provided).")
    args = p.parse_args()

    # Normalize and validate
    args.ghostparser_output = Path(args.ghostparser_output).expanduser().resolve()
    args.species_tree = Path(args.species_tree).expanduser().resolve()
    args.gene_trees = Path(args.gene_trees).expanduser().resolve()
    args.output_file = Path(args.output_file).expanduser().resolve()

    for label, path in [
        ("--ghostparser_output", args.ghostparser_output),
        ("--species_tree", args.species_tree),
        ("--gene_trees", args.gene_trees),
    ]:
        if not path.exists():
            p.error(f"{label} not found: {path}")

    return args


# -------------------- ghostparser parsing --------------------

def _split_header_and_choose_policy(header_line: str) -> Tuple[List[str], str]:
    if "\t" in header_line:
        fields = [f.strip() for f in header_line.rstrip("\n").split("\t")]
        return fields, "tab"
    fields = re.split(r"\s+", header_line.strip())
    return fields, "ws"


def _split_line(line: str, policy: str) -> List[str]:
    if policy == "tab":
        return [f.strip() for f in line.rstrip("\n").split("\t")]
    return re.split(r"\s+", line.strip())


def _iter_rows(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        header_fields: List[str] = []
        policy = "tab"
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            header_fields, policy = _split_header_and_choose_policy(line)
            break
        if not header_fields:
            raise ValueError(f"No header found in {path}")

        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cols = _split_line(line, policy)
            if len(cols) < len(header_fields):
                sys.stderr.write(
                    f"[warn] Skipping malformed line (expected {len(header_fields)} cols, got {len(cols)}): {line}"
                )
                continue
            row = {h: cols[i] if i < len(cols) else "" for i, h in enumerate(header_fields)}
            yield row


def extract_pairs_out(ghostparser_path: Path) -> Tuple[List[Tuple[str, str]], str]:
    """
    Returns:
      - ordered list of unique (Taxon_A, Taxon_B) where Test_conclusion == 'Evidence of unsampled introgression'
      - the single Taxon_out in the file (errors if multiple/none)
    """
    required = {"Taxon_A", "Taxon_B", "Taxon_out", "Test_conclusion"}
    seen: Set[Tuple[str, str]] = set()
    ordered_pairs: List[Tuple[str, str]] = []
    outs: Set[str] = set()
    rows_seen = 0

    for row in _iter_rows(ghostparser_path):
        rows_seen += 1
        if not required.issubset(row.keys()):
            missing = required - set(row.keys())
            raise ValueError(f"ghostparser_output missing columns: {', '.join(sorted(missing))}")

        outs.add(row["Taxon_out"].strip())

        if row["Test_conclusion"].strip() == "Evidence of unsampled introgression":
            a = row["Taxon_A"].strip()
            b = row["Taxon_B"].strip()
            pair = (a, b)
            if pair not in seen:
                seen.add(pair)
                ordered_pairs.append(pair)

    if rows_seen == 0:
        raise ValueError("ghostparser_output has no data rows.")
    if len(outs) != 1:
        raise ValueError(
            f"Expected exactly one Taxon_out in ghostparser_output; found {len(outs)}: {', '.join(sorted(outs))}"
        )

    return ordered_pairs, next(iter(outs))


# -------------------- species tree logic --------------------

def load_species_tree(species_tree_path: Path) -> Tree:
    try:
        return Tree(str(species_tree_path), format=1)
    except Exception as e:
        raise RuntimeError(f"Failed to load species tree: {e}")


def leaves_set(t: Tree) -> Set[str]:
    return set(t.get_leaf_names())


def ab_sister_wrt_c(species_tree: Tree, a: str, b: str, c: str) -> bool:
    """
    In the SPECIES tree, A and B are sisters w.r.t. C if, in the induced triplet {A,B,C},
    the unique child of MRCA(A,B,C) containing 2 leaves is {A,B}.
    """
    try:
        m = species_tree.get_common_ancestor(a, b, c)
    except Exception:
        return False
    child_sets = []
    for ch in m.get_children():
        s = set(ch.get_leaf_names()) & {a, b, c}
        if s:
            child_sets.append(s)
    pairs = [s for s in child_sets if len(s) == 2]
    if len(pairs) != 1:
        return False
    return pairs[0] == {a, b}


def build_triplets_from_species(species_tree: Tree,
                                pairs: List[Tuple[str, str]],
                                out_taxon: str) -> List[Tuple[str, str, str]]:
    """
    For each (A,B), find all taxa C (C != A,B,out_taxon) s.t. (A,B)|C in the species tree.
    """
    all_leaves = leaves_set(species_tree)
    if out_taxon not in all_leaves:
        raise RuntimeError(f"Taxon_out '{out_taxon}' not found in species tree.")

    triplets: List[Tuple[str, str, str]] = []
    for a, b in pairs:
        missing = [x for x in (a, b) if x not in all_leaves]
        if missing:
            sys.stderr.write(f"[warn] Skipping pair ({a},{b}): not in species tree: {', '.join(missing)}\n")
            continue
        for c in sorted(all_leaves - {a, b, out_taxon}):
            if ab_sister_wrt_c(species_tree, a, b, c):
                triplets.append((a, b, c))
    return triplets


# -------------------- gene tree topology counting (parallel) --------------------

def classify_triplet_in_rooted_tree(t: Tree, a: str, b: str, c: str) -> str:
    """
    Given a ROOTED gene tree (already rooted at out_taxon), return which pair is sister:
      'AB' if (A,B)|C
      'BC' if (B,C)|A
      'AC' if (A,C)|B
      ''   if unresolved.
    """
    try:
        m = t.get_common_ancestor(a, b, c)
    except Exception:
        return ""

    child_sets = []
    for ch in m.get_children():
        s = set(ch.get_leaf_names()) & {a, b, c}
        if s:
            child_sets.append(s)

    pairs = [s for s in child_sets if len(s) == 2]
    singles_count = sum(len(s) for s in child_sets if len(s) == 1)
    if len(pairs) != 1 or singles_count != 1:
        return ""

    pair = pairs[0]
    if pair == {a, b}:
        return "AB"
    if pair == {b, c}:
        return "BC"
    if pair == {a, c}:
        return "AC"
    return ""


def iter_gene_trees_lines(path: Path) -> Iterator[str]:
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            yield s


def chunked(it: Iterable[str], size: int) -> Iterator[List[str]]:
    chunk: List[str] = []
    for item in it:
        chunk.append(item)
        if len(chunk) >= size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def _worker_chunk(newicks: List[str],
                  out_taxon: str,
                  triplets: List[Tuple[str, str, str]],
                  tri_set_to_indices: Dict[frozenset, List[int]]) -> Tuple[Dict[int, Tuple[int, int, int]],
                                                                           Tuple[int, int, int, int]]:
    """
    Worker: process a chunk of gene trees and return:
      - counts: dict {triplet_index -> (ab_inc, bc_inc, ac_inc)} (only non-zero entries)
      - stats: (processed, parse_errors, missing_out, missing_triplet_taxa)
    """
    local_counts: Dict[int, List[int]] = {}
    processed = 0
    parse_errors = 0
    missing_out = 0
    missing_trip_taxa = 0

    for newick in newicks:
        try:
            t = Tree(newick, format=1)
        except Exception:
            parse_errors += 1
            continue

        leaf_names = set(t.get_leaf_names())
        if out_taxon not in leaf_names:
            missing_out += 1
            continue

        try:
            t.set_outgroup(t & out_taxon)
        except Exception:
            missing_out += 1
            continue

        processed += 1

        present_any = False
        for tri_set, idxs in tri_set_to_indices.items():
            if tri_set.issubset(leaf_names):
                present_any = True
                a, b, c = triplets[idxs[0]]
                which = classify_triplet_in_rooted_tree(t, a, b, c)
                if not which:
                    continue
                i = idxs[0]
                cnt = local_counts.setdefault(i, [0, 0, 0])
                if which == "AB":
                    cnt[0] += 1
                elif which == "BC":
                    cnt[1] += 1
                elif which == "AC":
                    cnt[2] += 1
        if not present_any:
            missing_trip_taxa += 1

    compact = {i: (v[0], v[1], v[2]) for i, v in local_counts.items()}
    return compact, (processed, parse_errors, missing_out, missing_trip_taxa)


def count_topologies_parallel(gene_trees_path: Path,
                              out_taxon: str,
                              triplets: List[Tuple[str, str, str]],
                              threads: int,
                              chunk_size: int = 128) -> Tuple[Dict[Tuple[str, str, str], Tuple[int, int, int]],
                                                              Tuple[int, int, int, int]]:
    """
    Parallel counting across gene trees.

    Returns:
      - counts: {(A,B,C): (AB_count, BC_count, AC_count)}
      - stats:  (processed, parse_errors, missing_out, missing_triplet_taxa)
    """
    n = len(triplets)
    if n == 0:
        return {}, (0, 0, 0, 0)

    tri_set_to_indices: Dict[frozenset, List[int]] = {}
    for i, (a, b, c) in enumerate(triplets):
        tri_set_to_indices.setdefault(frozenset({a, b, c}), []).append(i)

    counts_accum: List[List[int]] = [[0, 0, 0] for _ in range(n)]
    total_processed = total_parse_err = total_missing_out = total_missing_trip = 0

    chunks = list(chunked(iter_gene_trees_lines(gene_trees_path), chunk_size))
    if not chunks:
        return {t: (0, 0, 0) for t in triplets}, (0, 0, 0, 0)

    with ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [
            ex.submit(_worker_chunk, chunk, out_taxon, triplets, tri_set_to_indices)
            for chunk in chunks
        ]
        for fut in as_completed(futures):
            try:
                compact_counts, stats = fut.result()
            except Exception as e:
                sys.stderr.write(f"[warn] Worker failed: {e}\n")
                continue

            for i, (ab, bc, ac) in compact_counts.items():
                counts_accum[i][0] += ab
                counts_accum[i][1] += bc
                counts_accum[i][2] += ac

            p, pe, mo, mt = stats
            total_processed += p
            total_parse_err += pe
            total_missing_out += mo
            total_missing_trip += mt

    final_counts = {triplets[i]: (counts_accum[i][0], counts_accum[i][1], counts_accum[i][2])
                    for i in range(n)}
    return final_counts, (total_processed, total_parse_err, total_missing_out, total_missing_trip)


# -------------------- final CSV filter --------------------

def write_triplets_csv_from_intermediate(intermediate_path: Path, csv_path: Path) -> int:
    """
    Read ghostfinder_intermediate1.txt (TSV with header) and write ghostfinder_triplets.csv (no header),
    with lines Taxon_A,Taxon_B,Taxon_C for rows where AC_count > BC_count.
    Returns number of lines written.
    """
    if not intermediate_path.exists():
        raise FileNotFoundError(f"Intermediate file not found: {intermediate_path}")

    written = 0
    with intermediate_path.open("r", encoding="utf-8") as fh, csv_path.open("w", encoding="utf-8") as out:
        header = fh.readline()
        if not header:
            return 0
        hdr = [h.strip() for h in header.strip().split("\t")]
        col = {name: i for i, name in enumerate(hdr)}
        required = {"Taxon_A", "Taxon_B", "Taxon_C", "BC_count", "AC_count"}
        missing = required - set(col)
        if missing:
            raise ValueError(f"Intermediate missing columns: {', '.join(sorted(missing))}")

        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            try:
                ac = int(parts[col["AC_count"]]); bc = int(parts[col["BC_count"]])
            except (IndexError, ValueError):
                continue
            if ac > bc:
                a = parts[col["Taxon_A"]].strip()
                b = parts[col["Taxon_B"]].strip()
                c = parts[col["Taxon_C"]].strip()
                out.write(f"{a},{b},{c}\n")
                written += 1
    return written


# -------------------- main --------------------

def main() -> int:
    args = parse_args()
    threads_for_count = args.threads or (os.cpu_count() or 1)

    print("[info] Configuration:", file=sys.stderr)
    print(f"  ghostparser_output: {args.ghostparser_output}", file=sys.stderr)
    print(f"  species_tree      : {args.species_tree}", file=sys.stderr)
    print(f"  gene_trees        : {args.gene_trees}", file=sys.stderr)
    print(f"  output_file (for ghostparser): {args.output_file}", file=sys.stderr)
    print(f"  threads (counting): {threads_for_count}", file=sys.stderr)
    if args.threads is not None:
        print(f"  threads (forwarded to ghostparser): {args.threads}", file=sys.stderr)

    # Parse ghostparser output: get (A,B) pairs and single outgroup
    try:
        pairs, out_taxon = extract_pairs_out(args.ghostparser_output)
    except Exception as e:
        print(f"[error] Failed parsing {args.ghostparser_output}: {e}", file=sys.stderr)
        return 2
    print(f"[info] Taxon_out: {out_taxon}", file=sys.stderr)

    # Load species tree and build (A,B,C) triplets using relational logic
    try:
        stree = load_species_tree(args.species_tree)
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        return 3

    try:
        triplets = build_triplets_from_species(stree, pairs, out_taxon)
    except Exception as e:
        print(f"[error] Failed to build triplets from species tree: {e}", file=sys.stderr)
        return 4

    if not triplets:
        print("[warn] No (A,B,C) triplets found from species tree constraints.", file=sys.stderr)

    # Parallel count
    counts, stats = count_topologies_parallel(
        gene_trees_path=args.gene_trees,
        out_taxon=out_taxon,
        triplets=triplets,
        threads=threads_for_count,
        chunk_size=128
    )
    processed, parse_err, missing_out, missing_trip = stats
    print(f"[info] Gene trees processed: {processed}", file=sys.stderr)
    if parse_err:
        print(f"[info] Parse errors: {parse_err}", file=sys.stderr)
    if missing_out:
        print(f"[info] Skipped (outgroup missing/fail to root): {missing_out}", file=sys.stderr)
    if missing_trip:
        print(f"[info] Skipped (no complete triplet present): {missing_trip}", file=sys.stderr)

    # Write intermediate output
    intermediate_path = Path.cwd() / "ghostfinder_intermediate1.txt"
    with intermediate_path.open("w", encoding="utf-8") as outfh:
        outfh.write("Taxon_A\tTaxon_B\tTaxon_C\tAB_count\tBC_count\tAC_count\n")
        for a, b, c in triplets:
            ab, bc, ac = counts.get((a, b, c), (0, 0, 0))
            outfh.write(f"{a}\t{b}\t{c}\t{ab}\t{bc}\t{ac}\n")
    print(f"[info] Wrote counts for {len(triplets)} triplet(s) to {intermediate_path}", file=sys.stderr)

    # Create ghostfinder_triplets.csv with AC_count > BC_count
    csv_path = Path.cwd() / "ghostfinder_triplets.csv"
    try:
        n = write_triplets_csv_from_intermediate(intermediate_path, csv_path)
        print(f"[info] Wrote {n} triplet(s) to {csv_path} (AC_count > BC_count)", file=sys.stderr)
    except Exception as e:
        print(f"[error] Failed writing ghostfinder_triplets.csv: {e}", file=sys.stderr)
        return 5

    # 1) Remove ghostfinder_intermediate1.txt
    try:
        intermediate_path.unlink()
        print(f"[info] Removed {intermediate_path}", file=sys.stderr)
    except FileNotFoundError:
        pass
    except Exception as e:
        print(f"[warn] Could not remove {intermediate_path}: {e}", file=sys.stderr)

    # 3) Write ghostfinder_OUT.txt with the outgroup
    out_file = Path.cwd() / "ghostfinder_OUT.txt"
    try:
        out_file.write_text(f"{out_taxon}\n", encoding="utf-8")
        print(f"[info] Wrote outgroup to {out_file}", file=sys.stderr)
    except Exception as e:
        print(f"[error] Failed writing {out_file}: {e}", file=sys.stderr)
        return 6

    # 2 & 4) Build, print, and run ghostparser.py from the same directory as this script
    script_dir = Path(__file__).resolve().parent
    gp = script_dir / "ghostparser.py"
    if not gp.exists():
        print(f"[error] ghostparser.py not found at {gp}", file=sys.stderr)
        return 7

    # Build the command
    # Printed version must match the requested syntax: "python ghostparser.py ..."
    printed_cmd_parts = [
        "python",
        "ghostparser.py",
        "--out_taxa", "ghostfinder_OUT.txt",
        "--input_trees", str(args.gene_trees),
        "--species_tree", str(args.species_tree),
        "--output_file", str(args.output_file),
        "--triplets", "ghostfinder_triplets.csv",
    ]
    if args.threads is not None:
        printed_cmd_parts += ["--threads", str(args.threads)]
    printed_cmd = " ".join(shlex.quote(x) for x in printed_cmd_parts)
    print(f"[info] Running ghostparser:\n  {printed_cmd}", file=sys.stderr)

    # Executed version uses absolute paths and the current Python interpreter
    run_cmd = [
        sys.executable, str(gp),
        "--out_taxa", str(out_file),
        "--input_trees", str(args.gene_trees),
        "--species_tree", str(args.species_tree),
        "--output_file", str(args.output_file),
        "--triplets", str(csv_path),
    ]
    if args.threads is not None:
        run_cmd += ["--threads", str(args.threads)]

    try:
        result = subprocess.run(run_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[error] ghostparser.py exited with code {e.returncode}", file=sys.stderr)
        return e.returncode
    except Exception as e:
        print(f"[error] Failed to execute ghostparser.py: {e}", file=sys.stderr)
        return 8

    print("[info] ghostparser.py finished successfully.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

