"""Benchmark custom vs SciPy-backed statistical helpers used by GhostParser.

This script measures aggregate performance for:
- Pearson chi-square discordant count test
- Two-sample KS test

It avoids per-call logging and instead reports summary metrics over many runs.
"""

import argparse
import sys
import random
import statistics
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ghostparser import triplet_processor as tp


def _quantile(values: list[float], q: float) -> float:
    if not values:
        return 0.0
    if q <= 0:
        return min(values)
    if q >= 1:
        return max(values)
    sorted_values = sorted(values)
    index = q * (len(sorted_values) - 1)
    lower = int(index)
    upper = min(lower + 1, len(sorted_values) - 1)
    weight = index - lower
    return sorted_values[lower] * (1 - weight) + sorted_values[upper] * weight


def _time_batch(func, batch):
    start = time.perf_counter()
    for args in batch:
        func(*args)
    end = time.perf_counter()
    elapsed = end - start
    return elapsed, elapsed / len(batch)


def _summarize(label: str, values: list[float]) -> str:
    mean_val = statistics.mean(values)
    median_val = statistics.median(values)
    stdev_val = statistics.stdev(values) if len(values) > 1 else 0.0
    p05 = _quantile(values, 0.05)
    p95 = _quantile(values, 0.95)
    return (
        f"{label}: mean={mean_val:.9e}s median={median_val:.9e}s "
        f"std={stdev_val:.9e}s p05={p05:.9e}s p95={p95:.9e}s"
    )


def _generate_chi_square_batch(rng: random.Random, batch_size: int, max_count: int):
    batch = []
    for _ in range(batch_size):
        n_dis1 = rng.randint(0, max_count)
        n_dis2 = rng.randint(0, max_count)
        batch.append((n_dis1, n_dis2))
    return batch


def _generate_ks_batch(
    rng: random.Random,
    batch_size: int,
    min_sample_size: int,
    max_sample_size: int,
):
    batch = []
    for _ in range(batch_size):
        n1 = rng.randint(min_sample_size, max_sample_size)
        n2 = rng.randint(min_sample_size, max_sample_size)
        sample_a = [rng.random() for _ in range(n1)]
        sample_b = [rng.random() for _ in range(n2)]
        batch.append((sample_a, sample_b))
    return batch


def _run_comparison(custom_func, scipy_func, batch, repeats: int):
    custom_per_call = []
    scipy_per_call = []
    for _ in range(repeats):
        _, custom_avg = _time_batch(custom_func, batch)
        _, scipy_avg = _time_batch(scipy_func, batch)
        custom_per_call.append(custom_avg)
        scipy_per_call.append(scipy_avg)
    return custom_per_call, scipy_per_call


def main():
    parser = argparse.ArgumentParser(description="Profile custom vs SciPy GhostParser statistics helpers.")
    parser.add_argument("--seed", type=int, default=13, help="Random seed for reproducible synthetic data")
    parser.add_argument("--repeats", type=int, default=12, help="Number of repeated timing rounds")
    parser.add_argument(
        "--chi-batch-size",
        type=int,
        default=50000,
        help="Number of chi-square calls per timing round",
    )
    parser.add_argument(
        "--ks-batch-size",
        type=int,
        default=4000,
        help="Number of KS calls per timing round",
    )
    parser.add_argument(
        "--max-discordant-count",
        type=int,
        default=200,
        help="Max synthetic discordant count for chi-square inputs",
    )
    parser.add_argument(
        "--min-ks-sample-size",
        type=int,
        default=20,
        help="Minimum synthetic sample size per KS input",
    )
    parser.add_argument(
        "--max-ks-sample-size",
        type=int,
        default=80,
        help="Maximum synthetic sample size per KS input",
    )
    args = parser.parse_args()

    if args.repeats < 1:
        raise ValueError("--repeats must be >= 1")
    if args.chi_batch_size < 1 or args.ks_batch_size < 1:
        raise ValueError("--chi-batch-size and --ks-batch-size must be >= 1")
    if args.min_ks_sample_size < 1 or args.max_ks_sample_size < args.min_ks_sample_size:
        raise ValueError("Invalid KS sample size range")

    rng = random.Random(args.seed)

    chi_batch = _generate_chi_square_batch(rng, args.chi_batch_size, args.max_discordant_count)
    ks_batch = _generate_ks_batch(
        rng,
        args.ks_batch_size,
        args.min_ks_sample_size,
        args.max_ks_sample_size,
    )

    chi_custom, chi_scipy = _run_comparison(
        tp.pearson_discordant_chi_square_test,
        tp._pearson_discordant_chi_square_test_scipy,
        chi_batch,
        repeats=args.repeats,
    )

    ks_custom, ks_scipy = _run_comparison(
        tp.two_sample_ks_test,
        tp._two_sample_ks_test_scipy,
        ks_batch,
        repeats=args.repeats,
    )

    print("=== GhostParser Statistical Helper Profiling ===")
    print(f"seed={args.seed} repeats={args.repeats}")
    print(f"chi_batch_size={args.chi_batch_size} ks_batch_size={args.ks_batch_size}")
    print(
        "ks_sample_size_range="
        f"[{args.min_ks_sample_size}, {args.max_ks_sample_size}]"
    )
    print()

    print("[Chi-square per call]")
    print(_summarize("custom", chi_custom))
    print(_summarize("scipy_backup", chi_scipy))
    print(f"speedup(custom/scipy): {statistics.mean(chi_scipy) / statistics.mean(chi_custom):.3f}x")
    print()

    print("[KS per call]")
    print(_summarize("custom", ks_custom))
    print(_summarize("scipy_backup", ks_scipy))
    print(f"speedup(custom/scipy): {statistics.mean(ks_scipy) / statistics.mean(ks_custom):.3f}x")


if __name__ == "__main__":
    main()
