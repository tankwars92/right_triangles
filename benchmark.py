#!/usr/bin/env python3

import argparse
import csv
import json
import math
import os
import platform
import random
import statistics
import sys
import time
from math import isqrt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from rt_algorithms import OEIS, brute, quad_formula, R_fast
from nearlin import R_nearlinear


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check correctness and benchmark right-triangle counting algorithms."
    )
    parser.add_argument("--bfile", help="optional OEIS b-file for additional checks")
    parser.add_argument("--full", action="store_true", help="include n = 25600 (longer run)")
    parser.add_argument("--quick", action="store_true", help="run a smaller test set")
    parser.add_argument("--out", default=".", help="directory for result files (default: current directory)")
    return parser.parse_args()


args = parse_args()
OUT = args.out
os.makedirs(OUT, exist_ok=True)
JSON_PATH = os.path.join(OUT, "results.json")

results = {
    "environment": {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "processor": platform.processor() or "not reported",
        "logical_cpus": os.cpu_count(),
        "started_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    },
    "correctness": {},
    "timings": [],
    "analytic_pairs": [],
    "fits": {},
}


def save_progress():
    with open(JSON_PATH, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)


def heading(title):
    print()
    print(title)
    print("-" * len(title))


def status(ok, message):
    print(f"[{'PASS' if ok else 'FAIL'}] {message}", flush=True)


def format_time(value):
    if value is None:
        return "-"
    if value < 0.001:
        return f"{value * 1000:.3f} ms"
    return f"{value:.4f} s"


def timed(func, n, repeats):
    samples = []
    value = operations = None
    for _ in range(repeats):
        start = time.perf_counter()
        value, operations = func(n)
        samples.append(time.perf_counter() - start)
    return value, operations, min(samples), statistics.median(samples)


def call_base(n):
    value, operations = quad_formula(n)
    return value, operations


def call_fast(n):
    counter = []
    value = R_fast(n, counter)
    return value, counter[0]


def call_nearlinear(n):
    counter = [0]
    value = R_nearlinear(n, counter)
    return value, counter[0]


def pair_count(n):
    total = 0
    d = 1
    while d <= n - 1:
        quotient = (n - 1) // d
        last = (n - 1) // quotient
        root = isqrt(quotient)
        total += root * (root + 1) // 2
        d = last + 1
    return total


def power_fit(xs, ys):
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    mean_x = sum(lx) / len(lx)
    mean_y = sum(ly) / len(ly)
    denominator = sum((x - mean_x) ** 2 for x in lx)
    exponent = sum((x - mean_x) * (y - mean_y) for x, y in zip(lx, ly)) / denominator
    coefficient = math.exp(mean_y - exponent * mean_x)
    return exponent, coefficient


def check_bfile(path):
    values = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            fields = line.split()
            if len(fields) >= 2 and fields[0].isdigit() and fields[1].lstrip("-").isdigit():
                values[int(fields[0])] = int(fields[1])

    rng = random.Random(2026)
    limit = 2000 if args.quick else 10000
    fixed_fast = [n for n in (100, 200, 300, 500, 800, 1000, 2000, 3000, 5000, 8000, 10000)
                  if n <= limit and n in values]
    pool_fast = [n for n in range(36, min(1500, max(values, default=35)) + 1) if n in values]
    random_fast = rng.sample(pool_fast, min(40 if args.quick else 150, len(pool_fast)))
    tests_fast = sorted(set(fixed_fast + random_fast))
    bad_fast = [n for n in tests_fast if R_fast(n) != values[n]]
    results["correctness"]["bfile_n^1.5"] = {
        "checked": len(tests_fast), "mismatches": bad_fast, "source": path
    }
    status(not bad_fast, f"R_fast vs b-file: {len(tests_fast)} values checked; {len(bad_fast)} mismatches")

    fixed_nl = [n for n in (100, 200, 300, 500, 800, 1000, 1500, 2000, 3000)
                if n in values and (not args.quick or n <= 1000)]
    pool_nl = [n for n in range(36, min(1500, max(values, default=35)) + 1) if n in values]
    random_nl = rng.sample(pool_nl, min(15 if args.quick else 60, len(pool_nl)))
    tests_nl = sorted(set(fixed_nl + random_nl))
    bad_nl = [n for n in tests_nl if R_nearlinear(n) != values[n]]
    results["correctness"]["bfile_nearlinear"] = {
        "checked": len(tests_nl), "mismatches": bad_nl, "source": path
    }
    status(not bad_nl, f"R_nearlinear vs b-file: {len(tests_nl)} values checked; {len(bad_nl)} mismatches")


def run_correctness_checks():
    heading("1. Correctness checks")
    basic_ok = all(
        brute(n) == quad_formula(n)[0] == R_fast(n) == R_nearlinear(n) == OEIS[n - 1]
        for n in range(1, 8)
    )
    results["correctness"]["brute_n_1_to_7"] = "PASS" if basic_ok else "FAIL"
    status(basic_ok, "Brute force, baseline, R_fast, R_nearlinear and OEIS agree for n = 1..7")
    if not basic_ok:
        raise RuntimeError("The basic correctness check failed.")

    oeis_ok = all(R_fast(n) == R_nearlinear(n) == OEIS[n - 1] for n in range(1, 36))
    results["correctness"]["oeis_first_35"] = "PASS" if oeis_ok else "FAIL"
    status(oeis_ok, "R_fast and R_nearlinear match the first 35 OEIS values")
    if not oeis_ok:
        raise RuntimeError("The OEIS check failed.")

    upper = 60 if args.quick else 150
    agree = all(R_fast(n) == R_nearlinear(n) == quad_formula(n)[0] for n in range(1, upper + 1))
    results["correctness"][f"three_methods_through_{upper}"] = "PASS" if agree else "FAIL"
    status(agree, f"All three methods agree for every n <= {upper}")
    if not agree:
        raise RuntimeError("The three-method comparison failed.")

    if args.bfile:
        check_bfile(args.bfile)
    else:
        print("[SKIP] No b-file supplied; external OEIS checks above n = 35 were skipped.")
    save_progress()


def run_benchmarks():
    small = [50, 100, 200, 400] if args.quick else [50, 100, 200, 400, 800]
    large = [1600] if args.quick else [1600, 3200, 6400, 12800] + ([25600] if args.full else [])
    test_sizes = small + large

    heading("2. Runtime and operation counts")
    print("Times are the fastest sample; small cases use 3 runs, larger cases use 1.")
    header = (
        f"{'n':>7}  {'Baseline':>12}  {'Baseline ops':>14}  "
        f"{'R_fast':>12}  {'Fast nodes':>12}  {'R_nearlinear':>14}  {'Pair calls':>12}  {'Check':>7}"
    )
    print(header)
    print("-" * len(header))

    for n in test_sizes:
        row = {"n": n}
        repeats = 3 if n <= 800 else 1
        values = []

        if n in small:
            value, ops, best, median = timed(call_base, n, repeats)
            row.update(baseline_s=best, baseline_median_s=median, baseline_ops=ops)
            values.append(value)

        value, ops, best, median = timed(call_fast, n, repeats)
        row.update(fast_s=best, fast_median_s=median, fast_ops=ops)
        values.append(value)

        value, ops, best, median = timed(call_nearlinear, n, repeats)
        row.update(nearlinear_s=best, nearlinear_median_s=median, nearlinear_pairs=ops)
        values.append(value)

        row["agree"] = len(set(values)) == 1
        row["value"] = str(values[0])
        results["timings"].append(row)
        save_progress()

        base_time = format_time(row.get("baseline_s"))
        base_ops = f"{row['baseline_ops']:,}" if "baseline_ops" in row else "-"
        check = "PASS" if row["agree"] else "FAIL"
        print(
            f"{n:>7,}  {base_time:>12}  {base_ops:>14}  "
            f"{format_time(row['fast_s']):>12}  {row['fast_ops']:>12,}  "
            f"{format_time(row['nearlinear_s']):>14}  {row['nearlinear_pairs']:>12,}  {check:>7}"
        )
        if not row["agree"]:
            raise RuntimeError(f"Algorithms returned different results for n = {n}.")


def run_analysis():
    heading("3. Analytic pair counts")
    for exponent in range(2, 10):
        n = 10 ** exponent
        count = pair_count(n)
        ratio = count / (n * math.log(n))
        results["analytic_pairs"].append({"n": n, "pairs": count, "pairs_over_nlogn": ratio})
        print(f"n = 10^{exponent}: {count:>15,} pairs   pairs / (n ln n) = {ratio:.4f}")

    heading("4. Empirical growth rates")
    timings = results["timings"]

    def get_column(key, minimum_n=0):
        rows = [r for r in timings if key in r and r["n"] >= minimum_n]
        return [r["n"] for r in rows], [r[key] for r in rows]

    fits = {}
    fit_specs = [
        ("Baseline runtime", "baseline_s", 100),
        ("R_fast runtime", "fast_s", 1600),
        ("R_nearlinear runtime", "nearlinear_s", 1600),
        ("Baseline operation count", "baseline_ops", 100),
        ("R_fast node count", "fast_ops", 1600),
        ("R_nearlinear pair count", "nearlinear_pairs", 1600),
    ]
    for name, key, minimum_n in fit_specs:
        xs, ys = get_column(key, minimum_n)
        if len(xs) >= 2 and all(y > 0 for y in ys):
            exponent, coefficient = power_fit(xs, ys)
            fits[name] = {"exponent": exponent, "coefficient": coefficient, "n_range": [xs[0], xs[-1]]}

    analytic_xs = [10 ** e for e in range(5, 10)]
    analytic_ys = [pair_count(n) for n in analytic_xs]
    fits["Analytic pair count"] = {
        "exponent": power_fit(analytic_xs, analytic_ys)[0],
        "n_range": [analytic_xs[0], analytic_xs[-1]],
    }
    results["fits"] = fits

    for name, item in fits.items():
        n_range = item.get("n_range")
        range_text = f" (n = {n_range[0]:,}..{n_range[1]:,})" if n_range else ""
        print(f"{name:<29} exponent = {item['exponent']:.3f}{range_text}")

    crossover = None
    if "R_fast runtime" in fits and "R_nearlinear runtime" in fits:
        fast_fit = fits["R_fast runtime"]
        near_fit = fits["R_nearlinear runtime"]
        if fast_fit["exponent"] > near_fit["exponent"]:
            crossover = math.exp(
                math.log(near_fit["coefficient"] / fast_fit["coefficient"])
                / (fast_fit["exponent"] - near_fit["exponent"])
            )
            results["estimated_crossover_n"] = crossover
            print(f"\nEstimated crossover from fitted curves: n ~ {crossover:,.0f}")

    measured = [r["n"] for r in timings if r["nearlinear_s"] < r["fast_s"]]
    results["first_measured_crossover_n"] = measured[0] if measured else None
    if measured:
        print(f"First measured case where R_nearlinear is faster: n = {measured[0]:,}")
    else:
        print("R_fast was faster throughout the tested range.")
    save_progress()



def main():
    run_correctness_checks()
    run_benchmarks()
    run_analysis()
    print("\nBenchmark completed.")


if __name__ == "__main__":
    main()
