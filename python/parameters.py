#!/usr/bin/env python3
import os
import subprocess
import itertools
import csv
import re
import glob
from pathlib import Path


def run_mls():
    """
    Run ./mls for all parameter combinations.
    """
    dir_arg = "../afe_data_param_tuning/"
    position_tolerances = [2]#[0, 1, 3, 5, 7, 9]
    coverage_tolerances = [0.1] #[0.05, 0.1, 0.15, 0.2]
    min_lengths = [5]#[3, 5, 7, 10, 15]
    min_coverages =[0.0001] #[0.25, 0.01, 0.005, 0.0001, 0.0005, 0.00001, 0.000001]

    for pt, ct, ml, mc in itertools.product(
            position_tolerances,
            coverage_tolerances,
            min_lengths,
            min_coverages,
    ):
        cmd = [
            "./mls",
            "--dir", dir_arg,
            "--position-tolerance", str(pt),
            "--min-coverage", str(mc),
            "--min-length", str(ml),
            "--coverage-tolerance", str(ct),
        ]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)


def parse_stats_file(stats_path: Path) -> dict:
    """
    Parse COMP_*.stats file from gffcompare and return a dict of metrics.
    """
    metrics = {}

    # Default values (optional but can be handy)
    levels = ["Base", "Exon", "Intron", "Intron chain", "Transcript", "Locus"]
    for level in levels:
        key_base = level.replace(" ", "_")
        metrics[f"{key_base}_Sensitivity"] = ""
        metrics[f"{key_base}_Precision"] = ""

    for key in [
        "Matching_intron_chains",
        "Matching_transcripts",
        "Matching_loci",
    ]:
        metrics[key] = ""

    missed_novel_keys = [
        "Missed_exons",
        "Novel_exons",
        "Missed_introns",
        "Novel_introns",
        "Missed_loci",
        "Novel_loci",
    ]
    for key in missed_novel_keys:
        metrics[f"{key}_n"] = ""
        metrics[f"{key}_total"] = ""
        metrics[f"{key}_pct"] = ""

    with stats_path.open() as fh:
        for line in fh:
            # Sensitivity/Precision table
            m = re.match(r"^\s*(.+?) level:\s+([\d.]+)\s*\|\s*([\d.]+)", line)
            if m:
                level = m.group(1).strip()
                sens = m.group(2)
                prec = m.group(3)
                key_base = level.replace(" ", "_")
                metrics[f"{key_base}_Sensitivity"] = sens
                metrics[f"{key_base}_Precision"] = prec
                continue

            # Matching counts
            m = re.match(r"^\s*Matching (intron chains|transcripts|loci):\s+(\d+)", line)
            if m:
                what = m.group(1).replace(" ", "_")  # intron_chains / transcripts / loci
                value = m.group(2)
                metrics[f"Matching_{what}"] = value
                continue

            # Missed / Novel:
            m = re.match(
                r"^\s*(Missed exons|Novel exons|Missed introns|Novel introns|Missed loci|Novel loci):\s+(\d+)/(\d+)\s+\(\s*([\d.]+)%\)",
                line,
            )
            if m:
                key = m.group(1).replace(" ", "_")   # e.g. Missed_exons
                n = m.group(2)
                total = m.group(3)
                pct = m.group(4)
                metrics[f"{key}_n"] = n
                metrics[f"{key}_total"] = total
                metrics[f"{key}_pct"] = pct
                continue

    return metrics


def run_gffcompare_and_collect(outputs_dir: Path, csv_name: str = "gffcompare_summary.csv"):
    """
    In outputs_dir, run gffcompare for each result file and parse COMP_*.stats
    into a CSV.
    """
    # Change into outputs dir for gffcompare command
    os.chdir(outputs_dir)

    # Collect result files (e.g. *.gtf) but exclude the reference and any previous COMP_*.gtf
    all_gtf = glob.glob("*.gtf")
    result_files = [
        f for f in all_gtf
        if f != "splicing_variants.gtf" and not f.startswith("COMP_")
    ]

    rows = []
    index = 0
    for result_file in sorted(result_files):
        stem = Path(result_file).stem
        prefix = f"COMP_{index}"
        index += 1

        cmd = [
            "gffcompare",
            "-r", "splicing_variants.gtf",
            "-o", prefix,
            result_file,
            "splicing_variants.gtf",
        ]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)

        stats_file = outputs_dir / f"{prefix}.stats"
        if not stats_file.exists():
            print(f"WARNING: stats file {stats_file} not found, skipping.")
            continue

        metrics = parse_stats_file(stats_file)
        metrics["file_name"] = result_file
        rows.append(metrics)

    if not rows:
        print("No metrics collected; nothing to write to CSV.")
        return

    # Collect all keys for CSV header (union over all rows)
    fieldnames = sorted(set().union(*[row.keys() for row in rows]))
    # Ensure file_name is first
    fieldnames = ["file_name"] + [f for f in fieldnames if f != "file_name"]

    csv_path = outputs_dir / csv_name
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Written CSV summary to {csv_path}")


def main():
    # Step 1: run ./mls for all parameter combinations
    #run_mls()

    # Step 2: change into dir/outputs where dir is ../afe_data/
    outputs_dir = Path("../afe_data_param_tuning")
    outputs_dir.mkdir(parents=True, exist_ok=True)

    # Step 3: run gffcompare on each result file and collect metrics
    run_gffcompare_and_collect(outputs_dir)


if __name__ == "__main__":
    main()
