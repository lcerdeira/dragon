#!/usr/bin/env python3
"""Parse /usr/bin/time -v output to extract resource metrics."""

import argparse
import os
import re


def parse_time_output(path):
    """Parse GNU time -v output file."""
    metrics = {
        "wall_time_seconds": 0.0,
        "cpu_time_seconds": 0.0,
        "peak_ram_kb": 0,
        "exit_status": -1,
    }

    if not os.path.exists(path):
        return metrics

    with open(path) as f:
        for line in f:
            line = line.strip()

            # Wall clock time
            m = re.search(r"Elapsed \(wall clock\) time.*?(\d+):(\d+\.\d+)", line)
            if m:
                minutes = int(m.group(1))
                seconds = float(m.group(2))
                metrics["wall_time_seconds"] = minutes * 60 + seconds

            # User time
            m = re.search(r"User time \(seconds\):\s*(\d+\.\d+)", line)
            if m:
                metrics["cpu_time_seconds"] += float(m.group(1))

            # System time
            m = re.search(r"System time \(seconds\):\s*(\d+\.\d+)", line)
            if m:
                metrics["cpu_time_seconds"] += float(m.group(1))

            # Maximum resident set size
            m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", line)
            if m:
                metrics["peak_ram_kb"] = int(m.group(1))

            # Exit status
            m = re.search(r"Exit status:\s*(\d+)", line)
            if m:
                metrics["exit_status"] = int(m.group(1))

    return metrics


def get_dir_size(path):
    """Get total size of a directory in bytes."""
    total = 0
    if os.path.isdir(path):
        for dirpath, dirnames, filenames in os.walk(path):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                total += os.path.getsize(fp)
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resource-files", nargs="+", required=True)
    parser.add_argument("--index-resources", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tool", required=True)
    parser.add_argument("--dataset", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Parse index construction resources
    index_metrics = parse_time_output(args.index_resources)

    # Parse search resources for each divergence level
    search_metrics = []
    for rf in args.resource_files:
        m = parse_time_output(rf)
        div = "0.0"
        for part in os.path.basename(rf).split("_"):
            if part.startswith("div"):
                div = part.replace("div", "").replace("_resources.txt", "")
        m["divergence"] = div
        search_metrics.append(m)

    with open(args.output, "w") as out:
        out.write("tool\tdataset\tmetric_type\tdivergence\t"
                  "wall_time_s\tcpu_time_s\tpeak_ram_gb\tindex_size_gb\n")

        # Index construction row
        out.write(
            f"{args.tool}\t{args.dataset}\tindex\t-\t"
            f"{index_metrics['wall_time_seconds']:.2f}\t"
            f"{index_metrics['cpu_time_seconds']:.2f}\t"
            f"{index_metrics['peak_ram_kb'] / 1048576:.4f}\t"
            f"-\n"
        )

        # Search rows
        for m in search_metrics:
            out.write(
                f"{args.tool}\t{args.dataset}\tsearch\t{m['divergence']}\t"
                f"{m['wall_time_seconds']:.2f}\t"
                f"{m['cpu_time_seconds']:.2f}\t"
                f"{m['peak_ram_kb'] / 1048576:.4f}\t"
                f"-\n"
            )

    print(f"Resource metrics written to {args.output}")


if __name__ == "__main__":
    main()
