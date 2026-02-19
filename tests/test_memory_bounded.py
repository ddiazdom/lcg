#!/usr/bin/env python3
"""Integration test: prove LCG peak RSS stays bounded as input scales.

Creates N copies of scerevisiae8.fa.gz with renamed sequences, compresses
them with lcg cmp, and asserts that peak RSS for 4x input is less than
2x the peak RSS for 1x input.

Requirements: Python 3 stdlib only, /usr/bin/time -v (Linux).
"""

import gzip
import os
import re
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
LCG_BIN = os.path.join(REPO_ROOT, "build", "lcg")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")

FASTA_GZ = os.path.join(DATA_DIR, "scerevisiae8.fa.gz")


def generate_renamed_fasta(src_gz, dest_fa, prefix):
    """Read a gzipped FASTA and write plain FASTA with prefixed sequence names."""
    with gzip.open(src_gz, "rt") as fin, open(dest_fa, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                name = line[1:].rstrip("\n")
                fout.write(f">{prefix}_{name}\n")
            else:
                fout.write(line)


def run_lcg_cmp_with_memory(lcg_bin, file_list_path, output_path):
    """Run lcg cmp under /usr/bin/time -v, return peak RSS in kB."""
    cmd = [
        "/usr/bin/time", "-v",
        lcg_bin, "cmp",
        "-l", file_list_path,
        "-o", output_path,
        "-r",
        "-t", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"FAIL: lcg cmp returned {result.returncode}")
        print(f"  cmd: {' '.join(cmd)}")
        print(f"  stdout: {result.stdout.strip()}")
        print(f"  stderr: {result.stderr.strip()}")
        sys.exit(1)

    # Parse "Maximum resident set size (kbytes): NNNN" from /usr/bin/time stderr
    match = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", result.stderr)
    if not match:
        print("FAIL: Could not parse peak RSS from /usr/bin/time output")
        print(f"  stderr: {result.stderr.strip()}")
        sys.exit(1)

    return int(match.group(1))


def run_test(tmpdir, n_copies):
    """Generate n_copies of renamed FASTA, compress, return peak RSS in kB."""
    fasta_paths = []
    for i in range(n_copies):
        dest = os.path.join(tmpdir, f"copy{i}.fa")
        print(f"  Generating copy {i} -> {dest}")
        generate_renamed_fasta(FASTA_GZ, dest, f"copy{i}")
        fasta_paths.append(dest)

    list_path = os.path.join(tmpdir, f"fasta_{n_copies}.list")
    with open(list_path, "w") as f:
        for p in fasta_paths:
            f.write(p + "\n")

    output_path = os.path.join(tmpdir, f"out_{n_copies}.lcg")
    print(f"  Compressing {n_copies} copies with lcg cmp -t 1 ...")
    rss_kb = run_lcg_cmp_with_memory(LCG_BIN, list_path, output_path)
    print(f"  Peak RSS: {rss_kb} kB ({rss_kb / 1024:.1f} MB)")
    return rss_kb


def main():
    # Verify binary exists
    if not os.path.isfile(LCG_BIN):
        print(f"FAIL: lcg binary not found at {LCG_BIN}")
        sys.exit(1)

    # Verify /usr/bin/time exists
    if not os.path.isfile("/usr/bin/time"):
        print("FAIL: /usr/bin/time not found (needed for peak RSS measurement)")
        sys.exit(1)

    # Verify test data exists
    if not os.path.isfile(FASTA_GZ):
        print(f"FAIL: test data not found at {FASTA_GZ}")
        sys.exit(1)

    with tempfile.TemporaryDirectory() as tmpdir:
        print("=== N=1 (baseline) ===")
        rss_1 = run_test(tmpdir, 1)

        print()
        print("=== N=4 (4x input) ===")
        rss_4 = run_test(tmpdir, 4)

        print()
        print("=== Results ===")
        print(f"  N=1 peak RSS: {rss_1} kB ({rss_1 / 1024:.1f} MB)")
        print(f"  N=4 peak RSS: {rss_4} kB ({rss_4 / 1024:.1f} MB)")
        ratio = rss_4 / rss_1 if rss_1 > 0 else float("inf")
        print(f"  Ratio (N=4 / N=1): {ratio:.2f}x")
        print()

        if rss_4 < 2 * rss_1:
            print(f"PASS: Peak RSS ratio {ratio:.2f}x < 2.0x — memory is bounded")
            sys.exit(0)
        else:
            print(f"FAIL: Peak RSS ratio {ratio:.2f}x >= 2.0x — memory grew too much")
            sys.exit(1)


if __name__ == "__main__":
    main()
