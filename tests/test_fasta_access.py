#!/usr/bin/env python3
"""Integration test for LCG FASTA compression and random access.

Compresses two real FASTA files (scerevisiae8.fa.gz + B-3106.fa),
verifies metadata, then checks 50 random-range accesses byte-for-byte
against the originals.

Requirements: Python 3 stdlib only (gzip, subprocess, random).
"""

import gzip
import os
import random
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
LCG_BIN = os.path.join(REPO_ROOT, "build", "lcg")
DATA_DIR = os.path.join(SCRIPT_DIR, "data")

FASTA_GZ = os.path.join(DATA_DIR, "scerevisiae8.fa.gz")
FASTA_PLAIN = os.path.join(DATA_DIR, "B-3106.fa")

EXPECTED_NUM_SEQUENCES = 145
NUM_RANDOM_CHECKS = 50
RANDOM_SEED = 123


def parse_fasta(path):
    """Parse a FASTA file (plain or gzip) and return {name: sequence}."""
    seqs = {}
    current_name = None
    current_seq = []

    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_name is not None:
                    seqs[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name is not None:
            seqs[current_name] = "".join(current_seq)

    return seqs


def run_lcg(*args):
    """Run lcg with the given arguments, return CompletedProcess."""
    result = subprocess.run(
        [LCG_BIN] + list(args),
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"FAIL: lcg {' '.join(args)}")
        print(f"  stderr: {result.stderr.strip()}")
        sys.exit(1)
    return result


def main():
    # Verify binary exists
    if not os.path.isfile(LCG_BIN):
        print(f"FAIL: lcg binary not found at {LCG_BIN}")
        sys.exit(1)

    # Verify test data exists
    for path in (FASTA_GZ, FASTA_PLAIN):
        if not os.path.isfile(path):
            print(f"FAIL: test data not found at {path}")
            sys.exit(1)

    # Load all reference sequences
    print("Loading reference sequences...")
    seqs = {}
    seqs.update(parse_fasta(FASTA_GZ))
    seqs.update(parse_fasta(FASTA_PLAIN))
    print(f"  Loaded {len(seqs)} sequences from reference files")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create file list
        list_path = os.path.join(tmpdir, "fasta.list")
        with open(list_path, "w") as f:
            f.write(FASTA_GZ + "\n")
            f.write(FASTA_PLAIN + "\n")

        grammar_path = os.path.join(tmpdir, "test.lcg")

        # Compress
        print("Compressing with lcg cmp...")
        run_lcg("cmp", "-l", list_path, "-o", grammar_path, "-r")
        print("  Compression complete")

        # Check metadata
        print("Checking metadata...")
        result = run_lcg("met", grammar_path)
        met_output = result.stdout

        num_strings = None
        for line in met_output.splitlines():
            if "Number of compressed strings:" in line:
                # Format: "  Number of compressed strings:   145 (1.17 KB in pointers)"
                num_strings = int(line.split(":")[1].strip().split()[0])
                break

        if num_strings is None:
            print("FAIL: Could not parse number of sequences from metadata")
            print(f"  Metadata output:\n{met_output}")
            sys.exit(1)

        if num_strings != EXPECTED_NUM_SEQUENCES:
            print(f"FAIL: Expected {EXPECTED_NUM_SEQUENCES} sequences, got {num_strings}")
            sys.exit(1)

        print(f"  Metadata OK: {num_strings} sequences")

        # Random access checks
        print(f"Running {NUM_RANDOM_CHECKS} random access checks (seed={RANDOM_SEED})...")
        random.seed(RANDOM_SEED)
        seq_names = list(seqs.keys())
        errors = 0

        for i in range(NUM_RANDOM_CHECKS):
            name = random.choice(seq_names)
            seq_len = len(seqs[name])
            if seq_len < 100:
                continue

            start = random.randint(0, seq_len - 100)
            end = start + random.randint(10, min(500, seq_len - start - 1))

            # LCG end coordinate is inclusive
            expected = seqs[name][start:end + 1]

            result = subprocess.run(
                [LCG_BIN, "acc", grammar_path, f"{name}:{start}-{end}"],
                capture_output=True, text=True,
            )

            if result.returncode != 0:
                print(f"  [{i+1}] FAIL: lcg acc returned {result.returncode} for {name}:{start}-{end}")
                print(f"    stderr: {result.stderr.strip()}")
                errors += 1
                continue

            # Output is two lines: "index:start-end\nsequence_data"
            actual = result.stdout.strip().split("\n")[-1]

            if actual != expected:
                print(f"  [{i+1}] FAIL: {name}:{start}-{end}")
                print(f"    expected ({len(expected)} bp): {expected[:80]}...")
                print(f"    actual   ({len(actual)} bp): {actual[:80]}...")
                errors += 1
            else:
                print(f"  [{i+1}] OK: {name}:{start}-{end} ({len(expected)} bp)")

        print()
        if errors > 0:
            print(f"FAIL: {errors}/{NUM_RANDOM_CHECKS} random access checks failed")
            sys.exit(1)
        else:
            print(f"PASS: All {NUM_RANDOM_CHECKS} random access checks passed")
            sys.exit(0)


if __name__ == "__main__":
    main()
