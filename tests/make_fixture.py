#!/usr/bin/env python3
"""Generate a small read set with known truth for end-to-end alignment tests.

Every read is cut from the reference at a known position, optionally with
SNPs or a single 1bp indel injected, so the aligner's POS and CIGAR can be
checked exactly. The RNG is seeded, so the output is reproducible.

usage: make_fixture.py <reference.fa> <out_reads.fa> <out_truth.tsv>
"""
import random
import sys

READ_LEN = 100
SEED_LEN = READ_LEN // 5  # matches seedSkip() in utils/align.c
BASES = "ACGT"


def load_reference(path):
    with open(path) as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">")).upper()


def snp(read, offset, rng):
    alt = rng.choice([b for b in BASES if b != read[offset]])
    return read[:offset] + alt + read[offset + 1:]


def main(ref_path, reads_path, truth_path):
    ref = load_reference(ref_path)
    rng = random.Random(2020)
    cases = []  # (name, read, 0-based ref start, snps, expected indel op or "")

    def start():
        # stay clear of the genome edges; edge reads get their own cases
        return rng.randrange(READ_LEN, len(ref) - 2 * READ_LEN)

    for i in range(10):
        p = start()
        cases.append((f"exact_{i}", ref[p:p + READ_LEN], p, 0, ""))

    # a variant in the first seed used to drop the whole read
    for i in range(10):
        p = start()
        cases.append((f"snp_first_seed_{i}", snp(ref[p:p + READ_LEN], rng.randrange(SEED_LEN), rng), p, 1, ""))

    # variants in seeds 1 and 3: only seeds 2, 4 and 5 match exactly
    for i in range(10):
        p = start()
        read = ref[p:p + READ_LEN]
        read = snp(read, rng.randrange(0, SEED_LEN), rng)
        read = snp(read, rng.randrange(2 * SEED_LEN, 3 * SEED_LEN), rng)
        cases.append((f"snp_two_seeds_{i}", read, p, 2, ""))

    # 1bp insertion / deletion in the middle of the read
    for i in range(10):
        p = start()
        k = rng.randrange(40, 60)
        ins = rng.choice(BASES)
        read = ref[p:p + k] + ins + ref[p + k:p + READ_LEN - 1]
        cases.append((f"insertion_{i}", read, p, 0, "I"))
    for i in range(10):
        p = start()
        k = rng.randrange(40, 60)
        read = ref[p:p + k] + ref[p + k + 1:p + READ_LEN + 1]
        cases.append((f"deletion_{i}", read, p, 0, "D"))

    # reads touching both ends of the genome
    cases.append(("genome_start", ref[:READ_LEN], 0, 0, ""))
    cases.append(("genome_end", ref[-READ_LEN:], len(ref) - READ_LEN, 0, ""))

    with open(reads_path, "w", newline="\n") as reads, open(truth_path, "w", newline="\n") as truth:
        truth.write("name\tpos\tsnps\tindel\n")
        for name, read, p, snps, indel in cases:
            assert len(read) == READ_LEN, name
            reads.write(f">{name}\n{read}\n")
            truth.write(f"{name}\t{p + 1}\t{snps}\t{indel}\n")  # SAM POS is 1-based


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
