#!/usr/bin/env python3
"""Validate aligner SAM output against the reference (and optionally a truth table).

Checks on every mapped record:
  * CIGAR only uses M/I/D and consumes exactly len(SEQ) query bases
  * replaying the CIGAR at POS stays inside the reference
  * replayed mismatches are few (a wrong POS or reversed CIGAR shows up as many)
  * TLEN is 0 (single-end reads)
Checks across the file:
  * one primary record per read, secondary records flagged 0x100, no duplicate hits
  * mapped fraction >= --min-mapped
With --truth (from make_fixture.py):
  * every read is mapped at exactly its true POS
  * replayed mismatches == injected SNPs, and the expected indel op is present

usage: validate_sam.py --ref REF.fa --sam OUT.sam [--truth TRUTH.tsv] [--min-mapped F]
"""
import argparse
import collections
import re
import sys

FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
MAX_MISMATCH_FRACTION = 0.10  # simulated reads should be far below this


def load_reference(path):
    with open(path) as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">")).upper()


def parse_cigar(cigar):
    ops = [(int(n), op) for n, op in re.findall(r"(\d+)([^\d])", cigar)]
    if "".join(f"{n}{op}" for n, op in ops) != cigar:
        raise ValueError(f"malformed CIGAR {cigar!r}")
    return ops


def replay(ref, seq, pos0, ops):
    """Return mismatches when SEQ is laid onto REF at 0-based pos0 using ops."""
    r, q, mismatches = pos0, 0, 0
    for n, op in ops:
        if op == "M":
            if r + n > len(ref):
                raise ValueError("alignment runs past the end of the reference")
            mismatches += sum(1 for k in range(n) if ref[r + k] != seq[q + k] and seq[q + k] != "N")
            r += n
            q += n
        elif op == "I":
            q += n
        elif op == "D":
            r += n
        else:
            raise ValueError(f"unexpected CIGAR op {op!r}")
    return mismatches


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--sam", required=True)
    ap.add_argument("--truth")
    ap.add_argument("--min-mapped", type=float, default=0.0)
    args = ap.parse_args()

    ref = load_reference(args.ref)
    errors = []
    records = collections.defaultdict(list)  # qname -> [(flag, pos, ops, mismatches)] for valid records
    mapped_names = set()  # reads with at least one mapped record, valid or not

    with open(args.sam) as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                errors.append(f"line {lineno}: expected 11 fields, got {len(f)}")
                continue
            qname, flag, pos, cigar, tlen, seq = f[0], int(f[1]), int(f[3]), f[5], int(f[8]), f[9]
            where = f"line {lineno} ({qname})"
            if flag & FLAG_UNMAPPED:
                records[qname].append((flag, 0, None, None))
                continue
            mapped_names.add(qname)
            records[qname]  # count the read even if its record turns out invalid
            try:
                ops = parse_cigar(cigar)
                query_len = sum(n for n, op in ops if op in "MI")
                if query_len != len(seq):
                    raise ValueError(f"CIGAR {cigar} consumes {query_len} bases but SEQ has {len(seq)}")
                if pos < 1:
                    raise ValueError(f"POS {pos} is not 1-based")
                mismatches = replay(ref, seq, pos - 1, ops)
                if mismatches > MAX_MISMATCH_FRACTION * len(seq):
                    raise ValueError(f"{mismatches} mismatches replaying {cigar} at POS {pos}")
                if tlen != 0:
                    raise ValueError(f"TLEN {tlen} should be 0 for single-end reads")
            except ValueError as e:
                errors.append(f"{where}: {e}")
                continue
            records[qname].append((flag, pos, ops, mismatches))

    for qname, recs in records.items():
        primaries = [r for r in recs if not r[0] & FLAG_SECONDARY]
        if recs and len(primaries) != 1:
            errors.append(f"{qname}: {len(primaries)} primary records (expected 1)")
        hits = [(r[1], tuple(r[2])) for r in recs if r[2] is not None]
        if len(hits) != len(set(hits)):
            errors.append(f"{qname}: duplicate alignments")

    total = len(records)
    mapped = len(mapped_names)
    fraction = mapped / total if total else 0.0

    if args.truth:
        with open(args.truth) as fh:
            truth = [line.rstrip("\n").split("\t") for line in fh][1:]
        for name, pos, snps, indel in truth:
            recs = [r for r in records.get(name, []) if r[2] is not None]
            if not recs:
                errors.append(f"{name}: not mapped (truth POS {pos})")
                continue
            primary = next((r for r in recs if not r[0] & FLAG_SECONDARY), recs[0])
            _, got_pos, ops, mismatches = primary
            if got_pos != int(pos):
                errors.append(f"{name}: POS {got_pos}, expected {pos}")
            if mismatches != int(snps):
                errors.append(f"{name}: {mismatches} mismatches, expected {snps}")
            if indel and not any(op == indel for _, op in ops):
                errors.append(f"{name}: expected a {indel} in CIGAR {''.join(f'{n}{o}' for n, o in ops)}")

    print(f"{args.sam}: {total} reads, {mapped} mapped ({fraction:.1%}), {len(errors)} problems")
    if fraction < args.min_mapped:
        errors.append(f"mapped fraction {fraction:.1%} is below --min-mapped {args.min_mapped:.1%}")
    for e in errors[:25]:
        print("  FAIL", e)
    if len(errors) > 25:
        print(f"  ... and {len(errors) - 25} more")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
