#!/usr/bin/env python3
"""Check the command-line contract the desktop app relies on.

  * --progress-json: every stderr line is a JSON object; progress never goes
    backwards and ends at 100; "done" counts agree with the SAM file
  * exit codes and "error" events for usage, input and output problems
  * `-o -` writes exactly the bytes the file output does (no CRLF on Windows)
  * non-ASCII paths work (UTF-8 manifest on Windows)

usage: check_cli.py --run ./run --out tests/out
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

REF = "2019-nCoV.fa"
READS_1K = os.path.join("reads", "reads_1K.fa.gz")
FIXTURE_READS = os.path.join("tests", "fixtures", "reads.fa")

failures = []


def check(condition, message):
    if not condition:
        failures.append(message)
    return condition


def run(exe, *args):
    return subprocess.run([exe, *args], capture_output=True)


def events(proc, label):
    """Parse stderr as JSON lines; record a failure for any line that isn't one."""
    parsed = []
    for line in proc.stderr.decode("utf-8", "replace").splitlines():
        try:
            obj = json.loads(line)
        except json.JSONDecodeError:
            check(False, f"{label}: stderr line is not JSON: {line!r}")
            continue
        check(isinstance(obj, dict) and "event" in obj, f"{label}: JSON line without an event: {line!r}")
        parsed.append(obj)
    return parsed


def sam_length(path):
    """Reference length from the @SQ header line."""
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b"@SQ"):
                return int(line.rstrip().split(b"LN:")[1])
    return None


def sam_counts(path):
    primary = mapped = 0
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b"@"):
                continue
            flag = int(line.split(b"\t")[1])
            if flag & 0x100:
                continue
            primary += 1
            if not flag & 0x4:
                mapped += 1
    return primary, mapped


def expect_error(exe, label, args, code, name):
    proc = run(exe, "--progress-json", *args)
    check(proc.returncode == code, f"{label}: exit code {proc.returncode}, expected {code}")
    evs = events(proc, label)
    last = evs[-1] if evs else {}
    check(last.get("event") == "error" and last.get("code") == name and last.get("message"),
          f"{label}: expected a final error event with code {name!r}, got {last}")
    check(proc.stdout == b"", f"{label}: unexpected stdout output")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="aligner executable")
    ap.add_argument("--out", required=True, help="scratch directory")
    args = ap.parse_args()
    exe = os.path.abspath(args.run)
    os.makedirs(args.out, exist_ok=True)

    # progress events on the 1K reads
    sam = os.path.join(args.out, "cli_1K.sam")
    proc = run(exe, "--progress-json", "-r", REF, "-q", READS_1K, "-o", sam)
    check(proc.returncode == 0, f"progress run: exit code {proc.returncode}")
    evs = events(proc, "progress run")
    kinds = [e.get("event") for e in evs]
    check(kinds[:1] == ["index"] and evs[0].get("length") == sam_length(sam),
          f"progress run: first event {evs[:1]}, SAM @SQ LN {sam_length(sam)}")
    # progress follows the compressed bytes read, so a small gzip file moves in a few large steps
    percents = [e["percent"] for e in evs if e.get("event") == "progress"]
    check(len(percents) >= 2 and percents == sorted(percents) and percents[-1] == 100,
          f"progress run: percents not increasing to 100: {percents}")
    done = evs[-1] if evs else {}
    check(done.get("event") == "done", f"progress run: last event is {done}")
    if done.get("event") == "done":
        primary, mapped = sam_counts(sam)
        check(done["reads"] == 1000 == primary, f"done.reads {done['reads']}, SAM primaries {primary}")
        check(done["mapped"] == mapped, f"done.mapped {done['mapped']}, SAM mapped {mapped}")
        check(done["mapped"] + done["unmapped"] == done["reads"], f"done counts don't add up: {done}")
        check(isinstance(done.get("seconds"), (int, float)), f"done.seconds missing: {done}")

    # errors: exit codes and events
    expect_error(exe, "unknown option", ["--bogus"], 2, "usage")
    expect_error(exe, "missing reads file", ["-r", REF, "-q", os.path.join(args.out, "missing.fa"), "-o",
                                             os.path.join(args.out, "x.sam")], 3, "input")
    expect_error(exe, "unwritable output", ["-r", REF, "-q", FIXTURE_READS, "-o",
                                            os.path.join(args.out, "no-such-dir", "x.sam")], 4, "output")
    proc = run(exe, "--bogus")
    check(proc.returncode == 2 and b"unknown option" in proc.stderr, "unknown option without JSON: wrong exit or message")

    # stdout output is byte-identical to file output
    ascii_sam = os.path.join(args.out, "cli_fixture.sam")
    proc = run(exe, "--quiet", "-r", REF, "-q", FIXTURE_READS, "-o", ascii_sam)
    check(proc.returncode == 0, f"fixture run: exit code {proc.returncode}")
    with open(ascii_sam, "rb") as fh:
        expected = fh.read()
    check(b"\r" not in expected, "SAM output contains carriage returns")
    proc = run(exe, "-r", REF, "-q", FIXTURE_READS, "-o", "-")
    check(proc.returncode == 0 and proc.stdout == expected, "`-o -` output differs from file output")
    check(proc.stderr == b"", f"`-o -` wrote to stderr: {proc.stderr[:200]!r}")

    # non-ASCII paths
    udir = os.path.join(args.out, "ünï cødé")
    shutil.rmtree(udir, ignore_errors=True)
    os.makedirs(udir)
    shutil.copy(FIXTURE_READS, os.path.join(udir, "reads.fa"))
    usam = os.path.join(udir, "out.sam")
    proc = run(exe, "--quiet", "-r", REF, "-q", os.path.join(udir, "reads.fa"), "-o", usam)
    check(proc.returncode == 0, f"non-ASCII path: exit code {proc.returncode}: {proc.stderr[:300]!r}")
    if os.path.exists(usam):
        with open(usam, "rb") as fh:
            check(fh.read() == expected, "non-ASCII path: output differs from ASCII-path run")
    else:
        check(False, "non-ASCII path: no output file")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"check_cli: {len(failures)} problem(s)")
        return 1
    print("check_cli: progress events, exit codes, stdout output and non-ASCII paths OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
