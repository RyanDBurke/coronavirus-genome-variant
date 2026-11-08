// Tests for the SAM analysis behind the results viewer: depth, allele counts, variants, region queries, exports.
import assert from "node:assert/strict";
import { existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { test } from "node:test";
import { fileURLToPath } from "node:url";
import { createRequire } from "node:module";
import { gzipSync } from "node:zlib";

const require = createRequire(import.meta.url);
const a = require("../build/main/analysis.js");
const { runAligner } = require("../build/main/aligner.js");

const repo = fileURLToPath(new URL("../../", import.meta.url));
const exe = path.join(repo, process.platform === "win32" ? "run.exe" : "run");

//            1234567890123456789012345678901234567890
const REF = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

function sam(records) {
  const header = `@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:toy\tLN:${REF.length}\n`;
  return header + records.map(([name, flag, pos, cigar, seq]) =>
    `${name}\t${flag}\ttoy\t${pos}\t255\t${cigar}\t*\t0\t0\t${seq}\t*`).join("\n") + "\n";
}

async function withSam(records, fn) {
  const dir = mkdtempSync(path.join(tmpdir(), "gva-an-"));
  try {
    const file = path.join(dir, "x.sam");
    writeFileSync(file, sam(records));
    return await fn(file, dir);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
}

const ref = { name: "toy", seq: REF };

test("depth and allele counts follow the CIGAR; unmapped and secondary records are ignored", async () => {
  await withSam([
    ["exact", 0, 1, "10M", REF.slice(0, 10)],
    ["snp", 0, 1, "10M", "ACGTTCGTAC"],             // position 5: A -> T
    ["del", 0, 1, "4M2D4M", "ACGTGTAC"],             // positions 5-6 deleted
    ["ins", 0, 1, "4M1I4M", "ACGTGACGT"],            // insertion after position 4
    ["second", 256, 1, "10M", REF.slice(0, 10)],     // secondary: ignored
    ["none", 4, 0, "*", "ACGT"],                     // unmapped: ignored
  ], async (file) => {
    const an = await a.analyzeSam(file, ref);
    assert.deepEqual(an.totals, { primary: 5, mapped: 4 });
    assert.equal(an.reads.count, 4);
    const counts = (p) => Array.from(an.counts.subarray((p - 1) * 7, p * 7));   // 1-based position
    assert.equal(an.depth[0], 4);                    // all four reads cover position 1
    assert.deepEqual(counts(5), [2, 0, 0, 1, 0, 1, 0]);  // A x2 (exact, ins), T (snp), deleted (del)
    assert.deepEqual(counts(4), [0, 0, 0, 4, 0, 0, 1]);  // T x4, insertion after it
    assert.equal(an.depth[4], 4);
    assert.equal(an.depth[9], 3);                    // del read ends at 10, ins read at 8
    assert.equal(an.reads.maxSpan, 10);
  });
});

test("variants: SNVs, merged deletions and insertions above the thresholds", async () => {
  const reads = [];
  for (let i = 0; i < 6; i++) reads.push([`ref${i}`, 0, 1, "12M", REF.slice(0, 12)]);
  for (let i = 0; i < 4; i++) reads.push([`snp${i}`, 0, 1, "12M", "ACGTACGTAGGT"]);    // position 10: C -> G (40%)
  for (let i = 0; i < 3; i++) reads.push([`del${i}`, 0, 1, "2M3D7M", "ACCGTACGT"]);    // positions 3-5 deleted (3/14)
  for (let i = 0; i < 1; i++) reads.push([`ins${i}`, 0, 1, "6M2I6M", "ACGTACTTGTACGT"]); // 1/14: below 20%
  await withSam(reads, async (file) => {
    const an = await a.analyzeSam(file, ref);
    const variants = a.callVariants(an);
    assert.deepEqual(variants.map((v) => `${v.type}@${v.pos}:${v.ref}>${v.alt}/${v.length}`),
      ["deletion@3:GTA>-/3", "SNV@10:C>G/1"]);
    const snv = variants.find((v) => v.type === "SNV");
    assert.equal(snv.depth, 14);
    assert.equal(snv.altReads, 4);
    // below minimum depth nothing is called
    assert.deepEqual(a.callVariants(an, 100, 0.2), []);
  });
});

test("summary: bins, depth statistics and coverage", async () => {
  await withSam([["r", 0, 11, "20M", REF.slice(10, 30)]], async (file) => {
    const s = a.summarize(await a.analyzeSam(file, ref), 4);
    assert.equal(s.length, 40);
    assert.equal(s.binSize, 10);
    assert.deepEqual(s.binMax, [0, 1, 1, 0]);
    assert.equal(s.maxDepth, 1);
    assert.equal(s.covered1, 0.5);
    assert.equal(s.covered10, 0);
    assert.equal(s.meanDepth, 0.5);
  });
});

test("region queries return overlapping reads straight from the file, sampled to the limit", async () => {
  const reads = [];
  for (let s = 1; s <= 31; s += 2) reads.push([`r${s}`, 0, s, "10M", REF.slice(s - 1, s + 9)]);
  await withSam(reads, async (file) => {
    const an = await a.analyzeSam(file, ref);
    const r = await a.readRegion(an, file, 20, 22);
    // a read starting at 1-based s covers 0-based [s-1, s+9), so s = 13..21 overlap [20, 22)
    assert.deepEqual(r.reads.map((x) => x.name).sort(), ["r13", "r15", "r17", "r19", "r21"]);
    assert.equal(r.totalReads, 5);
    assert.equal(r.ref, REF.slice(20, 22));
    assert.equal(r.counts.length, 14);
    assert.equal(r.reads.find((x) => x.name === "r13").cigar, "10M");
    const sampled = await a.readRegion(an, file, 0, 40, 4);
    assert.equal(sampled.totalReads, 16);
    assert.equal(sampled.reads.length, 4);
    // the region is clamped to the genome
    const clamped = await a.readRegion(an, file, -50, 1000);
    assert.equal(clamped.start, 0);
    assert.equal(clamped.end, 40);
  });
});

test("exports and reference loading (plain and gzipped)", async () => {
  await withSam([["r", 0, 1, "4M", "ACGA"]], async (file, dir) => {
    const an = await a.analyzeSam(file, ref);
    const csv = a.coverageCsv(an).split("\n");
    assert.equal(csv[0], "position,reference,depth,A,C,G,T,N,deleted,insertion_after");
    assert.equal(csv[4], "4,T,1,1,0,0,0,0,0,0");
    assert.equal(csv.length, 42);
    const v = a.variantsCsv([{ type: "SNV", pos: 4, length: 1, ref: "T", alt: "A", depth: 10, altReads: 9, freq: 0.9 }], "toy");
    assert.equal(v, "reference,position,type,length,ref,alt,depth,alt_reads,frequency\ntoy,4,SNV,1,T,A,10,9,0.9000\n");

    const fasta = ">chr1 description\r\nacgt\r\nNNac\r\n>second\nTTTT\n";
    writeFileSync(path.join(dir, "r.fa"), fasta);
    writeFileSync(path.join(dir, "r.fa.gz"), gzipSync(fasta));
    for (const name of ["r.fa", "r.fa.gz"]) {
      assert.deepEqual(await a.loadReference(path.join(dir, name)), { name: "chr1", seq: "ACGTNNAC" });
    }
  });
});

test("on real aligner output, depth adds up to the reference bases each read covers", { skip: !existsSync(exe) && "aligner not built" }, async () => {
  const dir = mkdtempSync(path.join(tmpdir(), "gva-an-real-"));
  try {
    const samPath = path.join(dir, "out.sam");
    const run = runAligner({ exe, reference: path.join(repo, "2019-nCoV.fa"), reads: path.join(repo, "tests", "fixtures", "reads.fa"),
      samOut: samPath, indexOut: path.join(dir, "fm.txt"), cwd: dir, onEvent: () => {} });
    assert.equal((await run.done).code, 0);
    const reference = await a.loadReference(path.join(repo, "2019-nCoV.fa"));
    const an = await a.analyzeSam(samPath, reference);
    let expected = 0;
    for (const line of readFileSync(samPath, "latin1").split("\n")) {
      const f = line.split("\t");
      if (line.startsWith("@") || f.length < 10 || Number(f[1]) & 0x104) continue;
      for (const [n, op] of a.parseCigar(f[5])) if (op === "M" || op === "D") expected += n;
    }
    assert.equal(an.depth.reduce((s, d) => s + d, 0), expected);
    assert.equal(an.totals.mapped, 52);
    assert.equal(an.length, 29882);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});
