// Turns an aligner SAM file into what the results viewer needs, without loading the SAM into memory:
//   * per-base depth and allele counts (A, C, G, T, N, deleted, insertion-after) over the reference
//   * a position index of mapped reads (start, end, byte offset), so any region's reads can be read
//     straight from the file
// The SARS-CoV-2 genome is ~30 kb, so even 1M reads give a few MB of typed arrays.
// Kept free of Electron imports so it can be tested with plain Node.

import { createReadStream, promises as fs } from "node:fs";
import { gunzipSync } from "node:zlib";
import type { ReadRow, RegionData, Variant, ViewSummary } from "../shared/types";

/** columns of `counts`, per reference position */
export const COL = { A: 0, C: 1, G: 2, T: 3, N: 4, DEL: 5, INS: 6 } as const;
const NCOL = 7;
const BASE_COL: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };

const FLAG_UNMAPPED = 0x4;
const FLAG_SECONDARY = 0x100;

export interface Analysis {
  refName: string;
  ref: string;
  length: number;
  /** reads covering each position (bases plus deletions) */
  depth: Uint32Array;
  /** length * 7 counts, see COL */
  counts: Uint32Array;
  /** primary mapped reads, sorted by start */
  reads: { start: Int32Array; end: Int32Array; offset: Float64Array; bytes: Int32Array; count: number; maxSpan: number };
  totals: { primary: number; mapped: number };
}

// ---- reference ------------------------------------------------------------------

/** First sequence of a FASTA file (optionally gzipped), upper-cased. */
export async function loadReference(path: string): Promise<{ name: string; seq: string }> {
  let data = await fs.readFile(path);
  if (data[0] === 0x1f && data[1] === 0x8b) data = gunzipSync(data);
  const text = data.toString("latin1");
  const header = text.indexOf(">");
  if (header < 0) throw new Error("reference has no FASTA header");
  const headerEnd = text.indexOf("\n", header);
  const name = text.slice(header + 1, headerEnd < 0 ? undefined : headerEnd).trim().split(/\s+/)[0] ?? "ref";
  const next = headerEnd < 0 ? -1 : text.indexOf(">", headerEnd);
  const body = headerEnd < 0 ? "" : text.slice(headerEnd + 1, next < 0 ? undefined : next);
  return { name, seq: body.replace(/[^A-Za-z]/g, "").toUpperCase() };
}

// ---- CIGAR ------------------------------------------------------------------------

const CIGAR_RE = /(\d+)([MIDNSHP=X])/g;

export function parseCigar(cigar: string): Array<[number, string]> {
  const ops: Array<[number, string]> = [];
  for (const m of cigar.matchAll(CIGAR_RE)) ops.push([Number(m[1]), m[2]]);
  return ops;
}

// ---- growable typed arrays -------------------------------------------------------------

class IntList {
  data: Int32Array | Float64Array;
  length = 0;
  constructor(private readonly kind: "i32" | "f64", capacity = 1 << 16) {
    this.data = kind === "i32" ? new Int32Array(capacity) : new Float64Array(capacity);
  }
  push(value: number) {
    if (this.length === this.data.length) {
      const bigger = this.kind === "i32" ? new Int32Array(this.data.length * 2) : new Float64Array(this.data.length * 2);
      bigger.set(this.data);
      this.data = bigger;
    }
    this.data[this.length++] = value;
  }
}

// ---- the main pass --------------------------------------------------------------------

/** One streaming pass over the SAM file. onProgress gets the fraction of bytes read. */
export async function analyzeSam(samPath: string, reference: { name: string; seq: string },
                                 onProgress?: (fraction: number) => void): Promise<Analysis> {
  const ref = reference.seq;
  const length = ref.length;
  const depth = new Uint32Array(length);
  const counts = new Uint32Array(length * NCOL);
  const starts = new IntList("i32");
  const ends = new IntList("i32");
  const offsets = new IntList("f64");
  const sizes = new IntList("i32");
  let primary = 0;
  let mapped = 0;

  const total = (await fs.stat(samPath)).size;
  let fileOffset = 0;        // byte offset of `pending` in the file
  let pending: Buffer = Buffer.alloc(0);
  let lastReported = -1;

  const handleLine = (line: string, offset: number, bytes: number) => {
    if (line.length === 0 || line.charCodeAt(0) === 64 /* @ */) return;
    const f = line.split("\t", 10);
    if (f.length < 10) return;
    const flag = Number(f[1]);
    if (flag & FLAG_SECONDARY) return;
    primary++;
    if (flag & FLAG_UNMAPPED || f[5] === "*") return;
    mapped++;

    const seq = f[9];
    let r = Number(f[3]) - 1;
    const start = r;
    let q = 0;
    for (const [n, op] of parseCigar(f[5])) {
      switch (op) {
        case "M": case "=": case "X":
          for (let k = 0; k < n; k++, r++, q++) {
            if (r < 0 || r >= length) continue;
            depth[r]++;
            counts[r * NCOL + (BASE_COL[seq[q]] ?? COL.N)]++;
          }
          break;
        case "D":
          for (let k = 0; k < n; k++, r++) {
            if (r < 0 || r >= length) continue;
            depth[r]++;
            counts[r * NCOL + COL.DEL]++;
          }
          break;
        case "N":
          r += n;
          break;
        case "I":
          if (r > 0 && r <= length) counts[(r - 1) * NCOL + COL.INS]++;
          q += n;
          break;
        case "S":
          q += n;
          break;
        // H and P consume nothing
      }
    }
    starts.push(start);
    ends.push(r);
    offsets.push(offset);
    sizes.push(bytes);
  };

  for await (const chunk of createReadStream(samPath, { highWaterMark: 1 << 20 }) as AsyncIterable<Buffer>) {
    const buf = pending.length ? Buffer.concat([pending, chunk]) : chunk;
    let lineStart = 0;
    for (let nl = buf.indexOf(10, lineStart); nl >= 0; nl = buf.indexOf(10, lineStart)) {
      handleLine(buf.toString("latin1", lineStart, nl), fileOffset + lineStart, nl - lineStart);
      lineStart = nl + 1;
    }
    pending = buf.subarray(lineStart);
    fileOffset += lineStart;
    const fraction = fileOffset / Math.max(total, 1);
    const percent = Math.floor(fraction * 100);
    if (onProgress && percent !== lastReported) {
      lastReported = percent;
      onProgress(fraction);
    }
  }
  if (pending.length) handleLine(pending.toString("latin1"), fileOffset, pending.length);

  // counting sort by start position (positions are bounded by the genome length)
  const n = starts.length;
  const order = new Int32Array(n);
  const bucket = new Int32Array(length + 2);
  const clampStart = (s: number) => Math.min(Math.max(s, 0), length);
  for (let i = 0; i < n; i++) bucket[clampStart(starts.data[i]) + 1]++;
  for (let p = 1; p < bucket.length; p++) bucket[p] += bucket[p - 1];
  for (let i = 0; i < n; i++) order[bucket[clampStart(starts.data[i])]++] = i;

  const reads = {
    start: new Int32Array(n), end: new Int32Array(n), offset: new Float64Array(n), bytes: new Int32Array(n),
    count: n, maxSpan: 0,
  };
  for (let k = 0; k < n; k++) {
    const i = order[k];
    reads.start[k] = starts.data[i];
    reads.end[k] = ends.data[i];
    reads.offset[k] = offsets.data[i];
    reads.bytes[k] = sizes.data[i];
    reads.maxSpan = Math.max(reads.maxSpan, reads.end[k] - reads.start[k]);
  }

  onProgress?.(1);
  return { refName: reference.name, ref, length, depth, counts, reads, totals: { primary, mapped } };
}

// ---- variants -----------------------------------------------------------------------------

export const MIN_DEPTH = 10;
export const MIN_FREQ = 0.2;
/** Insertion sites this close together are one event. The aligner can split a multi-base insertion
 *  around a matching base (1I1M1I instead of 2I), which otherwise shows up as two neighbouring sites. */
export const INSERTION_MERGE_BP = 2;

/** Sites where at least MIN_FREQ of at least MIN_DEPTH reads differ from the reference.
 *  Consecutive deleted positions are merged into one deletion, and insertion sites within
 *  INSERTION_MERGE_BP of each other into one insertion spanning them (counts from its strongest site).
 *  Positions are 1-based. */
export function callVariants(a: Analysis, minDepth = MIN_DEPTH, minFreq = MIN_FREQ): Variant[] {
  const out: Variant[] = [];
  let openDeletion: Variant | null = null;
  let openInsertion: Variant | null = null;

  for (let p = 0; p < a.length; p++) {
    const d = a.depth[p];
    const row = p * NCOL;
    const refBase = a.ref[p];

    // deletions, merged across consecutive positions
    const del = a.counts[row + COL.DEL];
    if (d >= minDepth && del / d >= minFreq) {
      if (openDeletion && openDeletion.pos + openDeletion.length === p + 1) {
        openDeletion.length++;
        openDeletion.ref += refBase;
        openDeletion.depth = Math.min(openDeletion.depth, d);
        openDeletion.altReads = Math.min(openDeletion.altReads, del);
        openDeletion.freq = openDeletion.altReads / openDeletion.depth;
      } else {
        openDeletion = { type: "deletion", pos: p + 1, length: 1, ref: refBase, alt: "-", depth: d, altReads: del, freq: del / d };
        out.push(openDeletion);
      }
    } else {
      openDeletion = null;
    }
    if (d < minDepth) continue;

    for (const base of ["A", "C", "G", "T"] as const) {
      if (base === refBase) continue;
      const c = a.counts[row + COL[base]];
      if (c / d >= minFreq) {
        out.push({ type: "SNV", pos: p + 1, length: 1, ref: refBase, alt: base, depth: d, altReads: c, freq: c / d });
      }
    }
    const ins = a.counts[row + COL.INS];
    if (ins / d >= minFreq) {
      const pos = p + 1;
      if (openInsertion && pos - (openInsertion.pos + openInsertion.length - 1) <= INSERTION_MERGE_BP) {
        openInsertion.length = pos - openInsertion.pos + 1;
        openInsertion.ref = a.ref.slice(openInsertion.pos - 1, pos);
        if (ins > openInsertion.altReads) {
          openInsertion.altReads = ins;
          openInsertion.depth = d;
          openInsertion.freq = ins / d;
        }
      } else {
        openInsertion = { type: "insertion", pos, length: 1, ref: refBase, alt: "+", depth: d, altReads: ins, freq: ins / d };
        out.push(openInsertion);
      }
    }
  }
  return out;
}

// ---- summary for the window --------------------------------------------------------------

export const SUMMARY_BINS = 1200;

export function summarize(a: Analysis, bins = SUMMARY_BINS): ViewSummary {
  const binSize = Math.max(1, Math.ceil(a.length / bins));
  const binMean: number[] = [];
  const binMax: number[] = [];
  let sum = 0;
  let maxDepth = 0;
  let covered1 = 0;
  let covered10 = 0;
  for (let b = 0; b * binSize < a.length; b++) {
    let s = 0;
    let m = 0;
    const end = Math.min(a.length, (b + 1) * binSize);
    for (let p = b * binSize; p < end; p++) {
      const d = a.depth[p];
      s += d;
      if (d > m) m = d;
    }
    binMean.push(s / (end - b * binSize));
    binMax.push(m);
  }
  for (let p = 0; p < a.length; p++) {
    const d = a.depth[p];
    sum += d;
    if (d > maxDepth) maxDepth = d;
    if (d >= 1) covered1++;
    if (d >= 10) covered10++;
  }
  return {
    refName: a.refName,
    length: a.length,
    meanDepth: sum / Math.max(a.length, 1),
    maxDepth,
    covered1: covered1 / Math.max(a.length, 1),
    covered10: covered10 / Math.max(a.length, 1),
    binSize,
    binMean,
    binMax,
    variants: callVariants(a),
    minDepth: MIN_DEPTH,
    minFreq: MIN_FREQ,
    insertionMergeBp: INSERTION_MERGE_BP,
  };
}

// ---- regions -----------------------------------------------------------------------------

/** Largest region (bp) for which reads are returned; wider regions get depth only. */
export const MAX_READ_REGION = 2000;

/** first index with start >= value */
function lowerBound(arr: Int32Array, count: number, value: number) {
  let lo = 0;
  let hi = count;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid] < value) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

/** Reads overlapping [start, end) (0-based), evenly sampled down to `limit`, with depth and reference. */
export async function readRegion(a: Analysis, samPath: string, start: number, end: number, limit = 600): Promise<RegionData> {
  start = Math.max(0, Math.min(a.length - 1, Math.floor(start)));
  end = Math.max(start + 1, Math.min(a.length, Math.ceil(end)));
  const region: RegionData = {
    start, end,
    ref: a.ref.slice(start, end),
    depth: Array.from(a.depth.subarray(start, end)),
    counts: [],
    reads: [],
    totalReads: 0,
  };
  if (end - start > MAX_READ_REGION) return region;
  region.counts = Array.from(a.counts.subarray(start * NCOL, end * NCOL));

  const { reads } = a;
  const overlapping: number[] = [];
  for (let i = lowerBound(reads.start, reads.count, start - reads.maxSpan); i < reads.count && reads.start[i] < end; i++) {
    if (reads.end[i] > start) overlapping.push(i);
  }
  region.totalReads = overlapping.length;

  const picked = overlapping.length <= limit
    ? overlapping
    : Array.from({ length: limit }, (_, k) => overlapping[Math.floor((k * overlapping.length) / limit)]);

  const file = await fs.open(samPath, "r");
  try {
    region.reads = await Promise.all(picked.map(async (i): Promise<ReadRow> => {
      const buf = Buffer.alloc(reads.bytes[i]);
      await file.read(buf, 0, buf.length, reads.offset[i]);
      const f = buf.toString("latin1").split("\t", 10);
      return { name: f[0], flag: Number(f[1]), start: Number(f[3]) - 1, cigar: f[5], seq: f[9] };
    }));
  } finally {
    await file.close();
  }
  return region;
}

// ---- exports -------------------------------------------------------------------------------

export function coverageCsv(a: Analysis): string {
  const lines = ["position,reference,depth,A,C,G,T,N,deleted,insertion_after"];
  for (let p = 0; p < a.length; p++) {
    const c = a.counts.subarray(p * NCOL, p * NCOL + NCOL);
    lines.push(`${p + 1},${a.ref[p]},${a.depth[p]},${c[0]},${c[1]},${c[2]},${c[3]},${c[4]},${c[5]},${c[6]}`);
  }
  return lines.join("\n") + "\n";
}

export function variantsCsv(variants: Variant[], refName: string): string {
  const lines = ["reference,position,type,length,ref,alt,depth,alt_reads,frequency"];
  for (const v of variants) {
    lines.push(`${refName},${v.pos},${v.type},${v.length},${v.ref},${v.alt},${v.depth},${v.altReads},${v.freq.toFixed(4)}`);
  }
  return lines.join("\n") + "\n";
}
