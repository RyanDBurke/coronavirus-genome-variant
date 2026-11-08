// Types shared by the main process, the preload bridge and the window.

/** A built-in reference or read set shipped with the app. */
export interface Preset {
  id: string;
  label: string;
  detail: string;
}

export interface Presets {
  references: Preset[];
  reads: Preset[];
}

/** A file the user picked in a dialog. The window only ever sees the token, never sends paths. */
export interface PickedFile {
  token: string;
  name: string;
  bytes: number;
}

/** Where an input comes from: a built-in preset, or a file picked earlier in this session. */
export type InputChoice = { kind: "preset"; id: string } | { kind: "file"; token: string };

export interface RunRequest {
  reference: InputChoice;
  reads: InputChoice;
}

/** One line of the aligner's --progress-json output. */
export type AlignerEvent =
  | { event: "index"; length: number }
  | { event: "progress"; percent: number }
  | { event: "done"; reads: number; mapped: number; unmapped: number; seconds: number }
  | { event: "warning"; message: string }
  | { event: "error"; code: string; message: string };

/** Events the window receives during a run: the aligner's own, plus the app's results-preparation step. */
export type RunEvent = AlignerEvent | { event: "analyze"; percent: number };

/** A site where reads disagree with the reference. Positions are 1-based. */
export interface Variant {
  type: "SNV" | "deletion" | "insertion";
  pos: number;
  /** reference positions covered: deleted bases, or for insertions the span of merged insertion sites */
  length: number;
  ref: string;
  alt: string;
  depth: number;
  altReads: number;
  freq: number;
}

/** Everything the viewer shows before zooming in. */
export interface ViewSummary {
  refName: string;
  length: number;
  meanDepth: number;
  maxDepth: number;
  /** share of positions with depth >= 1 and >= 10 */
  covered1: number;
  covered10: number;
  /** genome-wide depth, binned */
  binSize: number;
  binMean: number[];
  binMax: number[];
  variants: Variant[];
  minDepth: number;
  minFreq: number;
  /** insertion sites at most this many bp apart are reported as one */
  insertionMergeBp: number;
}

/** One aligned read, as the pileup draws it. start is 0-based. */
export interface ReadRow {
  name: string;
  flag: number;
  start: number;
  cigar: string;
  seq: string;
}

/** A zoomed-in region: [start, end) 0-based, with per-base depth and (for small regions) reads. */
export interface RegionData {
  start: number;
  end: number;
  ref: string;
  depth: number[];
  /** per position: A, C, G, T, N, deleted, insertion-after (7 numbers each); empty for wide regions */
  counts: number[];
  reads: ReadRow[];
  /** reads overlapping the region; reads[] may be a sample of them */
  totalReads: number;
}

/** Plain-language description of a failed run. */
export interface FriendlyError {
  title: string;
  detail: string;
  hint: string;
}

export type RunResult =
  | {
      status: "done"; runId: string; reads: number; mapped: number; unmapped: number; seconds: number; samBytes: number;
      /** null if the results couldn't be prepared for viewing (summaryError says why) */
      summary: ViewSummary | null;
      summaryError?: string;
    }
  | { status: "cancelled"; runId: string }
  | { status: "failed"; runId: string; error: FriendlyError };

/** Shown in the About box. */
export interface AppInfo {
  name: string;
  version: string;
}

/** License files the About box can open: this app's, Electron's, Chromium's third-party notices,
 *  and the bundled fonts'. */
export type LicenseId = "app" | "electron" | "chromium" | "fonts";

/** The API the preload script exposes to the window as `window.aligner`. */
export interface AlignerApi {
  getAppInfo(): Promise<AppInfo>;
  /** opens a bundled license file in the default viewer */
  openLicense(id: LicenseId): Promise<void>;
  getPresets(): Promise<Presets>;
  chooseFile(kind: "reference" | "reads"): Promise<PickedFile | null>;
  start(request: RunRequest): Promise<{ runId: string } | { error: FriendlyError }>;
  cancel(): Promise<void>;
  saveSam(runId: string): Promise<string | null>;
  reveal(runId: string): Promise<void>;
  /** depth, reference and reads for [start, end) (0-based); null if the run is no longer available */
  region(runId: string, start: number, end: number): Promise<RegionData | null>;
  exportVariants(runId: string): Promise<string | null>;
  exportCoverage(runId: string): Promise<string | null>;
  onEvent(listener: (event: RunEvent) => void): () => void;
  onFinished(listener: (result: RunResult) => void): () => void;
}
