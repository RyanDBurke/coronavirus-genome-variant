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

/** Plain-language description of a failed run. */
export interface FriendlyError {
  title: string;
  detail: string;
  hint: string;
}

export type RunResult =
  | { status: "done"; runId: string; reads: number; mapped: number; unmapped: number; seconds: number; samBytes: number }
  | { status: "cancelled"; runId: string }
  | { status: "failed"; runId: string; error: FriendlyError };

/** The API the preload script exposes to the window as `window.aligner`. */
export interface AlignerApi {
  getPresets(): Promise<Presets>;
  chooseFile(kind: "reference" | "reads"): Promise<PickedFile | null>;
  start(request: RunRequest): Promise<{ runId: string } | { error: FriendlyError }>;
  cancel(): Promise<void>;
  saveSam(runId: string): Promise<string | null>;
  reveal(runId: string): Promise<void>;
  onEvent(listener: (event: AlignerEvent) => void): () => void;
  onFinished(listener: (result: RunResult) => void): () => void;
}
