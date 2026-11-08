// Runs the bundled aligner as a child process and turns its --progress-json output into events.
// Kept free of Electron imports so it can be tested with plain Node.

import { spawn } from "node:child_process";
import type { AlignerEvent } from "../shared/types";

export interface AlignerExit {
  /** process exit code, or null if it was killed */
  code: number | null;
  cancelled: boolean;
  /** the last "error" event, if any */
  error?: { code: string; message: string };
  /** stderr lines that weren't JSON (e.g. a crash message) */
  otherOutput: string[];
  /** set when the process couldn't be started at all (e.g. ENOENT) */
  spawnError?: string;
}

export interface AlignerRun {
  done: Promise<AlignerExit>;
  cancel(): void;
}

/** Splits a stream of chunks into complete lines, keeping any partial line for the next chunk. */
export class LineSplitter {
  private pending = "";

  push(chunk: string): string[] {
    const text = this.pending + chunk;
    const lines = text.split(/\r?\n/);
    this.pending = lines.pop() ?? "";
    return lines.filter((line) => line.length > 0);
  }

  flush(): string[] {
    const rest = this.pending;
    this.pending = "";
    return rest.length > 0 ? [rest] : [];
  }
}

/** Parses one stderr line; returns null for anything that isn't an aligner event. */
export function parseEvent(line: string): AlignerEvent | null {
  try {
    const value = JSON.parse(line) as unknown;
    if (value !== null && typeof value === "object" && typeof (value as { event?: unknown }).event === "string") {
      return value as AlignerEvent;
    }
  } catch {
    // not JSON
  }
  return null;
}

export interface RunOptions {
  exe: string;
  reference: string;
  reads: string;
  samOut: string;
  indexOut: string;
  cwd: string;
  onEvent(event: AlignerEvent): void;
}

export function runAligner(options: RunOptions): AlignerRun {
  const args = [
    "--progress-json",
    "-r", options.reference,
    "-q", options.reads,
    "-o", options.samOut,
    "-i", options.indexOut,
  ];

  // no shell: arguments go straight to the process, so file names can't be interpreted as commands
  const child = spawn(options.exe, args, { cwd: options.cwd, windowsHide: true, stdio: ["ignore", "ignore", "pipe"] });
  let cancelled = false;
  let error: AlignerExit["error"];
  const otherOutput: string[] = [];
  const splitter = new LineSplitter();

  const handleLine = (line: string) => {
    const event = parseEvent(line);
    if (event === null) {
      otherOutput.push(line);
      return;
    }
    if (event.event === "error") {
      error = { code: event.code, message: event.message };
    }
    options.onEvent(event);
  };

  child.stderr.setEncoding("utf8");
  child.stderr.on("data", (chunk: string) => splitter.push(chunk).forEach(handleLine));

  const done = new Promise<AlignerExit>((resolve) => {
    child.on("error", (err) => {
      resolve({ code: null, cancelled, otherOutput, spawnError: err.message });
    });
    child.on("close", (code) => {
      splitter.flush().forEach(handleLine);
      resolve({ code, cancelled, error, otherOutput });
    });
  });

  return {
    done,
    cancel() {
      if (child.exitCode === null && !cancelled) {
        cancelled = true;
        child.kill();
      }
    },
  };
}
