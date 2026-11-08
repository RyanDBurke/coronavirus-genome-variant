// Electron main process: owns the window, the file dialogs and the aligner process.
// The window can't touch the file system or start processes itself; it asks through the
// narrow IPC API below, and only ever refers to files by tokens handed out by our own dialogs.

import { app, BrowserWindow, dialog, ipcMain, Menu, nativeTheme, session, shell, type IpcMainInvokeEvent } from "electron";
import { randomUUID } from "node:crypto";
import { promises as fs } from "node:fs";
import path from "node:path";
import { pathToFileURL } from "node:url";
import type { InputChoice, PickedFile, Presets, RunRequest, RunResult, Variant, ViewSummary } from "../shared/types";
import { runAligner, type AlignerRun } from "./aligner";
import { analyzeSam, coverageCsv, loadReference, readRegion, summarize, variantsCsv, type Analysis } from "./analysis";
import { friendlyError } from "./messages";

// ---- locations ---------------------------------------------------------------

// Tests point the app at a throwaway profile; packaged builds always use the normal one.
if (!app.isPackaged && process.env.GVA_USER_DATA) {
  app.setPath("userData", process.env.GVA_USER_DATA);
}

/** Installed: everything is under resources/. From source: the repo root, one level above app/. */
const resources = app.isPackaged ? process.resourcesPath : path.resolve(app.getAppPath(), "..");
const dataDir = app.isPackaged ? path.join(resources, "data") : resources;
const alignerExe = app.isPackaged
  ? path.join(resources, "aligner", "run.exe")
  : path.join(resources, process.platform === "win32" ? "run.exe" : "run");
const runsDir = () => path.join(app.getPath("userData"), "runs");
const KEEP_RUNS = 3;

const indexHtml = path.join(__dirname, "..", "renderer", "index.html");
const indexUrl = pathToFileURL(indexHtml).href;

// ---- built-in inputs -----------------------------------------------------------

const referenceFiles: Record<string, string> = {
  "sars-cov-2": "2019-nCoV.fa",
};
const readFiles: Record<string, string> = {
  "1k": path.join("reads", "reads_1K.fa.gz"),
  "10k": path.join("reads", "reads_10K.fa.gz"),
  "1m": path.join("reads", "reads_1M.fa.gz"),
};

const presets: Presets = {
  references: [{ id: "sars-cov-2", label: "SARS-CoV-2", detail: "MN988713.1 · 29,882 bp" }],
  reads: [
    { id: "1k", label: "1,000 reads", detail: "Simulated · under a second" },
    { id: "10k", label: "10,000 reads", detail: "Simulated · a few seconds" },
    { id: "1m", label: "1,000,000 reads", detail: "Simulated · 2–3 minutes" },
  ],
};

/** Files picked in dialogs this session: token -> absolute path. */
const pickedFiles = new Map<string, string>();

function resolveInput(choice: InputChoice, table: Record<string, string>): string | null {
  if (choice?.kind === "preset" && typeof choice.id === "string" && Object.hasOwn(table, choice.id)) {
    return path.join(dataDir, table[choice.id]);
  }
  if (choice?.kind === "file" && typeof choice.token === "string") {
    return pickedFiles.get(choice.token) ?? null;
  }
  return null;
}

// ---- runs ----------------------------------------------------------------------

interface ActiveRun {
  id: string;
  aligner: AlignerRun;
  /** set by Cancel; also stops the results-preparation step */
  cancelled: boolean;
}
let activeRun: ActiveRun | null = null;

interface FinishedRun {
  samPath: string;
  label: string;
  analysis: Analysis | null;
  variants: Variant[];
}
const finishedRuns = new Map<string, FinishedRun>();

class Cancelled extends Error {}

function samPathFor(runId: string) {
  return path.join(runsDir(), runId, "alignments.sam");
}

/** Deletes all but the newest KEEP_RUNS run folders (the 1M result is ~150 MB each). */
async function pruneRuns() {
  const entries = await fs.readdir(runsDir(), { withFileTypes: true }).catch(() => []);
  const dirs = entries.filter((e) => e.isDirectory()).map((e) => e.name).sort();
  for (const name of dirs.slice(0, Math.max(0, dirs.length - KEEP_RUNS))) {
    finishedRuns.delete(name);
    await fs.rm(path.join(runsDir(), name), { recursive: true, force: true });
  }
}

async function startRun(window: BrowserWindow, request: RunRequest) {
  if (activeRun !== null) {
    return { error: { title: "A run is already in progress", detail: "", hint: "Wait for it to finish, or cancel it." } };
  }
  const reference = resolveInput(request?.reference, referenceFiles);
  const reads = resolveInput(request?.reads, readFiles);
  if (reference === null || reads === null) {
    return { error: { title: "Choose a reference and a reads file", detail: "", hint: "" } };
  }

  // run ids sort by start time, which pruneRuns relies on
  const runId = new Date().toISOString().replace(/[:.]/g, "-") + "-" + randomUUID().slice(0, 8);
  const runDir = path.join(runsDir(), runId);
  await fs.mkdir(runDir, { recursive: true });
  await pruneRuns();

  const label = request.reads.kind === "preset" ? `reads-${request.reads.id}` : path.parse(reads).name.replace(/\.(fa|fasta|fq|fastq)$/i, "");
  let counts = { reads: 0, mapped: 0, unmapped: 0, seconds: 0 };
  const aligner = runAligner({
    exe: alignerExe,
    reference,
    reads,
    samOut: samPathFor(runId),
    indexOut: path.join(runDir, "FMindex.txt"),
    cwd: runDir,
    onEvent: (event) => {
      if (event.event === "done") {
        counts = { reads: event.reads, mapped: event.mapped, unmapped: event.unmapped, seconds: event.seconds };
      }
      if (!window.isDestroyed()) window.webContents.send("run:event", event);
    },
  });
  const run: ActiveRun = { id: runId, aligner, cancelled: false };
  activeRun = run;
  const send = (channel: string, value: unknown) => {
    if (!window.isDestroyed()) window.webContents.send(channel, value);
  };

  void aligner.done.then(async (exit) => {
    let result: RunResult | undefined;
    if (exit.cancelled || run.cancelled) {
      result = { status: "cancelled", runId };
    } else if (exit.code === 0) {
      const samPath = samPathFor(runId);
      const stat = await fs.stat(samPath).catch(() => null);
      const finished: FinishedRun = { samPath, label, analysis: null, variants: [] };
      let summary: ViewSummary | null = null;
      let summaryError: string | undefined;

      // prepare the viewer: one pass over the SAM for depth, allele counts and a read index
      try {
        send("run:event", { event: "analyze", percent: 0 });
        const ref = await loadReference(reference);
        finished.analysis = await analyzeSam(samPath, ref, (fraction) => {
          if (run.cancelled) throw new Cancelled();
          send("run:event", { event: "analyze", percent: Math.floor(fraction * 100) });
        });
        summary = summarize(finished.analysis);
        finished.variants = summary.variants;
      } catch (err) {
        if (err instanceof Cancelled) {
          result = { status: "cancelled", runId };
        } else {
          summaryError = err instanceof Error ? err.message : String(err);
        }
      }
      if (result === undefined) {
        result = { status: "done", runId, ...counts, samBytes: stat?.size ?? 0, summary, summaryError };
        finishedRuns.set(runId, finished);
      }
    } else {
      result = { status: "failed", runId, error: friendlyError(exit) };
    }

    if (result.status === "cancelled") await fs.rm(runDir, { recursive: true, force: true });
    activeRun = null;
    send("run:finished", result);
  });

  return { runId };
}

// ---- IPC -----------------------------------------------------------------------

/** Only our own page may call the API. */
function fromOurPage(event: IpcMainInvokeEvent) {
  return BrowserWindow.fromWebContents(event.sender) !== null &&
    event.senderFrame === event.sender.mainFrame && event.senderFrame.url === indexUrl;
}

function handle<T extends unknown[], R>(channel: string, fn: (event: IpcMainInvokeEvent, ...args: T) => R) {
  ipcMain.handle(channel, (event, ...args) => {
    if (!fromOurPage(event)) throw new Error("blocked: request from an unexpected page");
    return fn(event, ...(args as T));
  });
}

function registerIpc() {
  handle("presets:get", () => presets);

  handle("file:choose", async (event, kind: "reference" | "reads"): Promise<PickedFile | null> => {
    const window = BrowserWindow.fromWebContents(event.sender);
    const options = {
      title: kind === "reference" ? "Choose a reference genome" : "Choose a reads file",
      properties: ["openFile" as const],
      filters: [
        { name: "FASTA or FASTQ", extensions: ["fa", "fasta", "fna", "fq", "fastq", "gz"] },
        { name: "All files", extensions: ["*"] },
      ],
    };
    const picked = window ? await dialog.showOpenDialog(window, options) : await dialog.showOpenDialog(options);
    if (picked.canceled || picked.filePaths.length === 0) return null;
    const filePath = picked.filePaths[0];
    const stat = await fs.stat(filePath);
    const token = randomUUID();
    pickedFiles.set(token, filePath);
    return { token, name: path.basename(filePath), bytes: stat.size };
  });

  handle("run:start", (event, request: RunRequest) => {
    const window = BrowserWindow.fromWebContents(event.sender);
    if (!window) throw new Error("no window");
    return startRun(window, request);
  });

  handle("run:cancel", () => {
    if (activeRun) {
      activeRun.cancelled = true;
      activeRun.aligner.cancel();
    }
  });

  handle("view:region", (_event, runId: string, start: number, end: number) => {
    const run = finishedRuns.get(runId);
    if (!run?.analysis || !Number.isFinite(start) || !Number.isFinite(end)) return null;
    return readRegion(run.analysis, run.samPath, start, end);
  });

  const saveText = async (event: IpcMainInvokeEvent, name: string, kind: string, text: string) => {
    const window = BrowserWindow.fromWebContents(event.sender);
    if (!window) return null;
    const target = await dialog.showSaveDialog(window, {
      title: `Save ${kind}`,
      defaultPath: path.join(app.getPath("documents"), name),
      filters: [{ name: "CSV", extensions: ["csv"] }],
    });
    if (target.canceled || !target.filePath) return null;
    await fs.writeFile(target.filePath, text);
    return target.filePath;
  };

  handle("export:variants", (event, runId: string) => {
    const run = finishedRuns.get(runId);
    if (!run?.analysis) return null;
    return saveText(event, `${run.label}-variants.csv`, "variants", variantsCsv(run.variants, run.analysis.refName));
  });

  handle("export:coverage", (event, runId: string) => {
    const run = finishedRuns.get(runId);
    if (!run?.analysis) return null;
    return saveText(event, `${run.label}-coverage.csv`, "coverage", coverageCsv(run.analysis));
  });

  handle("run:save", async (event, runId: string) => {
    const run = finishedRuns.get(runId);
    const window = BrowserWindow.fromWebContents(event.sender);
    if (!run || !window) return null;
    const target = await dialog.showSaveDialog(window, {
      title: "Save alignments",
      defaultPath: path.join(app.getPath("documents"), `${run.label}.sam`),
      filters: [{ name: "SAM alignments", extensions: ["sam"] }],
    });
    if (target.canceled || !target.filePath) return null;
    await fs.copyFile(run.samPath, target.filePath);
    return target.filePath;
  });

  handle("run:reveal", (_event, runId: string) => {
    const run = finishedRuns.get(runId);
    if (run) shell.showItemInFolder(run.samPath);
  });
}

// ---- window and hardening ------------------------------------------------------

function createWindow() {
  const window = new BrowserWindow({
    width: 1100,
    height: 860,
    minWidth: 720,
    minHeight: 600,
    show: false,
    title: "Genome Variant Aligner",
    backgroundColor: nativeTheme.shouldUseDarkColors ? "#11181b" : "#f5f7f7",
    autoHideMenuBar: true,
    webPreferences: {
      preload: path.join(__dirname, "..", "preload", "preload.js"),
      contextIsolation: true,
      sandbox: true,
      nodeIntegration: false,
      webSecurity: true,
      spellcheck: false,
      devTools: !app.isPackaged,
    },
  });
  window.once("ready-to-show", () => window.show());
  void window.loadFile(indexHtml);
  return window;
}

function harden() {
  // no network: the app is fully offline, so refuse every http(s)/ws request outright
  const allowedResources = new Set(["index.html", "styles.css", "renderer.js", "viewer.js"]
    .map((name) => pathToFileURL(path.join(path.dirname(indexHtml), name)).href));
  session.defaultSession.webRequest.onBeforeRequest((details, callback) => {
    callback({ cancel: !allowedResources.has(details.url) });
  });
  // no camera, microphone, notifications, etc.
  session.defaultSession.setPermissionRequestHandler((_wc, _permission, callback) => callback(false));
  session.defaultSession.setPermissionCheckHandler(() => false);

  app.on("web-contents-created", (_event, contents) => {
    contents.on("will-navigate", (event) => event.preventDefault());
    contents.on("will-attach-webview", (event) => event.preventDefault());
    contents.setWindowOpenHandler(() => ({ action: "deny" }));
  });
}

if (!app.requestSingleInstanceLock()) {
  app.quit();
} else {
  app.on("second-instance", () => {
    const [window] = BrowserWindow.getAllWindows();
    if (window) {
      if (window.isMinimized()) window.restore();
      window.focus();
    }
  });

  void app.whenReady().then(() => {
    if (app.isPackaged) Menu.setApplicationMenu(null);
    harden();
    registerIpc();
    createWindow();
  });

  app.on("before-quit", () => {
    if (activeRun) {
      activeRun.cancelled = true;
      activeRun.aligner.cancel();
    }
  });
  app.on("window-all-closed", () => app.quit());
}
