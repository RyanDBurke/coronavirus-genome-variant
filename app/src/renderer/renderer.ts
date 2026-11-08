// The window's behaviour: pick inputs, start a run, show progress, then the result or a clear error.
// Each step is its own page (pages.ts); a finished run unlocks the Explore and Variants pages.

import type { AlignerApi, FriendlyError, InputChoice, LicenseId, PickedFile, Preset, RunEvent, RunResult } from "../shared/types.js";
import { Pages } from "./pages.js";
import { Viewer } from "./viewer.js";

declare global {
  interface Window {
    aligner: AlignerApi;
  }
}

const api = window.aligner;
const pages = new Pages();
const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;

const FILE_CHOICE = "__file__";

type Kind = "reference" | "reads";
interface InputState {
  selected: string;      // preset id or FILE_CHOICE
  file: PickedFile | null;
}
const inputs: Record<Kind, InputState> = {
  reference: { selected: "sars-cov-2", file: null },
  reads: { selected: "10k", file: null },
};
const presetLabels: Record<Kind, Map<string, string>> = { reference: new Map(), reads: new Map() };

let running = false;
let currentRunId: string | null = null;
let viewer: Viewer | null = null;
let timer: number | undefined;
let startedAt = 0;

// ---- formatting -------------------------------------------------------------

const count = (n: number) => n.toLocaleString("en-US");

function bytes(n: number) {
  if (n >= 1e9) return `${(n / 1e9).toFixed(1)} GB`;
  if (n >= 1e6) return `${(n / 1e6).toFixed(n >= 1e8 ? 0 : 1)} MB`;
  if (n >= 1e3) return `${Math.round(n / 1e3)} KB`;
  return `${n} bytes`;
}

function duration(seconds: number) {
  if (seconds < 60) return `${seconds.toFixed(1)} s`;
  const m = Math.floor(seconds / 60);
  return `${m} min ${Math.round(seconds - m * 60)} s`;
}

// ---- input choices -----------------------------------------------------------

function choiceCard(kind: Kind, value: string, title: string, detail: string, extraClass = "") {
  const label = document.createElement("label");
  label.className = `choice ${extraClass}`.trim();
  const input = document.createElement("input");
  input.type = "radio";
  input.name = kind;
  input.value = value;
  input.id = `${kind}-${value}`;
  input.checked = inputs[kind].selected === value;
  const strong = document.createElement("strong");
  strong.textContent = title;
  const small = document.createElement("small");
  small.textContent = detail;
  label.append(input, strong, small);
  return { label, input, small };
}

function renderChoices(kind: Kind, presets: Preset[]) {
  const container = $(kind === "reference" ? "referenceChoices" : "readsChoices");
  container.replaceChildren();

  for (const preset of presets) {
    presetLabels[kind].set(preset.id, preset.label);
    const { label, input } = choiceCard(kind, preset.id, preset.label, preset.detail);
    input.addEventListener("change", () => {
      inputs[kind].selected = preset.id;
      updateSummary();
      inputsChanged();
    });
    container.append(label);
  }

  const picked = inputs[kind].file;
  const { label, input, small } = choiceCard(
    kind,
    FILE_CHOICE,
    kind === "reference" ? "Your own reference…" : "Your own reads…",
    picked ? `${picked.name} · ${bytes(picked.bytes)}` : "FASTA or FASTQ, optionally .gz",
    "file",
  );
  if (picked) small.classList.add("picked");

  // picking happens on click, so choosing the same card again lets you pick a different file
  input.addEventListener("click", async () => {
    const previous = inputs[kind].selected;
    const file = await api.chooseFile(kind);
    if (file) {
      inputs[kind] = { selected: FILE_CHOICE, file };
    } else if (!inputs[kind].file) {
      inputs[kind].selected = previous === FILE_CHOICE ? presets[0].id : previous;
    } else {
      inputs[kind].selected = FILE_CHOICE;
    }
    renderChoices(kind, presets);
    updateSummary();
    if (file) inputsChanged();
  });
  container.append(label);
}

function choiceFor(kind: Kind): InputChoice | null {
  const state = inputs[kind];
  if (state.selected === FILE_CHOICE) {
    return state.file ? { kind: "file", token: state.file.token } : null;
  }
  return { kind: "preset", id: state.selected };
}

function describe(kind: Kind) {
  const state = inputs[kind];
  return state.selected === FILE_CHOICE ? (state.file?.name ?? "no file") : (presetLabels[kind].get(state.selected) ?? "");
}

function updateSummary() {
  $("runSummary").textContent = `${describe("reads")} against ${describe("reference")}`;
}

/** A result belongs to the inputs it was made from: changing them clears it and locks the results pages. */
function inputsChanged() {
  if (running) return;
  closeViewer();
  show(null);
  currentRunId = null;
}

function setInputsEnabled(enabled: boolean) {
  document.querySelectorAll<HTMLInputElement>(".choices input").forEach((input) => (input.disabled = !enabled));
}

// ---- running -----------------------------------------------------------------

function show(panel: string | null) {
  for (const id of ["progressPanel", "resultPanel", "errorPanel", "cancelledNote"]) {
    $(id).hidden = id !== panel;
  }
}

let shownPercent: number | null = null;

/** percent null shows the indeterminate animation. The bar only animates forward: leaving the
 *  indeterminate state or starting a new stage at a lower percent jumps straight there. */
function setProgress(percent: number | null, stage: string) {
  const bar = $("progressBar");
  const fill = $("progressFill");
  $("stage").textContent = stage;
  const snap = percent === null || shownPercent === null || percent < shownPercent;
  if (snap) fill.style.transition = "none";
  if (percent === null) {
    bar.classList.add("indeterminate");
    bar.removeAttribute("aria-valuenow");
    fill.style.width = "";
  } else {
    bar.classList.remove("indeterminate");
    bar.setAttribute("aria-valuenow", String(percent));
    fill.style.width = `${percent}%`;
  }
  if (snap) {
    void fill.offsetWidth; // apply the jump before turning the transition back on
    fill.style.transition = "";
  }
  shownPercent = percent;
}

// ---- warnings ------------------------------------------------------------------

let runWarnings: string[] = [];

/** The aligner's warnings in plain words. */
function plainWarning(message: string) {
  const long = /reads longer than (\d+) bp are reported as unmapped/.exec(message);
  if (long) {
    return `Some reads are longer than ${long[1]} bases. The aligner can't place reads that long, so they are counted as unmapped.`;
  }
  return message.charAt(0).toUpperCase() + message.slice(1) + ".";
}

function showWarnings() {
  const box = $("resultWarning");
  box.replaceChildren();
  const unique = [...new Set(runWarnings)];
  if (unique.length > 0) {
    const title = document.createElement("strong");
    title.textContent = unique.length === 1 ? "Warning" : "Warnings";
    box.append(title);
    for (const message of unique) {
      const line = document.createElement("p");
      line.textContent = plainWarning(message);
      box.append(line);
    }
  }
  box.hidden = unique.length === 0;
}

function setRunning(isRunning: boolean) {
  running = isRunning;
  $<HTMLButtonElement>("runButton").disabled = isRunning;
  $("cancelButton").hidden = !isRunning;
  setInputsEnabled(!isRunning);
  pages.setBusy(isRunning);
  window.clearInterval(timer);
  if (isRunning) {
    startedAt = performance.now();
    $("elapsed").textContent = "0.0 s";
    timer = window.setInterval(() => {
      $("elapsed").textContent = duration((performance.now() - startedAt) / 1000);
    }, 100);
  }
}

function showError(error: FriendlyError) {
  $("errorTitle").textContent = error.title;
  $("errorHint").textContent = error.hint;
  $("errorDetail").textContent = error.detail || "No further details.";
  $<HTMLDetailsElement>("errorPanel").querySelector("details")!.hidden = !error.detail;
  show("errorPanel");
  $("errorPanel").scrollIntoView({ block: "nearest", behavior: "smooth" });
}

async function start() {
  const reference = choiceFor("reference");
  const reads = choiceFor("reads");
  if (!reference || !reads) {
    showError({ title: "Choose a file first", detail: "", hint: "Pick a file for each step, or use a built-in option." });
    return;
  }
  setRunning(true);
  closeViewer();
  runWarnings = [];
  setProgress(null, "Building the genome index…");
  show("progressPanel");
  $("progressPanel").scrollIntoView({ block: "nearest" });
  $("savedNote").textContent = "";

  const started = await api.start({ reference, reads });
  if ("error" in started) {
    setRunning(false);
    showError(started.error);
    return;
  }
  currentRunId = started.runId;
}

function closeViewer() {
  viewer?.destroy();
  viewer = null;
  $<HTMLButtonElement>("alignContinue").disabled = true;
  pages.lockAfter("align");
  $("viewerUnavailable").hidden = true;
  $("exportNote").textContent = "";
}

function onEvent(event: RunEvent) {
  if (!running) return;
  if (event.event === "index") {
    setProgress(0, "Aligning reads…");
  } else if (event.event === "progress") {
    setProgress(event.percent, `Aligning reads… ${event.percent}%`);
  } else if (event.event === "analyze") {
    setProgress(event.percent, `Preparing the results viewer… ${event.percent}%`);
  } else if (event.event === "warning") {
    runWarnings.push(event.message);
  }
}

function onFinished(result: RunResult) {
  if (result.runId !== currentRunId) return;
  setRunning(false);

  if (result.status === "cancelled") {
    show("cancelledNote");
    return;
  }
  if (result.status === "failed") {
    showError(result.error);
    return;
  }

  $("statReads").textContent = count(result.reads);
  $("statMapped").textContent = count(result.mapped);
  $("statUnmapped").textContent = count(result.unmapped);
  const share = result.reads > 0 ? (100 * result.mapped) / result.reads : 0;
  $("mapbarFill").style.width = `${share}%`;
  $("resultTime").textContent = `${share.toFixed(1)}% mapped · ${duration(result.seconds)}`;

  showWarnings();
  $("resultNote").textContent = `The SAM file is ${bytes(result.samBytes)}. The app keeps your last 3 runs; save a copy to keep this one.`;
  show("resultPanel");
  $("resultPanel").scrollIntoView({ block: "nearest", behavior: "smooth" });

  if (result.summary) {
    viewer = new Viewer(api, result.runId, result.summary, () => pages.goTo("explore"));
    $<HTMLButtonElement>("alignContinue").disabled = false;
    pages.unlock("explore");
  } else {
    const unavailable = $("viewerUnavailable");
    unavailable.textContent = `The results viewer isn't available for this run: ${result.summaryError ?? "unknown error"}. The SAM file is still complete.`;
    unavailable.hidden = false;
  }
}

// ---- section help --------------------------------------------------------------

/** Each "?" button shows or hides the note it controls (aria-controls). */
function initHelp() {
  document.querySelectorAll<HTMLButtonElement>(".help-toggle").forEach((button) => {
    const note = $(button.getAttribute("aria-controls")!);
    button.addEventListener("click", () => {
      const open = button.getAttribute("aria-expanded") !== "true";
      button.setAttribute("aria-expanded", String(open));
      note.hidden = !open;
    });
  });
}

// ---- About box -----------------------------------------------------------------

async function initAbout() {
  const dialog = $<HTMLDialogElement>("aboutDialog");
  $("aboutButton").addEventListener("click", () => {
    $("aboutError").textContent = "";
    dialog.showModal();
  });
  $("aboutClose").addEventListener("click", () => dialog.close());
  // clicking the backdrop closes it too
  dialog.addEventListener("click", (event) => {
    if (event.target === dialog) dialog.close();
  });
  dialog.querySelectorAll<HTMLButtonElement>("[data-license]").forEach((button) => {
    button.addEventListener("click", async () => {
      $("aboutError").textContent = "";
      try {
        await api.openLicense(button.dataset.license as LicenseId);
      } catch {
        $("aboutError").textContent = "That license file couldn't be opened.";
      }
    });
  });
  const info = await api.getAppInfo();
  $("aboutVersion").textContent = `Version ${info.version}`;
}

// ---- wiring ------------------------------------------------------------------

async function init() {
  initHelp();
  void initAbout();
  const presets = await api.getPresets();
  renderChoices("reference", presets.references);
  renderChoices("reads", presets.reads);
  updateSummary();

  api.onEvent(onEvent);
  api.onFinished(onFinished);

  $("runButton").addEventListener("click", () => void start());
  $("startOver").addEventListener("click", () => pages.goTo("reference"));
  $("cancelButton").addEventListener("click", () => {
    setProgress(null, "Cancelling…");
    void api.cancel();
  });
  $("saveButton").addEventListener("click", async () => {
    if (!currentRunId) return;
    const saved = await api.saveSam(currentRunId);
    if (saved) $("savedNote").textContent = `Saved to ${saved}`;
  });
  $("revealButton").addEventListener("click", () => {
    if (currentRunId) void api.reveal(currentRunId);
  });
  $("exportVariants").addEventListener("click", async () => {
    if (!currentRunId) return;
    const saved = await api.exportVariants(currentRunId);
    if (saved) $("exportNote").textContent = `Saved to ${saved}`;
  });
  $("exportCoverage").addEventListener("click", async () => {
    if (!currentRunId) return;
    const saved = await api.exportCoverage(currentRunId);
    if (saved) $("exportNote").textContent = `Saved to ${saved}`;
  });
}

void init();
