// The window's behaviour: pick inputs, start a run, show progress, then the result or a clear error.

import type { AlignerApi, AlignerEvent, FriendlyError, InputChoice, PickedFile, Preset, RunResult } from "../shared/types.js";

declare global {
  interface Window {
    aligner: AlignerApi;
  }
}

const api = window.aligner;
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

function setInputsEnabled(enabled: boolean) {
  document.querySelectorAll<HTMLInputElement>(".choices input").forEach((input) => (input.disabled = !enabled));
}

// ---- running -----------------------------------------------------------------

function show(panel: string | null) {
  for (const id of ["progressPanel", "resultPanel", "errorPanel", "cancelledNote"]) {
    $(id).hidden = id !== panel;
  }
}

function setProgress(percent: number | null, stage: string) {
  const bar = $("progressBar");
  $("stage").textContent = stage;
  if (percent === null) {
    bar.classList.add("indeterminate");
    bar.removeAttribute("aria-valuenow");
    $("progressFill").style.width = "";
  } else {
    bar.classList.remove("indeterminate");
    bar.setAttribute("aria-valuenow", String(percent));
    $("progressFill").style.width = `${percent}%`;
  }
}

function setRunning(isRunning: boolean) {
  running = isRunning;
  $<HTMLButtonElement>("runButton").disabled = isRunning;
  $("cancelButton").hidden = !isRunning;
  setInputsEnabled(!isRunning);
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
  // reset the bar without animating down from the previous run
  const fill = $("progressFill");
  fill.style.transition = "none";
  setProgress(0, "");
  void fill.offsetWidth;
  fill.style.transition = "";
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

function onEvent(event: AlignerEvent) {
  if (!running) return;
  if (event.event === "index") {
    setProgress(0, "Aligning reads…");
  } else if (event.event === "progress") {
    setProgress(event.percent, `Aligning reads… ${event.percent}%`);
  } else if (event.event === "warning") {
    $("resultNote").dataset.warning = event.message;
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

  const note = $("resultNote");
  const warning = note.dataset.warning ? ` Note: ${note.dataset.warning}.` : "";
  delete note.dataset.warning;
  note.textContent = `The SAM file is ${bytes(result.samBytes)}. The app keeps your last 3 runs; save a copy to keep this one.${warning}`;
  show("resultPanel");
  $("resultPanel").scrollIntoView({ block: "nearest", behavior: "smooth" });
}

// ---- wiring ------------------------------------------------------------------

async function init() {
  const presets = await api.getPresets();
  renderChoices("reference", presets.references);
  renderChoices("reads", presets.reads);
  updateSummary();

  api.onEvent(onEvent);
  api.onFinished(onFinished);

  $("runButton").addEventListener("click", () => void start());
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
}

void init();
