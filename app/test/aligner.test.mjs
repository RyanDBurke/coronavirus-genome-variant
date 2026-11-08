// Tests for the aligner wrapper and error messages, against the real aligner built in the repo root.
import assert from "node:assert/strict";
import { existsSync, mkdtempSync, rmSync, statSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { test } from "node:test";
import { fileURLToPath } from "node:url";
import { createRequire } from "node:module";

const require = createRequire(import.meta.url);
const { LineSplitter, parseEvent, runAligner } = require("../build/main/aligner.js");
const { friendlyError } = require("../build/main/messages.js");

const repo = fileURLToPath(new URL("../../", import.meta.url));
const exe = path.join(repo, process.platform === "win32" ? "run.exe" : "run");
const ref = path.join(repo, "2019-nCoV.fa");
const reads1k = path.join(repo, "reads", "reads_1K.fa.gz");
const reads1m = path.join(repo, "reads", "reads_1M.fa.gz");
const haveExe = existsSync(exe);

test("LineSplitter joins lines split across chunks and handles CRLF", () => {
  const s = new LineSplitter();
  assert.deepEqual(s.push('{"event":"pro'), []);
  assert.deepEqual(s.push('gress","percent":1}\r\n{"a":'), ['{"event":"progress","percent":1}']);
  assert.deepEqual(s.push("2}\n\n"), ['{"a":2}']);
  assert.deepEqual(s.push("tail"), []);
  assert.deepEqual(s.flush(), ["tail"]);
});

test("parseEvent accepts events and rejects everything else", () => {
  assert.deepEqual(parseEvent('{"event":"progress","percent":5}'), { event: "progress", percent: 5 });
  assert.equal(parseEvent("Segmentation fault"), null);
  assert.equal(parseEvent('{"percent":5}'), null);
  assert.equal(parseEvent("[1,2]"), null);
  assert.equal(parseEvent("null"), null);
});

test("friendlyError maps exit codes to actionable messages", () => {
  const base = { cancelled: false, otherOutput: [] };
  assert.match(friendlyError({ ...base, code: 3, error: { code: "input", message: "no reads" } }).title, /input file/);
  assert.equal(friendlyError({ ...base, code: 3, error: { code: "input", message: "no reads" } }).detail, "no reads");
  assert.match(friendlyError({ ...base, code: 4 }).title, /save the results/);
  assert.match(friendlyError({ ...base, code: 5 }).title, /memory/);
  assert.match(friendlyError({ ...base, code: null, spawnError: "spawn ENOENT" }).title, /couldn't start/);
  const crash = friendlyError({ ...base, code: 3221225477, otherOutput: ["boom"] });
  assert.match(crash.title, /stopped unexpectedly/);
  assert.equal(crash.detail, "boom");
});

function scratch() {
  return mkdtempSync(path.join(tmpdir(), "gva-test-"));
}

test("runs the aligner on 1,000 reads and reports progress and counts", { skip: !haveExe && "aligner not built" }, async () => {
  const dir = scratch();
  try {
    const events = [];
    const run = runAligner({ exe, reference: ref, reads: reads1k, samOut: path.join(dir, "out.sam"),
      indexOut: path.join(dir, "fm.txt"), cwd: dir, onEvent: (e) => events.push(e) });
    const exit = await run.done;
    assert.equal(exit.code, 0, exit.otherOutput.join("\n"));
    assert.equal(events[0].event, "index");
    assert.equal(events.at(-1).event, "done");
    assert.equal(events.at(-1).reads, 1000);
    assert.equal(events.at(-1).mapped + events.at(-1).unmapped, 1000);
    assert.ok(statSync(path.join(dir, "out.sam")).size > 100_000);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

test("reports a missing reads file as an input error", { skip: !haveExe && "aligner not built" }, async () => {
  const dir = scratch();
  try {
    const run = runAligner({ exe, reference: ref, reads: path.join(dir, "missing.fa"), samOut: path.join(dir, "out.sam"),
      indexOut: path.join(dir, "fm.txt"), cwd: dir, onEvent: () => {} });
    const exit = await run.done;
    assert.equal(exit.code, 3);
    assert.equal(exit.error.code, "input");
    assert.match(friendlyError(exit).title, /input file/);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

test("cancelling stops a long run quickly", { skip: !haveExe && "aligner not built" }, async () => {
  const dir = scratch();
  try {
    let run;
    const firstProgress = new Promise((resolve) => {
      run = runAligner({ exe, reference: ref, reads: reads1m, samOut: path.join(dir, "out.sam"),
        indexOut: path.join(dir, "fm.txt"), cwd: dir, onEvent: (e) => e.event === "progress" && resolve() });
    });
    await firstProgress;
    const t0 = Date.now();
    run.cancel();
    const exit = await run.done;
    assert.equal(exit.cancelled, true);
    assert.ok(Date.now() - t0 < 5000, "took too long to stop");
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

test("a missing aligner executable is reported, not thrown", async () => {
  const dir = scratch();
  try {
    const run = runAligner({ exe: path.join(dir, "nope.exe"), reference: ref, reads: reads1k, samOut: path.join(dir, "o.sam"),
      indexOut: path.join(dir, "fm.txt"), cwd: dir, onEvent: () => {} });
    const exit = await run.done;
    assert.ok(exit.spawnError);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});
