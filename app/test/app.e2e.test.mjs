// Drives the real app window through its five pages: choose inputs, align the 10K preset, explore the
// result and its variants, then try an over-long read and cancel a 1M run.
// Screenshots go to test-results/ for a visual check.
import assert from "node:assert/strict";
import { existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { test } from "node:test";
import { fileURLToPath, pathToFileURL } from "node:url";
import { _electron as electron } from "playwright-core";

const appDir = fileURLToPath(new URL("..", import.meta.url));
const exe = path.join(appDir, "..", process.platform === "win32" ? "run.exe" : "run");
const shots = path.join(appDir, "test-results");
const PAGES = ["reference", "reads", "align", "explore", "variants"];

test("walk through the pages: align 10,000 reads, explore them, then cancel a 1,000,000-read run", { skip: !existsSync(exe) && "aligner not built", timeout: 120_000 }, async () => {
  mkdirSync(shots, { recursive: true });
  const userData = mkdtempSync(path.join(tmpdir(), "gva-e2e-"));
  const electronApp = await electron.launch({ args: ["."], cwd: appDir, env: { ...process.env, GVA_USER_DATA: userData } });
  try {
    const page = await electronApp.firstWindow();
    const shot = (name) => page.screenshot({ path: path.join(shots, `${name}.png`) });
    const shownPage = () => page.$$eval(".page", (els) => els.filter((el) => !el.hidden).map((el) => el.dataset.page));
    const expectPage = async (name) => {
      await page.waitForSelector(`.page[data-page="${name}"]:not([hidden])`);
      assert.deepEqual(await shownPage(), [name], `expected only the ${name} page`);
      assert.equal(await page.getAttribute(`.steps [data-step="${name}"]`, "aria-current"), "step");
    };
    const stepEnabled = async () => Object.fromEntries(await Promise.all(
      PAGES.map(async (name) => [name, await page.isEnabled(`.steps [data-step="${name}"]`)])));
    const next = () => page.click(".page:not([hidden]) [data-next]");
    const back = () => page.click(".page:not([hidden]) [data-back]");
    // each page's "?" opens its note and closes it again
    const checkHelp = async (id) => {
      const toggle = `.help-toggle[aria-controls="${id}"]`;
      assert.equal(await page.isVisible(`#${id}`), false, `${id} open by default`);
      await page.click(toggle);
      assert.equal(await page.isVisible(`#${id}`), true, `${id} didn't open`);
      assert.equal(await page.getAttribute(toggle, "aria-expanded"), "true");
      assert.ok((await page.textContent(`#${id}`)).trim().length > 80, `${id} is nearly empty`);
      await page.click(toggle);
      assert.equal(await page.isVisible(`#${id}`), false, `${id} didn't close`);
    };

    // ---- opens on the reference page, with later steps locked ----
    await page.waitForSelector("#reference-sars-cov-2");
    await expectPage("reference");
    assert.deepEqual(await stepEnabled(), { reference: true, reads: false, align: false, explore: false, variants: false });
    for (const id of ["progressPanel", "resultPanel", "errorPanel", "cancelledNote"]) {
      assert.equal(await page.isVisible(`#${id}`), false, `${id} visible before any run`);
    }

    // both bundled fonts load (the window has no network access, so they must come from the app)
    const loadedFonts = await page.evaluate(async () => {
      await document.fonts.ready;
      return [...document.fonts].filter((f) => f.status === "loaded").map((f) => f.family.replace(/"/g, ""));
    });
    assert.ok(loadedFonts.includes("Playfair Display Variable"), `fonts loaded: ${loadedFonts}`);
    assert.ok(loadedFonts.includes("Oi"), `fonts loaded: ${loadedFonts}`);

    await page.click('.help-toggle[aria-controls="help-ref"]');
    await shot("1-reference");
    await page.click('.help-toggle[aria-controls="help-ref"]');
    await checkHelp("help-ref");

    // ---- security ----
    // the window has no Node.js access
    assert.equal(await page.evaluate(() => typeof window.require), "undefined");
    assert.equal(await page.evaluate(() => typeof window.process), "undefined");
    // the request policy allows the app's own font files and nothing else outside the page's files
    const fontUrl = pathToFileURL(path.join(appDir, "build", "renderer", "fonts", "oi-latin-400-normal.woff2")).href;
    const requests = await electronApp.evaluate(async ({ session }, urls) => {
      const result = [];
      for (const url of urls) result.push(await session.defaultSession.fetch(url).then(() => "allowed", () => "blocked"));
      return result;
    }, ["file:///C:/Windows/win.ini", "https://example.com", "file://server/share/private.txt", fontUrl]);
    assert.deepEqual(requests, ["blocked", "blocked", "blocked", "allowed"]);
    await page.evaluate(() => { window.open("https://example.com"); });
    assert.equal(electronApp.windows().length, 1);
    const foreignIpc = await electronApp.evaluate(async ({ BrowserWindow }) => {
      const foreign = new BrowserWindow({ show: false, webPreferences: { nodeIntegration: true, contextIsolation: false, sandbox: false } });
      try {
        await foreign.loadURL("about:blank");
        return await foreign.webContents.executeJavaScript(`require('electron').ipcRenderer.invoke('presets:get').then(() => 'allowed', () => 'blocked')`);
      } finally {
        foreign.destroy();
      }
    });
    assert.equal(foreignIpc, "blocked");

    // ---- About box: version from the main process, licenses listed, closes again ----
    await page.click("#aboutButton");
    await page.waitForSelector("#aboutDialog[open]");
    assert.match(await page.textContent("#aboutVersion"), /^Version \d+\.\d+\.\d+/);
    for (const id of ["app", "fonts", "electron", "chromium"]) {
      assert.equal(await page.locator(`#aboutDialog [data-license="${id}"]`).count(), 1, `no ${id} license link`);
    }
    assert.match(await page.textContent("#aboutDialog"), /kseq\.h[\s\S]*zlib[\s\S]*Playfair Display[\s\S]*Electron[\s\S]*Chromium/);
    await shot("8-about");
    await page.click("#aboutClose");
    await page.waitForSelector("#aboutDialog:not([open])", { state: "attached" });
    // the files those links open exist when running from source
    const electronDist = path.join(appDir, "node_modules", "electron", "dist");
    for (const file of [path.join(appDir, "..", "LICENSE"), path.join(appDir, "build", "fonts-OFL.txt"),
      path.join(electronDist, "LICENSE"), path.join(electronDist, "LICENSES.chromium.html")]) {
      assert.ok(existsSync(file), `missing ${file}`);
    }

    // ---- reads, with Back and the step indicator ----
    await next();
    await expectPage("reads");
    assert.equal(await page.isChecked("#reads-10k"), true, "10K preset isn't the default");
    await checkHelp("help-reads");
    await shot("2-reads");
    await back();
    await expectPage("reference");
    await page.click('.steps [data-step="reads"]');   // reached once, so it stays open
    await expectPage("reads");
    await next();

    // ---- align: Continue stays off until a run finishes ----
    await expectPage("align");
    await checkHelp("help-align");
    assert.equal(await page.isEnabled("#alignContinue"), false);
    assert.deepEqual(await stepEnabled(), { reference: true, reads: true, align: true, explore: false, variants: false });
    assert.match(await page.textContent("#runSummary"), /10,000 reads against SARS-CoV-2/);
    await page.click("#runButton");
    await page.waitForSelector("#resultPanel:not([hidden])", { timeout: 60_000 });
    assert.equal(await page.textContent("#statReads"), "10,000");
    const mapped = Number((await page.textContent("#statMapped")).replace(/,/g, ""));
    const unmapped = Number((await page.textContent("#statUnmapped")).replace(/,/g, ""));
    assert.equal(mapped + unmapped, 10000);
    assert.ok(mapped >= 9900, `only ${mapped} mapped`);
    assert.equal(await page.isVisible("#resultWarning"), false, "warning shown for the 10K run");
    await page.waitForSelector("#alignContinue:not([disabled])");
    assert.equal((await stepEnabled()).explore, true);
    await shot("3-align");

    // ---- explore ----
    await next();
    await expectPage("explore");
    await checkHelp("help-viewer");
    assert.match(await page.textContent("#covMean"), /^\d[\d,]*×$/);
    const variantRows = await page.locator("#variantRows tr").count();
    assert.ok(variantRows > 0, "no variants listed");
    assert.equal(await page.textContent("#varCount"), String(variantRows));
    // the tracks are drawn once their page is visible
    assert.ok(await page.$eval("#overviewCanvas", (c) => c.width > 0 && c.height > 0), "overview not drawn");
    await shot("4-explore");
    // go-to box and zoom controls
    await page.fill("#gotoInput", "27,950-28,050");
    await page.press("#gotoInput", "Enter");
    await page.waitForFunction(() => document.getElementById("regionLabel").textContent.includes("27,950–28,050"));
    await page.click("#showAll");
    await page.waitForFunction(() => /Zoom in to/.test(document.getElementById("detailNote").textContent));

    // ---- variants: a row click goes back to Explore, centred on that site ----
    await next();
    await expectPage("variants");
    await checkHelp("help-variants");
    await shot("5-variants");
    await page.click("#variantRows tr:first-child");
    await expectPage("explore");
    await page.waitForFunction(() => /^Showing [\d,]+ of [\d,]+ reads/.test(document.getElementById("detailNote").textContent));
    assert.match(await page.textContent("#regionLabel"), /bp$/);
    await page.locator("#detailCanvas").scrollIntoViewIfNeeded();
    await shot("6-explore-reads");
    await page.click('.steps [data-step="variants"]');
    await expectPage("variants");

    // exports (the save dialog is replaced so the test doesn't block on it)
    const exportDir = mkdtempSync(path.join(tmpdir(), "gva-export-"));
    await electronApp.evaluate(({ dialog }, dir) => {
      let n = 0;
      dialog.showSaveDialog = async () => ({ canceled: false, filePath: `${dir}/export-${++n}.csv` });
    }, exportDir.replace(/\\/g, "/"));
    await page.click("#exportVariants");
    await page.waitForFunction(() => /Saved to/.test(document.getElementById("exportNote").textContent));
    const variantsCsv = readFileSync(path.join(exportDir, "export-1.csv"), "utf8").trim().split("\n");
    assert.equal(variantsCsv.length, variantRows + 1);
    await page.click("#exportCoverage");
    await page.waitForFunction(() => /export-2/.test(document.getElementById("exportNote").textContent));
    assert.equal(readFileSync(path.join(exportDir, "export-2.csv"), "utf8").trim().split("\n").length, 29882 + 1);
    rmSync(exportDir, { recursive: true, force: true });

    // Start over returns to the first page and keeps the result reachable
    await page.click("#startOver");
    await expectPage("reference");
    assert.equal((await stepEnabled()).variants, true);

    // ---- changing the reads clears the result; an over-long read gets a plain warning ----
    const longDir = mkdtempSync(path.join(tmpdir(), "gva-long-"));
    const genome = readFileSync(path.join(appDir, "..", "2019-nCoV.fa"), "utf8").split("\n").slice(1).join("").replace(/\s/g, "");
    const longReads = path.join(longDir, "long.fa");
    writeFileSync(longReads, `>ok\n${genome.slice(1000, 1100)}\n>long\n${genome.slice(2000, 2200)}\n`);
    await electronApp.evaluate(({ dialog }, file) => {
      dialog.showOpenDialog = async () => ({ canceled: false, filePaths: [file] });
    }, longReads);
    await page.click('.steps [data-step="reads"]');
    await expectPage("reads");
    await page.click("label:has(#reads-__file__)");
    await page.waitForFunction(() => /long\.fa/.test(document.getElementById("runSummary").textContent));
    assert.deepEqual(await stepEnabled(), { reference: true, reads: true, align: true, explore: false, variants: false });
    await next();
    await expectPage("align");
    assert.equal(await page.isVisible("#resultPanel"), false, "old result still shown");
    assert.equal(await page.isEnabled("#alignContinue"), false);
    await page.click("#runButton");
    await page.waitForSelector("#resultPanel:not([hidden])", { timeout: 30_000 });
    assert.equal(await page.isVisible("#resultWarning"), true);
    assert.match(await page.textContent("#resultWarning"), /longer than 189 bases[\s\S]*counted as unmapped/);
    assert.equal(await page.textContent("#statUnmapped"), "1");
    await shot("7-long-read-warning");
    rmSync(longDir, { recursive: true, force: true });

    // ---- a long run can be cancelled; moving between pages pauses while it runs ----
    await back();
    await expectPage("reads");
    await page.click("label:has(#reads-1m)");
    await next();
    await page.click("#runButton");
    await page.waitForFunction(() => /\d+%/.test(document.getElementById("stage").textContent), null, { timeout: 30_000 });
    assert.equal(await page.isEnabled(".page:not([hidden]) [data-back]"), false, "Back enabled during a run");
    assert.equal((await stepEnabled()).reads, false, "steps enabled during a run");
    await shot("9-running");
    await page.click("#cancelButton");
    await page.waitForSelector("#cancelledNote:not([hidden])", { timeout: 10_000 });
    assert.equal(await page.isEnabled("#runButton"), true);
    assert.equal(await page.isEnabled(".page:not([hidden]) [data-back]"), true);
    assert.equal(await page.isVisible("#resultPanel"), false);
    assert.equal(await page.isEnabled("#alignContinue"), false);
    assert.equal((await stepEnabled()).explore, false);

    // the window can't reach the network
    const fetched = await page.evaluate(() => fetch("https://example.com").then(() => "reached", () => "blocked"));
    assert.equal(fetched, "blocked");
    // Keep navigation probing last: Playwright can retain a pending navigation
    // after Electron cancels it and then wait forever before subsequent clicks.
    const navigation = electronApp.evaluate(({ BrowserWindow }) => new Promise((resolve) => {
      const contents = BrowserWindow.getAllWindows()[0].webContents;
      contents.once("will-navigate", (event) => resolve(event.defaultPrevented));
    }));
    await page.evaluate(() => {
      const link = document.createElement("a");
      link.href = "https://example.com";
      document.body.append(link);
      link.click();
    });
    assert.equal(await navigation, true);
    assert.match(page.url(), /\/renderer\/index\.html$/);
  } finally {
    await electronApp.close();
    rmSync(userData, { recursive: true, force: true });
  }
});
