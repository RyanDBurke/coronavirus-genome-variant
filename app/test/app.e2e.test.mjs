// Drives the real app window: align the 10K preset, check the result, then cancel a 1M run.
// Screenshots go to test-results/ for a visual check.
import assert from "node:assert/strict";
import { existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { test } from "node:test";
import { fileURLToPath } from "node:url";
import { _electron as electron } from "playwright-core";

const appDir = fileURLToPath(new URL("..", import.meta.url));
const exe = path.join(appDir, "..", process.platform === "win32" ? "run.exe" : "run");
const shots = path.join(appDir, "test-results");

test("align 10,000 reads, then cancel a 1,000,000-read run", { skip: !existsSync(exe) && "aligner not built", timeout: 120_000 }, async () => {
  mkdirSync(shots, { recursive: true });
  const userData = mkdtempSync(path.join(tmpdir(), "gva-e2e-"));
  const electronApp = await electron.launch({ args: ["."], cwd: appDir, env: { ...process.env, GVA_USER_DATA: userData } });
  try {
    const page = await electronApp.firstWindow();
    await page.waitForSelector("#reads-10k");
    await page.screenshot({ path: path.join(shots, "1-start.png") });

    // only the inputs show before a run
    for (const id of ["progressPanel", "resultPanel", "errorPanel", "cancelledNote"]) {
      assert.equal(await page.isVisible(`#${id}`), false, `${id} visible before any run`);
    }

    // the window has no Node.js access
    assert.equal(await page.evaluate(() => typeof window.require), "undefined");
    assert.equal(await page.evaluate(() => typeof window.process), "undefined");

    // Exercise the request policy without relying on the page's CSP.
    const requests = await electronApp.evaluate(async ({ session }) => {
      const result = [];
      for (const url of ['file:///C:/Windows/win.ini', 'https://example.com', 'file://server/share/private.txt']) {
        result.push(await session.defaultSession.fetch(url).then(() => 'allowed', () => 'blocked'));
      }
      return result;
    });
    assert.deepEqual(requests, ['blocked', 'blocked', 'blocked']);
    await page.evaluate(() => { window.open('https://example.com'); });
    assert.equal(electronApp.windows().length, 1);
    const foreignIpc = await electronApp.evaluate(async ({ BrowserWindow }) => {
      const foreign = new BrowserWindow({ show: false, webPreferences: { nodeIntegration: true, contextIsolation: false, sandbox: false } });
      try {
        await foreign.loadURL('about:blank');
        return await foreign.webContents.executeJavaScript(`require('electron').ipcRenderer.invoke('presets:get').then(() => 'allowed', () => 'blocked')`);
      } finally {
        foreign.destroy();
      }
    });
    assert.equal(foreignIpc, 'blocked');

    // 10K preset is selected by default
    assert.equal(await page.isChecked("#reads-10k"), true);
    await page.click("#runButton");
    await page.waitForSelector("#resultPanel:not([hidden])", { timeout: 60_000 });
    assert.equal(await page.textContent("#statReads"), "10,000");
    const mapped = Number((await page.textContent("#statMapped")).replace(/,/g, ""));
    const unmapped = Number((await page.textContent("#statUnmapped")).replace(/,/g, ""));
    assert.equal(mapped + unmapped, 10000);
    assert.ok(mapped >= 9900, `only ${mapped} mapped`);
    await page.screenshot({ path: path.join(shots, "2-result.png") });

    // results viewer: stats, variants table, jumping to a variant shows reads
    await page.waitForSelector("#viewerSection:not([hidden])");
    assert.match(await page.textContent("#covMean"), /^\d[\d,]*×$/);
    const variantRows = await page.locator("#variantRows tr").count();
    assert.ok(variantRows > 0, "no variants listed");
    assert.equal(await page.textContent("#varCount"), String(variantRows));
    await page.locator("#viewerSection").scrollIntoViewIfNeeded();
    await page.screenshot({ path: path.join(shots, "3-viewer-genome.png") });
    await page.click("#variantRows tr:first-child");
    await page.waitForFunction(() => /^Showing [\d,]+ of [\d,]+ reads/.test(document.getElementById("detailNote").textContent));
    assert.match(await page.textContent("#regionLabel"), /bp$/);
    await page.locator("#detailCanvas").scrollIntoViewIfNeeded();
    await page.screenshot({ path: path.join(shots, "4-viewer-variant.png") });

    // go-to box and zoom controls
    await page.fill("#gotoInput", "27,950-28,050");
    await page.press("#gotoInput", "Enter");
    await page.waitForFunction(() => document.getElementById("regionLabel").textContent.includes("27,950–28,050"));
    await page.click("#showAll");
    await page.waitForFunction(() => /Zoom in to/.test(document.getElementById("detailNote").textContent));

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

    // a long run can be cancelled
    await page.click("label:has(#reads-1m)");
    await page.click("#runButton");
    await page.waitForFunction(() => /\d+%/.test(document.getElementById("stage").textContent), null, { timeout: 30_000 });
    await page.screenshot({ path: path.join(shots, "5-running.png") });
    await page.click("#cancelButton");
    await page.waitForSelector("#cancelledNote:not([hidden])", { timeout: 10_000 });
    assert.equal(await page.isEnabled("#runButton"), true);
    assert.equal(await page.isVisible("#resultPanel"), false);
    assert.equal(await page.isVisible("#viewerSection"), false);

    // the window can't reach the network
    const fetched = await page.evaluate(() => fetch("https://example.com").then(() => "reached", () => "blocked"));
    assert.equal(fetched, "blocked");
    // Keep navigation probing last: Playwright can retain a pending navigation
    // after Electron cancels it and then wait forever before subsequent clicks.
    const navigation = electronApp.evaluate(({ BrowserWindow }) => new Promise((resolve) => {
      const contents = BrowserWindow.getAllWindows()[0].webContents;
      contents.once('will-navigate', (event) => resolve(event.defaultPrevented));
    }));
    await page.evaluate(() => {
      const link = document.createElement('a');
      link.href = 'https://example.com';
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
