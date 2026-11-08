// Packaging bundles ../run.exe; fail early with a useful message if it hasn't been built.
import { existsSync } from "node:fs";
import { fileURLToPath } from "node:url";

const exe = fileURLToPath(new URL("../../run.exe", import.meta.url));
if (!existsSync(exe)) {
  console.error(`run.exe not found at ${exe}.\nBuild it first: run \`make\` in the repo root from the MSYS2 UCRT64 terminal.`);
  process.exit(1);
}
