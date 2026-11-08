// Copies the renderer's HTML and CSS next to its compiled script.
import { cpSync } from "node:fs";
import { fileURLToPath } from "node:url";

const root = fileURLToPath(new URL("..", import.meta.url));
for (const file of ["index.html", "styles.css"]) {
  cpSync(`${root}src/renderer/${file}`, `${root}build/renderer/${file}`);
}
