// Copies the renderer's HTML and CSS next to its compiled script, bundles the font
// (the app is offline, so they can't come from a font service), and writes its license.
import { cpSync, mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";

const root = fileURLToPath(new URL("..", import.meta.url));
for (const file of ["index.html", "styles.css"]) {
  cpSync(`${root}src/renderer/${file}`, `${root}build/renderer/${file}`);
}

/** Font files the window loads; main.ts allows exactly these. */
const fonts = [
  ["@fontsource-variable/playfair-display", "playfair-display-latin-wght-normal.woff2"],
];
mkdirSync(`${root}build/renderer/fonts`, { recursive: true });
for (const [pkg, file] of fonts) {
  cpSync(`${root}node_modules/${pkg}/files/${file}`, `${root}build/renderer/fonts/${file}`);
}

// The font license, opened from the About box (and shipped with the installer).
const licenses = fonts.map(([pkg]) => readFileSync(`${root}node_modules/${pkg}/LICENSE`, "utf8").trim());
writeFileSync(`${root}build/fonts-OFL.txt`,
  `Font bundled with Genome Variant Aligner: Playfair Display.\n\n${licenses.join("\n\n" + "-".repeat(72) + "\n\n")}\n`);
