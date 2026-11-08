// Draws the app icon and writes build-resources/icon.ico (16-256 px) and build-resources/icon.png (256 px).
// The design is a stylised read pileup: reads stacked over the genome, sharing one variant base.
// Colours follow the app (and rburke.dev): near-black tile, cream reads, red variant column.
// No dependencies: shapes are rasterised with 4x4 supersampling, then encoded as PNG and BMP by hand.
// Run with `node scripts/make-icon.mjs` from app/ after changing the design; the outputs are committed.

import { mkdirSync, writeFileSync } from "node:fs";
import { deflateSync } from "node:zlib";
import { fileURLToPath } from "node:url";

const outDir = fileURLToPath(new URL("../build-resources/", import.meta.url));

const TILE_TOP = [0x1c, 0x1b, 0x18];
const TILE_BOTTOM = [0x0c, 0x0c, 0x0c];
/** a thin edge so the dark tile still stands out on a dark taskbar */
const EDGE = [0x6e, 0x6c, 0x5c];
const READ = [0xfd, 0xfb, 0xd4];
const VARIANT = [0xe3, 0x17, 0x17];

// Shapes in unit coordinates (0..1). Small sizes use three thicker reads so they stay legible.
const layouts = {
  large: {
    barHeight: 0.11,
    rows: [[0.29, 0.16, 0.66], [0.43, 0.30, 0.84], [0.57, 0.16, 0.76], [0.71, 0.26, 0.84]],
    variant: [0.50, 0.60],
  },
  small: {
    barHeight: 0.16,
    rows: [[0.30, 0.14, 0.68], [0.50, 0.30, 0.86], [0.70, 0.18, 0.80]],
    variant: [0.47, 0.62],
  },
};

function inRoundedRect(x, y, x0, y0, x1, y1, r) {
  if (x < x0 || x > x1 || y < y0 || y > y1) return false;
  const cx = Math.min(Math.max(x, x0 + r), x1 - r);
  const cy = Math.min(Math.max(y, y0 + r), y1 - r);
  return (x - cx) ** 2 + (y - cy) ** 2 <= r * r;
}

/** Colour at a point, or null outside the icon. */
function colourAt(x, y, layout) {
  if (!inRoundedRect(x, y, 0.03, 0.03, 0.97, 0.97, 0.19)) return null;
  const h = layout.barHeight;
  for (const [cy, x0, x1] of layout.rows) {
    if (inRoundedRect(x, y, x0, cy - h / 2, x1, cy + h / 2, h / 2)) {
      return x >= layout.variant[0] && x <= layout.variant[1] ? VARIANT : READ;
    }
  }
  if (!inRoundedRect(x, y, 0.045, 0.045, 0.955, 0.955, 0.175)) return EDGE;
  const t = (y - 0.03) / 0.94;
  return TILE_TOP.map((c, i) => c + (TILE_BOTTOM[i] - c) * t);
}

/** RGBA pixels, top row first. */
function render(size) {
  const layout = size <= 32 ? layouts.small : layouts.large;
  const SS = 4;
  const rgba = Buffer.alloc(size * size * 4);
  for (let py = 0; py < size; py++) {
    for (let px = 0; px < size; px++) {
      let r = 0, g = 0, b = 0, a = 0;
      for (let sy = 0; sy < SS; sy++) {
        for (let sx = 0; sx < SS; sx++) {
          const c = colourAt((px + (sx + 0.5) / SS) / size, (py + (sy + 0.5) / SS) / size, layout);
          if (c) { r += c[0]; g += c[1]; b += c[2]; a++; }
        }
      }
      const o = (py * size + px) * 4;
      if (a > 0) {
        rgba[o] = Math.round(r / a);
        rgba[o + 1] = Math.round(g / a);
        rgba[o + 2] = Math.round(b / a);
        rgba[o + 3] = Math.round((255 * a) / (SS * SS));
      }
    }
  }
  return rgba;
}

// ---- PNG ----------------------------------------------------------------------

const CRC_TABLE = Array.from({ length: 256 }, (_, n) => {
  let c = n;
  for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
  return c >>> 0;
});
function crc32(buf) {
  let c = 0xffffffff;
  for (const byte of buf) c = CRC_TABLE[(c ^ byte) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}
function chunk(type, data) {
  const body = Buffer.concat([Buffer.from(type, "ascii"), data]);
  const out = Buffer.alloc(body.length + 8);
  out.writeUInt32BE(data.length, 0);
  body.copy(out, 4);
  out.writeUInt32BE(crc32(body), body.length + 4);
  return out;
}
function png(size, rgba) {
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(size, 0);
  ihdr.writeUInt32BE(size, 4);
  ihdr[8] = 8;  // bit depth
  ihdr[9] = 6;  // RGBA
  const raw = Buffer.alloc(size * (size * 4 + 1));
  for (let y = 0; y < size; y++) rgba.copy(raw, y * (size * 4 + 1) + 1, y * size * 4, (y + 1) * size * 4);
  return Buffer.concat([
    Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]),
    chunk("IHDR", ihdr),
    chunk("IDAT", deflateSync(raw, { level: 9 })),
    chunk("IEND", Buffer.alloc(0)),
  ]);
}

// ---- ICO ----------------------------------------------------------------------

/** 32-bit BMP entry (BGRA, bottom row first, plus an empty AND mask), the most compatible ICO format. */
function bmpEntry(size, rgba) {
  const header = Buffer.alloc(40);
  header.writeUInt32LE(40, 0);
  header.writeInt32LE(size, 4);
  header.writeInt32LE(size * 2, 8); // colour rows + mask rows
  header.writeUInt16LE(1, 12);
  header.writeUInt16LE(32, 14);
  const pixels = Buffer.alloc(size * size * 4);
  for (let y = 0; y < size; y++) {
    for (let x = 0; x < size; x++) {
      const s = (y * size + x) * 4;
      const d = ((size - 1 - y) * size + x) * 4;
      pixels[d] = rgba[s + 2];
      pixels[d + 1] = rgba[s + 1];
      pixels[d + 2] = rgba[s];
      pixels[d + 3] = rgba[s + 3];
    }
  }
  const maskRow = Math.ceil(size / 32) * 4;
  return Buffer.concat([header, pixels, Buffer.alloc(maskRow * size)]);
}

function ico(images) {
  const header = Buffer.alloc(6 + 16 * images.length);
  header.writeUInt16LE(1, 2);
  header.writeUInt16LE(images.length, 4);
  let offset = header.length;
  images.forEach(({ size, data }, i) => {
    const e = 6 + 16 * i;
    header[e] = size >= 256 ? 0 : size;
    header[e + 1] = size >= 256 ? 0 : size;
    header.writeUInt16LE(1, e + 4);   // planes
    header.writeUInt16LE(32, e + 6);  // bits per pixel
    header.writeUInt32LE(data.length, e + 8);
    header.writeUInt32LE(offset, e + 12);
    offset += data.length;
  });
  return Buffer.concat([header, ...images.map((image) => image.data)]);
}

mkdirSync(outDir, { recursive: true });
const sizes = [256, 128, 64, 48, 32, 24, 16];
const images = sizes.map((size) => {
  const rgba = render(size);
  // 256 px is stored as PNG to keep the file small; Windows Vista and later read both
  return { size, data: size >= 256 ? png(size, rgba) : bmpEntry(size, rgba) };
});
writeFileSync(`${outDir}icon.ico`, ico(images));
writeFileSync(`${outDir}icon.png`, images[0].data);
console.log(`wrote ${outDir}icon.ico (${sizes.join(", ")} px) and icon.png`);
