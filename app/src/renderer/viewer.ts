// Results viewer: whole-genome depth overview, a zoomable detail track with reads, and the variants table.
// Everything is drawn on canvases from data the main process prepares; nothing is loaded from the network.

import type { AlignerApi, RegionData, Variant, ViewSummary } from "../shared/types.js";

const MAX_READ_REGION = 2000;   // must match analysis.ts
const MIN_WIDTH = 40;           // narrowest region, bp
// tuned for the dark background; the legend swatches in styles.css must match
const BASE_COLOURS: Record<string, string> = { A: "#3ECF6E", C: "#5B93FF", G: "#F5A524", T: "#FF5A4E", N: "#8A8A80" };
const INSERTION_COLOUR = "#C77DFF";

const ROW_HEIGHT = 9;
const ROW_GAP = 2;
const MAX_ROWS = 40;
const RULER = 22;
const REF_ROW = 16;
const DEPTH_TRACK = 64;
const TRACK_GAP = 8;

const fmt = (n: number) => Math.round(n).toLocaleString("en-US");

function css(name: string) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

/** Sizes a canvas for its CSS box and the screen's pixel density; returns a context in CSS pixels. */
function prepare(canvas: HTMLCanvasElement, cssHeight: number) {
  const dpr = window.devicePixelRatio || 1;
  const width = canvas.clientWidth;
  canvas.style.height = `${cssHeight}px`;
  canvas.width = Math.round(width * dpr);
  canvas.height = Math.round(cssHeight * dpr);
  const ctx = canvas.getContext("2d")!;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, width, cssHeight);
  return { ctx, width };
}

/** A "nice" tick step (1, 2 or 5 x 10^k) giving roughly `target` ticks over `span`. */
function tickStep(span: number, target: number) {
  const raw = span / target;
  const pow = 10 ** Math.floor(Math.log10(raw));
  for (const m of [1, 2, 5, 10]) if (m * pow >= raw) return m * pow;
  return 10 * pow;
}

function parseCigar(cigar: string): Array<[number, string]> {
  return Array.from(cigar.matchAll(/(\d+)([MIDNSHP=X])/g), (m) => [Number(m[1]), m[2]] as [number, string]);
}

export class Viewer {
  private start = 0;
  private end: number;
  private region: RegionData | null = null;
  private requestId = 0;
  private drag: { x: number; start: number; end: number } | null = null;
  private resize = new ResizeObserver(() => this.drawAll());

  private readonly overview = document.getElementById("overviewCanvas") as HTMLCanvasElement;
  private readonly detail = document.getElementById("detailCanvas") as HTMLCanvasElement;

  /** onJump runs before a variant-table click moves the view, so the page holding the tracks can be shown */
  constructor(private readonly api: AlignerApi, private readonly runId: string, private readonly summary: ViewSummary,
              private readonly onJump: () => void = () => {}) {
    this.end = summary.length;
    this.renderStats();
    this.renderVariants();
    this.bind();
    this.resize.observe(this.overview);
    void this.load();
  }

  destroy() {
    this.resize.disconnect();
    this.unbind.forEach((fn) => fn());
  }

  // ---- navigation -----------------------------------------------------------------

  /** Shows [start, end), clamped to the genome and to the minimum width. */
  goTo(start: number, end: number) {
    const length = this.summary.length;
    let width = Math.max(MIN_WIDTH, Math.min(length, Math.round(end - start)));
    let s = Math.round(start);
    if (s < 0) s = 0;
    if (s + width > length) s = Math.max(0, length - width);
    width = Math.min(width, length - s);
    if (s === this.start && s + width === this.end && this.region) return;
    this.start = s;
    this.end = s + width;
    void this.load();
  }

  centerOn(pos0: number, width = this.end - this.start) {
    this.goTo(pos0 - width / 2, pos0 + width / 2);
  }

  zoom(factor: number) {
    const mid = (this.start + this.end) / 2;
    this.centerOn(mid, (this.end - this.start) * factor);
  }

  pan(fraction: number) {
    const shift = (this.end - this.start) * fraction;
    this.goTo(this.start + shift, this.end + shift);
  }

  private async load() {
    const id = ++this.requestId;
    this.updateLabel();
    this.drawOverview();
    const region = await this.api.region(this.runId, this.start, this.end);
    if (id !== this.requestId) return;   // a newer request superseded this one
    this.region = region;
    this.drawDetail();
  }

  // ---- wiring -----------------------------------------------------------------------

  private unbind: Array<() => void> = [];

  private on<K extends keyof HTMLElementEventMap>(el: HTMLElement, type: K, fn: (e: HTMLElementEventMap[K]) => void) {
    el.addEventListener(type, fn as EventListener);
    this.unbind.push(() => el.removeEventListener(type, fn as EventListener));
  }

  private bind() {
    const byId = (id: string) => document.getElementById(id)!;
    this.on(byId("zoomIn"), "click", () => this.zoom(0.5));
    this.on(byId("zoomOut"), "click", () => this.zoom(2));
    this.on(byId("panLeft"), "click", () => this.pan(-0.5));
    this.on(byId("panRight"), "click", () => this.pan(0.5));
    this.on(byId("showAll"), "click", () => this.goTo(0, this.summary.length));
    this.on(byId("gotoForm"), "submit", (e) => {
      e.preventDefault();
      this.goToText((byId("gotoInput") as HTMLInputElement).value);
    });

    // overview: click or drag to move the window
    const overviewPos = (e: MouseEvent) => {
      const rect = this.overview.getBoundingClientRect();
      return ((e.clientX - rect.left) / rect.width) * this.summary.length;
    };
    const moveTo = (e: MouseEvent) => {
      const width = this.end - this.start;
      this.centerOn(overviewPos(e), width >= this.summary.length ? 1000 : width);
    };
    this.on(this.overview, "mousedown", (e) => {
      moveTo(e);
      const move = (ev: MouseEvent) => moveTo(ev);
      const up = () => { window.removeEventListener("mousemove", move); window.removeEventListener("mouseup", up); };
      window.addEventListener("mousemove", move);
      window.addEventListener("mouseup", up);
    });

    // detail: drag to pan, hover for the base under the pointer
    this.on(this.detail, "mousedown", (e) => {
      this.drag = { x: e.clientX, start: this.start, end: this.end };
      const move = (ev: MouseEvent) => {
        if (!this.drag) return;
        const bpPerPx = (this.drag.end - this.drag.start) / this.detail.clientWidth;
        const shift = (this.drag.x - ev.clientX) * bpPerPx;
        this.goTo(this.drag.start + shift, this.drag.end + shift);
      };
      const up = () => {
        this.drag = null;
        window.removeEventListener("mousemove", move);
        window.removeEventListener("mouseup", up);
      };
      window.addEventListener("mousemove", move);
      window.addEventListener("mouseup", up);
    });
    this.on(this.detail, "mousemove", (e) => this.hover(e));
    this.on(this.detail, "mouseleave", () => this.setHover(""));
    this.on(this.detail, "dblclick", (e) => {
      const rect = this.detail.getBoundingClientRect();
      const pos = this.start + ((e.clientX - rect.left) / rect.width) * (this.end - this.start);
      this.centerOn(pos, (this.end - this.start) / 2);
    });
  }

  /** Accepts "28000", "28,000", "27950-28050" or "27,950–28,050" (1-based, inclusive). */
  private goToText(text: string) {
    const nums = text.replace(/,/g, "").match(/\d+/g)?.map(Number) ?? [];
    if (nums.length === 0) return;
    if (nums.length === 1) {
      this.centerOn(nums[0] - 1, Math.min(this.end - this.start, 200));
    } else {
      const [a, b] = nums[0] <= nums[1] ? [nums[0], nums[1]] : [nums[1], nums[0]];
      this.goTo(a - 1, b);
    }
  }

  // ---- stats and table -------------------------------------------------------------------

  private renderStats() {
    const s = this.summary;
    const set = (id: string, text: string) => (document.getElementById(id)!.textContent = text);
    set("covMean", `${fmt(s.meanDepth)}×`);
    set("covMax", `${fmt(s.maxDepth)}×`);
    set("cov1", `${(s.covered1 * 100).toFixed(1)}%`);
    set("cov10", `${(s.covered10 * 100).toFixed(1)}%`);
    set("varCount", fmt(s.variants.length));
    set("variantRule", `Sites where at least ${Math.round(s.minFreq * 100)}% of at least ${s.minDepth} reads differ from ${s.refName}. Deletions over neighbouring bases, and insertions within ${s.insertionMergeBp} bp of each other, count as one site.`);
  }

  private renderVariants() {
    const body = document.getElementById("variantRows")!;
    body.replaceChildren();
    document.getElementById("variantEmpty")!.hidden = this.summary.variants.length > 0;
    document.getElementById("variantTable")!.hidden = this.summary.variants.length === 0;

    for (const v of this.summary.variants) {
      const tr = document.createElement("tr");
      tr.tabIndex = 0;
      const cells = [
        fmt(v.pos) + (v.length > 1 ? `–${fmt(v.pos + v.length - 1)}` : ""),
        change(v),
        v.type === "deletion" && v.length > 1 ? `deletion (${v.length} bp)` : v.type,
        fmt(v.depth),
      ];
      for (const text of cells) {
        const td = document.createElement("td");
        td.textContent = text;
        tr.append(td);
      }
      const freq = document.createElement("td");
      freq.className = "freq";
      const bar = document.createElement("span");
      bar.className = "freqbar";
      bar.style.setProperty("--f", String(v.freq));
      const pct = document.createElement("span");
      pct.textContent = `${(v.freq * 100).toFixed(1)}%`;
      freq.append(bar, pct);
      tr.append(freq);

      const jump = () => {
        this.onJump();
        this.centerOn(v.pos - 1 + v.length / 2, Math.max(100, v.length + 60));
      };
      tr.addEventListener("click", jump);
      tr.addEventListener("keydown", (e) => {
        if (e.key === "Enter" || e.key === " ") {
          e.preventDefault();
          jump();
        }
      });
      body.append(tr);
    }
  }

  // ---- drawing --------------------------------------------------------------------------------

  private drawAll() {
    this.drawOverview();
    this.drawDetail();
  }

  private updateLabel() {
    const label = document.getElementById("regionLabel")!;
    label.textContent = `${this.summary.refName}:${fmt(this.start + 1)}–${fmt(this.end)} · ${fmt(this.end - this.start)} bp`;
  }

  private drawOverview() {
    const s = this.summary;
    const height = 86;
    const axis = 16;
    const { ctx, width } = prepare(this.overview, height);
    const plotH = height - axis - 6;
    const top = 6;
    const x = (pos: number) => (pos / s.length) * width;
    const y = (d: number) => top + plotH - (d / Math.max(s.maxDepth, 1)) * plotH;

    // depth: bin max as a light area, bin mean as the solid area
    for (const [values, colour] of [[s.binMax, css("--accent-soft")], [s.binMean, css("--accent")]] as const) {
      ctx.beginPath();
      ctx.moveTo(0, y(0));
      values.forEach((d, i) => {
        ctx.lineTo(x(i * s.binSize), y(d));
        ctx.lineTo(x(Math.min(s.length, (i + 1) * s.binSize)), y(d));
      });
      ctx.lineTo(width, y(0));
      ctx.closePath();
      ctx.fillStyle = colour;
      ctx.fill();
    }

    // variant ticks along the top
    ctx.fillStyle = css("--bad");
    for (const v of s.variants) ctx.fillRect(Math.round(x(v.pos - 1)) - 0.5, 0, 1.5, 5);

    // axis
    ctx.fillStyle = css("--muted");
    ctx.font = `11px ${css("--mono")}`;
    ctx.textBaseline = "bottom";
    const step = tickStep(s.length, 6);
    for (let p = 0; p <= s.length; p += step) {
      const px = x(p);
      ctx.fillRect(px, y(0), 1, 3);
      ctx.textAlign = p === 0 ? "left" : px > width - 30 ? "right" : "center";
      ctx.fillText(p === 0 ? "1" : fmt(p), px, height);
    }

    // current window
    if (this.end - this.start < s.length) {
      const x0 = x(this.start);
      const w = Math.max(3, x(this.end) - x0);
      ctx.strokeStyle = css("--fg");
      ctx.lineWidth = 1.5;
      ctx.strokeRect(x0, top - 2, w, plotH + 4);
    }
  }

  private drawDetail() {
    const r = this.region;
    const span = this.end - this.start;
    const showReads = r !== null && span <= MAX_READ_REGION;
    const rowsUsed = showReads ? Math.min(MAX_ROWS, packRows(r!.reads).rows) : 0;
    const height = RULER + REF_ROW + TRACK_GAP + DEPTH_TRACK + TRACK_GAP + (showReads ? rowsUsed * (ROW_HEIGHT + ROW_GAP) + 4 : 0);
    const { ctx, width } = prepare(this.detail, Math.max(height, RULER + REF_ROW + DEPTH_TRACK + 2 * TRACK_GAP));
    const scale = width / span;
    const x = (pos: number) => (pos - this.start) * scale;
    const fg = css("--fg");
    const muted = css("--muted");
    const line = css("--line");
    const mono = css("--mono");

    // ruler
    ctx.fillStyle = muted;
    ctx.font = `11px ${mono}`;
    ctx.textBaseline = "top";
    const step = tickStep(span, Math.max(2, Math.floor(width / 110)));
    const first = Math.ceil((this.start + 1) / step) * step;
    for (let p = first; p <= this.end; p += step) {
      const px = x(p - 0.5);
      ctx.fillRect(px, 14, 1, 6);
      ctx.textAlign = px < 30 ? "left" : px > width - 30 ? "right" : "center";
      ctx.fillText(fmt(p), px, 0);
    }

    if (!r) {
      this.setNote("Loading…");
      return;
    }

    // reference
    let yTop = RULER;
    if (scale >= 1) {
      for (let i = 0; i < r.ref.length; i++) {
        const base = r.ref[i];
        if (scale >= 8) {
          ctx.fillStyle = BASE_COLOURS[base] ?? muted;
          ctx.font = `600 ${Math.min(13, scale)}px ${mono}`;
          ctx.textAlign = "center";
          ctx.fillText(base, x(r.start + i + 0.5), yTop + 1);
        } else {
          ctx.fillStyle = BASE_COLOURS[base] ?? muted;
          ctx.fillRect(x(r.start + i), yTop + 3, Math.max(1, scale), REF_ROW - 6);
        }
      }
    } else {
      ctx.fillStyle = line;
      ctx.fillRect(0, yTop + 6, width, 4);
    }
    yTop += REF_ROW + TRACK_GAP;

    // depth: one bar per base, or the max per pixel when zoomed out
    const maxDepth = Math.max(1, ...r.depth);
    ctx.fillStyle = css("--accent");
    if (scale >= 1) {
      r.depth.forEach((d, i) => {
        const h = (d / maxDepth) * (DEPTH_TRACK - 12);
        ctx.fillRect(x(r.start + i), yTop + DEPTH_TRACK - h, Math.max(1, scale - (scale > 4 ? 1 : 0)), h);
      });
      // bases that differ from the reference, stacked in their colours
      if (r.counts.length) {
        for (let i = 0; i < r.depth.length; i++) {
          const d = r.depth[i];
          if (d === 0) continue;
          let yy = yTop + DEPTH_TRACK;
          for (const [col, base] of [[0, "A"], [1, "C"], [2, "G"], [3, "T"]] as const) {
            if (base === r.ref[i]) continue;
            const c = r.counts[i * 7 + col];
            if (c / d < 0.05) continue;
            const h = (c / maxDepth) * (DEPTH_TRACK - 12);
            yy -= h;
            ctx.fillStyle = BASE_COLOURS[base];
            ctx.fillRect(x(r.start + i), yy, Math.max(1, scale - (scale > 4 ? 1 : 0)), h);
          }
          ctx.fillStyle = css("--accent");
        }
      }
    } else {
      for (let px = 0; px < width; px++) {
        const a = Math.floor(px / scale);
        const b = Math.min(r.depth.length, Math.max(a + 1, Math.floor((px + 1) / scale)));
        let m = 0;
        for (let i = a; i < b; i++) if (r.depth[i] > m) m = r.depth[i];
        const h = (m / maxDepth) * (DEPTH_TRACK - 12);
        ctx.fillRect(px, yTop + DEPTH_TRACK - h, 1, h);
      }
    }
    ctx.fillStyle = muted;
    ctx.font = `11px ${mono}`;
    ctx.textAlign = "left";
    ctx.fillText(`depth (max ${fmt(maxDepth)}×)`, 2, yTop);
    yTop += DEPTH_TRACK + TRACK_GAP;

    if (!showReads) {
      this.setNote(`Zoom in to ${fmt(MAX_READ_REGION)} bp or less to see individual reads. Drag to pan, double-click to zoom in.`);
      return;
    }

    // reads
    const { placed } = packRows(r.reads);
    const readFill = css("--read");
    for (const { read, row } of placed) {
      if (row >= MAX_ROWS) continue;
      const y0 = yTop + row * (ROW_HEIGHT + ROW_GAP);
      let pos = read.start;
      let q = 0;
      for (const [n, op] of parseCigar(read.cigar)) {
        if (op === "M" || op === "=" || op === "X") {
          ctx.fillStyle = readFill;
          ctx.fillRect(x(pos), y0, n * scale, ROW_HEIGHT);
          for (let k = 0; k < n; k++) {
            const refBase = r.ref[pos + k - r.start];
            const base = read.seq[q + k];
            if (refBase !== undefined && base !== refBase) {
              ctx.fillStyle = BASE_COLOURS[base] ?? muted;
              ctx.fillRect(x(pos + k), y0, Math.max(1.5, scale), ROW_HEIGHT);
            }
          }
          pos += n;
          q += n;
        } else if (op === "D" || op === "N") {
          ctx.fillStyle = fg;
          ctx.fillRect(x(pos), y0 + ROW_HEIGHT / 2 - 0.5, n * scale, 1);
          pos += n;
        } else if (op === "I") {
          ctx.fillStyle = INSERTION_COLOUR;
          ctx.fillRect(x(pos) - 1, y0 - 1, 2, ROW_HEIGHT + 2);
          q += n;
        } else if (op === "S") {
          q += n;
        }
      }
    }

    const hidden = placed.filter((p) => p.row >= MAX_ROWS).length;
    const shown = r.reads.length - hidden;
    this.setNote(
      `Showing ${fmt(shown)} of ${fmt(r.totalReads)} reads overlapping this region` +
        (shown < r.totalReads ? " (an even sample across the region)" : "") +
        ". Coloured marks are bases that differ from the reference. Drag to pan, double-click to zoom in.",
    );
  }

  private setNote(text: string) {
    document.getElementById("detailNote")!.textContent = text;
  }

  private hover(e: MouseEvent) {
    const r = this.region;
    if (!r || this.drag) return;
    const rect = this.detail.getBoundingClientRect();
    const i = Math.floor(((e.clientX - rect.left) / rect.width) * (this.end - this.start));
    if (i < 0 || i >= r.depth.length) return;
    let text = `Position ${fmt(r.start + i + 1)} · reference ${r.ref[i]} · depth ${fmt(r.depth[i])}`;
    if (r.counts.length) {
      const c = r.counts.slice(i * 7, i * 7 + 7);
      text += ` · A ${fmt(c[0])} C ${fmt(c[1])} G ${fmt(c[2])} T ${fmt(c[3])}`;
      if (c[5]) text += ` · deleted ${fmt(c[5])}`;
      if (c[6]) text += ` · insertion after ${fmt(c[6])}`;
    }
    this.setHover(text);
  }

  private setHover(text: string) {
    document.getElementById("hoverInfo")!.textContent = text || "Hover over the track for base counts.";
  }
}

function change(v: Variant) {
  if (v.type === "SNV") return `${v.ref} → ${v.alt}`;
  if (v.type === "deletion") return v.ref.length <= 12 ? `${v.ref} → –` : `${v.ref.length} bases deleted`;
  return v.length > 1 ? `bases inserted within ${v.ref}` : `${v.ref} → ${v.ref}+…`;
}

/** Greedy row packing: each read goes in the first row whose last read ended before it starts. */
function packRows(reads: RegionData["reads"]) {
  const sorted = [...reads].sort((a, b) => a.start - b.start);
  const rowEnds: number[] = [];
  const placed = sorted.map((read) => {
    let span = 0;
    for (const [n, op] of parseCigar(read.cigar)) if ("MDN=X".includes(op)) span += n;
    let row = rowEnds.findIndex((end) => end < read.start);
    if (row < 0) {
      row = rowEnds.length;
      rowEnds.push(0);
    }
    rowEnds[row] = read.start + span;
    return { read, row };
  });
  return { placed, rows: rowEnds.length };
}
