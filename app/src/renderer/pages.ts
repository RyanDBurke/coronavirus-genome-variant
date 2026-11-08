// One page per step. Only the current page is shown; the step indicator at the top lets you go back
// to any step you've reached, and later steps unlock as you go (Continue, or a finished run).

export const PAGES = ["reference", "reads", "align", "explore", "variants"] as const;
export type Page = (typeof PAGES)[number];

export class Pages {
  private current = 0;
  /** index of the furthest page that can be opened */
  private unlocked = 0;
  /** while a run is in progress, the inputs can't change, so moving between pages is paused */
  private busy = false;

  private readonly sections = PAGES.map((page) => document.querySelector<HTMLElement>(`.page[data-page="${page}"]`)!);
  private readonly steps = PAGES.map((page) => document.querySelector<HTMLButtonElement>(`.steps [data-step="${page}"]`)!);

  constructor() {
    this.steps.forEach((step, i) => step.addEventListener("click", () => this.goTo(PAGES[i])));
    document.querySelectorAll<HTMLButtonElement>("[data-next]").forEach((button) =>
      button.addEventListener("click", () => {
        this.unlock(PAGES[this.current + 1]);
        this.goTo(PAGES[this.current + 1]);
      }));
    document.querySelectorAll<HTMLButtonElement>("[data-back]").forEach((button) =>
      button.addEventListener("click", () => this.goTo(PAGES[this.current - 1])));
    fillBanners();
    this.render(false);
  }

  get page(): Page {
    return PAGES[this.current];
  }

  goTo(page: Page) {
    const index = PAGES.indexOf(page);
    if (index < 0 || index > this.unlocked || this.busy || index === this.current) return;
    this.current = index;
    this.render(true);
  }

  unlock(page: Page) {
    this.unlocked = Math.max(this.unlocked, PAGES.indexOf(page));
    this.render(false);
  }

  /** Makes pages after `page` unreachable again, leaving them if one is open. */
  lockAfter(page: Page) {
    const index = PAGES.indexOf(page);
    this.unlocked = Math.min(this.unlocked, index);
    if (this.current > index) {
      this.current = index;
      this.render(true);
    } else {
      this.render(false);
    }
  }

  setBusy(busy: boolean) {
    this.busy = busy;
    this.render(false);
  }

  private render(moved: boolean) {
    this.sections.forEach((section, i) => (section.hidden = i !== this.current));
    this.steps.forEach((step, i) => {
      step.disabled = this.busy ? i !== this.current : i > this.unlocked;
      step.classList.toggle("done", i < this.current);
      if (i === this.current) step.setAttribute("aria-current", "step");
      else step.removeAttribute("aria-current");
    });
    document.querySelectorAll<HTMLButtonElement>("[data-back]").forEach((button) => (button.disabled = this.busy));
    if (moved) {
      window.scrollTo(0, 0);
      // move keyboard and screen-reader focus to the new page's title
      this.sections[this.current].querySelector<HTMLElement>("h2")?.focus({ preventScroll: true });
    }
  }
}

/** Fills each decorative banner with its word repeated, separated by stars (as on rburke.dev).
 *  Two identical rows let the CSS animation loop seamlessly. */
function fillBanners() {
  document.querySelectorAll<HTMLElement>(".banner[data-word]").forEach((banner) => {
    const track = document.createElement("div");
    track.className = "banner-track";
    for (let row = 0; row < 2; row++) {
      const line = document.createElement("div");
      line.className = "banner-row";
      for (let i = 0; i < 8; i++) {
        const item = document.createElement("span");
        item.className = "banner-item";
        const star = document.createElement("span");
        star.className = "star";
        star.textContent = "*";
        item.append(banner.dataset.word!, " ", star);
        line.append(item);
      }
      track.append(line);
    }
    banner.replaceChildren(track);
  });
}
