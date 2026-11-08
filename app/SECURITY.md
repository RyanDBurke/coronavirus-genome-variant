# Desktop security and release checks

The app is offline. Inputs remain on this computer. No automatic updates or update checks are enabled; users install newer releases manually.

## Electron checklist

Reviewed against https://www.electronjs.org/docs/latest/tutorial/security on 2026-11-07.

| # | Recommendation | Implementation / outcome |
|---|---|---|
| 1 | Secure content | Only bundled HTML, CSS and JavaScript; network requests denied. |
| 2 | Disable Node integration | Explicitly disabled in the renderer. |
| 3 | Context isolation | Explicitly enabled. |
| 4 | Sandbox | Explicitly enabled. |
| 5 | Permissions | Requests and checks denied. |
| 6 | Web security | Enabled. |
| 7 | CSP | No inline scripts, network connections, forms or base overrides. |
| 8 | Insecure content | Not enabled. |
| 9 | Experimental features | Not enabled. |
| 10 | Blink features | No additional features enabled. |
| 11 | Webview popups | No webviews; attachment blocked. |
| 12 | Webview options | Attachment blocked entirely. |
| 13 | Navigation | Renderer navigation blocked. |
| 14 | New windows | All requests denied. |
| 15 | External URLs | No openExternal handler; reveal opens only a known run in Explorer. The About box opens one of three fixed, bundled license files (app LICENSE, Electron LICENSE, LICENSES.chromium.html) with `shell.openPath`, chosen by a fixed id; the page cannot supply a path. |
| 16 | Current Electron | Release maintainer must check supported versions and patch availability. |
| 17 | IPC sender | Exact index URL, owning BrowserWindow and main-frame identity required. |
| 18 | Avoid file protocol | Deliberate bounded exception: exact URL allowlist for index.html, styles.css, renderer.js, viewer.js, pages.js and the two bundled font files (Playfair Display and Oi, under fonts/). All other request URLs denied, including other local files and UNC URLs. No general directory access. |
| 19 | Fuses | runAsNode, NODE_OPTIONS and CLI inspector disabled; asar integrity and asar-only loading enabled; cookie encryption enabled. |
| 20 | API exposure | Fixed contextBridge API; no raw ipcRenderer or IPC event objects exposed. |

Packaged windows disable DevTools and remove the application menu. Development builds retain debugging for tests. File-picker tokens are random, held only in main-process memory, and valid for the current process session. They do not grant arbitrary path access. `app:info` returns only the app name and version. The aligner is spawned without a shell. Renderer isolation does not protect against every possible bug in native parsing of untrusted input.

## Dependency and update policy

Before each release, verify Electron is on the latest patch of a supported major (Electron supports its latest three majors). Prefer dependency releases at least two weeks old. If the latest security patch is newer, resolve that conflict explicitly before release: take the security patch with owner agreement or delay publication; do not silently ship an outdated patch. Re-run all app and packaging checks after upgrades.

Run `npm audit --omit=dev` before packaging and record its result. All current dependencies are declared as development dependencies; this audit alone does not establish the safety of the bundled Electron/Chromium runtime. Check Electron advisories separately. Automatic updates require a future owner decision because they add network activity.

## Installer behavior

NSIS installs Windows x64 for the current user, with one click and asInvoker execution. Apps & features registration is provided by NSIS. `deleteAppDataOnUninstall` removes the app profile, including saved runs and settings, on ordinary uninstall; builder's update uninstall preserves data. Exported files outside the profile remain.

Packaged Electron uses the product name for userData: `%APPDATA%\Genome Variant Aligner`, with results in `runs`. Confirm this path in the installed acceptance test. Source tests use an isolated temporary profile. The app retains three run folders (roughly 450 MB for three 1M-read runs).

Signing is not configured. Future Azure Trusted Signing must sign the standalone run.exe before it is bundled and hashed, and sign the app and installer through builder. See workflow comments and package metadata for configuration and credentials. Checksums detect file changes but do not replace publisher signatures.

## Clean Windows acceptance (owner, pending)

Use a standard account on clean Windows 10/11 x64, without Node, MSYS2 or WSL:

1. Download the installer, run.exe and SHA256SUMS.txt from the authorized test release. Compare both executables using `Get-FileHash -Algorithm SHA256`.
2. Run the installer. Expect unsigned SmartScreen behavior (More info → Run anyway), with no administrator prompt.
3. Launch the installed app, align the bundled 10K preset and inspect coverage, variant sites and the read viewer. Export a SAM and CSV outside the app profile.
4. Try F12 and Ctrl+Shift+I/J; DevTools must stay closed. Confirm no application menu appears.
5. Confirm the Apps & features entry and the profile/runs location above. Uninstall through Apps & features.
6. Confirm the app, uninstall entry and profile/runs are removed, while exported files remain.

Tag-triggered publication and this clean-machine test remain acceptance gates; local builds alone do not complete phase 7.

## Local verification record (2026-11-07)

- Linux sanitizer suite and native Windows suite passed (25 unit checks each, SAM known-answer/1K validation and CLI checks).
- All 14 desktop tests passed, including direct network/local/UNC request rejection, popup denial, foreign-page IPC rejection and prevented navigation. The navigation probe runs last because Playwright retains the cancelled navigation before subsequent clicks.
- `npm audit --omit=dev`: zero vulnerabilities. Installer build passed; local unsigned x64 installer is 146,940,536 bytes. Bundled aligner completed a 10K run; packaged fuses inspected.
- Workflow PowerShell blocks parsed locally; stable/prerelease/build-metadata/invalid tag cases, exact version metadata restoration, checksum contents and release-note formatting checked with disposable fixtures. Content tests used Windows PowerShell 5's UTF-8 output because pwsh 7 is unavailable locally; CI's utf8NoBOM encoding remains to be exercised by an authorized release run.
- An installed app was already running. Its child-process command lines confirm `%APPDATA%\Genome Variant Aligner` as userData. The attempted second packaged launch exited successfully under the single-instance lock; this does not establish that the newly built packaged UI works. Existing owner app processes were left running.
- [Electron releases](https://releases.electronjs.org/releases/stable) and [44.10.2 notes](https://releases.electronjs.org/release/v44.10.2) checked: pinned 44.10.2 is the current patch of supported major 44 (latest three stable majors: 43–45), released November 3. Notes include upstream backports and Node updates; this is not an exhaustive Chromium advisory audit. The two-week maturity rule would permit that version from November 17, absent a separately agreed exception; recheck patch/support/advisory status when publishing.
- Pending: actual GitHub release execution, clean-machine install/uninstall and new packaged UI/DevTools acceptance. No release published and no tag pushed.

## Phase 7 release and phase 8 verification (2026-11-07)

- CI for `20b6601` passed (Linux and Windows jobs). The authorized test tag `v0.1.0-rc.1` ran the release workflow successfully and published a prerelease with the installer, run.exe and SHA256SUMS.txt. Downloaded files matched the published checksums. Release notes rendered correctly. Tags with `+build` metadata remain unexercised in CI.
- Partial owner acceptance in Windows Sandbox with the rc.1 installer: download, checksum verification, install, launch and a 10K run passed. Edge held both unsigned executables until **Keep** was chosen; this is now documented in the README and release notes. Sandbox became unstable before the viewer, export, DevTools-shortcut and uninstall checks, and its account is an administrator, so the no-UAC check is not established. The full checklist above will be repeated once on the final installer, using a standard Windows account.
- Phase 8 added two IPC calls (`app:info`, `app:license`), both behind the existing sender check. The click-through test covers the About box and the long-read warning; all 15 desktop tests pass.

## Redesign (redesign branch, 2026-11-08)

- The window now shows one page per step, styled after rburke.dev. The renderer allowlist gained exactly three files: `pages.js` and two `.woff2` fonts. The fonts are copied from pinned devDependencies at build time; nothing is fetched at runtime. The click-through test checks that a bundled font URL is allowed while other local, UNC and HTTPS URLs stay blocked.
- The About box gained a fourth fixed license id (`fonts`), which opens `fonts-OFL.txt` shipped beside the app resources.
