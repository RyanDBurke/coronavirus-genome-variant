// Exposes a small, fixed API to the window as `window.aligner`. The window gets no Node.js access.

import { contextBridge, ipcRenderer, type IpcRendererEvent } from "electron";
import type { AlignerApi, RunEvent, RunResult } from "../shared/types";

function subscribe<T>(channel: string, listener: (value: T) => void) {
  const wrapped = (_event: IpcRendererEvent, value: T) => listener(value);
  ipcRenderer.on(channel, wrapped);
  return () => {
    ipcRenderer.removeListener(channel, wrapped);
  };
}

const api: AlignerApi = {
  getAppInfo: () => ipcRenderer.invoke("app:info"),
  openLicense: (id) => ipcRenderer.invoke("app:license", id),
  getPresets: () => ipcRenderer.invoke("presets:get"),
  chooseFile: (kind) => ipcRenderer.invoke("file:choose", kind),
  start: (request) => ipcRenderer.invoke("run:start", request),
  cancel: () => ipcRenderer.invoke("run:cancel"),
  saveSam: (runId) => ipcRenderer.invoke("run:save", runId),
  reveal: (runId) => ipcRenderer.invoke("run:reveal", runId),
  region: (runId, start, end) => ipcRenderer.invoke("view:region", runId, start, end),
  exportVariants: (runId) => ipcRenderer.invoke("export:variants", runId),
  exportCoverage: (runId) => ipcRenderer.invoke("export:coverage", runId),
  onEvent: (listener) => subscribe<RunEvent>("run:event", listener),
  onFinished: (listener) => subscribe<RunResult>("run:finished", listener),
};

contextBridge.exposeInMainWorld("aligner", api);
