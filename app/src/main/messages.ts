// Turns aligner exit codes and error events into messages a non-developer can act on.
// Exit codes are documented in the repo README ("For programs calling the aligner").

import type { FriendlyError } from "../shared/types";
import type { AlignerExit } from "./aligner";

export function friendlyError(exit: AlignerExit): FriendlyError {
  if (exit.spawnError !== undefined) {
    return {
      title: "The aligner couldn't start",
      detail: exit.spawnError,
      hint: "Reinstall the app. If it keeps happening, your antivirus may be blocking run.exe.",
    };
  }

  const message = exit.error?.message ?? exit.otherOutput.slice(-3).join("\n");

  switch (exit.code) {
    case 3:
      return {
        title: "Couldn't read an input file",
        detail: message,
        hint: "Check that the file is FASTA or FASTQ, optionally gzipped (.gz), and that it isn't empty.",
      };
    case 4:
      return {
        title: "Couldn't save the results",
        detail: message,
        hint: "Check that your disk has free space. The 1,000,000-read set needs about 160 MB.",
      };
    case 5:
      return {
        title: "Ran out of memory",
        detail: message,
        hint: "Close other programs, or try a smaller reads file.",
      };
    default:
      return {
        title: "The aligner stopped unexpectedly",
        detail: message || `Exit code ${exit.code ?? "unknown"}.`,
        hint: "Try again. If it keeps happening, please report it on GitHub with the details above.",
      };
  }
}
