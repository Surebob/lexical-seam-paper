"""Record canonical blocked state and preserve legacy coupled-logistic outputs.

This script does not run fresh seam-theory computation. It records that the
current decoupled-erf canonical experiment is blocked and delegates historical
coupled-logistic reconstruction to the archive producer.
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from source_config import (
    CANONICAL_MISSING_SOURCES,
    CANONICAL_MODEL,
    CANONICAL_STATUS,
    EXPERIMENT_DIR,
    EXPERIMENT_ID,
    LEGACY_ARCHIVE_DIR,
    LEGACY_MODEL,
)


def write_canonical_blocked_files() -> None:
    outputs = EXPERIMENT_DIR / "outputs"
    outputs.mkdir(parents=True, exist_ok=True)

    blocked = EXPERIMENT_DIR / "BLOCKED.md"
    blocked.write_text(
        "\n".join(
            [
                "# Experiment 4 Canonical Status: BLOCKED",
                "",
                "Canonical Phase 2 model: `decoupled_erf`.",
                "",
                "The historical seam-theory bundles available in `results/` are legitimate",
                "coupled-logistic experiments from the project's prior canonical model, but",
                "they do not supply the decoupled-erf seam-theory validation requested by the",
                "current Phase 2 model canon.",
                "",
                "Missing canonical sources:",
                *[f"- {item}" for item in CANONICAL_MISSING_SOURCES],
                "",
                "Per the approved Option C migration decision, the coupled-logistic work is",
                "preserved as first-class research record under",
                "`archive/legacy_coupled_logistic/`. The BLOCKED status applies only to the",
                "current canonical decoupled-erf output, not to the scientific validity of",
                "the legacy archive.",
                "",
            ]
        )
    )

    manifest = {
        "experiment_id": EXPERIMENT_ID,
        "canonical_model": CANONICAL_MODEL,
        "canonical_status": CANONICAL_STATUS,
        "blocked_reason": "No decoupled-erf seam sign/projection/second-order source outputs exist in the historical result bundles.",
        "missing_sources": CANONICAL_MISSING_SOURCES,
        "legacy_archive": str(LEGACY_ARCHIVE_DIR.relative_to(EXPERIMENT_DIR)),
        "legacy_model": LEGACY_MODEL,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "canonical_outputs": [],
    }
    (outputs / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def run_legacy_archive_producer() -> None:
    script = LEGACY_ARCHIVE_DIR / "scripts" / "run_legacy_coupled_logistic.py"
    subprocess.run([sys.executable, str(script)], check=True, cwd=EXPERIMENT_DIR.parents[1])


def main() -> None:
    write_canonical_blocked_files()
    run_legacy_archive_producer()
    print(f"{EXPERIMENT_ID}: canonical {CANONICAL_STATUS}; legacy archive regenerated.")


if __name__ == "__main__":
    main()

