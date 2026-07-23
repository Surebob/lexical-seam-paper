"""Record canonical blocked state and preserve legacy coupled-logistic 5c output."""

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


def write_blocked() -> None:
    outputs = EXPERIMENT_DIR / "outputs"
    outputs.mkdir(parents=True, exist_ok=True)
    (EXPERIMENT_DIR / "BLOCKED.md").write_text(
        "\n".join(
            [
                "# Experiment 5c Canonical Status: BLOCKED",
                "",
                "Canonical Phase 2 model: `decoupled_erf`.",
                "",
                "The available smooth-parameter-sweep source bundle varies the historical",
                "coupled-logistic/sigmoid-width parameterization. It is legitimate prior",
                "research, but it is not the decoupled-erf parameter-control experiment.",
                "",
                "Missing canonical sources:",
                *[f"- {item}" for item in CANONICAL_MISSING_SOURCES],
                "",
                "Per the approved legacy policy, the historical coupled-logistic sweep is",
                "preserved under `archive/legacy_coupled_logistic/` as first-class research",
                "record. The BLOCKED status applies only to the current decoupled-erf canon.",
                "",
            ]
        )
    )
    manifest = {
        "experiment_id": EXPERIMENT_ID,
        "canonical_model": CANONICAL_MODEL,
        "canonical_status": CANONICAL_STATUS,
        "missing_sources": CANONICAL_MISSING_SOURCES,
        "legacy_archive": str(LEGACY_ARCHIVE_DIR.relative_to(EXPERIMENT_DIR)),
        "legacy_model": LEGACY_MODEL,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "canonical_outputs": [],
    }
    (outputs / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    write_blocked()
    script = LEGACY_ARCHIVE_DIR / "scripts" / "run_legacy_coupled_logistic.py"
    subprocess.run([sys.executable, str(script)], check=True, cwd=EXPERIMENT_DIR.parents[1])
    print(f"{EXPERIMENT_ID}: canonical {CANONICAL_STATUS}; legacy archive regenerated.")


if __name__ == "__main__":
    main()

