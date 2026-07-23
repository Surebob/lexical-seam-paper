from __future__ import annotations

from pathlib import Path

import run_s2_v3_windows as runner


def main() -> int:
    output_root = Path(__file__).resolve().parent / "outputs"
    preflight = runner.preflight_check()
    if not preflight["ok"]:
        print("Pre-flight check failed:")
        for message in preflight["messages"]:
            print(f"- {message}")
        return 2
    smoke = runner.run_smoke_test(base_seed=20260415, n_starts=100, max_nfev=12000, output_root=output_root)
    max_diff = max(float(row["abs_diff"]) for row in smoke["comparisons"])
    print(f"Smoke test pass: {smoke['passes_smoke']} (max abs diff {max_diff:.12g})")
    return 0 if smoke["passes_smoke"] else 3


if __name__ == "__main__":
    raise SystemExit(main())

