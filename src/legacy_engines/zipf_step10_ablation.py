import json
from pathlib import Path

import numpy as np

import eml_zipf_enriched_search as enriched


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "results" / "zipf_step10_ablation"
GRID = np.linspace(0.05, 0.95, 1000, dtype=np.float64)
CORPORA = [
    {
        "name": "Shakespeare",
        "slug": "shakespeare",
        "summary_path": ROOT / "results" / "zipf_enriched_search_full" / "summary.json",
    },
    {
        "name": "War and Peace",
        "slug": "war_and_peace",
        "summary_path": ROOT / "results" / "zipf_enriched_war_and_peace_full_seq" / "summary.json",
    },
]


def parse_expr(expr: str):
    def parse_at(index: int):
        if expr[index].isalpha():
            end = index
            while end < len(expr) and expr[end].isalpha():
                end += 1
            token = expr[index:end]
            if end < len(expr) and expr[end] == "[":
                end += 1
                children = []
                while True:
                    child, end = parse_at(end)
                    children.append(child)
                    if expr[end] == ",":
                        end += 1
                        continue
                    if expr[end] == "]":
                        end += 1
                        break
                return (token, *children), end
            return token, end
        if expr[index].isdigit():
            end = index
            while end < len(expr) and expr[end].isdigit():
                end += 1
            return expr[index:end], end
        raise ValueError(f"Bad token near {expr[index:index+12]!r}")

    node, end = parse_at(0)
    if end != len(expr):
        raise ValueError(f"Trailing content in expression: {expr[end:]!r}")
    return node


def eval_node(node, x_values: np.ndarray):
    if isinstance(node, str):
        if node == "x":
            return np.asarray(x_values, dtype=np.float64)
        if node == "1":
            return np.ones_like(x_values, dtype=np.float64)
        raise ValueError(f"Unsupported terminal {node!r}")

    op = node[0]
    if op in enriched.UNARY_FUNCS:
        child = eval_node(node[1], x_values)
        result = enriched.UNARY_FUNCS[op](child)
    elif op in enriched.BINARY_FUNCS:
        left = eval_node(node[1], x_values)
        right = eval_node(node[2], x_values)
        result = enriched.BINARY_FUNCS[op](left, right)
    else:
        raise ValueError(f"Unsupported op {op!r}")

    if result is None:
        raise ValueError(f"Operation {op!r} returned None during evaluation")
    return result


def polyfit_r2(x: np.ndarray, y: np.ndarray, degree: int):
    coeffs = np.polyfit(x, y, degree)
    pred = np.polyval(coeffs, x)
    residual = float(np.sum((y - pred) ** 2))
    total = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if total == 0.0 else 1.0 - residual / total
    return coeffs, pred, r2


def _scale_points(x, y, left, top, width, height):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    xmin, xmax = float(np.min(x)), float(np.max(x))
    ymin, ymax = float(np.min(y)), float(np.max(y))
    if xmax <= xmin:
        xmax = xmin + 1.0
    if ymax <= ymin:
        ymax = ymin + 1.0
    xs = left + (x - xmin) / (xmax - xmin) * width
    ys = top + height - (y - ymin) / (ymax - ymin) * height
    return " ".join(f"{px:.2f},{py:.2f}" for px, py in zip(xs, ys)), ymin, ymax


def write_svg_plot(plot_path: Path, title: str, x, step2, step10, diff, poly_pred, poly_degree, poly_r2):
    width = 1000
    height = 1200
    margin_x = 80
    panel_w = width - 2 * margin_x
    panel_h = 260
    panels = [
        (70, step2, step10, "Step 2", "Step 10", "#1f77b4", "#ff7f0e", "Step-2 vs Step-10 ZM Corrections"),
        (430, diff, None, "Difference", None, "#d62728", None, "Difference Curve"),
        (790, diff, poly_pred, "Difference", f"Poly deg {poly_degree}, R2={poly_r2:.6f}", "#d62728", "#2ca02c", "Polynomial Approximation to Difference"),
    ]

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2:.0f}" y="35" text-anchor="middle" font-size="24" font-family="Helvetica, Arial, sans-serif">{title}</text>',
    ]

    for top, y1, y2, label1, label2, color1, color2, panel_title in panels:
        parts.append(f'<rect x="{margin_x}" y="{top}" width="{panel_w}" height="{panel_h}" fill="none" stroke="#cccccc" stroke-width="1"/>')
        pts1, ymin, ymax = _scale_points(x, y1, margin_x, top, panel_w, panel_h)
        parts.append(f'<polyline fill="none" stroke="{color1}" stroke-width="2" points="{pts1}"/>')
        if y2 is not None:
            pts2, ymin2, ymax2 = _scale_points(x, y2, margin_x, top, panel_w, panel_h)
            ymin = min(ymin, ymin2)
            ymax = max(ymax, ymax2)
            parts.append(f'<polyline fill="none" stroke="{color2}" stroke-width="2" points="{pts2}"/>')
        parts.append(f'<text x="{width/2:.0f}" y="{top-12}" text-anchor="middle" font-size="18" font-family="Helvetica, Arial, sans-serif">{panel_title}</text>')
        parts.append(f'<text x="{margin_x}" y="{top+panel_h+24}" text-anchor="start" font-size="12" font-family="Helvetica, Arial, sans-serif">x in [0.05, 0.95]</text>')
        parts.append(f'<text x="{margin_x+panel_w}" y="{top+panel_h+24}" text-anchor="end" font-size="12" font-family="Helvetica, Arial, sans-serif">y in [{ymin:.3e}, {ymax:.3e}]</text>')
        legend_y = top + 18
        parts.append(f'<line x1="{margin_x+10}" y1="{legend_y}" x2="{margin_x+35}" y2="{legend_y}" stroke="{color1}" stroke-width="2"/>')
        parts.append(f'<text x="{margin_x+42}" y="{legend_y+4}" font-size="12" font-family="Helvetica, Arial, sans-serif">{label1}</text>')
        if y2 is not None:
            parts.append(f'<line x1="{margin_x+170}" y1="{legend_y}" x2="{margin_x+195}" y2="{legend_y}" stroke="{color2}" stroke-width="2"/>')
            parts.append(f'<text x="{margin_x+202}" y="{legend_y+4}" font-size="12" font-family="Helvetica, Arial, sans-serif">{label2}</text>')

    parts.append("</svg>")
    plot_path.write_text("\n".join(parts), encoding="utf-8")


def run_one(entry):
    summary = json.loads(entry["summary_path"].read_text())
    step2_expr = summary["zm_search"]["step_summary"][1]["top_candidates"][0]["expr"]
    step10_expr = summary["zm_search"]["best"]["expr"]

    step2_values = eval_node(parse_expr(step2_expr), GRID)
    step10_values = eval_node(parse_expr(step10_expr), GRID)
    diff = step10_values - step2_values

    poly_results = []
    best_poly = None
    for degree in (3, 4, 5):
        coeffs, pred, r2 = polyfit_r2(GRID, diff, degree)
        item = {
            "degree": degree,
            "r2": float(r2),
            "coefficients": [float(c) for c in coeffs],
            "prediction": pred,
        }
        poly_results.append(item)
        if best_poly is None or item["r2"] > best_poly["r2"]:
            best_poly = item

    plot_path = OUTDIR / f"{entry['slug']}_step10_ablation.svg"
    write_svg_plot(
        plot_path,
        f"{entry['name']}: Step-10 Ablation",
        GRID,
        step2_values,
        step10_values,
        diff,
        best_poly["prediction"],
        best_poly["degree"],
        best_poly["r2"],
    )

    return {
        "name": entry["name"],
        "slug": entry["slug"],
        "step2_expr": step2_expr,
        "step10_expr": step10_expr,
        "grid_min": float(GRID.min()),
        "grid_max": float(GRID.max()),
        "difference_min": float(np.min(diff)),
        "difference_max": float(np.max(diff)),
        "difference_mean_abs": float(np.mean(np.abs(diff))),
        "difference_std": float(np.std(diff)),
        "polyfits": [
            {
                "degree": item["degree"],
                "r2": item["r2"],
                "coefficients": item["coefficients"],
            }
            for item in poly_results
        ],
        "best_poly_degree": int(best_poly["degree"]),
        "best_poly_r2": float(best_poly["r2"]),
        "plot_path": str(plot_path),
    }


def write_report(results):
    lines = [
        "# Zipf Step-10 Ablation",
        "",
        "- Evaluated the step-2 and step-10 Zipf-Mandelbrot corrections on 1000 evenly spaced normalized x points in `[0.05, 0.95]`.",
        "- Fitted degree 3, 4, and 5 polynomials to `step10 - step2` and report R².",
        "",
    ]
    for result in results:
        lines.extend(
            [
                f"## {result['name']}",
                "",
                f"- step-2 formula: `{result['step2_expr']}`",
                f"- step-10 formula: `{result['step10_expr']}`",
                f"- difference range: [`{result['difference_min']:.6e}`, `{result['difference_max']:.6e}`]",
                f"- mean abs difference: `{result['difference_mean_abs']:.6e}`",
                f"- std of difference: `{result['difference_std']:.6e}`",
            ]
        )
        for item in result["polyfits"]:
            lines.append(f"- polynomial degree {item['degree']} R²: `{item['r2']:.12f}`")
        lines.extend(
            [
                f"- best polynomial degree: `{result['best_poly_degree']}`",
                f"- best polynomial R²: `{result['best_poly_r2']:.12f}`",
                f"- plot: `{result['plot_path']}`",
                "",
            ]
        )
    (OUTDIR / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    results = [run_one(entry) for entry in CORPORA]
    (OUTDIR / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    write_report(results)
    print(f"Saved {OUTDIR / 'summary.json'}")
    print(f"Saved {OUTDIR / 'report.md'}")


if __name__ == "__main__":
    main()
