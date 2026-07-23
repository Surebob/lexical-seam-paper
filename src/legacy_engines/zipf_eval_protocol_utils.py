from __future__ import annotations

import json
from pathlib import Path


def load_rows_map(path: Path) -> dict[str, dict]:
    rows = json.loads(path.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def best_bundle(row: dict, key: str = "best_lambda") -> dict:
    return row[key]


def has_splitfit_schema(bundle: dict) -> bool:
    return "selection" in bundle and "full_refit" in bundle


def selection_record(bundle: dict) -> dict:
    if has_splitfit_schema(bundle):
        return bundle["selection"]
    return bundle


def full_refit_record(bundle: dict) -> dict:
    if has_splitfit_schema(bundle):
        return bundle["full_refit"]
    return bundle


def split_fit_params(row: dict, key: str = "best_lambda") -> dict:
    return selection_record(best_bundle(row, key))["params"]


def split_fit_test_avg_nll(row: dict, key: str = "best_lambda") -> float:
    return float(selection_record(best_bundle(row, key))["test_avg_nll"])


def split_fit_test_loglike(row: dict, key: str = "best_lambda") -> float | None:
    record = selection_record(best_bundle(row, key))
    value = record.get("test_loglike")
    return None if value is None else float(value)


def split_fit_train_avg_nll(row: dict, key: str = "best_lambda") -> float | None:
    record = selection_record(best_bundle(row, key))
    value = record.get("train_avg_nll")
    return None if value is None else float(value)


def selection_lambda(row: dict, key: str = "best_lambda", lambda_field: str = "lambda") -> float:
    return float(best_bundle(row, key)[lambda_field])


def full_refit_params(row: dict, key: str = "best_lambda") -> dict:
    return full_refit_record(best_bundle(row, key))["params"]


def full_refit_bic(row: dict, key: str = "best_lambda") -> float | None:
    record = full_refit_record(best_bundle(row, key))
    value = record.get("bic")
    return None if value is None else float(value)


def full_refit_rmse(row: dict, key: str = "best_lambda") -> float | None:
    record = full_refit_record(best_bundle(row, key))
    value = record.get("rmse")
    return None if value is None else float(value)


def full_refit_step2_gain(row: dict, key: str = "best_lambda") -> float | None:
    record = full_refit_record(best_bundle(row, key))
    value = record.get("step2_gain")
    return None if value is None else float(value)


def full_refit_step2_helps(row: dict, key: str = "best_lambda") -> bool | None:
    record = full_refit_record(best_bundle(row, key))
    value = record.get("step2_helps")
    return None if value is None else bool(value)


def transition_fraction_from_params(params: dict) -> float | None:
    value = params.get("transition_fraction")
    return None if value is None else float(value)
