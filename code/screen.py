#!/usr/bin/env python3
"""Genome-wide genotype-dependent trajectory screen for GSE248369.

The screen uses the trajectory and raw counts already produced by the original
notebook. Genes are filtered only by detection (>=20 trajectory cells), then
fit with the same six-knot condition-aware tradeSeq NB-GAM. Statistical tests
are lineage-specific and pairwise between genotypes. Effect magnitude is the
integrated absolute log2 fold-change over a uniform pseudotime grid.

Cells are not treated as biological replicates: this model is a candidate
screen, and embryo-level summaries are retained for robustness assessment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import random
import socket
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[1]
INPUT = Path(
    os.environ.get(
        "BIOBABEL_TRAJECTORY_INPUT",
        REPO_ROOT / "data" / "traj_ab.h5ad",
    )
).expanduser().resolve()
RAW_INPUT = Path(
    os.environ.get(
        "BIOBABEL_ANALYSIS_INPUT",
        REPO_ROOT / "data" / "analysis_input.h5ad",
    )
).expanduser().resolve()
ORIGINAL_NOTEBOOK = (
    REPO_ROOT / "agent.ipynb"
)
BEAM_RESULTS = REPO_ROOT / "tables" / "beam.csv.gz"

SEED = 2024
N_KNOTS = 6
N_POINTS = 100
MIN_CELLS_EXPRESSED = 20
L2FC_THRESHOLD = 0.5
DEFAULT_CHUNK_SIZE = 500

CONDITION_LEVELS = ("control", "dEndo", "DKO")
GENOTYPE_MAP = {
    "control": "control",
    "Brg1Δendo;Brm+/-": "dEndo",
    "DKO": "DKO",
}
LINEAGES = {1: "alpha", 2: "beta"}
CONTRASTS = (
    (1, 2, "dEndo-control", "dEndo", "control"),
    (1, 3, "DKO-control", "DKO", "control"),
    (2, 3, "DKO-dEndo", "DKO", "dEndo"),
)
ORIGINAL_MARKERS = (
    "Neurog3", "Sox4", "Fev", "Chga", "Chgb", "Neurod1", "Pax6", "Insm1",
    "Gcg", "Arx", "Irx1", "Irx2", "Ttr", "Pou3f4", "Peg3", "Isl1",
    "Ins1", "Ins2", "Pdx1", "Nkx6-1", "Pax4", "Iapp", "Mafa", "Slc2a2",
    "Ucn3", "Nnat",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv_gz(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    frame.to_csv(tmp, sep="\t", index=False, compression="gzip")
    os.replace(tmp, path)


def atomic_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, path)


def technical_class(symbol: str) -> str:
    symbol = str(symbol)
    if symbol.startswith(("mt-", "Mt-")):
        return "mitochondrial"
    if symbol.startswith(("Rpl", "Rps")):
        return "ribosomal"
    if symbol.startswith(("Hba", "Hbb")):
        return "haemoglobin"
    if symbol.startswith("Gm"):
        return "predicted_gene"
    return ""


def prepare(
    output_root: Path,
    chunk_size: int,
    universe_mode: str,
    adjust_sex: bool,
) -> None:
    sys.path.insert(0, str(SCRIPT_PATH.parent))
    import common
    from tradeseq._assign_cells import assign_cells
    from tradeseq._offset import compute_offset

    np.random.seed(SEED)
    random.seed(SEED)
    data = ad.read_h5ad(INPUT)
    info = common.assign_arms(data)
    # Freeze both quantities that otherwise depend on the matrix presented to
    # each fit. The offset is calculated once from the complete raw-count
    # matrix, exactly as in an unchunked fit. The lineage draw uses an explicit
    # Generator because fit_gam's default Generator is intentionally unseeded.
    full_offset = compute_offset(None, data.layers["counts"].T)
    frozen_assignment = assign_cells(
        data.obsm["cell_weights"], rng=np.random.default_rng(SEED)
    )

    detected = pd.to_numeric(
        data.var["num_cells_expressed"], errors="coerce"
    ).fillna(0).astype(int)
    all_genes = pd.DataFrame(
        {
            "gene_id": data.var_names.astype(str),
            "gene_symbol": data.var["gene_short_name"].astype(str).to_numpy(),
            "num_cells_expressed": detected.to_numpy(),
            "feature_type": data.var.get(
                "feature_type", pd.Series("", index=data.var_names)
            ).astype(str).to_numpy(),
            "source_index": np.arange(data.n_vars, dtype=int),
        }
    )
    if universe_mode == "all-expressed":
        keep = detected >= MIN_CELLS_EXPRESSED
        universe = all_genes.loc[keep.to_numpy()].reset_index(drop=True)
        selection_definition = (
            "All genes detected in at least 20 trajectory cells; no biological "
            "or technical class was removed before fitting."
        )
    elif universe_mode == "lineage-panel":
        symbols = data.var["gene_short_name"].astype(str).to_numpy()
        technical = np.array(
            [
                symbol.startswith(
                    ("mt-", "Rps", "Rpl", "Gm", "Rik", "Hb", "Hba", "Hbb")
                )
                for symbol in symbols
            ],
            dtype=bool,
        )
        expressed = (detected.to_numpy() >= MIN_CELLS_EXPRESSED) & (~technical)
        counts = data.layers["counts"]
        counts_dense = (
            counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
        )
        libraries = counts_dense.sum(axis=1, keepdims=True).astype(float)
        libraries[libraries == 0] = 1.0
        variance = np.where(
            expressed,
            np.log1p(counts_dense / libraries * 1e4).var(axis=0),
            -np.inf,
        )
        top_variable_ids = data.var_names[
            np.argsort(variance)[::-1][:400]
        ].astype(str).tolist()

        beam = pd.read_csv(BEAM_RESULTS, index_col=0)
        beam = beam.sort_values("qval", kind="mergesort")
        expressed_ids = set(data.var_names[expressed].astype(str))
        beam_ids = [
            str(gene_id)
            for gene_id in beam.index[:1000]
            if str(gene_id) in expressed_ids
        ]
        marker_set = set(ORIGINAL_MARKERS)
        marker_ids = data.var_names[
            data.var["gene_short_name"].astype(str).isin(marker_set)
        ].astype(str).tolist()
        selected_ids = list(
            dict.fromkeys(beam_ids + top_variable_ids + marker_ids)
        )[:1500]
        universe = (
            all_genes.set_index("gene_id")
            .loc[selected_ids]
            .reset_index()
        )
        selection_definition = (
            "The original notebook's genotype-blind main tradeSeq panel: "
            "BEAM top 1000 (expression/technical filtered), top 400 variable "
            "genes under the same filter, plus the predeclared marker panel."
        )
    else:
        raise ValueError("Unsupported universe mode: {}".format(universe_mode))
    universe["technical_class"] = universe["gene_symbol"].map(technical_class)
    universe["original_marker"] = universe["gene_symbol"].isin(ORIGINAL_MARKERS)

    y_markers = ("Ddx3y", "Eif2s3y", "Kdm5d", "Uty")
    symbol_values = data.var["gene_short_name"].astype(str).to_numpy()
    y_indices = np.flatnonzero(np.isin(symbol_values, y_markers))
    y_counts = data.layers["counts"][:, y_indices]
    y_detected = np.asarray((y_counts > 0).sum(axis=1)).ravel() > 0
    sample_values = data.obs["sample"].astype(str).to_numpy()
    sex_by_sample = {}
    for sample in pd.unique(sample_values):
        detection_fraction = float(y_detected[sample_values == sample].mean())
        sex_by_sample[str(sample)] = (
            "male" if detection_fraction >= 0.10 else "female"
        )

    lineage_input = pd.DataFrame(
        {
            "cell_id": data.obs_names.astype(str),
            "pseudotime": data.obs["Pseudotime"].astype(float).to_numpy(),
            "alpha_weight": data.obsm["cell_weights"][:, 0],
            "beta_weight": data.obsm["cell_weights"][:, 1],
            "alpha_assignment": frozen_assignment[:, 0],
            "beta_assignment": frozen_assignment[:, 1],
            "tmm_log_offset": full_offset,
            "genotype": data.obs["genotype"].astype(str).to_numpy(),
            "geno_short": data.obs["genotype"].map(GENOTYPE_MAP).astype(str).to_numpy(),
            "sample": data.obs["sample"].astype(str).to_numpy(),
            "sex": [
                sex_by_sample[str(sample)]
                for sample in data.obs["sample"].astype(str)
            ],
            "annotation": data.obs["broad_annotation"].astype(str).to_numpy(),
        }
    )

    atomic_tsv_gz(universe, output_root / "gene_universe.tsv.gz")
    atomic_tsv_gz(lineage_input, output_root / "lineage_input.tsv.gz")
    n_chunks = int(math.ceil(len(universe) / chunk_size))
    atomic_json(
        {
            "seed": SEED,
            "n_knots": N_KNOTS,
            "n_points": N_POINTS,
            "min_cells_expressed": MIN_CELLS_EXPRESSED,
            "l2fc_threshold": L2FC_THRESHOLD,
            "chunk_size": chunk_size,
            "n_chunks": n_chunks,
            "n_genes": int(len(universe)),
            "n_cells": int(data.n_obs),
            "universe_mode": universe_mode,
            "selection_definition": selection_definition,
            "adjust_sex": bool(adjust_sex),
            "sex_definition": (
                "Embryo classified male when >=10% of trajectory cells "
                "detected at least one of Ddx3y/Eif2s3y/Kdm5d/Uty"
            ),
            "sex_by_sample": sex_by_sample,
            "condition_levels": list(CONDITION_LEVELS),
            "lineages": {str(k): v for k, v in LINEAGES.items()},
            "offset": "TMM log offset frozen from all 23,308 raw-count genes",
            "lineage_assignment": (
                "Single multinomial draw frozen with numpy Generator seed 2024"
            ),
            "trajectory_assignment": info,
            "input": str(INPUT),
        },
        output_root / "screen_plan.json",
    )
    print(
        "PREPARED",
        json.dumps(
            {
                "n_cells": data.n_obs,
                "n_genes": len(universe),
                "n_chunks": n_chunks,
                "chunk_size": chunk_size,
            },
            sort_keys=True,
        ),
        flush=True,
    )


def load_chunk_data(
    output_root: Path, chunk_index: int
) -> tuple[ad.AnnData, pd.DataFrame, dict, np.ndarray, np.ndarray]:
    plan = json.loads((output_root / "screen_plan.json").read_text(encoding="utf-8"))
    universe = pd.read_csv(output_root / "gene_universe.tsv.gz", sep="\t")
    lineage = pd.read_csv(output_root / "lineage_input.tsv.gz", sep="\t")

    start = chunk_index * int(plan["chunk_size"])
    stop = min(start + int(plan["chunk_size"]), len(universe))
    if start >= stop:
        raise ValueError(
            "Chunk {} is empty for {} genes and chunk size {}".format(
                chunk_index, len(universe), plan["chunk_size"]
            )
        )
    chunk_meta = universe.iloc[start:stop].copy()

    data = ad.read_h5ad(INPUT)
    if list(data.obs_names.astype(str)) != lineage["cell_id"].astype(str).tolist():
        raise RuntimeError("Cell order differs from the frozen lineage input")
    data.obsm["pseudotime"] = np.column_stack(
        [lineage["pseudotime"].to_numpy(), lineage["pseudotime"].to_numpy()]
    )
    data.obsm["cell_weights"] = lineage[
        ["alpha_weight", "beta_weight"]
    ].to_numpy(dtype=float)
    data.obs["geno_short"] = pd.Categorical(
        lineage["geno_short"], categories=CONDITION_LEVELS, ordered=True
    )
    offset = lineage["tmm_log_offset"].to_numpy(dtype=float)
    frozen_assignment = lineage[
        ["alpha_assignment", "beta_assignment"]
    ].to_numpy(dtype=int)

    gene_ids = chunk_meta["gene_id"].astype(str).tolist()
    data = data[:, gene_ids].copy()
    return data, chunk_meta, plan, offset, frozen_assignment


def effect_table(
    fitted: ad.AnnData,
    fitted_gene_ids: list[str],
    gene_meta: pd.DataFrame,
) -> pd.DataFrame:
    import tradeseq as ts

    if not fitted_gene_ids:
        return pd.DataFrame()
    prediction = ts.predict_smooth(
        fitted, gene=fitted_gene_ids, n_points=N_POINTS, tidy=False
    )
    metadata = gene_meta.set_index("gene_id")
    records: list[pd.DataFrame] = []

    for lineage_id, lineage_name in LINEAGES.items():
        condition_values: dict[str, np.ndarray] = {}
        for condition in CONDITION_LEVELS:
            columns = [
                "lineage{}_condition{}_point{}".format(
                    lineage_id, condition, point
                )
                for point in range(1, N_POINTS + 1)
            ]
            condition_values[condition] = prediction.loc[
                fitted_gene_ids, columns
            ].to_numpy(dtype=float)

        for _, _, contrast_name, numerator, denominator in CONTRASTS:
            numerator_values = np.clip(condition_values[numerator], 1e-12, None)
            denominator_values = np.clip(condition_values[denominator], 1e-12, None)
            delta = np.log2(numerator_values / denominator_values)
            abs_delta = np.abs(delta)
            peak_index = np.argmax(abs_delta, axis=1)
            row_index = np.arange(len(fitted_gene_ids))
            peak_l2fc = delta[row_index, peak_index]
            signed_mean = delta.mean(axis=1)
            direction = np.where(
                signed_mean > 0, "up", np.where(signed_mean < 0, "down", "flat")
            )
            records.append(
                pd.DataFrame(
                    {
                        "gene_id": fitted_gene_ids,
                        "gene_symbol": metadata.loc[
                            fitted_gene_ids, "gene_symbol"
                        ].to_numpy(),
                        "lineage": lineage_name,
                        "contrast": contrast_name,
                        "iae_log2fc": abs_delta.mean(axis=1),
                        "rms_log2fc": np.sqrt(np.mean(delta ** 2, axis=1)),
                        "max_abs_log2fc": abs_delta.max(axis=1),
                        "mean_signed_log2fc": signed_mean,
                        "endpoint_log2fc": delta[:, -1],
                        "fraction_abs_log2fc_ge_0_5": (
                            abs_delta >= L2FC_THRESHOLD
                        ).mean(axis=1),
                        "peak_log2fc": peak_l2fc,
                        "peak_time_fraction": peak_index / float(N_POINTS - 1),
                        "direction_by_integral": direction,
                    }
                )
            )
    return pd.concat(records, ignore_index=True)


def fit_chunk(output_root: Path, chunk_index: int, n_jobs: int) -> None:
    import tradeseq as ts

    started = time.time()
    np.random.seed(SEED)
    random.seed(SEED)
    data, gene_meta, plan, offset, frozen_assignment = load_chunk_data(
        output_root, chunk_index
    )
    gene_ids = data.var_names.astype(str).tolist()

    fixed_effects = None
    if bool(plan.get("adjust_sex", False)):
        lineage = pd.read_csv(output_root / "lineage_input.tsv.gz", sep="\t")
        male = lineage["sex"].astype(str).eq("male").to_numpy(dtype=float)
        fixed_effects = np.column_stack(
            [np.ones(len(male), dtype=float), male]
        )

    ts.fit_gam(
        data,
        genes=gene_ids,
        n_knots=N_KNOTS,
        conditions_key="geno_short",
        parallel=n_jobs > 1,
        n_jobs=n_jobs,
        verbose=False,
        offset=offset,
        U=fixed_effects,
        _w_samp=frozen_assignment,
    )
    converged = data.var["tradeseq_converged"].fillna(False).astype(bool)
    fitted_gene_ids = data.var_names[converged].astype(str).tolist()

    test_zero = ts.condition_test(
        data,
        global_=False,
        pairwise=True,
        lineages=True,
        l2fc=0.0,
    )
    test_threshold = ts.condition_test(
        data,
        global_=False,
        pairwise=True,
        lineages=True,
        l2fc=L2FC_THRESHOLD,
    )
    meta = gene_meta.set_index("gene_id")
    test_records: list[pd.DataFrame] = []
    for lineage_id, lineage_name in LINEAGES.items():
        for condition_a, condition_b, contrast_name, _, _ in CONTRASTS:
            tag = "lineage{}_conds{}vs{}".format(
                lineage_id, condition_a, condition_b
            )
            test_records.append(
                pd.DataFrame(
                    {
                        "gene_id": gene_ids,
                        "gene_symbol": meta.loc[gene_ids, "gene_symbol"].to_numpy(),
                        "num_cells_expressed": meta.loc[
                            gene_ids, "num_cells_expressed"
                        ].to_numpy(),
                        "technical_class": meta.loc[
                            gene_ids, "technical_class"
                        ].fillna("").to_numpy(),
                        "original_marker": meta.loc[
                            gene_ids, "original_marker"
                        ].astype(bool).to_numpy(),
                        "lineage": lineage_name,
                        "contrast": contrast_name,
                        "converged": converged.reindex(gene_ids).to_numpy(),
                        "wald_stat": test_zero.loc[
                            gene_ids, "waldStat_" + tag
                        ].to_numpy(),
                        "df": test_zero.loc[
                            gene_ids, "df_" + tag
                        ].to_numpy(),
                        "pvalue": test_zero.loc[
                            gene_ids, "pvalue_" + tag
                        ].to_numpy(),
                        "wald_stat_l2fc_0_5": test_threshold.loc[
                            gene_ids, "waldStat_" + tag
                        ].to_numpy(),
                        "df_l2fc_0_5": test_threshold.loc[
                            gene_ids, "df_" + tag
                        ].to_numpy(),
                        "pvalue_l2fc_0_5": test_threshold.loc[
                            gene_ids, "pvalue_" + tag
                        ].to_numpy(),
                    }
                )
            )
    tests = pd.concat(test_records, ignore_index=True)
    effects = effect_table(data, fitted_gene_ids, gene_meta)
    combined = tests.merge(
        effects,
        how="left",
        on=["gene_id", "gene_symbol", "lineage", "contrast"],
        validate="one_to_one",
    )

    chunk_path = output_root / "chunks" / "chunk_{:04d}.tsv.gz".format(chunk_index)
    atomic_tsv_gz(combined, chunk_path)
    atomic_json(
        {
            "chunk_index": chunk_index,
            "chunk_size": int(len(gene_ids)),
            "first_gene_id": gene_ids[0],
            "last_gene_id": gene_ids[-1],
            "n_converged": int(converged.sum()),
            "n_failed": int((~converged).sum()),
            "elapsed_seconds": time.time() - started,
            "host": socket.gethostname(),
            "pid": os.getpid(),
            "n_jobs": n_jobs,
            "seed": SEED,
            "script_sha256": sha256_file(SCRIPT_PATH),
            "screen_plan": plan,
        },
        output_root / "chunks" / "chunk_{:04d}.json".format(chunk_index),
    )
    print(
        "CHUNK_DONE",
        json.dumps(
            {
                "chunk": chunk_index,
                "genes": len(gene_ids),
                "converged": int(converged.sum()),
                "seconds": round(time.time() - started, 2),
            },
            sort_keys=True,
        ),
        flush=True,
    )


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    finite_index = np.flatnonzero(np.isfinite(p))
    if finite_index.size == 0:
        return q
    finite_p = p[finite_index]
    order = np.argsort(finite_p)
    ranked = finite_p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    q[finite_index[order]] = adjusted
    return q


def git_metadata(path: Path) -> dict:
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
        ).strip()
        dirty = bool(
            subprocess.check_output(
                ["git", "-C", str(path), "status", "--porcelain"], text=True
            ).strip()
        )
        return {"path": str(path), "commit": commit, "dirty": dirty}
    except Exception as exc:
        return {"path": str(path), "error": repr(exc)}


def combine(output_root: Path) -> None:
    import monocle2py
    import tradeseq

    plan = json.loads((output_root / "screen_plan.json").read_text(encoding="utf-8"))
    expected = int(plan["n_chunks"])
    chunk_paths = [
        output_root / "chunks" / "chunk_{:04d}.tsv.gz".format(index)
        for index in range(expected)
    ]
    missing = [str(path) for path in chunk_paths if not path.exists()]
    if missing:
        raise RuntimeError("Missing {} chunk files: {}".format(len(missing), missing[:8]))

    table = pd.concat(
        [pd.read_csv(path, sep="\t") for path in chunk_paths],
        ignore_index=True,
    )
    expected_rows = int(plan["n_genes"]) * len(LINEAGES) * len(CONTRASTS)
    if len(table) != expected_rows:
        raise RuntimeError(
            "Expected {} combined rows, observed {}".format(expected_rows, len(table))
        )
    if table.duplicated(["gene_id", "lineage", "contrast"]).any():
        raise RuntimeError("Duplicate gene-lineage-contrast rows detected")

    table["qvalue"] = np.nan
    table["qvalue_l2fc_0_5"] = np.nan
    for (_, _), index in table.groupby(["lineage", "contrast"]).groups.items():
        table.loc[index, "qvalue"] = bh_adjust(table.loc[index, "pvalue"])
        table.loc[index, "qvalue_l2fc_0_5"] = bh_adjust(
            table.loc[index, "pvalue_l2fc_0_5"]
        )

    table["significant_q05"] = table["qvalue"] < 0.05
    table["significant_l2fc_0_5_q05"] = table["qvalue_l2fc_0_5"] < 0.05
    table["evidence_rank"] = table.groupby(
        ["lineage", "contrast"]
    )["pvalue"].rank(method="min", ascending=True)
    table["integrated_effect_rank"] = table.groupby(
        ["lineage", "contrast"]
    )["iae_log2fc"].rank(method="min", ascending=False)
    table["integrated_effect_rank_among_q05"] = np.nan
    significant = table["significant_q05"] & table["iae_log2fc"].notna()
    table.loc[significant, "integrated_effect_rank_among_q05"] = (
        table.loc[significant]
        .groupby(["lineage", "contrast"])["iae_log2fc"]
        .rank(method="min", ascending=False)
    )

    table = table.sort_values(
        ["lineage", "contrast", "qvalue", "iae_log2fc"],
        ascending=[True, True, True, False],
    ).reset_index(drop=True)
    atomic_tsv_gz(table, output_root / "all_genes_lineage_contrasts.tsv.gz")

    control = table[table["contrast"].isin(["dEndo-control", "DKO-control"])].copy()
    atomic_tsv_gz(control, output_root / "control_contrasts.tsv.gz")

    selected = control[
        control["significant_q05"]
        & control["iae_log2fc"].notna()
        & control["technical_class"].fillna("").eq("")
    ].copy()
    selected = (
        selected.sort_values(
            ["lineage", "contrast", "iae_log2fc", "qvalue"],
            ascending=[True, True, False, True],
        )
        .groupby(["lineage", "contrast"], group_keys=False)
        .head(50)
        .reset_index(drop=True)
    )
    atomic_tsv_gz(
        selected, output_root / "top50_integrated_effect_non_technical.tsv.gz"
    )

    markers = control[control["gene_symbol"].isin(ORIGINAL_MARKERS)].copy()
    atomic_tsv_gz(markers, output_root / "original_marker_results.tsv.gz")

    summary = (
        table.groupby(["lineage", "contrast"], as_index=False)
        .agg(
            n_genes=("gene_id", "size"),
            n_converged=("converged", "sum"),
            n_finite_p=("pvalue", lambda x: int(np.isfinite(x).sum())),
            n_q05=("significant_q05", "sum"),
            n_l2fc_0_5_q05=("significant_l2fc_0_5_q05", "sum"),
            median_iae_q05=(
                "iae_log2fc",
                lambda x: float(
                    np.nanmedian(
                        x[
                            table.loc[x.index, "significant_q05"].to_numpy()
                        ].to_numpy()
                    )
                )
                if table.loc[x.index, "significant_q05"].any()
                else np.nan,
            ),
        )
    )
    atomic_tsv_gz(summary, output_root / "screen_summary.tsv.gz")

    tradeseq_repo = Path(tradeseq.__file__).resolve().parents[1]
    monocle_repo = Path(monocle2py.__file__).resolve().parents[1]
    manifest = {
        "completed_at_unix": time.time(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version,
        "parameters": plan,
        "multiple_testing": (
            "Benjamini-Hochberg within each lineage x genotype contrast "
            "across all finite fitted-gene p-values"
        ),
        "primary_effect": (
            "mean absolute fitted log2 fold-change across 100 uniformly spaced "
            "pseudotime points (normalized integrated absolute effect)"
        ),
        "input_files": {
            str(INPUT): sha256_file(INPUT),
            str(RAW_INPUT): sha256_file(RAW_INPUT),
            str(ORIGINAL_NOTEBOOK): sha256_file(ORIGINAL_NOTEBOOK),
            str(SCRIPT_PATH): sha256_file(SCRIPT_PATH),
        },
        "packages": {
            "tradeseq_version": getattr(tradeseq, "__version__", None),
            "tradeseq_source": str(Path(tradeseq.__file__).resolve()),
            "monocle2py_version": getattr(monocle2py, "__version__", None),
            "monocle2py_source": str(Path(monocle2py.__file__).resolve()),
            "anndata_version": getattr(ad, "__version__", None),
            "numpy_version": np.__version__,
            "pandas_version": pd.__version__,
        },
        "source_repositories": {
            "tradeseq": git_metadata(tradeseq_repo),
            "monocle2py": git_metadata(monocle_repo),
        },
        "outputs": {
            "all_genes": "all_genes_lineage_contrasts.tsv.gz",
            "control_contrasts": "control_contrasts.tsv.gz",
            "top50": "top50_integrated_effect_non_technical.tsv.gz",
            "markers": "original_marker_results.tsv.gz",
            "summary": "screen_summary.tsv.gz",
        },
    }
    atomic_json(manifest, output_root / "manifest.json")
    print(
        "COMBINE_DONE",
        json.dumps(
            {
                "rows": len(table),
                "genes": int(plan["n_genes"]),
                "chunks": expected,
                "selected_top_rows": len(selected),
            },
            sort_keys=True,
        ),
        flush=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "fit-chunk", "combine"))
    parser.add_argument(
        "--output-root",
        type=Path,
        default=REPO_ROOT.parent / "bio-babel-pancreas-screen",
    )
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    parser.add_argument(
        "--universe-mode",
        choices=("all-expressed", "lineage-panel"),
        default="lineage-panel",
    )
    parser.add_argument(
        "--adjust-sex",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include an embryo-sex fixed effect in the condition-aware GAM.",
    )
    parser.add_argument("--chunk-index", type=int)
    parser.add_argument("--n-jobs", type=int, default=1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_root = args.output_root.resolve()
    if args.mode == "prepare":
        prepare(
            args.output_root,
            args.chunk_size,
            args.universe_mode,
            args.adjust_sex,
        )
    elif args.mode == "fit-chunk":
        if args.chunk_index is None:
            raise ValueError("--chunk-index is required for fit-chunk")
        fit_chunk(args.output_root, args.chunk_index, args.n_jobs)
    else:
        combine(args.output_root)


if __name__ == "__main__":
    main()
