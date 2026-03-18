#!/usr/bin/env python3
"""Step 4: RNA velocity analysis (stub for test pipeline)."""
import argparse
import json
import os

import anndata as ad
import numpy as np
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(description="scVelo RNA velocity")
    parser.add_argument("--adata", required=True, help="Input integrated h5ad")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--n_jobs", type=int, default=1, help="Number of parallel jobs")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading {args.adata}")
    adata = sc.read_h5ad(args.adata)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Add mock velocity fields for downstream compatibility
    n_obs, n_vars = adata.n_obs, adata.n_vars
    rng = np.random.default_rng(42)
    adata.layers["velocity"] = rng.normal(0, 0.1, (n_obs, n_vars)).astype(np.float32)
    adata.layers["spliced"] = np.abs(rng.normal(5, 2, (n_obs, n_vars))).astype(np.float32)
    adata.layers["unspliced"] = np.abs(rng.normal(1, 0.5, (n_obs, n_vars))).astype(np.float32)

    # Add velocity confidence
    adata.obs["velocity_confidence"] = rng.uniform(0.5, 1.0, n_obs).astype(np.float32)

    # Save velocity adata
    out_path = os.path.join(args.outdir, "velocity_adata.h5ad")
    adata.write_h5ad(out_path)
    print(f"Saved {out_path}")

    # Save QC report
    qc = {
        "script": "04_scvelo_velocity.py",
        "mode": "stub",
        "n_cells": int(n_obs),
        "n_genes": int(n_vars),
        "mean_velocity_confidence": float(adata.obs["velocity_confidence"].mean()),
    }
    qc_path = os.path.join(args.outdir, "velocity_qc.json")
    with open(qc_path, "w") as f:
        json.dump(qc, f, indent=2)
    print(f"QC report: {qc_path}")


if __name__ == "__main__":
    main()
