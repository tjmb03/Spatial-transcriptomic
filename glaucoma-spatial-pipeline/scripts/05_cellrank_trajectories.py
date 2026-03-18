#!/usr/bin/env python3
"""Step 5: CellRank fate probabilities (stub for test pipeline)."""
import argparse
import json
import os

import numpy as np
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(description="CellRank trajectory inference")
    parser.add_argument("--adata", required=True, help="Input velocity h5ad")
    parser.add_argument("--outdir", required=True, help="Output directory")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading {args.adata}")
    adata = sc.read_h5ad(args.adata)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Add mock fate probabilities
    rng = np.random.default_rng(42)
    n_fates = 3
    fate_probs = rng.dirichlet(np.ones(n_fates), size=adata.n_obs).astype(np.float32)
    for i in range(n_fates):
        adata.obs[f"fate_{i}"] = fate_probs[:, i]

    # Add mock pseudotime
    adata.obs["pseudotime"] = rng.uniform(0, 1, adata.n_obs).astype(np.float32)

    # Save fate adata
    out_path = os.path.join(args.outdir, "fate_adata.h5ad")
    adata.write_h5ad(out_path)
    print(f"Saved {out_path}")

    # Save kernel weights
    weights = {
        "script": "05_cellrank_trajectories.py",
        "mode": "stub",
        "n_cells": int(adata.n_obs),
        "n_fates": n_fates,
        "kernel_weights": {"VelocityKernel": 0.8, "ConnectivityKernel": 0.2},
    }
    weights_path = os.path.join(args.outdir, "kernel_weights.json")
    with open(weights_path, "w") as f:
        json.dump(weights, f, indent=2)
    print(f"Kernel weights: {weights_path}")


if __name__ == "__main__":
    main()
