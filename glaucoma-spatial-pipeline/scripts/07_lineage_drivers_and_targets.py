#!/usr/bin/env python3
"""Step 7: Lineage drivers and therapeutic targets (stub for test pipeline)."""
import argparse
import os

import numpy as np
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(description="Lineage drivers and therapeutic targets")
    parser.add_argument("--adata", required=True, help="Input fate h5ad")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--no_opentargets", action="store_true",
                        help="Skip OpenTargets lookup")
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading {args.adata}")
    adata = sc.read_h5ad(args.adata)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    if args.no_opentargets:
        print("Skipping OpenTargets lookup")

    # Generate mock therapeutic targets CSV
    rng = np.random.default_rng(42)
    n_targets = min(20, adata.n_vars)
    gene_names = adata.var_names[:n_targets].tolist()

    lines = ["gene,score,fate,druggable,mechanism"]
    for gene in gene_names:
        score = round(float(rng.uniform(0.5, 1.0)), 3)
        fate = rng.choice(["RGC_death", "RGC_survival", "glia_activation"])
        druggable = rng.choice(["yes", "no"])
        mechanism = rng.choice(["kinase", "receptor", "transcription_factor", "channel"])
        lines.append(f"{gene},{score},{fate},{druggable},{mechanism}")

    targets_path = os.path.join(args.outdir, "therapeutic_targets.csv")
    with open(targets_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Therapeutic targets: {targets_path} ({n_targets} targets)")


if __name__ == "__main__":
    main()
