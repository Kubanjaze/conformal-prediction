import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse, os, warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
RDLogger.DisableLog("rdApp.*")

FAMILY_COLORS = {"benz": "#4C72B0", "naph": "#DD8452", "ind": "#55A868",
                 "quin": "#C44E52", "pyr": "#8172B2", "bzim": "#937860", "other": "#808080"}

def load_compounds(path):
    df = pd.read_csv(path)
    records, n_bad = [], 0
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        if mol is None: n_bad += 1; continue
        try:
            pic50 = float(row["pic50"])
        except (KeyError, ValueError):
            continue
        if np.isnan(pic50): continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, useChirality=True)
        fam = str(row["compound_name"]).split("_")[0]
        records.append({"compound_name": str(row["compound_name"]),
                        "family": fam if fam in FAMILY_COLORS else "other",
                        "pic50": pic50, "fp": list(fp)})
    print(f"  {len(records)} valid ({n_bad} skipped)")
    return pd.DataFrame(records)

def conformal_coverage(X, y, coverage, n_reps=100, seed=42):
    """Monte Carlo estimate of empirical coverage across many random splits."""
    rng = np.random.default_rng(seed)
    n = len(y)
    n_cal = max(5, int(n * 0.30))
    intervals, coverages = [], []

    for _ in range(n_reps):
        idx = rng.permutation(n)
        cal_idx = idx[:n_cal]
        tr_idx = idx[n_cal:]
        te_idx = rng.choice(n, size=10, replace=False)  # random test set for each rep

        rf = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
        rf.fit(X[tr_idx], y[tr_idx])
        resids = np.abs(y[cal_idx] - rf.predict(X[cal_idx]))
        q_level = min(1.0, (n_cal * coverage + 1) / n_cal)
        q = np.quantile(resids, q_level)

        y_hat = rf.predict(X[te_idx])
        covered = np.mean(np.abs(y[te_idx] - y_hat) <= q)
        coverages.append(covered)
        intervals.append(2 * q)

    return np.array(coverages), np.array(intervals)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--coverage", type=float, default=0.90, help="Target coverage level (0-1)")
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nLoading: {args.input}")
    df = load_compounds(args.input)
    X = np.array(df["fp"].tolist(), dtype=float)
    y = df["pic50"].values

    coverages_list = []
    for cov in [0.80, 0.90, 0.95]:
        print(f"Estimating coverage at target={cov}...")
        cov_arr, int_arr = conformal_coverage(X, y, cov)
        coverages_list.append({
            "target_coverage": cov,
            "empirical_coverage_mean": round(cov_arr.mean(), 4),
            "empirical_coverage_std": round(cov_arr.std(), 4),
            "mean_interval_width": round(int_arr.mean(), 4),
        })

    res_df = pd.DataFrame(coverages_list)
    res_df.to_csv(os.path.join(args.output_dir, "conformal_results.csv"), index=False)
    print(f"Saved: {args.output_dir}/conformal_results.csv")

    # Plot: target vs empirical coverage
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    targets = res_df["target_coverage"].values
    empirical = res_df["empirical_coverage_mean"].values
    err = res_df["empirical_coverage_std"].values
    widths = res_df["mean_interval_width"].values

    ax1.errorbar(targets, empirical, yerr=err, fmt="o-", color="#4C72B0", lw=2, capsize=5, markersize=8)
    ax1.plot([0.7, 1.0], [0.7, 1.0], "k--", lw=1, label="Perfect calibration")
    ax1.set_xlabel("Target Coverage", fontsize=11); ax1.set_ylabel("Empirical Coverage", fontsize=11)
    ax1.set_title("Conformal Coverage Calibration", fontsize=12, fontweight="bold")
    ax1.legend(fontsize=9); ax1.set_xlim(0.75, 1.0); ax1.set_ylim(0.75, 1.0)
    ax1.spines["top"].set_visible(False); ax1.spines["right"].set_visible(False)

    ax2.bar([f"{t:.0%}" for t in targets], widths, color="#DD8452", edgecolor="white")
    for i, w in enumerate(widths):
        ax2.text(i, w+0.01, f"±{w/2:.2f}", ha="center", fontsize=9)
    ax2.set_xlabel("Target Coverage", fontsize=11); ax2.set_ylabel("Mean Interval Width (pIC50)", fontsize=11)
    ax2.set_title("Prediction Interval Width", fontsize=12, fontweight="bold")
    ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
    plt.suptitle("Split Conformal Prediction", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "coverage_plot.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {args.output_dir}/coverage_plot.png")

    print(f"\n--- Conformal Prediction Results (100 random splits) ---")
    print(f"  {'Target':>7}  {'Empirical':>10}  {'Interval Width':>14}")
    for _, row in res_df.iterrows():
        print(f"  {row['target_coverage']:>7.0%}  {row['empirical_coverage_mean']:>10.3f}±{row['empirical_coverage_std']:.3f}  {row['mean_interval_width']:>14.3f}")
    print("\nDone.")

if __name__ == "__main__":
    main()
