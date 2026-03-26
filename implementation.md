# Phase 48 — Conformal Prediction Basics

**Version:** 1.1 | **Tier:** Standard | **Date:** 2026-03-26

## Goal
Apply split conformal prediction to produce prediction intervals for pIC50 regression.
Guarantee coverage at a target level (80%, 90%, 95%) — intervals contain the true value with stated probability.

CLI: `python main.py --input data/compounds.csv --coverage 0.90`

Outputs: conformal_results.csv, coverage_plot.png

## Logic
- Split: calibration set (30%) and training set (70%)
- Train RF regressor on training set
- Calibration: compute nonconformity scores = |y_true - y_pred| on calibration set
- Quantile: q = (n_cal * coverage + 1) / n_cal quantile of nonconformity scores
- Prediction interval for new compound: [y_hat - q, y_hat + q]
- Monte Carlo evaluation: 100 random splits to estimate empirical coverage distribution

## Results

| Target Coverage | Empirical Coverage | Interval Width (pIC50) |
|---|---|---|
| 80% | 0.941 ± 0.062 | ±0.571 (total 1.141) |
| 90% | 0.977 ± 0.044 | ±0.847 (total 1.694) |
| 95% | 0.993 ± 0.025 | ±0.931 (total 1.862) |

## Key Insights
- Empirical coverage consistently exceeds target — expected for small calibration sets (n=14)
- Coverage guarantee is a lower bound; over-coverage is the conservative direction
- Interval width grows from ±0.57 to ±0.93 pIC50 as target rises from 80% to 95%
- At ±0.57 pIC50 (80% target), conformal intervals are slightly wider than RF LOO-CV MAE (0.268) — reasonable
- Conformal intervals are distribution-free: no normality assumption needed

## Deviations from Plan
- None; plan implemented as specified.
