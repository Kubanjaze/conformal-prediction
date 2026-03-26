# conformal-prediction — Phase 48

Split conformal prediction for pIC50 regression. Produces guaranteed prediction intervals at 80%, 90%, and 95% coverage targets using 100 Monte Carlo random splits.

## Usage

```bash
PYTHONUTF8=1 python main.py --input data/compounds.csv
```

## Outputs

- `output/conformal_results.csv` — empirical coverage and interval widths per target level
- `output/coverage_plot.png` — calibration curve + interval width bar chart
