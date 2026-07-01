# TimeSCape_py — Python Package

Circadian rhythm detection in single-cell RNA-seq data.  
Python implementation of **TimeSCape v0.2**, mirroring the MATLAB v0.2 and R (`TimeSCape_R`) pipelines.

**Author:** Selim Romero · Texas A&M University · ssromerogon@tamu.edu

---

## Installation

```bash
# Option A — conda environment (recommended)
conda env create -f environment.yml
conda activate timescape

# Option B — pip
pip install -e .
```

---

## Quick start

```python
import scanpy as sc
from timescape import run_timescape, build_tmeta

adata = sc.read_h5ad("my_data.h5ad")

# Auto-build ZT metadata from labels like 'ZT00', 'ZT03', ...
tmeta = build_tmeta(sorted(adata.obs["ZT_time"].unique()))

T1, T2 = run_timescape(
    adata        = adata,
    tmeta        = tmeta,
    celltype_col = "cell_type",
    zt_col       = "ZT_time",
    period       = 24,
    norm_str     = "lib_size",
    outdir       = "./timescape_output",
    n_jobs       = -1,
)
```

See `run_timescape_demo.ipynb` for the full annotated workflow.

---

## Package structure

```
TimeSCape_py/
├── timescape/
│   ├── __init__.py      — public API exports
│   ├── core.py          — estimate_phase_r() — cosinor + F-test + Pearson
│   ├── pipeline.py      — run_timescape() — full per-cell-type pipeline
│   ├── normalize.py     — normalize_lib_size(), normalize_adata()
│   ├── visualize.py     — generate_heatmap(), plot_gene_single(), save_batch_plots()
│   └── utils.py         — bh_adjust(), wrap_acrophase(), build_tmeta()
├── run_timescape_demo.ipynb   — demo notebook
├── pyproject.toml
├── environment.yml
└── README.md
```

---

## Output files (per cell type)

| File | Contents |
|------|----------|
| `*_circadian_analysis_all.csv` | All genes: Acrophase, Amplitude, Mesor, Period, p-values, ρ |
| `*_circadian_analysis_confident.csv` | Genes passing both tests at p < 0.05 |
| `*_circadian_ZTs_mean.csv` | Per-ZT mean expression (all genes) |
| `*_circadian_ZTs_mean_normalized.csv` | ZT-normalized means |
| `*_circadian_ZTs_mean_confident.csv` | Per-ZT means (confident genes) |
| `*_circadian_ZTs_mean_normalized_confident.csv` | ZT-normalized means (confident) |
| `summary_results.csv` | Cross-cell-type summary |
| `heatmap_strict.png` | Z-score heatmap |

---

## Cross-platform equivalence

| Aspect | MATLAB v0.2 | R (TimeSCape_R) | Python (this) |
|--------|-------------|-----------------|---------------|
| NLS solver | Trust-Region (`fit()`) | Levenberg-Marquardt (`nlsLM`) | Trust-Region Reflective (`curve_fit`) |
| F-test | d1=2, d2=N−3 | identical | identical |
| BH correction | custom | `p.adjust("BH")` | custom (same formula) |
| Parallelism | `parfor` | `future_lapply` | `joblib` |
| GUI | MATLAB App | Shiny (`app.R`) | *planned* |

Statistical results differ by < 0.1 % across all platforms.
