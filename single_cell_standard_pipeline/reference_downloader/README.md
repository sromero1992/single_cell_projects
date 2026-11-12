# reference_downloader

Standalone downloader for the reference databases the scRNA-seq pipeline uses.
**No pipeline or agent code is required** — copy this folder anywhere and run it.

It fetches two independent things:

| What | Output | Size | Used by |
|---|---|---|---|
| Enrichr cell-type gene-set libraries | `enrichr_gene_sets.json` | ~tens of MB | annotation + pathway enrichment |
| SCENIC / cisTarget databases | `scenic_refs/` (`*.feather` + motif `.tbl` + `allTFs_*.txt`) | ~1–2 GB | SCENIC / GRN stage |

## Install

```bash
pip install -r requirements.txt
```

(`requests` is the only third-party dependency; SCENIC downloads use the standard library.)

## Usage

```bash
# everything: Enrichr + mouse SCENIC (10kb database)
python download_references.py --all

# just the Enrichr libraries (fast)
python download_references.py --enrichr --out-dir .

# just SCENIC, smaller 500bp DB, into a chosen folder, with checksum verification
python download_references.py --scenic --scenic-db 500bp --scenic-dir ~/scenic_refs --verify

# human genome instead of mouse
python download_references.py --scenic --species human
```

### Flags
- `--all` — Enrichr + SCENIC
- `--enrichr` / `--scenic` — pick one or both
- `--out-dir` — where `enrichr_gene_sets.json` is written (default: current folder)
- `--scenic-dir` — SCENIC output folder (default: `scenic_refs`)
- `--species` — `mouse` (default) or `human`
- `--scenic-db` — `10kb` (standard) or `500bp` (smaller/faster)
- `--verify` — SHA-256 check SCENIC files against the Aerts lab checksum list

Already-downloaded files of the correct size are skipped, so re-running resumes safely.

## Wiring the downloads into the pipeline

- **Enrichr:** place `enrichr_gene_sets.json` next to `discover.py` (the repo root).
  `LibraryScoreAgent` and `DownstreamAgent` load it automatically.
- **SCENIC:** point the pipeline at the folder you downloaded into, e.g.
  `~/scenic_refs`, via `config.extra['scenic_refs_dir']`, and run **without**
  `--skip-scenic`. The SCENIC stage looks for `*.genes_vs_motifs.rankings.feather`
  there, plus the motif `.tbl` and the TF list.

## Source

SCENIC/cisTarget databases come from the Aerts lab:
<https://resources.aertslab.org/cistarget/databases/> (v10 `mc_v10_clust`).
If a URL 404s (they occasionally reorganise paths), browse that directory and update
`SCENIC_URLS` at the top of `download_references.py`.

Enrichr libraries come from <https://maayanlab.cloud/Enrichr/>.
