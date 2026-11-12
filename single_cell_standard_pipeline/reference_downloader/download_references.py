#!/usr/bin/env python3
"""
download_references.py  --  standalone reference-database downloader
===================================================================

NO pipeline / agent dependencies. Copy this folder anywhere and run it to fetch the
reference data the scRNA-seq pipeline needs:

  1. Enrichr cell-type gene-set libraries  ->  enrichr_gene_sets.json
     (used for annotation + pathway enrichment)

  2. SCENIC / cisTarget databases          ->  scenic_refs/
       <genome>_..._v10_clust.genes_vs_motifs.rankings.feather   (ranking DB, ~1 GB)
       motifs-v10nr_clust-nr.<org>-m0.001-o0.0.tbl               (motif -> TF annotation)
       allTFs_<genome>.txt                                       (transcription-factor list)
     (used by the SCENIC / GRN stage; several GB total)

Examples
--------
    pip install -r requirements.txt

    # everything (Enrichr + mouse SCENIC 10kb database)
    python download_references.py --all

    # just the Enrichr libraries
    python download_references.py --enrichr

    # just SCENIC, smaller 500bp database, into a chosen folder
    python download_references.py --scenic --scenic-db 500bp --scenic-dir ~/scenic_refs

    # human instead of mouse
    python download_references.py --scenic --species human

Notes
-----
* Downloads are RESUMABLE-ish: a file already present with the right size is skipped.
* Pass --verify to check SHA-256 against the Aerts lab checksum list.
* The Aerts lab occasionally reorganises paths; if a SCENIC URL 404s, browse
  https://resources.aertslab.org/cistarget/databases/ and update SCENIC_URLS below.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen, Request


# ===========================================================================
# Enrichr gene-set libraries  ->  enrichr_gene_sets.json
# ===========================================================================
ENRICHR_BASE = "https://maayanlab.cloud/Enrichr/geneSetLibrary"
ENRICHR_LIBRARIES = {
    "CellMarker_Augmented_2021": "CellMarker2",
    "PanglaoDB_Augmented_2021":  "PanglaoDB",
    "Tabula_Muris":              "Tabula_Muris",
    "Mouse_Gene_Atlas":          "Mouse_Gene_Atlas",
}
MOUSE_HINTS = ["mus musculus", " mm ", "murine", "mouse", "_mm_"]
HUMAN_ONLY  = ["homo sapiens", " hs ", "_hs_"]
MIN_GENES_PER_TERM = 3


def _to_mouse_gene(g: str) -> str:
    g = g.strip()
    if not g:
        return g
    if len(g) > 1 and g.isupper():
        return g[0].upper() + g[1:].lower()
    return g


def _is_mouse_term(term: str) -> bool:
    tl = term.lower()
    if any(h in tl for h in HUMAN_ONLY):
        return False
    return True


def _clean_term(term: str) -> str:
    term = re.sub(r"_Mus musculus.*$", "", term, flags=re.IGNORECASE)
    term = re.sub(r"_Homo sapiens.*$", "", term, flags=re.IGNORECASE)
    term = re.sub(r"\s+Mm\s+CL:\d+.*$", "", term, flags=re.IGNORECASE)
    term = re.sub(r"\s+Hs\s+CL:\d+.*$", "", term, flags=re.IGNORECASE)
    term = re.sub(r"\s+\(.*?\)$", "", term)
    return term.strip()


def download_enrichr(out_dir: Path) -> None:
    import requests  # only needed for this path
    out_path = out_dir / "enrichr_gene_sets.json"
    per_library: dict[str, dict[str, list[str]]] = {}
    sources: list[str] = []
    for lib_name, label in ENRICHR_LIBRARIES.items():
        try:
            print(f"  Fetching {lib_name} ...", flush=True)
            r = requests.get(ENRICHR_BASE,
                             params={"libraryName": lib_name, "mode": "text"},
                             timeout=180)
            r.raise_for_status()
            terms: dict[str, list[str]] = {}
            for line in r.text.splitlines():
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                term  = parts[0].strip()
                genes = [g.strip() for g in parts[2:] if g.strip()]
                if not genes or not _is_mouse_term(term):
                    continue
                clean = _clean_term(term)
                if not clean:
                    continue
                norm = sorted({_to_mouse_gene(g) for g in genes if len(g) >= 2})
                if len(norm) >= MIN_GENES_PER_TERM:
                    terms[clean] = sorted(set(terms.get(clean, [])) | set(norm))
            per_library[label] = terms
            sources.append(label)
            print(f"    {len(terms):,} mouse terms kept from {label}", flush=True)
        except Exception as exc:
            print(f"  WARNING: could not fetch {lib_name}: {exc}", flush=True)
        time.sleep(1)

    if not per_library:
        raise RuntimeError("No Enrichr gene sets downloaded — check internet connection.")

    combined: dict[str, set] = defaultdict(set)
    for terms in per_library.values():
        for term, genes in terms.items():
            combined[term].update(genes)
    merged = {t: sorted(g) for t, g in sorted(combined.items()) if len(g) >= MIN_GENES_PER_TERM}

    output = {
        "metadata": {
            "downloaded": datetime.date.today().isoformat(),
            "libraries":  sources,
            "species":    "mouse",
            "n_terms_merged": len(merged),
            "gene_convention": "Mouse title-case (Lgr5, Muc2). Human ALL-CAPS converted.",
        },
        "libraries":  {lbl: dict(sorted(t.items())) for lbl, t in per_library.items()},
        "gene_sets":  merged,
    }
    out_path.write_text(json.dumps(output, indent=2))
    total = sum(len(v) for v in per_library.values())
    print(f"\n[OK] Enrichr: {total:,} per-library terms "
          f"({len(merged):,} merged) -> {out_path}")


# ===========================================================================
# SCENIC / cisTarget databases  ->  scenic_refs/
# ===========================================================================
# v10 (mc_v10_clust, "SCENIC+") databases from resources.aertslab.org.
SCENIC_URLS = {
    "mouse": {
        "10kb":  "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/"
                 "refseq_r80/mc_v10_clust/gene_based/"
                 "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "500bp": "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/"
                 "refseq_r80/mc_v10_clust/gene_based/"
                 "mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "motif": "https://resources.aertslab.org/cistarget/motif2tf/"
                 "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
        "tfs":   "https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt",
    },
    "human": {
        "10kb":  "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/"
                 "refseq_r80/mc_v10_clust/gene_based/"
                 "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "500bp": "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/"
                 "refseq_r80/mc_v10_clust/gene_based/"
                 "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "motif": "https://resources.aertslab.org/cistarget/motif2tf/"
                 "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "tfs":   "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt",
    },
}
CHECKSUMS_URL = "https://resources.aertslab.org/cistarget/databases/sha256sum.txt"


def _human(n: int) -> str:
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024 or u == "GB":
            return f"{n:.1f}{u}"
        n /= 1024.0


def download_file(url: str, dest: Path, expected_size: int | None = None) -> Path:
    """Stream a URL to dest with a progress line. Skips if already fully present."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        head = urlopen(Request(url, method="HEAD"), timeout=60)
        total = int(head.headers.get("Content-Length", 0)) or (expected_size or 0)
    except Exception:
        total = expected_size or 0

    if dest.exists() and total and dest.stat().st_size == total:
        print(f"  [skip] {dest.name} already present ({_human(total)})")
        return dest

    print(f"  [get ] {dest.name}  ({_human(total) if total else 'size unknown'})")
    req = Request(url, headers={"User-Agent": "reference_downloader/1.0"})
    with urlopen(req, timeout=120) as resp, open(dest, "wb") as fh:
        got = 0
        last = time.time()
        while True:
            chunk = resp.read(1 << 20)  # 1 MB
            if not chunk:
                break
            fh.write(chunk)
            got += len(chunk)
            if time.time() - last > 2:
                pct = f"{100*got/total:.0f}%" if total else _human(got)
                print(f"        {dest.name}: {pct} ({_human(got)})", flush=True)
                last = time.time()
    print(f"  [done] {dest.name} -> {_human(dest.stat().st_size)}")
    return dest


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def download_scenic(out_dir: Path, species: str, db: str, verify: bool) -> None:
    urls = SCENIC_URLS[species]
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nSCENIC refs ({species}, {db} database) -> {out_dir}")
    targets = [
        (urls[db],   out_dir / Path(urls[db]).name),
        (urls["motif"], out_dir / Path(urls["motif"]).name),
        (urls["tfs"],   out_dir / Path(urls["tfs"]).name),
    ]
    for url, dest in targets:
        try:
            download_file(url, dest)
        except Exception as exc:
            print(f"  WARNING: failed {url}\n           {exc}")
            print("           (check the path at https://resources.aertslab.org/cistarget/databases/)")

    if verify:
        try:
            import requests
            print("\n  Verifying SHA-256 against Aerts lab checksums ...")
            sums = {}
            for line in requests.get(CHECKSUMS_URL, timeout=120).text.splitlines():
                p = line.split()
                if len(p) == 2:
                    sums[Path(p[1]).name] = p[0]
            for _, dest in targets:
                if dest.exists() and dest.name in sums:
                    ok = _sha256(dest) == sums[dest.name]
                    print(f"    {dest.name}: {'OK' if ok else 'MISMATCH'}")
        except Exception as exc:
            print(f"  [verify skipped]: {exc}")

    print(f"\n[OK] SCENIC refs in {out_dir}. Point the pipeline at it with "
          f"extra['scenic_refs_dir']='{out_dir}' (and drop --skip-scenic).")


# ===========================================================================
def main() -> None:
    ap = argparse.ArgumentParser(description="Download reference databases (standalone).")
    ap.add_argument("--all", action="store_true", help="Enrichr + SCENIC (defaults)")
    ap.add_argument("--enrichr", action="store_true", help="download Enrichr gene sets")
    ap.add_argument("--scenic",  action="store_true", help="download SCENIC/cisTarget DBs")
    ap.add_argument("--out-dir", default=".", help="where enrichr_gene_sets.json goes (default: here)")
    ap.add_argument("--scenic-dir", default="scenic_refs", help="SCENIC output folder")
    ap.add_argument("--species", choices=["mouse", "human"], default="mouse")
    ap.add_argument("--scenic-db", choices=["10kb", "500bp"], default="10kb",
                    help="ranking DB window (500bp is smaller/faster, 10kb is standard)")
    ap.add_argument("--verify", action="store_true", help="SHA-256 verify SCENIC files")
    args = ap.parse_args()

    do_enrichr = args.enrichr or args.all
    do_scenic  = args.scenic or args.all
    if not (do_enrichr or do_scenic):
        ap.error("nothing to do — pass --all, --enrichr, and/or --scenic")

    if do_enrichr:
        print("=== Enrichr gene-set libraries ===")
        download_enrichr(Path(args.out_dir).expanduser())
    if do_scenic:
        download_scenic(Path(args.scenic_dir).expanduser(), args.species,
                        args.scenic_db, args.verify)


if __name__ == "__main__":
    main()
