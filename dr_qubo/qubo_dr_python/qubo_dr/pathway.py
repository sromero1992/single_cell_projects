# Author: Selim Romero, Texas A&M University
"""
Pathway gene fetching and browsing — Enrichr unified backend.

Uses gseapy to access Enrichr gene-set collections, which cover:
  KEGG, GO Biological Process, GO Molecular Function, Reactome,
  WikiPathways, MSigDB Hallmarks, and hundreds more.

No KEGG REST API or QuickGO calls.  A single keyword search spans all
default libraries simultaneously.  Specify ``library`` only when you
want to restrict the search for speed.
"""

# Optional import — gseapy is the only external dependency
try:
    import gseapy
    HAS_GSEAPY = True
except ImportError:
    HAS_GSEAPY = False


# ---------------------------------------------------------------------------
# Default Enrichr libraries searched when no library is specified
# Priority order: KEGG first (most curated pathways), then GO/Reactome/Wiki
# ---------------------------------------------------------------------------
_DEFAULT_LIBRARIES = {
    'human': [
        'KEGG_2021_Human',
        'GO_Biological_Process_2025',
        'GO_Molecular_Function_2025',
        'Reactome_2024',
        'WikiPathway_2024_Human',
        'MSigDB_Hallmark_2020',
    ],
    'mouse': [
        'KEGG_2026',
        'GO_Biological_Process_2025',
        'GO_Molecular_Function_2025',
        'Reactome_2024',
    ],
}


def _get_default_libraries(organism='human'):
    """Return the default library list for an organism."""
    return _DEFAULT_LIBRARIES.get(organism.lower(), _DEFAULT_LIBRARIES['human'])


def _require_gseapy():
    if not HAS_GSEAPY:
        raise ImportError(
            "gseapy is required: pip install gseapy\n"
            "It provides access to KEGG, GO, Reactome, WikiPathways, and more via Enrichr."
        )


# ===========================================================================
#  Library discovery
# ===========================================================================

def list_enrichr_libraries(filter_keyword=None):
    """
    List available Enrichr libraries via gseapy.

    Parameters
    ----------
    filter_keyword : str, optional
        Case-insensitive substring filter.

    Returns
    -------
    libraries : list of str
        Sorted list of library names matching the filter.

    Examples
    --------
    >>> list_enrichr_libraries('KEGG')
    >>> list_enrichr_libraries('GO_Biological')
    >>> list_enrichr_libraries('Reactome')
    >>> list_enrichr_libraries()          # all libraries
    """
    _require_gseapy()
    libs = gseapy.get_library_name()
    if filter_keyword:
        kw = filter_keyword.lower()
        libs = [l for l in libs if kw in l.lower()]
    return sorted(libs)


# ===========================================================================
#  Core Enrichr fetch (internal)
# ===========================================================================

def _fetch_from_library(keyword, library):
    """
    Search one Enrichr library for terms matching keyword (substring, case-insensitive).

    Returns
    -------
    matches : dict  {term_name: [UPPERCASE gene symbols]}
    """
    _require_gseapy()
    try:
        gene_sets = gseapy.get_library(name=library)
    except Exception as e:
        raise ValueError(
            f"Failed to load Enrichr library '{library}': {e}\n"
            f"Run list_enrichr_libraries() to see valid names."
        )
    kw = keyword.lower()
    return {
        k: [g.upper() for g in v]
        for k, v in gene_sets.items()
        if kw in k.lower()
    }


# ===========================================================================
#  Public API
# ===========================================================================

def browse_pathways(keyword, organism='human', libraries=None,
                    show_genes=False, max_results=20):
    """
    Search all default Enrichr libraries for gene sets matching a keyword.

    Searches KEGG, GO Biological Process, GO Molecular Function, Reactome,
    WikiPathways, and MSigDB Hallmarks simultaneously (for human).
    No API key or special setup is required.

    Parameters
    ----------
    keyword : str
        Search keyword, e.g. 'Wnt', 'TGF-beta', 'NF-kB', 'apoptosis'.
        Case-insensitive substring match.
    organism : str, optional
        'human' or 'mouse' (default 'human').  Selects the organism-specific
        library set when libraries=None.
    libraries : list of str, optional
        Override the default library set.  Useful to restrict to specific
        databases.  E.g. ['KEGG_2021_Human', 'Reactome_2022'].
        Use list_enrichr_libraries() to explore all options.
    show_genes : bool, optional
        Print the first 10 genes for each matching term (default False).
    max_results : int, optional
        Max number of matches to display per library (default 20).

    Returns
    -------
    results : dict
        {library_name: {term_name: [gene_list]}}  for all matches found.

    Examples
    --------
    >>> browse_pathways('Wnt')
    >>> browse_pathways('TGF', organism='mouse')
    >>> browse_pathways('apoptosis', show_genes=True)
    >>> browse_pathways('NF-kB', libraries=['Reactome_2022', 'KEGG_2021_Human'])
    """
    _require_gseapy()
    libs = libraries or _get_default_libraries(organism)
    results = {}

    for lib in libs:
        try:
            matches = _fetch_from_library(keyword, lib)
        except ValueError as e:
            print(f"  [warning] Skipping '{lib}': {e}")
            continue

        if not matches:
            continue

        results[lib] = matches
        print(f"\n[{lib}  —  matches for '{keyword}']")
        for term, genes in list(matches.items())[:max_results]:
            print(f"  [{len(genes):4d} genes]  {term}")
            if show_genes:
                print(f"             → {', '.join(genes[:10])} ...")

    if not results:
        print(f"\nNo matches for '{keyword}' across {len(libs)} libraries.")
        print(f"  Libraries searched: {libs}")
        print(f"  Tips:")
        print(f"    • try a shorter keyword or alternate spelling")
        print(f"    • list_enrichr_libraries('{keyword[:4]}') to find relevant libraries")
        print(f"    • browse_pathways('{keyword}', libraries=['<specific_library>'])")
    else:
        total = sum(len(m) for m in results.values())
        print(f"\n→ {total} total match(es) across {len(results)} librar(ies).")
        print(f"  Copy the exact term name and pass it to get_pathway_genes() or preview_pathway_genes()")
        print(f"  Optionally add library='<library name>' once you know which database to use.")

    return results


def get_pathway_genes(pathway_name, organism='human', library=None, **kwargs):
    """
    Fetch gene list from Enrichr.

    Searches across all default Enrichr libraries (KEGG, GO BP, GO MF,
    Reactome, WikiPathways, …) and returns the first substring match found.
    Specify ``library`` to restrict the search and speed up the lookup once
    you know which database contains your pathway of interest.

    Parameters
    ----------
    pathway_name : str
        Pathway / gene-set name or substring (case-insensitive).
        E.g. 'Wnt signaling pathway', 'TGF-beta', 'HALLMARK_APOPTOSIS',
             'Intestinal immune network', 'IBD'.
    organism : str, optional
        'human' or 'mouse' (default 'human').
    library : str, optional
        Enrichr library to restrict the search to (default: auto-search
        all default libraries in order).
        E.g. 'KEGG_2021_Human', 'Reactome_2022',
             'GO_Biological_Process_2023', 'WikiPathway_2023_Human'.
        Use list_enrichr_libraries() to browse all ~200 options.

    Returns
    -------
    genes : list of str
        Uppercase gene symbols for the first matching pathway found.

    Raises
    ------
    ValueError
        If no matching pathway is found in any of the searched libraries.

    Examples
    --------
    >>> get_pathway_genes('Wnt signaling pathway')
    >>> get_pathway_genes('Wnt signaling', library='Reactome_2022')
    >>> get_pathway_genes('apoptosis', organism='mouse')
    >>> get_pathway_genes('intestinal', library='KEGG_2021_Human')
    """
    _require_gseapy()
    # Ignore legacy 'source' kwarg silently (backward compat with old pipeline calls)
    kwargs.pop('source', None)

    libs = [library] if library else _get_default_libraries(organism)

    for lib in libs:
        try:
            matches = _fetch_from_library(pathway_name, lib)
        except ValueError:
            continue  # library load error — skip silently

        if matches:
            term, genes = next(iter(matches.items()))
            print(f"  Found '{term}'\n  Library: [{lib}]  |  Genes: {len(genes)}")
            return genes

    raise ValueError(
        f"No pathway matching '{pathway_name}' found.\n"
        f"  Libraries searched: {libs}\n"
        f"  Tips:\n"
        f"    • browse_pathways('{pathway_name}') to see all matches\n"
        f"    • try a shorter or alternate keyword\n"
        f"    • list_enrichr_libraries('{pathway_name.split()[0]}') to discover libraries"
    )


def preview_pathway_genes(pathway_name, organism='human', library=None):
    """
    Fetch and display the full gene list for a specific pathway.

    Useful for inspecting a pathway before running the pipeline.

    Parameters
    ----------
    pathway_name : str
        Pathway name or substring (case-insensitive).
    organism : str, optional
        'human' or 'mouse' (default 'human').
    library : str, optional
        Enrichr library to search.  If None, all default libraries are tried.

    Returns
    -------
    genes : list of str
        Uppercase gene symbols.

    Examples
    --------
    >>> preview_pathway_genes('Wnt signaling pathway')
    >>> preview_pathway_genes('Wnt signaling', library='KEGG_2021_Human')
    >>> preview_pathway_genes('T cell receptor signaling', organism='mouse')
    """
    genes = get_pathway_genes(pathway_name, organism=organism, library=library)
    lib_tag = library or f'{organism} defaults'
    print(f"\n[{lib_tag}  —  '{pathway_name}']  {len(genes)} genes:")
    for i in range(0, len(genes), 6):
        print('  ' + '  '.join(f"{g:<10}" for g in genes[i:i+6]))
    return genes
