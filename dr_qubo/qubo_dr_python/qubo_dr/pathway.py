# Author: Selim Romero, Texas A&M University
"""
Pathway gene fetching module.

Retrieves gene lists from KEGG and Gene Ontology (GO Biological Process) databases.
Supports both REST API queries and local library access via gseapy.
"""

import time
import requests
import numpy as np

# Optional imports
try:
    import gseapy
    HAS_GSEAPY = True
except ImportError:
    HAS_GSEAPY = False


KEGG_REST_BASE = "https://rest.kegg.jp"
QUICKGO_REST_BASE = "https://www.ebi.ac.uk/QuickGO/services"


def fetch_kegg_pathway(pathway_id, organism='hsa'):
    """
    Fetch gene list from KEGG REST API.

    Queries the KEGG REST API for a given pathway and organism, extracts gene symbols.

    Parameters
    ----------
    pathway_id : str
        KEGG pathway ID. Can be bare ID (e.g., '04310') or prefixed (e.g., 'hsa04310').
        If bare, organism prefix is added.
    organism : str, optional
        KEGG organism code (default 'hsa' for human).
        Other options: 'mmu' (mouse), 'dme' (fly), 'cel' (worm), etc.

    Returns
    -------
    genes : list of str
        Gene symbols (uppercase) from the pathway.

    Raises
    ------
    ValueError
        If API request fails or pathway not found.
    requests.RequestException
        If network error occurs.

    Notes
    -----
    Uses KEGG REST API: https://rest.kegg.jp/get/path:{pathway_id}
    Extracts gene symbols from the GENE section of the response.
    Rate-limited with small delay to respect server.

    Examples
    --------
    >>> genes = fetch_kegg_pathway('04310', organism='hsa')  # Wnt signaling
    >>> len(genes) > 0
    True
    """
    # Normalize pathway_id (add organism prefix if missing)
    if not pathway_id.startswith(organism):
        if pathway_id.startswith('path:'):
            pathway_id = pathway_id[5:]  # Remove 'path:' prefix if present
        pathway_id = f"path:{organism}{pathway_id}"
    else:
        pathway_id = f"path:{pathway_id}"

    url = f"{KEGG_REST_BASE}/get/{pathway_id}"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
    except requests.RequestException as e:
        raise ValueError(f"Failed to fetch KEGG pathway {pathway_id}: {e}")

    text = response.text
    genes = []

    # Parse GENE section
    lines = text.split('\n')
    in_gene_section = False

    for line in lines:
        if line.startswith('GENE'):
            in_gene_section = True
            continue
        if in_gene_section:
            if line and not line[0].isspace():
                # End of GENE section
                break
            if line.strip():
                # Parse gene line: "        1234  SYMBOL_NAME; full name ..."
                parts = line.split()
                if len(parts) >= 2:
                    gene_symbol = parts[1].split(';')[0].strip()
                    genes.append(gene_symbol.upper())

    if not genes:
        raise ValueError(f"No genes found in KEGG pathway {pathway_id}")

    # Rate limiting
    time.sleep(0.5)

    return genes


def fetch_kegg_pathway_list(organism='hsa'):
    """
    Fetch list of all KEGG pathways for a given organism.

    Parameters
    ----------
    organism : str, optional
        KEGG organism code (default 'hsa').

    Returns
    -------
    pathways : dict
        Dictionary mapping pathway ID to pathway name.
        Keys are bare IDs (e.g., '04310'), values are descriptions.

    Raises
    ------
    ValueError
        If API request fails.
    """
    url = f"{KEGG_REST_BASE}/list/pathway/{organism}"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
    except requests.RequestException as e:
        raise ValueError(f"Failed to fetch KEGG pathway list: {e}")

    pathways = {}
    lines = response.text.strip().split('\n')

    for line in lines:
        parts = line.split('\t')
        if len(parts) >= 2:
            path_id = parts[0].replace(f'path:{organism}', '')
            path_name = parts[1]
            pathways[path_id] = path_name

    time.sleep(0.5)
    return pathways


def search_kegg_pathways(keyword, organism='hsa'):
    """
    Search KEGG pathways by keyword.

    Parameters
    ----------
    keyword : str
        Search term (e.g., 'Wnt', 'signaling').
    organism : str, optional
        KEGG organism code (default 'hsa').

    Returns
    -------
    matching_pathways : dict
        Dictionary of {pathway_id: pathway_name} for pathways matching keyword.
    """
    all_pathways = fetch_kegg_pathway_list(organism=organism)
    keyword_lower = keyword.lower()

    matching = {
        path_id: path_name
        for path_id, path_name in all_pathways.items()
        if keyword_lower in path_name.lower()
    }

    return matching


def fetch_gobp_geneset(go_term_name=None, go_term_id=None, organism='human', source='gseapy'):
    """
    Fetch gene list for a GO Biological Process term.

    Supports two sources:
    - 'gseapy': Uses local gseapy library (GO_Biological_Process_2021).
    - 'api': Queries QuickGO REST API directly.

    Parameters
    ----------
    go_term_name : str, optional
        GO term name (e.g., 'Wnt signaling pathway').
    go_term_id : str, optional
        GO term ID (e.g., 'GO:0016055'). Takes precedence if both provided.
    organism : str, optional
        Organism for mapping gene IDs ('human' or 'mouse', default 'human').
    source : str, {'gseapy', 'api'}, optional
        Data source (default 'gseapy').

    Returns
    -------
    genes : list of str
        Gene symbols (uppercase) annotated to the GO term.

    Raises
    ------
    ValueError
        If term not found or neither go_term_name nor go_term_id provided.
    ImportError
        If gseapy source selected but gseapy not installed.

    Notes
    -----
    gseapy source uses the 'GO_Biological_Process_2021' library.
    api source queries QuickGO with organism mapping (human→9606, mouse→10090).
    """
    if go_term_id is None and go_term_name is None:
        raise ValueError("Must provide either go_term_id or go_term_name")

    if source == 'gseapy':
        if not HAS_GSEAPY:
            raise ImportError("gseapy not installed. Install via: pip install gseapy")
        return _fetch_gobp_gseapy(go_term_name=go_term_name, go_term_id=go_term_id)

    elif source == 'api':
        return _fetch_gobp_quickgo(go_term_id=go_term_id, organism=organism)

    else:
        raise ValueError(f"Unknown source: {source}. Choose 'gseapy' or 'api'")


def _fetch_gobp_gseapy(go_term_name=None, go_term_id=None):
    """Fetch from gseapy GO library."""
    try:
        gene_sets = gseapy.get_library(name='GO_Biological_Process_2021')
    except Exception as e:
        raise ValueError(f"Failed to load GO library via gseapy: {e}")

    if go_term_id:
        # Search by exact ID
        go_term_id = go_term_id.replace('GO:', '')
        for term_key in gene_sets.keys():
            if go_term_id in term_key:
                return list(gene_sets[term_key])
        raise ValueError(f"GO term {go_term_id} not found in gseapy library")

    if go_term_name:
        # Search by name (substring match)
        name_lower = go_term_name.lower()
        for term_key in gene_sets.keys():
            if name_lower in term_key.lower():
                return list(gene_sets[term_key])
        raise ValueError(f"GO term '{go_term_name}' not found in gseapy library")


def _fetch_gobp_quickgo(go_term_id, organism='human'):
    """Fetch from QuickGO REST API."""
    if not go_term_id:
        raise ValueError("go_term_id required for QuickGO API fetch")

    # Clean up GO ID
    go_term_id = go_term_id.replace('GO:', '')

    # Map organism to taxon ID
    taxon_map = {'human': '9606', 'mouse': '10090'}
    taxon_id = taxon_map.get(organism.lower(), '9606')

    url = f"{QUICKGO_REST_BASE}/annotation/downloadSearch"
    params = {
        'goId': f'GO:{go_term_id}',
        'taxonId': taxon_id,
        'aspect': 'biological_process',
        'geneProductType': 'protein',
        'pageSize': 10000,
    }

    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
    except requests.RequestException as e:
        raise ValueError(f"Failed to fetch from QuickGO: {e}")

    # Parse TSV/CSV response
    genes = set()
    lines = response.text.strip().split('\n')

    for line in lines[1:]:  # Skip header
        parts = line.split('\t')
        if len(parts) > 2:
            gene_symbol = parts[2].strip().upper()
            if gene_symbol:
                genes.add(gene_symbol)

    time.sleep(0.5)
    return list(genes)


def list_gobp_terms(keyword, organism='human'):
    """
    Search GO Biological Process terms by keyword.

    Parameters
    ----------
    keyword : str
        Search term (e.g., 'Wnt', 'signaling').
    organism : str, optional
        Organism for mapping (default 'human').

    Returns
    -------
    matching_terms : dict
        Dictionary mapping term names to gene lists.

    Raises
    ------
    ImportError
        If gseapy not installed.
    """
    if not HAS_GSEAPY:
        raise ImportError("gseapy required for GO term listing. Install via: pip install gseapy")

    try:
        gene_sets = gseapy.get_library(name='GO_Biological_Process_2021')
    except Exception as e:
        raise ValueError(f"Failed to load GO library: {e}")

    keyword_lower = keyword.lower()
    matching = {}

    for term_key, genes in gene_sets.items():
        if keyword_lower in term_key.lower():
            matching[term_key] = list(genes)

    return matching


def get_pathway_genes(source, pathway_id_or_name, organism='human', **kwargs):
    """
    Unified interface to fetch pathway genes from KEGG or GO databases.

    Parameters
    ----------
    source : str, {'kegg', 'gobp'}
        Pathway database source.
    pathway_id_or_name : str
        KEGG pathway ID (e.g., '04310') or GO term name/ID (e.g., 'GO:0016055').
    organism : str, optional
        Organism for pathway mapping (default 'human').
        KEGG: organism code ('hsa', 'mmu', etc.).
        GOBP: 'human' or 'mouse'.
    **kwargs : dict
        Additional arguments passed to fetch functions
        (e.g., 'source' for GOBP API choice).

    Returns
    -------
    genes : list of str
        Gene symbols (uppercase) from the pathway.

    Raises
    ------
    ValueError
        If source unknown or pathway not found.

    Examples
    --------
    >>> genes = get_pathway_genes('kegg', '04310', organism='hsa')
    >>> genes = get_pathway_genes('gobp', 'Wnt signaling pathway', organism='human')
    """
    source_lower = source.lower()

    if source_lower == 'kegg':
        return fetch_kegg_pathway(pathway_id_or_name, organism=organism)

    elif source_lower == 'gobp':
        # Infer whether input is ID or name
        if pathway_id_or_name.startswith('GO:'):
            return fetch_gobp_geneset(go_term_id=pathway_id_or_name, organism=organism, **kwargs)
        else:
            return fetch_gobp_geneset(go_term_name=pathway_id_or_name, organism=organism, **kwargs)

    else:
        raise ValueError(f"Unknown pathway source: {source}. Choose 'kegg' or 'gobp'")
