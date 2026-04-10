# Author: Selim Romero, Texas A&M University
# Pathway gene fetching from Enrichr.

#' Fetch pathway gene list from Enrichr
#'
#' Queries the Enrichr API for gene sets matching a pathway name keyword
#' and returns the gene list from the best match.
#'
#' @param pathway_name Character. Pathway name or keyword for case-insensitive
#'   substring search (e.g., \code{"Wnt signaling"}, \code{"intestinal immune"}).
#' @param organism Character. Either \code{"human"} (default) or
#'   \code{"mouse"}.
#' @param library Character or NULL. Enrichr library to restrict to
#'   (e.g., \code{"KEGG_2021_Human"}, \code{"Reactome_2022"}). If NULL
#'   (default), searches across several major libraries automatically.
#'
#' @return Character vector of gene names from the best-matching pathway.
#'   Returns an empty character vector if no match is found.
#'
#' @details
#'   Requires the \pkg{httr} and \pkg{jsonlite} packages.
#'   Searches are performed by downloading the gene set library and doing
#'   local substring matching on pathway names.
#'
#' @importFrom httr GET content status_code
#' @importFrom jsonlite fromJSON
#' @export
fetch_pathway_genes <- function(pathway_name, organism = "human",
                                library = NULL) {
  if (!requireNamespace("httr",     quietly = TRUE)) stop("Package 'httr' required.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' required.")

  # Default libraries by organism
  if (is.null(library)) {
    if (organism == "mouse") {
      libraries <- c("KEGG_2019_Mouse", "GO_Biological_Process_2023",
                     "Reactome_2022")
    } else {
      libraries <- c("KEGG_2021_Human", "GO_Biological_Process_2023",
                     "Reactome_2022", "WikiPathway_2023_Human")
    }
  } else {
    libraries <- library
  }

  base_url <- "https://maayanlab.cloud/Enrichr"

  for (lib in libraries) {
    tryCatch({
      # Download full gene-set library
      url <- sprintf("%s/geneSetLibrary?mode=text&libraryName=%s", base_url, lib)
      resp <- httr::GET(url, httr::timeout(30))

      if (httr::status_code(resp) != 200) next

      raw_text <- httr::content(resp, as = "text", encoding = "UTF-8")
      lines    <- strsplit(raw_text, "\n")[[1]]
      lines    <- lines[nzchar(lines)]

      # Each line: "PathwayName\tDescription\tGene1\tGene2\t..."
      for (line in lines) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) < 3) next

        pw_name <- parts[1]
        if (grepl(pathway_name, pw_name, ignore.case = TRUE)) {
          genes <- parts[-(1:2)]   # skip name + description
          genes <- genes[nzchar(genes)]
          message(sprintf("Found pathway: '%s' in library '%s' (%d genes)",
                           pw_name, lib, length(genes)))
          return(genes)
        }
      }
    }, error = function(e) {
      message(sprintf("  Warning: library '%s' could not be fetched: %s", lib, e$message))
    })
  }

  warning(sprintf("No pathway matching '%s' found across searched libraries.", pathway_name))
  return(character(0))
}
