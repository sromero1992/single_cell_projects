# Lightweight tests that do not require Seurat/Bioconductor to be installed.
# They lock the two contracts most likely to break silently: the metadata schema
# and the marker CSV parser (including pipe-delimited gene explosion).

test_that("TamuScDSC_schema exposes the expected column names", {
  s <- TamuScDSC_schema()
  expect_true(all(c("df_score", "sc_score", "consensus", "cluster_flag") %in% names(s)))
  # Must match the names 01b_cluster_doublet_review.R expects, verbatim.
  expect_identical(s$df_score, "DF_score")
  expect_identical(s$sc_score, "scDblFinder_score")
  expect_identical(s$consensus, "doublet_consensus")
})

test_that("norm_key and clean_header normalise as expected", {
  expect_identical(norm_key("Sample.ID"), norm_key("sample_id"))
  expect_identical(clean_header("  Diet <x/> "), "Diet")
})

test_that("get_markers explodes pipe-delimited genes to one row each", {
  m <- get_markers()
  skip_if(nrow(m) == 0, "marker CSV not installed in this build")
  expect_true(all(c("cell_type", "gene", "weight") %in% colnames(m)))
  # No pipe characters should survive in the exploded gene column.
  expect_false(any(grepl("\\|", m$gene)))
  # A known cell type from the shipped CSV should map to multiple genes.
  bcells <- m$gene[m$cell_type == "B cells"]
  expect_true(length(bcells) >= 2)
})

test_that("add_markers and markers_as_list round-trip", {
  m <- get_markers()
  skip_if(nrow(m) == 0, "marker CSV not installed in this build")
  m2 <- add_markers(m, "Tuft", c("Dclk1", "Trpm5"))
  expect_true(all(c("Dclk1", "Trpm5") %in% m2$gene[m2$cell_type == "Tuft"]))
  lst <- markers_as_list(m2)
  expect_true(is.list(lst))
  expect_true(all(c("Dclk1", "Trpm5") %in% lst[["Tuft"]]))
})

test_that("consensus rule combines calls correctly", {
  df <- c(c1 = "Doublet", c2 = "Singlet", c3 = "Doublet")
  sc <- c(c1 = "Doublet", c2 = "Doublet", c3 = "Singlet")
  u  <- TamuScDSC:::.apply_consensus(df, sc, "union")
  i  <- TamuScDSC:::.apply_consensus(df, sc, "intersect")
  expect_identical(unname(u[c("c1", "c2", "c3")]), c("Doublet", "Doublet", "Doublet"))
  expect_identical(unname(i[c("c1", "c2", "c3")]), c("Doublet", "Singlet", "Singlet"))
})
