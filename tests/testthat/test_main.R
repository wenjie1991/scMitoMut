test_that("test_main", {
  # load the data
  ## Use the example data
  f <- system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")


  ## Load the data with parse_table function
  f_h5_tmp <- tempfile(fileext = ".h5")
  f_h5 <- parse_table(f, sep = "\t", h5_file = f_h5_tmp)

  # open the h5f file
  x <- open_h5_file(f_h5)
  # run the model fit
  run_model_fit(x)
  # expect_snapshot(format(x))
  expect_output(format(x), "Available loci: 16")
  expect_output(format(x), "Available cells: 1359")
  # Filter the loci based on the model fit results
  x <- filter_loc(x, min_cell = 10, model = "bb", p_threshold = 0.01, p_adj_method = "fdr")
  # expect_snapshot(format(x))
  expect_output(format(x), "Loci passed the filter: 10")
})
