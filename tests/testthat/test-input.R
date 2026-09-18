test_that("svp_input handles matrices and validates coordinates", {
  counts <- matrix(rpois(60, 5), nrow = 6,
                   dimnames = list(paste0("g", 1:6), paste0("s", 1:10)))
  coords <- data.frame(x = rnorm(10), y = rnorm(10))
  inp <- svp_input(counts, coords)
  expect_equal(dim(inp$counts), c(6, 10))
  expect_equal(nrow(inp$coords), 10)
  expect_named(inp$coords, c("x", "y"))

  # unnamed columns fall back to the first two
  coords2 <- data.frame(px = rnorm(10), py = rnorm(10))
  expect_equal(nrow(svp_input(counts, coords2)$coords), 10)

  # x/y column names are picked even when not first
  coords3 <- data.frame(junk = 1:10, x = rnorm(10), y = rnorm(10))
  expect_named(svp_input(counts, coords3)$coords, c("x", "y"))

  expect_error(svp_input(counts, NULL), "coords")
  expect_error(svp_input(counts, coords[1:5, ]), "match")
  expect_error(svp_input(counts, data.frame(a = 1)), "at least two")
  expect_error(svp_input("nope", coords), "Unsupported")
})

test_that("svp_input reads SpatialExperiment", {
  skip_if_not_installed("SpatialExperiment")
  counts <- matrix(rpois(60, 5), nrow = 6,
                   dimnames = list(paste0("g", 1:6), paste0("s", 1:10)))
  coords <- cbind(x = rnorm(10), y = rnorm(10))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts), spatialCoords = coords)
  inp <- svp_input(spe)
  expect_equal(dim(inp$counts), c(6, 10))
  expect_equal(nrow(inp$coords), 10)
})

test_that("svp_input gives informative errors for Seurat without images", {
  skip_if_not_installed("SeuratObject")
  counts <- matrix(rpois(60, 5), nrow = 6,
                   dimnames = list(paste0("g", 1:6), paste0("s", 1:10)))
  so <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  expect_error(svp_input(so), "tissue coordinates|coords")
})
