#' Normalize spatial transcriptomics input
#'
#' Converts matrices, \code{SpatialExperiment}/\code{SingleCellExperiment},
#' or \code{Seurat} objects into a common internal representation: a
#' genes x spots count matrix plus a two-column coordinate data.frame.
#'
#' @param x A genes x spots count/expression matrix or data.frame, a
#'   \code{SpatialExperiment}, a \code{SingleCellExperiment} (coordinates are
#'   then taken from \code{colData} columns \code{x}/\code{y}), or a
#'   \code{Seurat} object (coordinates from tissue images).
#' @param coords Optional data.frame/matrix with spot coordinates; required
#'   for plain matrix input. Ignored when coordinates are extracted from
#'   \code{x} itself, unless given, in which case it takes precedence.
#' @param assay Assay name or index to use for container classes.
#'
#' @return A list with elements \code{counts} (genes x spots matrix),
#'   \code{coords} (data.frame with columns \code{x}, \code{y}) and
#'   \code{meta} (per-spot metadata data.frame, possibly empty).
#'
#' @examples
#' counts <- matrix(rpois(60, 5), nrow = 6,
#'                  dimnames = list(paste0("g", 1:6), paste0("s", 1:10)))
#' coords <- data.frame(x = rnorm(10), y = rnorm(10))
#' inp <- svp_input(counts, coords)
#' dim(inp$counts)
#'
#' @export
svp_input <- function(x, coords = NULL, assay = NULL) {
  if (is.matrix(x) || is.data.frame(x)) {
    counts <- as.matrix(x)
    if (is.null(coords)) {
      stop("For matrix/data.frame input, 'coords' must be supplied ",
           "(data.frame or matrix with two numeric columns).", call. = FALSE)
    }
    cd <- .coords_xy(coords)
    meta <- data.frame(row.names = colnames(counts) %||%
                         paste0("spot", seq_len(ncol(counts))))
    if (is.null(colnames(counts))) {
      colnames(counts) <- rownames(meta)
    }
  } else if (inherits(x, "SpatialExperiment")) {
    if (!requireNamespace("SpatialExperiment", quietly = TRUE) ||
        !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Packages 'SpatialExperiment' and 'SummarizedExperiment' are ",
           "required to read SpatialExperiment objects.", call. = FALSE)
    }
    counts <- SummarizedExperiment::assay(x, i = assay %||% 1L)
    cd <- if (!is.null(coords)) {
      .coords_xy(coords)
    } else {
      .coords_xy(SpatialExperiment::spatialCoords(x))
    }
    meta <- as.data.frame(SummarizedExperiment::colData(x))
  } else if (inherits(x, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Package 'SummarizedExperiment' is required to read ",
           "SingleCellExperiment objects.", call. = FALSE)
    }
    counts <- SummarizedExperiment::assay(x, i = assay %||% 1L)
    if (is.null(coords)) {
      cd <- SummarizedExperiment::colData(x)
      nms <- tolower(colnames(cd))
      if (!("x" %in% nms && "y" %in% nms)) {
        stop("SingleCellExperiment input needs coordinates: supply 'coords' ",
             "or store colData columns named 'x' and 'y'.", call. = FALSE)
      }
      cd <- .coords_xy(as.data.frame(cd))
    } else {
      cd <- .coords_xy(coords)
    }
    meta <- as.data.frame(SummarizedExperiment::colData(x))
  } else if (inherits(x, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Package 'SeuratObject' is required to read Seurat objects.",
           call. = FALSE)
    }
    an <- assay %||% SeuratObject::DefaultAssay(x)
    counts <- as.matrix(SeuratObject::GetAssayData(x, assay = an, layer = "counts"))
    if (is.null(coords)) {
      tc <- tryCatch(SeuratObject::GetTissueCoordinates(x),
                     error = function(e) NULL)
      if (is.null(tc)) {
        stop("Seurat object carries no tissue coordinates: supply 'coords' ",
             "explicitly.", call. = FALSE)
      }
      coords <- tc
    }
    cd <- .coords_xy(coords)
    meta <- data.frame(row.names = colnames(counts))
  } else {
    stop("Unsupported input class: ", paste(class(x), collapse = "/"),
         ". Provide a matrix, SpatialExperiment, SingleCellExperiment, or ",
         "Seurat object.", call. = FALSE)
  }

  counts <- as.matrix(counts)
  if (is.null(rownames(counts))) {
    rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  }
  if (nrow(cd) != ncol(counts)) {
    stop("Coordinate rows (", nrow(cd), ") do not match spots (", ncol(counts),
         "). Check the orientation: counts must be genes x spots.",
         call. = FALSE)
  }
  if (is.null(rownames(meta)) || length(rownames(meta)) != ncol(counts)) {
    meta <- data.frame(row.names = colnames(counts))
  }
  list(counts = counts, coords = cd, meta = meta)
}
