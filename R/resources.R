#' Load PROGENy pathway footprint gene sets
#'
#' Fetches the PROGENy perturbation-response footprint resource at runtime
#' through decoupleR and returns it as a tidy data.frame
#' (\code{source}/\code{target}/\code{mor}) ready for
#' \code{\link{compute_pathway_activity}}. The \code{mor} column carries
#' PROGENy's perturbation-response weights, which the \code{mean_z} and
#' decoupleR back ends honour (negative weights = genes repressed by the
#' pathway). Resources are never redistributed with the package and are
#' fetched on demand; an internet connection is required for
#' the first call.
#'
#' PROGENy's 14 cancer-signalling footprints are a useful finer-grained
#' alternative to MSigDB Hallmark sets, whose coarse granularity saturates
#' the permutation floor on real Visium sections.
#'
#' @param organism \code{"human"} (default) or \code{"mouse"}.
#' @param top Number of top-weighted genes kept per pathway (default 500).
#'
#' @return A tidy data.frame with columns \code{source}, \code{target},
#'   \code{mor}.
#'
#' @examples
#' \dontrun{
#' pw <- get_progeny_sets("human", top = 250)
#' act <- compute_pathway_activity(counts, pw, coords = coords)
#' }
#'
#' @export
get_progeny_sets <- function(organism = c("human", "mouse"), top = 500L) {
  organism <- match.arg(organism)
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package 'decoupleR' is required to fetch the PROGENy resource.",
         call. = FALSE)
  }
  pw <- as.data.frame(decoupleR::get_progeny(organism = organism, top = top))
  wcol <- intersect(c("mor", "weight"), colnames(pw))[1]
  if (is.na(wcol) || !all(c("source", "target") %in% colnames(pw))) {
    stop("Unexpected PROGENy resource layout from decoupleR ",
         "(expected source/target/weight columns).", call. = FALSE)
  }
  data.frame(source = as.character(pw$source),
             target = as.character(pw$target),
             mor = as.numeric(pw[[wcol]]),
             stringsAsFactors = FALSE)
}
