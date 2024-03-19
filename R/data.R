
#' Load the example datasets for BPCellsExperiment
#'
#' Two subsetted versions of 10x Genomics' 3k PBMC multi-omics datasets
#'
#' @param dataset 'sorted' or 'unsorted'
#'
#' @return
#' \itemize{
#' \item \code{load_example_sce()} returns a \code{\link{SingleCellExperiment}}.
#' \item \code{load_example_csce()} returns a \code{\link{ChromExperiment}}.
#' \item \code{load_example_scme()} combines the two of above and returns a
#' \code{\link{SingleCellMultiExperiment}}
#' }
#'
#' @name load-examples
NULL

#' @examples
#' sce <- load_example_sce()
#' sce
#'
#' @export
#' @rdname load-examples
load_example_sce <- function(dataset = c("sorted", "unsorted")) {
  dataset <- match.arg(arg = dataset)
  mat <- .read_example_10x_h5(dataset = dataset, type = "rna")
  cdata <- .read_example_metadata(dataset = dataset)
  sce <- SingleCellExperiment(
    assays = list(counts = mat$mat),
    rowData = mat$rowData,
    colData = cdata
  )
  return(sce)
}

#' @examples
#' csce <- load_example_csce()
#' csce
#'
#' @importFrom Signac StringToGRanges
#' @export
#' @rdname load-examples
load_example_csce <- function(dataset = c("sorted", "unsorted")) {
  dataset <- match.arg(arg = dataset)
  mat <- .read_example_10x_h5(dataset = dataset, type = "peaks")
  cdata <- .read_example_metadata(dataset = dataset)

  annotations <- .read_example_annot()
  ranges <- StringToGRanges(
    regions = rownames(x = mat$mat),
    sep = c(":", "-")
  )
  frags <- .read_example_frags(dataset = dataset)
  cbpce <- ChromExperiment(
    assays = list(counts = mat$mat),
    rowData = mat$rowData,
    colData = cdata,
    ranges = ranges,
    fragments = frags,
    genome = NULL,
    annotations = annotations
  )
}

#' @examples
#' scme <- load_example_scme()
#' scme
#'
#' @export
#' @rdname load-examples
load_example_scme <- function(dataset = c("sorted", "unsorted")) {
  dataset <- match.arg(arg = dataset)
  bpce <- load_example_sce(dataset = dataset)
  cbpce <- load_example_csce(dataset = dataset)
  scme <- SingleCellMultiExperiment(
    experiments = list(RNA = bpce, ATAC = cbpce)
  )
  return(scme)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom utils read.csv
.read_example_metadata <- function(dataset = c("sorted", "unsorted")) {
  dataset <- match.arg(arg = dataset)
  cdata_file <- paste0("pbmc_", dataset, "_metadata.csv")
  cdata_file <- system.file("extdata", cdata_file, package = .packageName)
  cdata <- read.csv(file = cdata_file, row.names = 1, header = TRUE)
  return(cdata)
}

#' @importFrom BiocGenerics counts
#' @importFrom hdf5r.Extra h5Read
.read_example_10x_h5 <- function(
    dataset = c("sorted", "unsorted"),
    type = c("rna", "peaks")
) {
  dataset <- match.arg(arg = dataset)
  type <- match.arg(arg = type)

  mat_file <- paste0("pbmc_", dataset, "_", type, ".h5")
  mat_file <- system.file("extdata", mat_file, package = .packageName)
  mat <- read10xH5(path = mat_file, use.names = FALSE, use.BPCells = TRUE) %>%
    counts()

  rdata <- h5Read(x = mat_file, name = "matrix/features")[-1] %>%
    as.data.frame()
  rdata <- rdata[, c("id", "name", "feature_type", "genome")]
  rownames(x = rdata) <- rdata$id
  return(list(mat = mat, rowData = rdata))
}

#' @importFrom rtracklayer import
.read_example_annot <- function() {
  gtf_file <- system.file("extdata", "genes.gtf.gz", package = .packageName)
  annot <- import(con = gtf_file)
  return(annot)
}

.read_example_frags <- function(dataset = c("sorted", "unsorted")) {
  dataset <- match.arg(arg = dataset)

  frags_file <- paste0("pbmc_", dataset, "_fragments.h5")
  frags_file <- system.file("extdata", frags_file, package = .packageName)

  frags <- openFragments(x = frags_file, group = "fragments")
  return(frags)
}
