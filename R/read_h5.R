
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions ####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Transform a list to an S4 object
#'
#' Helper function to transform a list to an object with given S4 class.
#'
#' @param robj A list
#' @param S4Class A character specifying the target S4 class.
#'
#' @details
#' Currently supported S4 classes:
#' \itemize{
#' \item \code{\link[SingleCellExperiment]{LinearEmbeddingMatrix}}
#' \item \code{\link[GenomicRanges]{GRanges}}
#' \item \code{\link[S4Vectors]{SimpleList}}
#' \item \code{\link[GenomeInfoDb]{Seqinfo}}
#' }
#'
#' @return If \code{S4Class} is \code{NULL} or an unsupported class, the
#' original \code{robj} will return. Otherwise, a target S4 object will return.
#'
#' @importFrom S4Vectors List
#' @export
BPCE_listToS4 <- function(robj, S4Class = NULL) {
  if (is.null(x = S4Class)) {
    return(robj)
  }
  return(switch(
    EXPR = S4Class,
    "LinearEmbeddingMatrix" = .to_LinearEmbeddingMatrix(robj = robj),
    "GRanges" = .to_GRanges(robj = robj),
    "SimpleList" = List(robj),
    "Seqinfo" = .df_to_seqinfo(df = robj),
    .to_not_support_S4(robj = robj, S4Class = S4Class)
  ))
}

#' Read a SingleCellExperiment from HDF5 file
#'
#' @param file An HDF5 file.
#' @param name An HDF5 group where the \code{SingleCellExperiment} is stored.
#' @param verbose Display progress.
#' @param ... Arguments passed to other methods
#'
#' @return
#' \itemize{
#' \item \code{readH5SCE} returns an object inheriting from
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' \item \code{readH5SCME} returns an object inheriting from
#' \code{\link{SingleCellMultiExperiment}}
#' }
#' The concrete class is automatically deduced.
#'
#' @importFrom hdf5r.Extra h5AbsLinkName h5Attr h5Exists
#' @importFrom SingleCellExperiment int_metadata
#' @export
#' @rdname readH5SCE
readH5SCE <- function(file, name = "/", verbose = TRUE, ...) {
  file <- file_path_as_absolute(x = file)
  name <- h5AbsLinkName(name = name)
  S4Class <- h5Attr(x = file, which = "S4Class", name = name)
  verboseMsg(
    "Reading ", S4Class, ":",
    "\n  File: ", file,
    "\n  Group: ", name
  )
  out <- .read_h5_sce(file = file, name = name, verbose = verbose)
  if (h5Exists(x = file, name = file.path(name, .frag_key))) {
    out$fragments <- .read_h5_fragments(
      file = file,
      chrs = seqlevels(x = out$rowRanges),
      cells = rownames(x = out$colData),
      name = file.path(name, .frag_key),
      verbose = verbose
    )
  }
  out <- switch(
    EXPR = S4Class,
    "SingleCellExperiment" = .list_to_SCE(rlist = out),
    "ChromExperiment" = .list_to_CSCE(rlist = out),
    stop("Unsupported S4Class: ", S4Class)
  )
  int_metadata(x = out)[[.path_key]] <- c(file, name)
  return(out)
}

#' @importFrom hdf5r.Extra h5AbsLinkName h5Attr
#' @importFrom SingleCellExperiment int_metadata
#' @export
#' @rdname readH5SCE
readH5SCME <- function(file, name = "/", verbose = TRUE, ...) {
  name <- h5AbsLinkName(name = name)
  file <- file_path_as_absolute(x = file)
  S4Class <- h5Attr(x = file, which = "S4Class", name = name)
  verboseMsg(
    "Reading ", S4Class, ":",
    "\n  File: ", file,
    "\n  Group: ", name
  )
  out <- .read_h5_scme(file = file, name = name, verbose = verbose)
  out <- switch(
    EXPR = S4Class,
    "SingleCellMultiExperiment" = .list_to_SCME(rlist = out),
    stop("Unsupported S4Class: ", S4Class)
  )
  int_metadata(x = out)[[.path_key]] <- c(file, name)
  return(out)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Transform to S4 object ######################################################

#' @importFrom GenomeInfoDb Seqinfo
.df_to_seqinfo <- function(df) {
  seqlengths <- df[['seqlengths']] %||% NA
  isCircular <- df[['isCircular']] %||% NA
  genome <- df[['genome']] %||% NA
  return(Seqinfo(
    seqnames = rownames(x = df),
    seqlengths = seqlengths,
    isCircular = isCircular,
    genome = genome
  ))
}

#' @importFrom GenomicRanges makeGRangesFromDataFrame
.to_GRanges <- function(robj) {
  sinfo <- .df_to_seqinfo(robj[['seqinfo']])
  gr <- makeGRangesFromDataFrame(
    df = robj[['granges']],
    keep.extra.columns = TRUE,
    seqinfo = sinfo
  )
  return(gr)
}

#' @importFrom SingleCellExperiment LinearEmbeddingMatrix
.to_LinearEmbeddingMatrix <- function(robj) {
  return(LinearEmbeddingMatrix(
    sampleFactors = as.matrix(x = robj[['sampleFactors']]),
    featureLoadings = as.matrix(x = robj[['featureLoadings']]),
    factorData = robj[['factorData']],
    metadata = robj[['metadata']]
  ))
}

.to_not_support_S4 <- function(robj, S4Class) {
  warning(
    "Reading ", S4Class, " from HDF5 file is not supported.",
    call. = FALSE, immediate. = FALSE
  )
  return(robj)
}

## SingleCellExperiment ########################################################

#' @importFrom hdf5r.Extra h5AbsLinkName h5Attr h5List h5Read
.read_h5_sce_assays <- function(
    file,
    rownames,
    colnames,
    name = "/assays",
    verbose = TRUE,
    ...
) {
  name <- h5AbsLinkName(name = name)
  assays <- list()
  assay.names <- h5List(x = file, name = name)
  for (i in assay.names) {
    verboseMsg("Reading assay '",  i, "'")
    h5a <- h5Attr(x = file, which = "version", name = file.path(name, i))
    if (any(h5a %in% .bpcells_matrix_attr)) {
      assays[[i]] <- open_matrix_hdf5(path = file, group = file.path(name, i))
    } else {
      assays[[i]] <- h5Read(x = file, name = file.path(name, i))
    }
    rownames(x = assays[[i]]) <- rownames
    colnames(x = assays[[i]]) <- colnames
  }
  return(assays)
}

.read_h5_sce_Exps <- function(
    file,
    name = "/altExps",
    key = "altExp",
    verbose = TRUE,
    ...
) {
  name <- h5AbsLinkName(name = name)
  altExps <- list()
  altExp.names <- h5List(x = file, name = name)
  for (i in altExp.names) {
    verboseMsg("Reading ", key, " '",  i, "'")
    altExps[[i]] <- readH5SCE(
      file = file,
      name = file.path(name, i),
      verbose = verbose
    )
  }
  return(altExps)
}

.read_h5_sce <- function(file, name = "/", verbose = TRUE, ...) {
  name <- h5AbsLinkName(name = name)
  slot.names <- h5List(x = file, name = name)
  out <- list()
  to.read <- c("rowData", "colData")
  for (i in to.read) {
    verboseMsg("Reading '", i, "'")
    out[[i]] <- h5Read(x = file, name = file.path(name, i))
  }

  out$assays <- .read_h5_sce_assays(
    file = file,
    rownames = rownames(x = out$rowData),
    colnames = rownames(x = out$colData),
    name = file.path(name, "assays"),
    verbose = verbose
  )

  out$altExps <- .read_h5_sce_Exps(
    file = file,
    name = file.path(name, .alt_key),
    verbose = verbose
  )

  to.read <- c(
    "rowRanges",
    .red_key,
    .colp_key,
    .rowp_key,
    .sinfo_key,
    .annot_key,
    "metadata"
  )
  to.read <- intersect(x = to.read, y = slot.names)
  for (i in to.read) {
    verboseMsg("Reading ", i)
    out[[i]] <- h5Read(
      x = file,
      name = file.path(name, i),
      toS4.func = BPCE_listToS4
    )
  }
  return(out)
}

.list_to_SCE <- function(rlist) {
  out <- SingleCellExperiment(
    assays = rlist$assays,
    rowData = rlist$rowData,
    colData = rlist$colData,
    reducedDims = rlist[[.red_key]],
    colPairs = rlist[[.colp_key]],
    rowPairs = rlist[[.rowp_key]],
    altExps = rlist[[.alt_key]],
    metadata = rlist$metadata
  )
  if (length(x = rlist$rowRanges) == 0) {
    return(out)
  }
  rowRanges(x = out) <- rlist$rowRanges
  rowData(x = out) <- DataFrame(rlist$rowData)
  return(out)
}

## ChromExperiment #############################################################

#' @importFrom hdf5r.Extra h5AbsLinkName h5List
.read_h5_fragments <- function(
    file,
    chrs,
    cells,
    name = "/fragments",
    verbose = TRUE,
    ...
) {
  name <- h5AbsLinkName(name = name)
  slot.names <- h5List(x = file, name = name)
  if (length(x = slot.names) == 0) {
    return(NULL)
  }
  verboseMsg("Reading fragments")
  open_fragments_hdf5(path = file, group = name) %>%
    select_cells(cell_selection = cells) %>%
    select_chromosomes(chromosome_selection = chrs)
}

.list_to_CSCE <- function(rlist) {
  out <- ChromExperiment(
    assays = rlist$assays,
    rowData = rlist$rowData,
    colData = rlist$colData,
    reducedDims = rlist[[.red_key]],
    colPairs = rlist[[.colp_key]],
    rowPairs = rlist[[.rowp_key]],
    altExps = rlist[[.alt_key]],
    metadata = rlist$metadata,
    fragments = rlist[[.frag_key]],
    annotations = rlist[[.annot_key]],
    ranges = rlist$rowRanges
  )
  seqinfo(x = out) <- rlist[[.sinfo_key]]
  return(out)
}

## SingleCellMultiExperiment ###################################################

.read_h5_scme <- function(file, name = "/", verbose = TRUE, ...) {
  name <- h5AbsLinkName(name = name)
  slot.names <- h5List(x = file, name = name)
  out <- list()
  out$experiments <- .read_h5_sce_Exps(
    file = file,
    name = file.path(name, "experiments"),
    key = "experiments",
    verbose = verbose
  )
  to.read <- c(
    "sampleMap",
    "colData",
    "defaultExp",
    .colp_key,
    .red_key,
    "metadata"
  )
  for (i in to.read) {
    verboseMsg("Reading '", i, "'")
    out[[i]] <- h5Read(x = file, name = file.path(name, i))
  }
  to.read <- intersect(x = to.read, y = slot.names)
  for (i in to.read) {
    verboseMsg("Reading ", i)
    out[[i]] <- h5Read(
      x = file,
      name = file.path(name, i),
      toS4.func = BPCE_listToS4
    )
  }
  return(out)
}

#' @importFrom MultiAssayExperiment ExperimentList
.list_to_SCME <- function(rlist) {
  if (is_true(rlist$defaultExp == 0)) {
    rlist$defaultExp <- NA
  }
  out <- SingleCellMultiExperiment(
    experiments = ExperimentList(rlist$experiments),
    sampleMap = listToMap(listmap = rlist$sampleMap),
    colData = rlist$colData,
    defaultExp = rlist$defaultExp,
    reducedDims = rlist[[.red_key]],
    colPairs = rlist[[.colp_key]],
    metadata = rlist$metadata
  )
  return(out)
}
