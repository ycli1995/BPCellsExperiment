#' @importFrom utils getS3method
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 Methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## h5Prep ######################################################################

# Re-export hdf5r.Extra generic
#' @importFrom hdf5r.Extra h5Prep
#' @export
hdf5r.Extra::h5Prep

#' @importClassesFrom S4Vectors DataFrame
#' @export
#' @method h5Prep DataFrame
h5Prep.DataFrame <- function(x, ...) {
  return(as.data.frame(x = x, ...))
}

#' @importClassesFrom S4Vectors SimpleList
#' @export
#' @method h5Prep SimpleList
h5Prep.SimpleList <- function(x, ...) {
  return(as.list(x = x, ...))
}

#' @importFrom SingleCellExperiment factorData featureLoadings
#' LinearEmbeddingMatrix sampleFactors
#' @importFrom S4Vectors metadata
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
#' @export
#' @method h5Prep LinearEmbeddingMatrix
h5Prep.LinearEmbeddingMatrix <- function(x, ...) {
  return(list(
    sampleFactors = sampleFactors(x = x),
    featureLoadings = featureLoadings(x = x),
    factorData = as.data.frame(x = factorData(x = x)),
    metadata = metadata(x = x)
  ))
}

#' @importFrom GenomeInfoDb seqinfo
#' @importClassesFrom GenomicRanges GRanges
#' @export
#' @method h5Prep GRanges
h5Prep.GRanges <- function(x, ...) {
  df <- as.data.frame(x = x)
  sinfo.df <- seqinfo(x = x) %>%
    h5Prep()
  return(list(granges = df, seqinfo = sinfo.df))
}

#' @importClassesFrom GenomeInfoDb Seqinfo
#' @export
#' @method h5Prep Seqinfo
h5Prep.Seqinfo <- function(x, ...) {
  sinfo.df <- as.data.frame(x = x)
  sinfo.df <- sinfo.df[, colSums(x = is.na(x = sinfo.df)) == 0, drop = FALSE]
  return(sinfo.df)
}

#' @importFrom S4Vectors mcols<- metadata
#' @importFrom SummarizedExperiment assays colData rowData rowRanges
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @export
#' @method h5Prep RangedSummarizedExperiment
h5Prep.RangedSummarizedExperiment <- function(x, ...) {
  rr <- rowRanges(x = x)
  mcols(x = rr) <- NULL
  if (!inherits(x = rr, what = "GRanges")) {
    rr <- list()
  }
  return(list(
    rowRanges = rr,
    rowData = as.data.frame(x = rowData(x = x)),
    colData = as.data.frame(x = colData(x = x)),
    assays = assays(x = x),
    metadata = metadata(x = x)
  ))
}

#' @importFrom SingleCellExperiment altExps colPairs reducedDims rowPairs
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @method h5Prep SingleCellExperiment
h5Prep.SingleCellExperiment <- function(x, ...) {
  old_func <- getS3method(f = "h5Prep", class = "RangedSummarizedExperiment")
  out <- old_func(x = x, ...)
  out[[.alt_key]] <- altExps(x = x)
  out[[.red_key]] <- reducedDims(x = x)
  out[[.colp_key]] <- colPairs(x = x, asSparse = FALSE)
  for (i in seq_along(out[[.colp_key]])) {
    out[[.colp_key]][[i]] <- out[[.colp_key]][[i]] %>%
      as(Class = "dgCMatrix")
    gc(verbose = FALSE)
  }
  out[[.rowp_key]] <- rowPairs(x = x, asSparse = FALSE)
  for (i in seq_along(out[[.rowp_key]])) {
    out[[.rowp_key]][[i]] <- out[[.rowp_key]][[i]] %>%
      as(Class = "dgCMatrix")
    gc(verbose = FALSE)
  }
  return(out)
}

#' @importFrom GenomeInfoDb seqinfo
#' @export
#' @method h5Prep ChromExperiment
h5Prep.ChromExperiment <- function(x, ...) {
  old_func <- getS3method(f = "h5Prep", class = "SingleCellExperiment")
  out <- old_func(x = x, ...)
  out[[.frag_key]] <- fragments(x = x)
  out[[.annot_key]] <- annotations(x = x)
  out[[.sinfo_key]] <- seqinfo(x = x)
  return(out)
}

#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment colPairs reducedDims rowPairs
#' @importFrom MultiAssayExperiment experiments mapToList sampleMap
#' @export
#' @method h5Prep SingleCellMultiExperiment
h5Prep.SingleCellMultiExperiment <- function(x, ...) {
  smap <- as.list(x = mapToList(dfmap = sampleMap(x = x)))
  exps <- experiments(x = x)
  for (i in names(x = exps)) {
    rownames(x = smap[[i]]) <- smap[[i]][["colname"]]
    smap[[i]] <- smap[[i]][colnames(x = exps[[i]]), ]
  }
  defaultExp <- defaultExp(x = x)
  if (is.na(x = defaultExp)) {
    defaultExp <- 0
  }
  out <- list(
    experiments = as.list(x = exps),
    colData = as.data.frame(x = colData(x = x)),
    sampleMap = smap,
    defaultExp = defaultExp,
    metadata =  metadata(x = x)
  )
  out[[.red_key]] <- reducedDims(x = int_SCE(x = x))
  out[[.colp_key]] <- colPairs(x = int_SCE(x = x), asSparse = FALSE)
  for (i in seq_along(out[[.colp_key]])) {
    out[[.colp_key]][[i]] <- out[[.colp_key]][[i]] %>%
      as(Class = "dgCMatrix")
    gc(verbose = FALSE)
  }
  return(out)
}

## h5Write #####################################################################

# Re-export hdf5r.Extra generic
#' @importFrom hdf5r.Extra h5Copy h5Write
#' @export
hdf5r.Extra::h5Write

#' @importFrom hdf5r h5garbage_collect
#' @importFrom hdf5r.Extra h5AbsLinkName h5CreateFile h5CreateGroup h5Copy
#' h5Overwrite h5Write
#' @export
#' @method h5Write IterableMatrix
h5Write.IterableMatrix <- function(
    x,
    file,
    name,
    overwrite = FALSE,
    gzip_level = 0L,
    ...
) {
  name <- h5AbsLinkName(name = name)
  file <- normalizePath(path = file, mustWork = FALSE)
  .check_bpcells_file_occupy(x = x, path = file, name = name)
  file <- h5Overwrite(file = file, name = name, overwrite = overwrite)
  message(
    "Writing a IterableMatrix:",
    "\n Class: ", class(x = x),
    "\n File: ", file,
    "\n Group: ", name
  )
  x <- write_matrix_hdf5(
    mat = x,
    path = file,
    group = name,
    gzip_level = gzip_level,
    ...
  )
  gc(verbose = FALSE)
  return(invisible(x = NULL))
}

#' @importFrom hdf5r.Extra h5AbsLinkName h5CreateGroup h5Copy h5Overwrite
#' h5Write
#' @importFrom BPCells write_fragments_hdf5
#' @export
#' @method h5Write IterableFragments
h5Write.IterableFragments <- function(
    x,
    file,
    name,
    overwrite = FALSE,
    gzip_level = 0L,
    ...
) {
  name <- h5AbsLinkName(name = name)
  .check_bpcells_file_occupy(x = x, path = file, name = name)
  file <- h5Overwrite(file = file, name = name, overwrite = overwrite)
  message(
    "Writing a IterableFragments:",
    "\n Class: ", class(x = x),
    "\n File: ", file,
    "\n Group: ", name
  )
  x <- write_fragments_hdf5(
    fragments = x,
    path = file,
    group = name,
    gzip_level = gzip_level,
    ...
  )
  gc(verbose = FALSE)
  return(invisible(x = NULL))
}

## writeH5SCE ##################################################################

#' @export
#' @rdname writeH5SCE
#' @method writeH5SCE default
writeH5SCE.default <- function(object, file, ...) {
  stop("Writing ", class(x = object), " to HDF5 file is not supported.")
}

#' @param name Name of the HDF5 group where to save \code{object}.
#' @param overwrite Whether or not to overwrite the existing HDF5 link.
#' @param gzip_level Enable zipping at the level given here.
#' @param verbose Display progress.
#'
#' @importFrom SingleCellExperiment altExps altExps<- int_metadata<-
#' @importFrom SummarizedExperiment assays<-
#' @importFrom hdf5r.Extra h5AbsLinkName
#' @export
#' @rdname writeH5SCE
#' @method writeH5SCE SingleCellExperiment
writeH5SCE.SingleCellExperiment <- function(
    object,
    file,
    name = "/",
    overwrite = FALSE,
    gzip_level = 0L,
    verbose = TRUE,
    ...
) {
  name <- h5AbsLinkName(name = name)
  out.file <- normalizePath(path = file, mustWork = FALSE)
  file <- .h5_overwrite_before(
    file = file,
    name = name,
    overwrite = overwrite
  )
  if (!identical(x = file, y = out.file)) {
    on.exit(expr = unlink(x = file, force = TRUE))
  }
  .write_h5_sce(
    object = object,
    file = file,
    name = name,
    gzip_level = gzip_level,
    verbose = verbose
  )
  .set_h5attr_sce(object = object, file = file, name = name, verbose = verbose)
  file <- file_path_as_absolute(x = file)
  out.file <- out.file %>%
    file_path_as_absolute() %>%
    .h5_overwrite_after(file = file, name = name)

  verboseMsg("Updating altExps")
  altExps(x = object) <- .read_h5_sce_Exps(
    file = out.file,
    name = file.path(name, .alt_key),
    verbose = verbose
  )

  verboseMsg("Updating assays")
  assays(x = object) <- .read_h5_sce_assays(
    file = out.file,
    name = file.path(name, "assays"),
    rownames = rownames(x = object),
    colnames = colnames(x = object),
    verbose = verbose
  )
  int_metadata(x = object)[[.path_key]] <- c(out.file, name)
  return(object)
}

#' @export
#' @rdname writeH5SCE
#' @method writeH5SCE ChromExperiment
writeH5SCE.ChromExperiment <- function(
    object,
    file,
    name = "/",
    overwrite = FALSE,
    gzip_level = 0L,
    verbose = TRUE,
    ...
) {
  old_func <- getS3method(f = "writeH5SCE", class = "SingleCellExperiment")
  object <- old_func(
    object = object,
    file = file,
    name = name,
    overwrite = overwrite,
    gzip_level = gzip_level,
    verbose = verbose,
    ...
  )
  verboseMsg("Updating fragments")
  fragments(x = object) <- .read_h5_fragments(
    file = file,
    cells = colnames(x = object),
    chrs = seqlevels(x = rowRanges(x = object)),
    name = file.path(name, "fragments"),
    verbose = verbose
  )
  return(object)
}

## writeH5SCME #################################################################

#' @export
#' @rdname writeH5SCME
#' @method writeH5SCME default
writeH5SCME.default <- function(object, file, ...) {
  stop("Writing ", class(x = object), " to HDF5 file is not supported.")
}

#' @param name Name of the HDF5 group where to save \code{object}.
#' @param overwrite Whether or not to overwrite the existing HDF5 link.
#' @param gzip_level Enable zipping at the level given here.
#' @param verbose Display progress.
#'
#' @importFrom hdf5r.Extra h5AbsLinkName
#' @importFrom MultiAssayExperiment ExperimentList
#' @export
#' @rdname writeH5SCME
#' @method writeH5SCME SingleCellMultiExperiment
writeH5SCME.SingleCellMultiExperiment <- function(
    object,
    file,
    name = "/",
    overwrite = FALSE,
    gzip_level = 0L,
    verbose = TRUE,
    ...
) {
  name <- h5AbsLinkName(name = name)
  out.file <- file
  file <- .h5_overwrite_before(
    file = file,
    name = name,
    overwrite = overwrite
  )
  if (!identical(x = file, y = out.file)) {
    on.exit(expr = unlink(x = file, force = TRUE))
  }
  .write_h5_scme(
    object = object,
    file = file,
    name = name,
    gzip_level = gzip_level,
    verbose = verbose
  )
  .set_h5attr_sce(object = object, file = file, name = name, verbose = verbose)
  file <- file_path_as_absolute(x = file)
  out.file <- out.file %>%
    file_path_as_absolute() %>%
    .h5_overwrite_after(file = file, name = name)

  verboseMsg("Updating experiments")
  exps <- .read_h5_sce_Exps(
    file = out.file,
    name = file.path(name, "experiments"),
    key = "experiment",
    verbose = verbose,
    ...
  )
  experiments(object = object) <- ExperimentList(exps)
  int_metadata(x = object)[[.path_key]] <- c(out.file, name)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Attributes ##################################################################

.get_sce_class <- function(object) {
  if (inherits(x = object, what = "SingleCellExperiment")) {
    return("SCE")
  }
  if (inherits(x = object, what = "SingleCellMultiExperiment")) {
    return("SCME")
  }
  stop(
    "Invalid class(object): ", class(x = object),
    "\n 'object' must inherits from either ",
    "SingleCellExperiment or SingleCellMultiExperiment"
  )
}

## SingleCellExperiment ########################################################

#' Check 'overwrite' before writing SCE into HDF5 file
#'
#' @return A character specifying the target file where to be actually written.
#'
#' @details
#' First check whether \code{file} already exists. If not, just return the
#' original \code{file}.
#' When \code{file} exists, check whether \code{name} exists. If not, also just
#' return the orignial \code{file}
#' When both \code{file} and \code{name} exist, check \code{overwrite}. If not,
#' raise an error.
#' When actually do \code{overwrite}, generate a temporary HDF5 file name
#' (without creating file) to hold the written data.
#'
#' @noRd
#' @importFrom hdf5r.Extra h5AbsLinkName h5Backup h5CreateFile h5Exists
.h5_overwrite_before <- function(file, name, overwrite) {
  name <- h5AbsLinkName(name = name)
  if (!file.exists(file)) {
    return(normalizePath(path = file, mustWork = FALSE))
  }
  file <- normalizePath(path = file)
  tmp.file <- tempfile(tmpdir = dirname(path = file), fileext = ".h5")
  if (!h5Exists(x = file, name = name)) {
    message(
      "HDF5 file already exists. Use temporary file to write:",
      "\n  File: ", file,
      "\n  Object: ", name,
      "\n  Temporary file: ", tmp.file
    )
    return(tmp.file)
  }
  if (!overwrite) {
    stop(
      "\nFound object that already exists: ",
      "\n  File: ", file,
      "\n  Object: ", name,
      "\nSet 'overwrite = TRUE' to remove it."
    )
  }
  message(
    "HDF5 object already exists. Use temporary file to overwrite:",
    "\n  File: ", file,
    "\n  Object: ", name,
    "\n  Temporary file: ", tmp.file
  )
  return(normalizePath(path = tmp.file, mustWork = FALSE))
}

#' @importFrom hdf5r.Extra h5WriteAttr
.set_h5attr_sce <- function(object, file, name = "/", verbose = TRUE) {
  attrs <- list(
    S4Class = class(x = object)[1],
    SCE_Class = .get_sce_class(object = object),
    Version = "v0.0.1"
  )
  verboseMsg("Add H5 attribute: ")
  for (i in seq_along(along.with = attrs)) {
    which <- names(x = attrs)[i]
    verboseMsg(" ", which, ": ", attrs[[i]])
    h5WriteAttr(x = file, which = which, robj = attrs[[i]], name = name)
  }
  return(invisible(x = NULL))
}

#' @importFrom hdf5r.Extra h5CreateGroup h5Prep h5Write
.write_h5_sce <- function(
    object,
    file,
    name = "/",
    gzip_level = 0L,
    verbose = TRUE,
    ...
) {
  object <- h5Prep(x = object)

  altExps <- object$altExps
  altExpNames <- names(x = altExps)
  h5CreateGroup(
    x = file,
    name = file.path(name, .alt_key),
    show.warnings = FALSE
  )
  for (i in seq_along(altExps)) {
    name_i <- altExpNames[i]
    verboseMsg("Writing altExp: ", name_i)
    altExps[[i]] <- writeH5SCE(
      object = altExps[[i]],
      file = file,
      name = file.path(name, .alt_key, name_i),
      gzip_level = gzip_level,
      verbose = verbose,
      overwrite = FALSE,
      ...
    )
    verboseMsg("\n")
  }
  object$altExps <- NULL

  assays <- object$assays
  assay.names <- names(x = assays)
  h5CreateGroup(
    x = file,
    name = file.path(name, "assays"),
    show.warnings = FALSE
  )
  for (i in seq_along(assays)) {
    verboseMsg("Writing assay '", assay.names[i], "'")
    h5Write(
      x = assays[[i]],
      file = file,
      name = file.path(name, "assays", assay.names[i]),
      overwrite = FALSE,
      gzip_level = gzip_level
    )
  }
  object$assays <- NULL

  slot.names <- names(x = object)
  for (i in seq_along(object)) {
    verboseMsg("Writing ", slot.names[i])
    h5Write(
      x = object[[i]],
      file = file,
      name = file.path(name, slot.names[i]),
      overwrite = FALSE,
      gzip_level = gzip_level
    )
  }
  gc(verbose = FALSE)
  return(invisible(x = NULL))
}

#' @importFrom hdf5r.Extra h5Copy h5Overwrite
.h5_overwrite_after <- function(file, out.file, name) {
  if (identical(x = file, y = out.file)) {
    return(out.file)
  }
  out.file <- h5Overwrite(file = out.file, name = name, overwrite = TRUE)
  h5Copy(
    from.file = file,
    from.name = name,
    to.file = out.file,
    to.name = name,
    overwrite = TRUE
  )
  return(out.file)
}

## SingleCellMultiExperiment ###################################################

#' @importFrom hdf5r.Extra h5List h5Write
.write_h5_scme <- function(
    object,
    file,
    name = "/",
    gzip_level = 0L,
    verbose = TRUE,
    ...
) {
  object <- h5Prep(x = object)
  exp.names <- names(x = object$experiments)
  for (i in exp.names) {
    verboseMsg("Writing experiment: ", i)
    object$experiments[[i]] <- writeH5SCE(
      object = object$experiments[[i]],
      file = file,
      name = file.path(name, "experiments", i),
      gzip_level = gzip_level,
      verbose = verbose
    )
    verboseMsg("Done\n")
  }
  if (!h5Exists(x = file, name = file.path(name, "experiments"))) {
    h5Write(x = list(), file = file, name = file.path(name, "experiments"))
  }
  object$experiments <- NULL
  slot.names <- names(x = object)
  for (i in seq_along(object)) {
    verboseMsg("Writing ", slot.names[i])
    h5Write(
      x = object[[i]],
      file = file,
      name = file.path(name, slot.names[i]),
      overwrite = FALSE,
      gzip_level = gzip_level
    )
  }
  gc(verbose = FALSE)
  return(invisible(x = NULL))
}
