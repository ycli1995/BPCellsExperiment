
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 Methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param name Name of the HDF5 group where to save \code{object}.
#' @param X Which \code{\link[SummarizedExperiment]{assay}} to be used as
#' `adata.X`. Default is the first one.
#' @param layers Which \code{\link[SummarizedExperiment]{assays}} to be used in
#' `adata.layers`. Default will be all assays except the `X` one.
#' @param obsm Which \code{\link[SingleCellExperiment]{reducedDim}} to be used
#' in `adata.obsm`. Default will be all reduction data.
#' @param obsp Which \code{\link[SingleCellExperiment]{colPair}} to be used in
#' `adata.obsp`.
#' @param varp Which \code{\link[SingleCellExperiment]{rowPair}} to be used in
#' `adata.varp`.
#' @param raw Which \code{\link[SingleCellExperiment]{altExp}} to be used in
#' `adata.raw`. The `nrows()` of `raw` must be larger than the `mainExp`.
#' @param overwrite Whether or not to overwrite the existing HDF5 link.
#' @param gzip_level Enable zipping at the level given here.
#' @param verbose Display progress.
#'
#' @export
#' @rdname exportH5AD
#' @method exportH5AD SingleCellExperiment
exportH5AD.SingleCellExperiment <- function(
    object,
    file,
    name = "/",
    X = NULL,
    layers = NULL,
    obsm = NULL,
    obsp = NULL,
    varp = NULL,
    uns = NULL,
    raw = NULL,
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
  object %>%
    .prep_h5ad(
      X = X,
      layers = layers,
      obsm = obsm,
      obsp = obsp,
      varp = varp,
      uns = uns,
      raw = raw
    ) %>%
    .write_h5ad_sce(
      file = file,
      name = name,
      gzip_level = gzip_level,
      verbose = verbose
    )
  file <- file_path_as_absolute(x = file)
  out.file <- out.file %>%
    file_path_as_absolute() %>%
    .h5_overwrite_after(file = file, name = name)
  return(invisible(x = NULL))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SummarizedExperiment assayNames colData
#' @importFrom SingleCellExperiment altExp altExpNames colPair colPairNames
#' reducedDim rowPair rowPairNames
.prep_h5ad <- function(
    object,
    X = NULL,
    layers = NULL,
    obsm = NULL,
    obsp = NULL,
    varp = NULL,
    uns = NULL,
    raw = NULL,
    raw_X = NULL
) {
  X <- X %||% assayNames(x = object)[1]
  out <- .prep_h5ad_raw(object = object, X = X)

  layers <- layers %||% assayNames(x = object)
  layers <- layers %>%
    intersect(y = assayNames(x = object)) %>%
    setdiff(y = X)
  out$layers <- .prep_h5ad_layers(object = object, layers = layers)

  out$obs <- object %>%
    colData() %>%
    as.data.frame()

  obsm <- obsm %||% reducedDimNames(x = object)
  obsm <- intersect(x = obsm, y = reducedDimNames(x = object))
  out$obsm <- list()
  for (i in obsm) {
    out$obsm[[paste0("X_", i)]] <- object %>%
      reducedDim(type = i, withDimnames = TRUE) %>%
      as.matrix()
  }

  obsp <- obsp %||% colPairNames(x = object)
  obsp <- intersect(x = obsp, y = colPairNames(x = object))
  out$obsp <- list()
  for (i in obsp) {
    out$obsm[[i]] <- colPair(x = object, type = i, asSparse = TRUE)
  }

  varp <- varp %||% rowPairNames(x = object)
  varp <- intersect(x = varp, y = rowPairNames(x = object))
  out$varp <- list()
  for (i in varp) {
    out$varp[[i]] <- rowPair(x = object, type = i, asSparse = TRUE)
  }

  raw <- raw[1] %||% altExpNames(x = object)[1]
  raw <- intersect(x = raw, y = altExpNames(x = object))
  if (length(x = raw) > 0) {
    raw_exp <- altExp(x = object, e = raw)
    if (nrows(x = raw_exp) >= nrows(x = object)) {
      out$raw <- .prep_h5ad_raw(object = raw_exp, X = raw_X)
    }
  }

  uns <- uns %||% names(x = metadata(x = object))
  uns <- intersect(x = uns, y = names(x = metadata(x = object)))
  out$uns <- metadata(x = object)[uns]

  return(out)
}

.prep_h5ad_layers <- function(object, layers = NULL) {
  out <- list()
  layers <- intersect(x = layers, y = assayNames(x = object))
  if (length(x = layers) == 0) {
    return(out)
  }
  for (i in layers) {
    out[[i]] <- assay(x = object, i = i, withDimnames = TRUE)
    if (inherits(x = out[[i]], what = "dgCMatrix")) {
      out[[i]] <- as(object = out[[i]], Class = "IterableMatrix")
    }
  }
  return(out)
}

.prep_h5ad_raw <- function(object, X = NULL) {
  out <- list(varm = list())
  out$var <- object %>%
    rowData() %>%
    as.data.frame()
  X <- X %||% assayNames(x = object)[1]
  if (length(x = X) == 0) {
    return(out)
  }
  out$X <- .prep_h5ad_layers(object = object, layers = X)[[1]]
  return(out)
}

.write_h5ad_sce <- function(
    object,
    file,
    name = "/",
    gzip_level = 0L,
    verbose = TRUE,
    ...
) {
  .write_raw_h5ad(
    object = object,
    file = file,
    name = name,
    gzip_level = gzip_level,
    verbose = verbose
  )
  for (i in names(x = object$layers)) {
    verboseMsg("Writing layer '", i, "'")
    .write_mat_h5ad(
      mat = object$layers[[i]],
      path = file,
      name = file.path(name, "layers", i),
      gzip_level = gzip_level
    )
  }
  if ("raw" %in% names(x = object)) {
    verboseMsg("Writing 'raw'")
    .write_raw_h5ad(
      object = object$raw,
      file = file,
      name = file.path(name, "raw"),
      gzip_level = gzip_level,
      verbose = verbose
    )
  }
  verboseMsg("Writing 'obs'")
  h5Write(
    x = object$obs,
    file = file,
    name = file.path(name, "obs"),
    overwrite = TRUE,
    gzip_level = gzip_level
  )
  for (i in c("obsm", "obsp", "varm", "varp", "uns")) {
    verboseMsg("Writing ", i)
    h5Write(
      x = object[[i]],
      file = file,
      name = file.path(name, i),
      gzip_level = gzip_level
    )
  }
  return(invisible(x = NULL))
}

.write_mat_h5ad <- function(
    mat,
    file,
    name = "X",
    gzip_level = 0L,
    verbose = FALSE
) {
  name <- h5AbsLinkName(name = name)
  if (inherits(x = mat, what = "IterableMatrix")) {
    mat <- write_matrix_anndata_hdf5(
      mat = mat,
      path = file,
      group = name,
      gzip_level = gzip_level
    )
    # h5Delete(x = file, name = file.path(dirname(path = name), "obs"))
    # h5Delete(x = file, name = file.path(dirname(path = name), "var"))
    return(invisible(x = NULL))
  }
  if (is_sparse(x = mat)) {
    mat <- writeTENxMatrix(
      x = mat,
      filepath = file,
      group = name,
      level = gzip_level,
      verbose = verbose
    )
    # h5Delete(x = file, name = file.path(name, "shape"))
    # h5Delete(x = file, name = file.path(name, "barcodes"))
    # h5Delete(x = file, name = file.path(name, "genes"))
    return(invisible(x = NULL))
  }
  mat <- writeHDF5Array(
    x = mat,
    filepath = file,
    name = name,
    level = gzip_level,
    with.dimnames = FALSE,
    verbose = verbose
  )
  return(invisible(x = NULL))
}

.write_raw_h5ad <- function(
    object,
    file,
    name = "raw",
    gzip_level = 0L,
    verbose = FALSE
) {
  name <- h5AbsLinkName(name = name)
  if ("X" %in% names(x = object)) {
    verboseMsg("Writing 'X'")
    .write_mat_h5ad(
      mat = object$X,
      file = file,
      name = file.path(name, "X"),
      gzip_level = gzip_level,
      verbose = verbose
    )
  }
  verboseMsg("Writing 'var'")
  h5Write(
    x = object$var,
    file = file,
    name = file.path(name, "var"),
    overwrite = TRUE,
    gzip_level = gzip_level
  )
  return(invisible(x = NULL))
}
