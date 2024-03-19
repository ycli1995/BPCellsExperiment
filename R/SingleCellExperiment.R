
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## dietSCE #####################################################################

#' @param assays Keep a subset of \code{\link[SummarizedExperiment]{assays}}
#' specified here. If is \code{NA}, remove all assays.
#' @param reducedDims Keep a subset of dimension reductions specified here. If
#' is \code{NA}, remove all \code{\link[SingleCellExperiment]{reducedDims}}.
#' @param colPairs Keep a subset of \code{\link[SingleCellExperiment]{colPairs}}
#' specified here. If is \code{NA}, remove all \code{colPairs}.
#' @param rowPairs Similar as \code{colPairs}.
#' @param altExps Keep a subset of \code{\link[SingleCellExperiment]{altExps}}
#' specified here. If is \code{NA}, remove all \code{altExps}.
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @rdname dietSCE
#' @export
#' @method dietSCE SingleCellExperiment
dietSCE.SingleCellExperiment <- function(
    object,
    assays = NULL,
    reducedDims = NULL,
    colPairs = NULL,
    rowPairs = NULL,
    altExps = NULL,
    ...
) {
  object <- .diet_sce_field(
    object = object,
    func = "assays",
    name.func = "assayNames",
    ns = "SummarizedExperiment",
    names = assays
  )
  fields <- list(
    reducedDims = reducedDims,
    colPairs = colPairs,
    rowPairs = rowPairs,
    altExps = altExps
  )
  for (i in names(x = fields)) {
    name.func <- paste0(sub(pattern = "s$", replacement = "", x = i), "Names")
    object <- .diet_sce_field(
      object = object,
      func = i,
      name.func = name.func,
      ns = "SingleCellExperiment",
      names = fields[[i]]
    )
  }
  return(object)
}

#' @param fragments Logical scalar specifying whether to keep
#' \code{\link{fragments}}.
#' @param annotations Logical scalar specifying whether to keep
#' \code{\link{annotations}}.
#' @param seqinfo Logical scalar specifying whether to keep
#' \code{\link{seqinfo}}.
#'
#' @rdname dietSCE
#' @export
#' @method dietSCE ChromExperiment
dietSCE.ChromExperiment <- function(
    object,
    assays = NULL,
    reducedDims = NULL,
    colPairs = NULL,
    rowPairs = NULL,
    altExps = NULL,
    fragments = NULL,
    annotations = NULL,
    seqinfo = NULL,
    ...
) {
  old_func <- getS3method(f = "dietSCE", class = "SingleCellExperiment")
  out <- old_func(
    object = object,
    assays = assays,
    reducedDims = reducedDims,
    colPairs = colPairs,
    rowPairs = rowPairs,
    altExps = altExps,
    ...
  )
  if (is_false(x = fragments)) {
    fragments(x = object) <- NULL
  }
  if (is_false(x = annotations)) {
    annotations(x = object) <- NULL
  }
  if (is_false(x = seqinfo)) {
    seqinfo(x = object) <- NULL
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## fetchData ###################################################################

#' Fetch cellular data
#'
#' Fetch data for a set of observations (columns) in an object
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @param ... Arguments to be passed to other methods
#'
#' @return
#' A \code{\link[base]{data.frame}} with cells as rows and \code{vars} data as
#' columns.
#'
#' @name fetchData
NULL

#' @param vars Vector of all variables to fetch
#' @param cells Cells to collect data for (default is all cells)
#' @param assay Assay to collect data from
#' @param use.Exp Which \code{\link[SingleCellExperiment]{altExp}} to collect
#' data from. If the selected \code{altExp} doesn't exist, will use the main
#' Experiment.
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom SingleCellExperiment altExp altExpNames reducedDims
#' @importFrom SeuratObject EmptyDF
#' @importFrom easy.utils fastIntersect fetchColnames
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname fetchData
setMethod(
  f = "fetchData",
  signature = "SingleCellExperiment",
  definition = function(
    object,
    vars = NULL,
    cells = NULL,
    assay = NULL,
    use.Exp = NULL
  ) {
    cells <- cells %||% colnames(x = object)
    if (is.numeric(x = cells)) {
      cells <- colnames(x = object)[cells]
    }
    cells <- fetchColnames(object = object, query = cells)
    data.fetched <- data.frame(row.names = cells)
    if (length(x = vars) == 0) {
      return(data.fetched)
    }
    # Check altExp in object
    if (length(x = use.Exp) > 1) {
      warning(
        "Only use the first one of 'use.Exp': ", use.Exp[1],
        call. = FALSE, immediate. = TRUE
      )
      use.Exp <- use.Exp[1]
    }
    if (is.null(x = use.Exp)) {
      assay <- assay %||% .select_default_assay(x = object)
    }
    if (any(use.Exp %in% altExpNames(x = object))) {
      alt.exp <- altExp(x = object, e = use.Exp)
      assay <- assay %||% .select_default_assay(x = alt.exp)
      assay.data <- assay(x = alt.exp, i = assay)
    } else {
      if (length(x = use.Exp) > 0) {
        warning(
          "Ignore the non-existing Exp '", use.Exp, "'",
          call. = FALSE, immediate. = TRUE
        )
        use.Exp <- NULL
      }
      assay <- assay %||% .select_default_assay(x = object)
      assay.data <- assay(x = object, i = assay)
    }

    raw.vars <- vars
    # Find vars in assays
    assay.vars <- fastIntersect(x = vars, y = rownames(x = assay.data))
    if (length(x = assay.vars) > 0) {
      assay.fetched <- assay.data[assay.vars, cells, drop = FALSE] %>%
        as.matrix()
      data.fetched <- cbind(data.fetched, t(x = assay.fetched))
      vars <- setdiff(x = vars, y = assay.vars)
    }
    # Find vars in colData
    colData.vars <- fastIntersect(x = vars, y = names(x = colData(x = object)))
    if (length(x = colData.vars) > 0) {
      colData.fetched <- colData(x = object)[cells, colData.vars, drop = FALSE]
      colData.fetched <- as.data.frame(x = colData.fetched, optional = TRUE)
      data.fetched <- cbind(data.fetched, colData.fetched)
      vars <- setdiff(x = vars, colData.vars)
    }
    # Find all vars in reducedDims(object)
    # (The same as Keys in SeuratObject)
    all.reducs <- reducedDims(x = object)
    for (i in names(x = all.reducs)) {
      use.reduc <- all.reducs[[i]]
      # First find vars already in colnames
      rd.vars <- fastIntersect(x = vars, y = colnames(x = use.reduc))
      if (length(x = rd.vars) > 0) {
        rd.fetched <- use.reduc[cells, rd.vars, drop = FALSE] %>%
          as.matrix()
        data.fetched <- cbind(data.fetched, rd.fetched)
      }
      vars <- setdiff(x = vars, y = rd.vars)
      # Then find vars in format 'Key_1', 'Key_2'...
      colnames(x = use.reduc) <- paste0(i, "_", 1:ncol(x = use.reduc))
      rd.vars <- fastIntersect(x = vars, y = colnames(x = use.reduc))
      if (length(x = rd.vars) > 0) {
        rd.fetched <- use.reduc[cells, rd.vars, drop = FALSE] %>%
          as.matrix()
        data.fetched <- cbind(data.fetched, rd.fetched)
      }
      vars <- setdiff(x = vars, y = rd.vars)
    }
    if (length(x = vars) > 0) {
      warning(
        "The following query variables were not found: \n  ",
        paste(vars, collapse = ", "),
        immediate. = TRUE
      )
    }
    final.vars <- fastIntersect(x = raw.vars, y = colnames(x = data.fetched))
    data.fetched <- data.fetched[, final.vars, drop = FALSE]
    return(data.fetched)
  }
)

#' @importFrom SummarizedExperiment rowData
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname fetchData
setMethod(
  f = "variableFeatures",
  signature = "SingleCellExperiment",
  definition = function(object, ...) {
    return(rownames(x = object)[rowData(x = object)[['variable']]])
  }
)

#' @param value An object of a class specified in the S4 method signature.
#'
#' @importFrom SummarizedExperiment rowData<-
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname fetchData
setMethod(
  f = "variableFeatures<-",
  signature = c("SingleCellExperiment", "character"),
  definition = function(object, ..., value) {
    rowData(x = object)[['variable']] <- rownames(x = object) %in% value
    return(object)
  }
)

#' @importFrom SummarizedExperiment rowData<-
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname fetchData
setMethod(
  f = "variableFeatures<-",
  signature = c("SingleCellExperiment", "NULL"),
  definition = function(object, ..., value) {
    rowData(x = object)[['variable']] <- rownames(x = object) %in% value
    return(object)
  }
)

## Backend path ################################################################

#' @importFrom BiocGenerics path
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod(
  f = "path",
  signature = "SingleCellExperiment",
  definition = function(object, ...) {
    return(int_metadata(x = object)[[.path_key]])
  }
)

## show ########################################################################

#' @importFrom methods getMethod show
setMethod(
  f = "show",
  signature = "SingleCellExperiment",
  definition = function(object) {
    show_func <- getMethod(
      f = "show",
      signature = "SingleCellExperiment",
      where = asNamespace(ns = "SingleCellExperiment")
    )
    show_func(object)
    .show_h5backed_path(object = object)
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom BiocGenerics path
.show_h5backed_path <- function(object) {
  p <- path(object = object)
  cat("Backend HDF5 file:\n")
  cat(" file: ", p[1], "\n")
  cat(" group: ", p[2], "\n")
  return(invisible(x = NULL))
}

#' Select default assay from SummarizedExperiment-like object
#'
#' Default order of priority is: logcounts, normcounts, counts, scaled
#'
#' @importFrom SummarizedExperiment assayNames
#' @noRd
.select_default_assay <- function(x) {
  default_assays <- c("logcounts", "normcounts", "counts", "scaled")
  all_assays <- assayNames(x = x)
  ret_assay <- intersect(x = default_assays, y = all_assays)
  if (length(x = ret_assay) > 0) {
    return(ret_assay[1])
  }
  stop(
    "\n  The SCE didn't contain any of the following assay(s): ",
    paste(default_assays, collapse = ", "), "\n",
    "\n  Cannot automatically select default assay."
  )
}

.diet_sce_field <- function(object, func, name.func, ns, names = NULL) {
  get.func <- getExportedValue(ns = ns, name = func)
  set.func <- getExportedValue(ns = ns, name = paste0(func, "<-"))
  getname.func <- getExportedValue(ns = ns, name = name.func)
  curr.names <- getname.func(object)
  if (length(x = curr.names) == 0) {
    return(object)
  }
  if (is.null(x = names)) {
    return(object)
  }
  if (is.na(x = names)) {
    message("Remove all ", func)
    return(set.func(object, value = list()))
  }
  old.names <- names
  names <- intersect(x = old.names, y = curr.names)
  if (length(x = names) == 0) {
    warning(
      "None of the following ", name.func, " is found: \n  ", old.names,
      call. = FALSE, immediate. = TRUE
    )
    return(object)
  }
  return(set.func(object, value = get.func(object)[[names]]))
}
