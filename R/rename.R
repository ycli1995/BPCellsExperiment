
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods ######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Rename the dimensions
#'
#' Methods to set new dimension names to an object. These helpers are mostly
#' designed for \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} and its
#' sub-classes.
#'
#' @param x A matrix-like object or object inheriting from
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#' @param new.names Strings to set as the new dimension names. Can be a single
#' vector of characters or a list containing multiple vectors.
#' @param i For \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} object,
#' should be integers or characters specifying the experiment(s) to be renamed.
#' \itemize{
#' \item When \code{new.names} is a list, the length of \code{i} should be the
#' same as \code{new.names}. If \code{i} is \code{NA} or missing, all
#' experiments will be used by default.
#' \item When \code{new.names} is a single vector, \code{i} should also be a
#' single integer or character. If \code{i} is \code{NA} or missing, will change
#' the primary dimension names of \code{x}.
#' }
#' @param ... Arguments passed to other methods
#'
#' @name set-dimnames
NULL

## setRownames =================================================================

#' @export
#' @rdname set-dimnames
setMethod(
  f = "setRownames",
  definition = function(x, new.names, ...) {
    rownames(x = x) <- new.names
    return(x)
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom S4Vectors list_OR_List
#' @export
#' @rdname set-dimnames
setMethod(
  f = "setRownames",
  signature = c("MultiAssayExperiment", "list_OR_List"),
  definition = function(x, new.names, i, ...) {
    if (missing(x = i)) {
      x <- .set_dimname_mae(
        x = x,
        i = names(x = x),
        new.names = new.names,
        type = "feature"
      )
      return(x)
    }
    if (inherits(x = i, what = "Character_OR_Numeric")) {
      x <- .set_dimname_mae(
        x = x,
        i = i,
        new.names = new.names,
        type = "feature"
      )
    }
    return(x)
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @rdname set-dimnames
setMethod(
  f = "setRownames",
  signature = c("MultiAssayExperiment", "character"),
  definition = function(x, new.names, i, ...) {
    if (missing(x = i)) {
      i <- NA
    }
    if (inherits(x = i, what = "Character_OR_Numeric") || is.na(x = i)) {
      x <- .set_dimname_mae(
        x = x,
        i = i,
        new.names = new.names,
        type = "feature"
      )
    }
    return(x)
  }
)

## setColnames =================================================================

#' @export
#' @rdname set-dimnames
setMethod(
  f = "setColnames",
  definition = function(x, new.names, ...) {
    colnames(x = x) <- new.names
    return(x)
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom S4Vectors list_OR_List
#' @export
#' @rdname set-dimnames
setMethod(
  f = "setColnames",
  signature = c("MultiAssayExperiment", "list_OR_List"),
  definition = function(x, new.names, i, ...) {
    if (missing(x = i)) {
      x <- .set_dimname_mae(
        x = x,
        i = names(x = x),
        new.names = new.names,
        type = "sample"
      )
      return(x)
    }
    if (inherits(x = i, what = "Character_OR_Numeric")) {
      x <- .set_dimname_mae(
        x = x,
        i = i,
        new.names = new.names,
        type = "sample"
      )
    }
    return(x)
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @rdname set-dimnames
setMethod(
  f = "setColnames",
  signature = c("MultiAssayExperiment", "character"),
  definition = function(x, new.names, i, ...) {
    if (missing(x = i)) {
      i <- NA
    }
    if (inherits(x = i, what = "Character_OR_Numeric") || is.na(x = i)) {
      x <- .set_dimname_mae(
        x = x,
        i = i,
        new.names = new.names,
        type = "sample"
      )
    }
    return(x)
  }
)

#' @importFrom BiocBaseUtils setSlots
#' @export
#' @rdname set-dimnames
setMethod(
  f = "setColnames",
  signature = c("SingleCellMultiExperiment", "character"),
  definition = function(x, new.names, i, ...) {
    if (missing(x = i)) {
      i <- NA
    }
    if (is.na(x = i)) {
      int_sce <- int_SCE(x = x)
      default.exp <- defaultExp(x = x)
      x <- .set_dimname_mae(
        x = as(object = x, Class = "MultiAssayExperiment"),
        i = i,
        new.names = new.names,
        type = "sample"
      )
      colnames(x = int_sce) <- rownames(x = colData(x = x))
      return(new(
        Class = "SingleCellMultiExperiment",
        x,
        int_SCE = int_sce,
        defaultExp = default.exp
      ))
    }
    callNextMethod(x = x, new.names = new.names, i = i)
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom fastmatch fmatch
#' @importFrom BiocBaseUtils setSlots
#' @importFrom MultiAssayExperiment experiments listToMap mapToList
#' renamePrimary
.set_dimname_mae <- function(x, i, new.names, type = c("sample", "feature")) {
  type <- match.arg(arg = type)
  if (any(is.na(x = i))) {
    # set primary names
    if (type == "sample") {
      stopifnot(is.character(x = new.names))
      x <- renamePrimary(x = x, value = new.names)
    }
    return(x)
  }
  stopifnot(length(x = i) == length(x = new.names))
  margin <- switch(EXPR = type, "feature" = 1L, "sample" = 2L)
  exps <- experiments(x = x)
  dmap <- mapToList(dfmap = sampleMap(x = x))
  for (j in seq_along(i)) {
    if (type == "sample") {
      idx <- fmatch(
        x = dmap[[i[j]]][['colname']],
        table = dimnames(x = exps[[i[j]]])[[margin]]
      )
      dmap[[i[j]]][['colname']] <- new.names[[j]][idx]
    }
    dimnames(x = exps[[i[j]]])[[margin]] <- new.names[[j]]
  }
  dmap <- listToMap(listmap = dmap)
  x <- setSlots(object = x, sampleMap = dmap, ExperimentList = exps)
  return(x)
}
