
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constructor ##################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' A single-cell MultiAssayExperiment class
#'
#' The \code{SingleCellMultiExperiment} class is designed to represent
#' single-cell multi-omics (or samples) data. It inherits from the
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} class and is used in
#' the same manner. In addition, the class supports storage of single-cell
#' specific data fields, such as \code{\link[SingleCellExperiment]{reducedDims}}
#' and \code{\link[SingleCellExperiment]{colPairs}}
#'
#' @param ... Arguments passed to the constructor function of
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}, to fill the slots
#' of the base class.
#' @param reducedDims,colPairs Arguments passed to the constructor of
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}}. The manners of
#' dimension alignments are also the same.
#' @param defaultExp String defining the default \code{\link{experiment}}. If
#' \code{NULL}, will set the first one by default.
#'
#' @details
#' In this class, the single-cell specific fields (\code{reducedDims} and
#' \code{colPairs}) are actually stored in \code{\link{int_SCE}} slot, which is
#' an internal \code{\link[SingleCellExperiment]{SingleCellExperiment}}. The
#' column names of \code{int_SCE} should be aligned with the row names of
#' \code{\link[SummarizedExperiment]{colData}} in the base
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} object.
#'
#' @return
#' A \code{SingleCellMultiExperiment} object.
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment colPairs reducedDims reducedDims<-
#' SingleCellExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @docType class
SingleCellMultiExperiment <- function(
    ...,
    reducedDims = list(),
    colPairs = list(),
    defaultExp = NULL
) {
  mae <- MultiAssayExperiment(...)
  defaultExp <- defaultExp %||% names(x = mae)[1]
  cdata <- DataFrame(row.names = rownames(x = colData(x = mae)))
  sce <- SingleCellExperiment(colData = cdata)
  reducedDims(x = sce) <- reducedDims
  colPairs(x = sce) <- colPairs
  return(.mae_to_scme(mae = mae, int_SCE = sce, defaultExp = defaultExp))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Accessing and modifying information in SingleCellMultiExperiment
#'
#' A set of getter and setter functions to extract the single-cell specific data
#' fields for a \code{\link{SingleCellMultiExperiment}} object.
#'
#' @param x A \code{\link{SingleCellMultiExperiment}} object
#' @param type,withDimnames,asSparse Arguments passed to methods for
#' \code{SingleCellExperiment} (\code{reducedDims}, \code{colPairs} and
#' \code{rowPairs})
#' @param e which \code{experiment} to get or set the data. If NA, will choose
#' the global fields. If NULL, will choose the \code{\link{defaultExp}}.
#'
#' @seealso
#' \itemize{
#' \item \code{\link[SingleCellExperiment]{reducedDims}}
#' \item \code{\link[SingleCellExperiment]{colPairs}}
#' }
#'
#' @name SingleCellMultiExperiment-methods
NULL

## reducedDims #################################################################

#' @param ... Arguments passed to other methods
#' @param value An object of a class specified in the S4 method signature.
#'
#' @importFrom SingleCellExperiment reducedDims
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "reducedDims",
  signature = "SingleCellMultiExperiment",
  definition = function(x, withDimnames = TRUE, e = NULL, ...) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      return(reducedDims(x = int_SCE(x = x), withDimnames = withDimnames))
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    return(reducedDims(
      x = experiment(x = x, e = e, ...),
      withDimnames = withDimnames
    ))
  }
)

#' @importFrom SingleCellExperiment reducedDims<-
#' @importFrom MultiAssayExperiment experiments<- mapToList
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "reducedDims<-",
  signature = "SingleCellMultiExperiment",
  definition = function(x, withDimnames = TRUE, e = NULL, ..., value) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      reducedDims(x = int_SCE(x = x), withDimnames = withDimnames) <- value
      return(x)
    }
    if (withDimnames) {
      sampMap <- mapToList(dfmap = sampleMap(x = x))[[e]]
      for (i in seq_along(value)) {
        value[[i]] <- .primary2name(
          x = value[[i]],
          primary = sampMap[['primary']],
          name = sampMap[['colname']],
          dim.use = 1
        )
      }
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    reducedDims(x = experiments(x)[[e]], withDimnames = withDimnames) <- value
    return(x)
  }
)

#' @importFrom SingleCellExperiment reducedDim
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "reducedDim",
  signature = c("SingleCellMultiExperiment", "Character_OR_Numeric"),
  definition = function(x, type, withDimnames = TRUE, e = NULL, ...) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      return(reducedDim(
        x = int_SCE(x = x),
        type = type,
        withDimnames = withDimnames
      ))
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    return(reducedDim(
      x = experiment(x = x, e = e, ...),
      type = type,
      withDimnames = withDimnames
    ))
  }
)

#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom MultiAssayExperiment experiments<- mapToList
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "reducedDim<-",
  signature = c("SingleCellMultiExperiment", "Character_OR_Numeric"),
  definition = function(x, type, withDimnames = TRUE, e = NULL, ..., value) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      reducedDim(
        x = int_SCE(x = x),
        type = type,
        withDimnames = withDimnames
      ) <- value
      return(x)
    }
    if (withDimnames) {
      sampMap <- mapToList(dfmap = sampleMap(x = x))[[e]]
      value <- .primary2name(
        x = value,
        primary = sampMap[['primary']],
        name = sampMap[['colname']],
        dim.use = 1
      )
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    reducedDim(
      x = experiments(x)[[e]],
      type = type,
      withDimnames = withDimnames
    ) <- value
    return(x)
  }
)

#' @importFrom SingleCellExperiment reducedDim
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "reducedDim",
  signature = c("SingleCellMultiExperiment", "missing"),
  definition = function(x, type, withDimnames = TRUE, e = NULL, ...) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      return(reducedDim(x = int_SCE(x = x), withDimnames = withDimnames))
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    return(reducedDim(
      x = experiment(x = x, e = e, ...),
      withDimnames = withDimnames
    ))
  }
)

#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom MultiAssayExperiment experiments<- mapToList
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "reducedDim<-",
  signature = c("SingleCellMultiExperiment", "missing"),
  definition = function(x, type, withDimnames = TRUE, e = NULL, ..., value) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      reducedDim(x = int_SCE(x = x), withDimnames = withDimnames) <- value
      return(x)
    }
    if (withDimnames) {
      sampMap <- mapToList(dfmap = sampleMap(x = x))[[e]]
      value <- .primary2name(
        x = value,
        primary = sampMap[['primary']],
        name = sampMap[['colname']],
        dim.use = 1
      )
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    reducedDim(x = experiments(x)[[e]], withDimnames = withDimnames) <- value
    return(x)
  }
)

#' @importFrom SingleCellExperiment reducedDimNames
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "reducedDimNames",
  signature = "SingleCellMultiExperiment",
  definition = function(x) {
    return(reducedDimNames(x = int_SCE(x = x)))
  }
)

#' @importFrom SingleCellExperiment reducedDimNames<-
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "reducedDimNames<-",
  signature = c("SingleCellMultiExperiment", "character"),
  definition = function(x, value) {
    reducedDimNames(x = int_SCE(x = x)) <- value
    return(x)
  }
)

## colPairs ####################################################################

#' @importFrom SingleCellExperiment colPairs
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "colPairs",
  signature = "SingleCellMultiExperiment",
  definition = function(x, asSparse = FALSE, e = NULL, ...) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      return(colPairs(x = int_SCE(x = x), asSparse = asSparse))
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    return(colPairs(x = experiment(x = x, e = e, ...), asSparse = asSparse))
  }
)

#' @importFrom SingleCellExperiment colPairs<-
#' @importFrom MultiAssayExperiment experiments<-
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "colPairs<-",
  signature = "SingleCellMultiExperiment",
  definition = function(x, e = NULL, ..., value) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      colPairs(x = int_SCE(x = x)) <- value
      return(x)
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    colPairs(x = experiments(x)[[e]]) <- value
    return(x)
  }
)

#' @importFrom SingleCellExperiment colPair
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "colPair",
  signature = c("SingleCellMultiExperiment", "Character_OR_Numeric"),
  definition = function(x, type, asSparse = FALSE, e = NULL, ...) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      return(colPair(x = int_SCE(x = x), type = type, asSparse = asSparse))
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    return(colPair(
      x = experiment(x = x, e = e, ...),
      type = type,
      asSparse = asSparse
    ))
  }
)

#' @importFrom SingleCellExperiment colPair<-
#' @importFrom MultiAssayExperiment experiments<-
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "colPair<-",
  signature = c("SingleCellMultiExperiment", "Character_OR_Numeric"),
  definition = function(x, type, e = NULL, ..., value) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      colPair(x = int_SCE(x = x), type = type) <- value
      return(x)
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    colPair(x = experiments(x)[[e]], type = type) <- value
    return(x)
  }
)

#' @importFrom SingleCellExperiment colPair
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "colPair",
  signature = c("SingleCellMultiExperiment", "missing"),
  definition = function(x, type, asSparse = FALSE, e = NULL, ...) {
    e <- e %||% defaultExp(x = x)
    if (is.na(x = e)) {
      return(colPair(x = int_SCE(x = x), asSparse = asSparse))
    }
    e <- .get_valid_exp_name(exps = x, exp.name = e)
    return(colPair(x = experiment(x = x, e = e, ...), asSparse = asSparse))
  }
)

#' @importFrom SingleCellExperiment colPairNames
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "colPairNames",
  signature = "SingleCellMultiExperiment",
  definition = function(x, ...) {
    return(colPairNames(x = int_SCE(x = x), ...))
  }
)

#' @importFrom SingleCellExperiment colPairNames<-
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "colPairNames<-",
  signature = c("SingleCellMultiExperiment", "character"),
  definition = function(x, ..., value) {
    colPairNames(x = int_SCE(x = x), ...) <- value
    return(x)
  }
)

## SingleCellMultiExperiment-internal ##########################################

#' Internal SingleCellExperiment object for SingleCellMultiExperiment class
#'
#' Methods to get or set the internal \code{SingleCellExperiment} slot in a
#' \code{\link{SingleCellMultiExperiment}}. These functions are intended for
#' package development, and should not be used by ordinary users.
#'
#' @param x A \code{\link{SingleCellMultiExperiment}} object
#' @param ... Arguments passed to other methods.
#'
#' @details
#' The internal \code{\link[SingleCellExperiment]{SingleCellExperiment}} is only
#' used for storage of \code{\link[SingleCellExperiment]{int_metadata}} and
#' \code{\link[SingleCellExperiment]{int_colData}}.
#'
#' @name SingleCellMultiExperiment-internal
NULL

#' @export
#' @rdname SingleCellMultiExperiment-internal
setMethod(
  f = "int_SCE",
  signature = "SingleCellMultiExperiment",
  definition = function(x) {
    return(getElement(object = x, name = "int_SCE"))
  }
)

#' @param value An object of a class specified in the S4 method signature.
#'
#' @importFrom BiocBaseUtils setSlots
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname SingleCellMultiExperiment-internal
setMethod(
  f = "int_SCE<-",
  signature = c("SingleCellMultiExperiment", "SingleCellExperiment"),
  definition = function(x, ..., value) {
    altExps(x = value) <- NULL
    assays(x = value, withDimnames = FALSE) <- list()
    colData(x = value) <- NULL
    return(setSlots(object = x, int_SCE = value))
  }
)

#' @importFrom SingleCellExperiment int_metadata
#' @export
#' @rdname SingleCellMultiExperiment-internal
setMethod(
  f = "int_metadata",
  signature = "SingleCellMultiExperiment",
  definition = function(x) {
    return(int_metadata(x = int_SCE(x = x)))
  }
)

#' @importFrom SingleCellExperiment int_metadata<-
#' @export
#' @rdname SingleCellMultiExperiment-internal
setMethod(
  f = "int_metadata<-",
  signature = "SingleCellMultiExperiment",
  definition = function(x, value) {
    int_metadata(x = int_SCE(x = x)) <- value
    return(x)
  }
)

## Backend path ################################################################

#' @param object A \code{\link{SingleCellMultiExperiment}} object
#'
#' @importFrom BiocGenerics path
#' @export
#' @rdname SingleCellMultiExperiment-methods
setMethod(
  f = "path",
  signature = "SingleCellMultiExperiment",
  definition = function(object, ...) {
    return(int_metadata(x = object)[[.path_key]])
  }
)

## defaultExp ##################################################################

#' Default experiment in SingleCellMultiExperiment
#'
#' Get or set the default experiment. This is used to control and switch the
#' experiment to be manipulated.
#'
#' @param x A \code{\link{SingleCellMultiExperiment}}
#' @param ... Arguments passed to other methods.
#' @param return.exp If \code{TRUE}, will return the default experiment object,
#' else will only return the strings specifying the default experiment. If
#' \code{defaultExp} is \code{NA}, this argument will be ignored.
#' @param value Update the default experiment strings. Use \code{NA} to
#' manipulate the \code{\link[SingleCellExperiment]{reducedDims}} and
#' \code{\link[SingleCellExperiment]{colPairs}} stored in \code{\link{int_SCE}}.
#'
#' @name defaultExp
NULL

#' @export
#' @rdname defaultExp
setMethod(
  f = "defaultExp",
  signature = "SingleCellMultiExperiment",
  definition = function(x, return.exp = FALSE, ...) {
    out <- getElement(object = x, name = "defaultExp")
    if (!is.na(x = out) && return.exp) {
      return(experiment(x = x, e = out, ...))
    }
    return(out)
  }
)

#' @importFrom BiocBaseUtils setSlots
#' @export
#' @rdname defaultExp
setMethod(
  f = "defaultExp<-",
  signature = c("SingleCellMultiExperiment", "character"),
  definition = function(x, ..., value) {
    value <- .get_valid_exp_name(exps = x, exp.name = value)
    return(setSlots(object = x, defaultExp = value, check = FALSE))
  }
)

## Subset ######################################################################

#' Subsetting a SingleCellMultiExperiment object
#'
#' Extracting and dividing a \code{SingleCellMultiExperiment}. These methods
#' actually inherit from \code{\link[MultiAssayExperiment]{subsetBy}} in
#' \pkg{MultiAssayExperiment} package.
#'
#' @param x A \code{\link{SingleCellMultiExperiment}} object.
#' @param ... Additional arguments passed on to lower level functions.
#' @param y Any argument used for subsetting, can be a \code{character},
#' \code{logical}, \code{integer}, \code{list} or \code{List} vector
#' @param i Either a \code{character}, \code{integer}, \code{logical} or
#' \code{GRanges} object for subsetting by rows
#' @param j Either a \code{character}, \code{logical}, or \code{numeric} vector
#' for subsetting by \code{colData} rows. See details for more information.
#' @param k Either a \code{character}, \code{logical}, or \code{numeric} vector
#' for subsetting by assays
#' @param drop logical (default FALSE) whether to drop all empty assay elements
#' in the \code{ExperimentList}
#'
#' @seealso
#' \code{\link[MultiAssayExperiment]{subsetBy}}
#'
#' @name SingleCellMultiExperiment-subset
NULL

#' @export
#' @rdname SingleCellMultiExperiment-subset
setMethod(
  f = "[",
  signature = c("SingleCellMultiExperiment", "ANY", "ANY", "ANY"),
  definition = function(x, i, j, k, ..., drop = FALSE) {
    sce <- int_SCE(x = x)
    default.exp <- defaultExp(x = x)
    x <- callNextMethod()
    default.exp <- .drop_defaultExp(exps = x, default.exp = default.exp)
    return(.mae_to_scme(mae = x, int_SCE = sce, defaultExp = default.exp))
  }
)

#' @param value An object of a class specified in the S4 method signature.
#'
#' @importFrom S4Vectors DataFrame
#' @export
#' @rdname SingleCellMultiExperiment-subset
setMethod(
  f = "[[<-",
  signature = "SingleCellMultiExperiment",
  definition = function(x, i, j, ..., value) {
    default.exp <- defaultExp(x = x)
    x <- callNextMethod()
    warning(
      "Directly modifying experiments in ", class(x = x), " by '[[<-' ",
      "will clean the internal SingleCellExperiment.",
      immediate. = TRUE, call. = FALSE
    )
    cdata <- DataFrame(row.names = rownames(x = colData(x = x)))
    sce <- SingleCellExperiment(colData = cdata)
    default.exp <- .drop_defaultExp(exps = x, default.exp = default.exp)
    return(.mae_to_scme(mae = x, int_SCE = sce, defaultExp = default.exp))
  }
)

#' @importFrom S4Vectors DataFrame
#' @export
#' @rdname SingleCellMultiExperiment-subset
setMethod(
  f = "[<-",
  signature = "SingleCellMultiExperiment",
  definition = function(x, i, j, ..., value) {
    default.exp <- defaultExp(x = x)
    x <- callNextMethod()
    warning(
      "Directly modifying ", class(x = x), " by '[<-' ",
      "will clean the internal SingleCellExperiment.",
      immediate. = TRUE, call. = FALSE
    )
    cdata <- DataFrame(row.names = rownames(x = colData(x = x)))
    sce <- SingleCellExperiment(colData = cdata)
    default.exp <- .drop_defaultExp(exps = x, default.exp = default.exp)
    return(.mae_to_scme(mae = x, int_SCE = sce, defaultExp = default.exp))
  }
)

#' @importFrom MultiAssayExperiment subsetByColData
#' @export
#' @rdname SingleCellMultiExperiment-subset
setMethod(
  f = "subsetByColData",
  signature = c("SingleCellMultiExperiment", "ANY"),
  definition = function(x, y) {
    sce <- int_SCE(x = x)
    default.exp <- defaultExp(x = x)
    x <- callNextMethod()
    default.exp <- .drop_defaultExp(exps = x, default.exp = default.exp)
    return(.mae_to_scme(mae = x, int_SCE = sce, defaultExp = default.exp))
  }
)

#' @importFrom MultiAssayExperiment subsetByColData
#' @export
#' @rdname SingleCellMultiExperiment-subset
setMethod(
  f = "subsetByColData",
  signature = c("SingleCellMultiExperiment", "character"),
  definition = function(x, y) {
    sce <- int_SCE(x = x)
    default.exp <- defaultExp(x = x)
    x <- callNextMethod()
    default.exp <- .drop_defaultExp(exps = x, default.exp = default.exp)
    return(.mae_to_scme(mae = x, int_SCE = sce, defaultExp = default.exp))
  }
)

#' @importFrom MultiAssayExperiment subsetByColumn
#' @export
#' @rdname SingleCellMultiExperiment-subset
setMethod(
  f = "subsetByColumn",
  signature = c("MultiAssayExperiment", "ANY"),
  definition = function(x, y) {
    sce <- int_SCE(x = x)
    default.exp <- defaultExp(x = x)
    x <- callNextMethod()
    default.exp <- .drop_defaultExp(exps = x, default.exp = default.exp)
    return(.mae_to_scme(mae = x, int_SCE = sce, defaultExp = default.exp))
  }
)

#' @importFrom MultiAssayExperiment subsetByAssay
#' @export
#' @rdname SingleCellMultiExperiment-subset
setMethod(
  f = "subsetByAssay",
  signature = c("SingleCellMultiExperiment", "ANY"),
  definition = function(x, y) {
    sce <- int_SCE(x = x)
    default.exp <- defaultExp(x = x)
    x <- callNextMethod()
    default.exp <- .drop_defaultExp(exps = x, default.exp = default.exp)
    return(.mae_to_scme(mae = x, int_SCE = sce, defaultExp = default.exp))
  }
)

## Show ########################################################################

#' @param object A \code{SingleCellMultiExperiment}
#'
#' @importFrom SingleCellExperiment colPairNames reducedDimNames
#' @importFrom S4Vectors coolcat
#' @importFrom methods show
#' @export
#' @rdname SingleCellMultiExperiment
setMethod(
  f = "show",
  signature = "SingleCellMultiExperiment",
  definition = function(object) {
    callNextMethod()
    default.exp <- defaultExp(x = object)
    cat("Default experiment: ", default.exp, "\n")
    coolcat("reducedDimNames(%d): %s\n", reducedDimNames(x = object))
    coolcat("colPairNames(%d): %s\n", colPairNames(x = object))
    .show_h5backed_path(object = object)
    return(invisible(x = NULL))
  }
)

## Validity ####################################################################

#' @importFrom S4Vectors isEmpty
#' @importFrom SummarizedExperiment rowRanges
.valid_scme <- function(object) {
  msg <- NULL
  default.exp <- defaultExp(x = object)
  if (!is.na(x = default.exp) && !default.exp %in% names(x = object)) {
    msg <- msg %>%
      c(paste0("'", default.exp, "' was not found in experiment names"))
  }
  int.sce <- int_SCE(x = object)
  cells <- rownames(x = colData(x = object))
  if (!identical(x = colnames(x = int.sce), y = cells)) {
    msg <- msg %>%
      c("'colnames(int_SCE(x))' must be identical to 'rownames(colData(x))'")
  }
  if (length(x = msg) > 0) {
    return(msg)
  }
  return(TRUE)
}

#' @importFrom S4Vectors setValidity2
setValidity2(Class = "SingleCellMultiExperiment", method = .valid_scme)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SummarizedExperiment colData
.mae_to_scme <- function(mae, int_SCE, defaultExp) {
  defaultExp <- defaultExp %||% names(x = mae)[1]
  cells <- rownames(x = colData(x = mae))
  int_SCE <- int_SCE[, cells, drop = FALSE]
  if (is.na(x = defaultExp)) {
    defaultExp <- NA_character_
  }
  return(new(
    Class = "SingleCellMultiExperiment",
    mae,
    int_SCE = int_SCE,
    defaultExp = defaultExp
  ))
}

#' @importFrom easy.utils fastIntersect
.primary2name <- function(x, primary, name, dim.use) {
  prim2name <- setNames(object = name, nm = primary)
  raw.dimn <- dimnames(x = x)[[dim.use]]
  any.primary <- fastIntersect(x = raw.dimn, y = primary)
  if (identical(x = any.primary, y = raw.dimn)) {
    dimnames(x = x)[[dim.use]] <- prim2name[raw.dimn]
  }
  return(x)
}

#' @importFrom easy.utils isValidCharacters
.get_valid_exp_name <- function(exps, exp.name) {
  exp.name <- exp.name[1]
  if (is.null(x = exp.name)) {
    return(names(x = exps)[1])
  }
  if (is.na(x = exp.name)) {
    return(NA_character_)
  }
  if (!any(isValidCharacters(x = exp.name))) {
    stop("Invalid experiment name was input.")
  }
  if (!any(exp.name %in% names(x = exps))) {
    stop("'", exp.name, "' was not found in experiment names")
  }
  return(exp.name)
}

.drop_defaultExp <- function(exps, default.exp) {
  if (default.exp %in% names(x = exps)) {
    return(default.exp)
  }
  default.exp2 <- names(x = exps)[1]
  warning(
    "The orignal defaultExp '", default.exp, "' has been removed.\n",
    "Set the 1st experiment in the output object to defaultExp: ",
    default.exp2,
    immediate. = TRUE, call. = FALSE
  )
  return(default.exp2)
}
