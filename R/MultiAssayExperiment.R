
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions ####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Extracting an assay with primaries as column names
#'
#' Change the column names of extracted experiment with primaries stored in the
#' \code{\link[SummarizedExperiment]{colData}}.
#'
#' @param x A \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} object.
#' @param i Which experiment to be fetched. Can be an integer or character.
#'
#' @return The extracted experiment object.
#'
#' @importFrom MultiAssayExperiment experiments mapToList sampleMap
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
getWithPrimaries <- function(x, i) {
  if (!inherits(x = x, what = "MultiAssayExperiment")) {
    stop("Provide a MultiAssayExperiment as input")
  }
  stopifnot(
    is.numeric(x = i) || is.character(x = i),
    identical(x = length(x = i), y = 1L),
    !is.na(x = i),
    !is.logical(x = i)
  )
  exp <- experiments(x = x)[[i]]
  sampMap <- mapToList(dfmap = sampleMap(x = x))[[i]]
  nameMap <- setNames(object = sampMap[['primary']], nm = sampMap[['colname']])
  colnames(x = exp) <- nameMap[colnames(x = exp)]
  return(exp)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Accessing single experiment data in MultiAssayExperiment
#'
#' Supplemental methods for accessing single experiment from
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} objects.
#'
#' @param x A \code{MultiAssayExperiment} object
#' @param e Which experiment to fetch. Can be a integer index or an experiment
#' name. If missing, will fetch the first experiment by default.
#' @param ... Arguments passed to other methods.
#'
#' @return
#' The extracted experiment object
#'
#' @seealso \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#'
#' @name experiment
NULL

#' @param withPrimaries Whether or not to use primary names to replace the
#' original column names
#' @param withColData Whether or not to also fetch the associated
#' \code{\link[SummarizedExperiment]{colData}} of \code{x}
#' @param mode String indicating how \code{MultiAssayExperiment} metadata should
#' be added or just replaced to the extracted experiment. Passed to
#' \code{\link[MultiAssayExperiment]{getWithColData}}.
#'
#' @importFrom MultiAssayExperiment experiments getWithColData mapToList
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @rdname experiment
setMethod(
  f = "experiment",
  signature = c("MultiAssayExperiment", "character"),
  definition = function(
    x,
    e,
    withPrimaries = FALSE,
    withColData = FALSE,
    mode = c("append", "replace"),
    ...
  ) {
    mode <- match.arg(arg = mode)
    if (!withColData && !withPrimaries) {
      return(experiments(x = x)[[e]])
    }
    if (!withColData && withPrimaries) {
      return(getWithPrimaries(x = x, i = e))
    }
    # MultiAssayExperiment::getWithColData will automatically set colnames to
    # primaries.
    return(tryCatch(
      expr = getWithColData(x = x, i = e, mode = mode),
      error = function(err) {
        message(
          "MultiAssayExperiment::getWithColData:\n ", err, "\n",
          "Force to set 'withColData = FALSE'"
        )
        if (withPrimaries) {
          return(getWithPrimaries(x = x, i = e))
        }
        return(experiments(x = x)[[e]])
      }
    ))
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @rdname experiment
setMethod(
  f = "experiment",
  signature = c("MultiAssayExperiment", "numeric"),
  definition = function(
    x,
    e,
    withPrimaries = FALSE,
    withColData = FALSE,
    mode = c("append", "replace"),
    ...
  ) {
    mode <- match.arg(arg = mode)
    exp.name <- names(x = x)[[e]]
    return(experiment(
      x = x,
      e = exp.name,
      withPrimaries = withPrimaries,
      withColData = withColData,
      mode = mode,
      ...
    ))
  }
)

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
#' @rdname experiment
setMethod(
  f = "experiment",
  signature = c("MultiAssayExperiment", "missing"),
  definition = function(
    x,
    e,
    withPrimaries = FALSE,
    withColData = FALSE,
    mode = c("append", "replace"),
    ...
  ) {
    mode <- match.arg(arg = mode)
    exp.name <- names(x = x)[[1]]
    return(experiment(
      x = x,
      e = exp.name,
      withPrimaries = withPrimaries,
      withColData = withColData,
      mode = mode,
      ...
    ))
  }
)
