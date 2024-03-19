#' @importFrom methods coerce setAs as
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods ######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## SelfHits and sparse matrix ##################################################

#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors from nLnode nRnode SelfHits to values values<-
#' @importClassesFrom Matrix dgCMatrix dgTMatrix
#' @importClassesFrom S4Vectors SelfHits
NULL

setAs(
  from = "SelfHits",
  to = "dgCMatrix",
  def = function(from) {
    return(sparseMatrix(
      i = from(x = from) - 1L,
      j = to(x = from) - 1L,
      x = values(x = from)[, 1],
      dims = c(nLnode(x = from), nRnode(x = from)),
      index1 = FALSE
    ))
  }
)

setAs(
  from = "SelfHits",
  to = "dgTMatrix",
  def = function(from) {
    return(sparseMatrix(
      i = from(x = from) - 1L,
      j = to(x = from) - 1L,
      x = values(x = from)[, 1],
      dims = c(nLnode(x = from), nRnode(x = from)),
      repr = "T",
      index1 = FALSE
    ))
  }
)

setAs(
  from = "dgTMatrix",
  to = "SelfHits",
  def = function(from) {
    if (!identical(x = nrow(x = from), y = ncol(x = from))) {
      stop(
        "The input matrix must have identical row numbers and column numbers ",
        "to be coerced to 'SelfHits'."
      )
    }
    hits <- SelfHits(
      from = slot(object = from, name = "i") + 1L,
      to = slot(object = from, name = "j") + 1L,
      nnode = unique(x = dim(x = from))
    )
    values(x = hits) <- slot(object = from, name = "x")
    return(hits)
  }
)

setAs(
  from = "dgCMatrix",
  to = "SelfHits",
  def = function(from) {
    from %>%
      as(Class = "dgTMatrix") %>%
      as(Class = "SelfHits")
  }
)
