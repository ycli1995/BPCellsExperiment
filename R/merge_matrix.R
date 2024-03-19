
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
#' @rdname mergeMatrices
#' @method mergeMatrices default
mergeMatrices.default <- function(x, y, ...) {
  stop("Merging unsupported matrices: ", class(x = x))
}

#' @examples
#' bpce1 <- load_example_sce()
#' bpce2 <- load_example_sce("unsorted")
#'
#' colnames(bpce1) <- paste0("sorted_", colnames(bpce1))
#' colnames(bpce2) <- paste0("unsorted_", colnames(bpce2))
#'
#' ## matrix
#' mat <- mergeMatrices(as.matrix(counts(bpce1)), as.matrix(counts(bpce2)))
#'
#' @importFrom rlang is_bare_list
#' @export
#' @rdname mergeMatrices
#' @method mergeMatrices matrix
mergeMatrices.matrix <- function(x, y, ...) {
  on.exit(expr = gc(verbose = FALSE))
  if (!is_bare_list(x = y)) {
    y <- list(y)
  }
  y <- c(list(x), y)
  all.colnames <- .check_merge_matrices_dimnames(
    mats = y,
    func = colnames
  )
  all.rownames <- .check_merge_matrices_dimnames(
    mats = y,
    func = rownames,
    allow.dup = TRUE
  )
  m <- matrix(
    data = 0,
    nrow = length(x = all.rownames),
    ncol = length(x = all.colnames),
    dimnames = list(all.rownames, all.colnames)
  )
  for (i in seq_along(y)) {
    m[rownames(x = y[[i]]), colnames(x = y[[i]])] <- as.matrix(x = y[[i]])
  }
  return(m)
}

#' @examples
#' ## dgCMatrix
#' ## The matrix 'y' will be coerced to dgCMatrix
#' mat <- mergeMatrices(as(counts(bpce1), "dgCMatrix"), counts(bpce2))
#'
#' @importFrom rlang is_bare_list
#' @importFrom SeuratObject RowMergeSparseMatrices
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @rdname mergeMatrices
#' @method mergeMatrices dgCMatrix
mergeMatrices.dgCMatrix <- function(x, y, ...) {
  on.exit(expr = gc(verbose = FALSE))
  if (!is_bare_list(x = y)) {
    y <- list(y)
  }
  y <- c(list(x), y)
  .check_merge_matrices_dimnames(
    mats = y,
    func = colnames,
    return.dnames = FALSE
  )
  .check_merge_matrices_dimnames(
    mats = y,
    func = rownames,
    allow.dup = TRUE,
    return.dnames = FALSE
  )
  for (i in seq_along(y)) {
    y[[i]] <- as(object = y[[i]], Class = "dgCMatrix")
  }
  return(RowMergeSparseMatrices(mat1 = y[[1]], mat2 = y[2:length(x = y)]))
}

#' @examples
#' ## IterableMatrix
#' ## The matrix 'y' will be coerced to IterableMatrix
#' mat <- mergeMatrices(counts(bpce1), as(counts(bpce2), "dgCMatrix"))
#'
#' @importFrom rlang is_bare_list
#' @export
#' @rdname mergeMatrices
#' @method mergeMatrices IterableMatrix
mergeMatrices.IterableMatrix <- function(x, y, ...) {
  on.exit(expr = gc(verbose = FALSE))
  if (!is_bare_list(x = y)) {
    y <- list(y)
  }
  y <- c(list(x), y)
  .check_merge_matrices_dimnames(
    mats = y,
    func = colnames,
    return.dnames = FALSE
  )
  all.rownames <- .check_merge_matrices_dimnames(
    mats = y,
    func = rownames,
    allow.dup = TRUE
  )
  for (i in seq_along(y)) {
    y[[i]] <- .prep_bpcells_mat_merge(mat = y[[i]], all.rownames = all.rownames)
  }
  m <- Reduce(f = cbind, x = y)
  return(m)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check_merge_matrices_dimnames <- function(
    mats,
    func,
    allow.dup = FALSE,
    return.dnames = TRUE
) {
  func.name <- func %>%
    substitute() %>%
    as.character()
  all.dnames <- lapply(X = mats, FUN = func)
  n.dnames <- lengths(x = all.dnames)
  if (any(n.dnames == 0)) {
    stop(func.name, "(input[[", which(x = n.dnames == 0)[1], "]]) is empty")
  }
  all.dnames <- unlist(x = all.dnames)
  if (anyDuplicated(x = all.dnames) && !allow.dup) {
    stop(func.name, " across all matrices must be unique")
  }
  if (return.dnames) {
    return(unique(x = all.dnames))
  }
  return(invisible(x = NULL))
}

#' @importFrom MatrixExtra emptySparse
.prep_bpcells_mat_merge <- function(mat, all.rownames) {
  if (!inherits(x = mat, what = "IterableMatrix")) {
    mat <- as(object = mat, Class = "IterableMatrix")
  }
  missing.rows <- setdiff(x = all.rownames, y = rownames(x = mat))
  if (length(x = missing.rows) == 0) {
    return(mat)
  }
  zero_mat <- emptySparse(
    nrow = length(x = missing.rows),
    ncol = ncol(x = mat),
    format = "C"
  )
  dimnames(x = zero_mat) <- list(missing.rows, colnames(x = mat))
  zero_mat <- zero_mat %>%
    as(Class = "IterableMatrix") %>%
    convert_matrix_type(type = matrix_type(x = mat))
  mat <- rbind(mat, zero_mat)[all.rownames, ]
  return(mat)
}
