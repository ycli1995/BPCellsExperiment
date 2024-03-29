
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Search common names for assays, reductions or altExps in SCEs
#'
#' @param SCEs A list of \code{SingleCellExperiment}
#' @param func Function to fetch the query names.
#' @param names Names to be searched. If \code{NULL}, search all common names.
#' If \code{FALSE}, no search will be done.
#' @param verbose Print progress.
#'
#' @noRd
.search_common_names <- function(SCEs, func, names = NULL, verbose = TRUE) {
  func.name <- as.character(x = substitute(expr = func))
  if (is_false(x = names)) {
    verboseMsg("Skip searching common ", func.name)
    return(NULL)
  }
  verboseMsg("Search common ", func.name, " ", paste(names, collapse = ", "))
  common.names <- SCEs %>%
    lapply(FUN = func) %>%
    Reduce(f = intersect)
  if (is.null(x = names)) {
    names <- common.names
  } else {
    names <- intersect(x = names, y = common.names)
  }
  if (length(x = names) == 0) {
    warning(
      "No common ", func.name, " was found.",
      call. = FALSE, immediate. = TRUE
    )
    return(NULL)
  }
  verboseMsg("Found common ", func.name, ": ", paste(names, collapse = ", "))
  return(names)
}

#' De-duplicate column names for input SCEs
#'
#' @param SCEs A list of \code{SingleCellExperiment}.
#' @param collapse Collapse character to link column names and the suffixes.
#' @param verbose Print progress.
#' @param stop Whether or not to just raise an error for duplicated column
#' names.
#'
#' @noRd
.check_duplicated_colnames <- function(SCEs, collapse = "_", stop = FALSE) {
  all.dimns <- SCEs %>%
    lapply(FUN = colnames) %>%
    unlist()
  if (!any(duplicated(x = all.dimns))) {
    return(SCEs)
  }
  if (stop) {
    stop(
      "Duplicated column names present across ",
      class(x = SCEs[[1]]), " objects provided."
    )
  }
  warning(
    "Some column names are duplicated across ", class(x = SCEs[[1]]),
    " objects provided. Adding suffixes to enforce unique.",
    call. = FALSE, immediate. = TRUE
  )
  for (i in seq_along(SCEs)) {
    colnames(x = SCEs[[i]]) <- paste0(colnames(x = SCEs[[i]]), collapse, i)
  }
  return(SCEs)
}

#' @importFrom S4Vectors merge
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment altExp altExpNames altExps<- reducedDims<-
.merge_sce_altExps <- function(
    SCEs,
    altExps = NULL,
    altExps.assays = NULL,
    verbose = TRUE
) {
  altExps <- .search_common_names(
    SCEs = SCEs,
    func = altExpNames,
    names = altExps,
    verbose = verbose
  )
  new.altExps <- list()
  for (i in altExps) {
    verboseMsg("Merging altExp: ", i)
    old.altExps <- lapply(X = SCEs, FUN = altExp, e = i, withDimnames = TRUE)
    for (j in seq_len(length.out = length(x = old.altExps))) {
      colData(x = old.altExps[[j]]) <- NULL
      altExps(x = old.altExps[[j]]) <- NULL
      reducedDims(x = old.altExps[[j]]) <- NULL
    }
    new.altExps[[i]] <- merge(
      x = old.altExps[[1]],
      y = old.altExps[2:length(x = old.altExps)],
      assays = altExps.assays,
      reductions = FALSE,
      altExps = FALSE,
      label = NULL,
      add.idx = NULL,
      verbose = verbose
    )
  }
  return(new.altExps)
}

#' @importFrom SingleCellExperiment reducedDim reducedDimNames
.merge_sce_reducedDims <- function(SCEs, reducedDims = NULL, verbose = TRUE) {
  reducs <- .search_common_names(
    SCEs = SCEs,
    func = reducedDimNames,
    names = reducedDims,
    verbose = verbose
  )
  new.reducs <- list()
  for (i in reducs) {
    verboseMsg("Merging reducedDim: ", i)
    old.reducs <- lapply(
      X = SCEs,
      FUN = function(x) x %>%
        reducedDim(type = i, withDimnames = TRUE) %>%
        as.matrix()
    )
    min.ncol <- vapply(
      X = old.reducs,
      FUN = ncol,
      FUN.VALUE = integer(length = 1L),
      USE.NAMES = FALSE
    )
    if (length(x = unique(x = min.ncol)) > 1) {
      min.ncol <- min(min.ncol)
      warning(
        "reducedDim '", i, "' contain differing numbers of dimensions, ",
        "merging first ", min.ncol,
        call. = FALSE, immediate. = TRUE
      )
      old.reducs <- lapply(X = old.reducs, FUN = function(x) x[, 1:min.ncol])
    }
    new.reducs[[i]] <- Reduce(f = rbind, x = old.reducs)
  }
  return(new.reducs)
}

#' @importFrom SummarizedExperiment assay assayNames
.merge_sce_assays <- function(SCEs, assays = NULL, verbose = TRUE) {
  assays <- .search_common_names(
    SCEs = SCEs,
    func = assayNames,
    names = assays,
    verbose = verbose
  )
  new.assays <- list()
  for (i in assays) {
    verboseMsg("Merging assay: ", i)
    old.assays <- lapply(X = SCEs, FUN = assay, i = i, withDimnames = TRUE)
    new.assays[[i]] <- mergeMatrices(
      x = old.assays[[1]],
      y = old.assays[2:length(x = old.assays)]
    )
  }
  return(new.assays)
}

#' @importFrom SummarizedExperiment colData
.prep_merge_SCEs <- function(
    SCEs,
    label = NULL,
    add.idx = NULL,
    idx.collapse = "_",
    verbose = TRUE
) {
  verboseMsg("Preparing SCEs to be merged:")
  if (!is.null(x = add.idx)) {
    stopifnot(length(x = add.idx) == length(x = SCEs))
    verboseMsg("  Adding idx prefixes: ", paste(add.idx, collapse = ", "))
    for (i in seq_along(SCEs)) {
      colnames(x = SCEs[[i]]) <- paste0(
        add.idx[i],
        idx.collapse,
        colnames(x = SCEs[[i]])
      )
    }
  }
  if (!is.null(x = label)) {
    label <- label[1]
    if (is.null(x = names(x = SCEs))) {
      names(x = SCEs) <- as.character(x = seq_along(SCEs))
    }
    verboseMsg(
      "  Adding column '", label, "' to specify batch info: ",
      paste(names(x = SCEs), collapse = ", ")
    )
    for (i in names(x = SCEs)) {
      colData(x = SCEs[[i]])[, label] <- i
    }
  }
  SCEs <- .check_duplicated_colnames(SCEs = SCEs)
  return(SCEs)
}

#' @importFrom SummarizedExperiment colData rowData
#' @importFrom SingleCellExperiment SingleCellExperiment
.merge_SCEs <- function(
    SCEs,
    assays = NULL,
    reducedDims = NULL,
    altExps = NULL,
    altExps.assays = NULL,
    label = NULL,
    add.idx = NULL,
    idx.collapse = "_",
    verbose = TRUE,
    ...
) {
  SCEs <- .prep_merge_SCEs(
    SCEs = SCEs,
    label = label,
    add.idx = add.idx,
    idx.collapse = idx.collapse,
    verbose = verbose
  )
  # colData and rowData
  verboseMsg("Merging rowData...")
  rdata <- SCEs %>%
    lapply(FUN = rowData) %>%
    .concat_DFs(r.join = "outer")

  verboseMsg("Merging colData...")
  cdata <- SCEs %>%
    lapply(FUN = colData) %>%
    .rbind_DFs(join = "outer")

  # merge altExps
  new.altExps <- .merge_sce_altExps(
    SCEs = SCEs,
    altExps = altExps,
    altExps.assays = altExps.assays,
    verbose = verbose
  )
  # merge assays
  new.assays <- .merge_sce_assays(
    SCEs = SCEs,
    assays = assays,
    verbose = verbose
  )
  # merge reductions
  new.reducs <- .merge_sce_reducedDims(
    SCEs = SCEs,
    reducedDims = reducedDims,
    verbose = verbose
  )
  new.SCE <- SingleCellExperiment(
    assays = new.assays,
    rowData = rdata,
    colData = cdata,
    reducedDims = new.reducs,
    altExps = new.altExps
  )
  return(new.SCE)
}
