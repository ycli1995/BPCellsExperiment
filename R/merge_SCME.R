
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SummarizedExperiment colData
.check_duplicated_colnames_SCMEs <- function(
    SCMEs,
    collapse = "_",
    stop = FALSE
) {
  all.dimns <- SCMEs %>%
    lapply(FUN = function(x) rownames(x = colData(x = x))) %>%
    unlist()
  if (!any(duplicated(x = all.dimns))) {
    return(SCMEs)
  }
  if (stop) {
    stop("Duplicate colnames present across objects provided.")
  }
  warning(
    "Some colnames are duplicated across ", class(x = SCMEs[[1]]), " objects ",
    "provided. Renaming to enforce unique.",
    call. = FALSE, immediate. = TRUE
  )
  for (i in seq_along(SCMEs)) {
    SCMEs[[i]] <- setColnames(
      x = SCMEs[[i]],
      new.names = paste0(rownames(x = colData(x = SCMEs[[i]])), collapse, i)
    )
    SCMEs[[i]] <- setColnames(
      x = SCMEs[[i]],
      new.names = lapply(
        X = colnames(x = SCMEs[[i]]),
        FUN = function(x) paste0(x, collapse, i)
      )
    )
  }
  return(SCMEs)
}

#' @importFrom SummarizedExperiment colData
.prep_merge_SCMEs <- function(
    SCMEs,
    label = NULL,
    add.idx = NULL,
    idx.collapse = "_",
    verbose = TRUE
) {
  verboseMsg("Preparing SCMEs to be merged:")
  if (!is.null(x = add.idx)) {
    stopifnot(length(x = add.idx) == length(x = SCMEs))
    verboseMsg("  Adding idx prefixes: ", paste(add.idx, collapse = ", "))
    for (i in seq_along(SCMEs)) {
      SCMEs[[i]] <- setColnames(
        x  = SCMEs[[i]],
        new.names = paste0(
          add.idx[i],
          idx.collapse,
          rownames(x = colData(x = SCMEs[[i]]))
        )
      )
      SCMEs[[i]] <- setColnames(
        x = SCMEs[[i]],
        new.names = lapply(
          X = colnames(x = SCMEs[[i]]),
          FUN = function(x) paste0(add.idx, idx.collapse, x)
        )
      )
    }
  }
  if (!is.null(x = label)) {
    if (is.null(x = names(x = SCMEs))) {
      names(x = SCMEs) <- as.character(x = seq_along(SCMEs))
    }
    verboseMsg(
      "  Adding column '", label, "' to specify batch info: ",
      paste(names(x = SCMEs), collapse = ", ")
    )
    for (i in names(x = SCMEs)) {
      colData(x = SCMEs[[i]])[, label] <- i
    }
  }
  SCMEs <- .check_duplicated_colnames_SCMEs(SCMEs = SCMEs)
  return(SCMEs)
}

#' @importFrom S4Vectors merge
#' @importFrom SummarizedExperiment assayNames assays assays<- colData rowData
#' @importFrom SingleCellExperiment reducedDimNames reducedDims
#' @importFrom MultiAssayExperiment listToMap mapToList
.merge_SCMEs <- function(
    SCMEs,
    experiments = NULL,
    global.reducedDims = NULL,
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
  SCMEs <- .prep_merge_SCMEs(
    SCMEs = SCMEs,
    label = label,
    add.idx = add.idx,
    idx.collapse = idx.collapse,
    verbose = verbose
  )
  # colData
  verboseMsg("Merging colData...")
  cdata <- SCMEs %>%
    lapply(FUN = colData) %>%
    .rbind_DFs(join = "outer")

  # merge experiments
  all.exps <- SCMEs %>%
    lapply(FUN = names) %>%
    unlist() %>%
    table()
  e.use <- experiments %||% names(x = all.exps)
  e.use <- intersect(x = e.use, y = names(x = all.exps))
  if (length(x = e.use) == 0) {
    stop("'experiments' not found: ", paste(experiments, collapse = "', '"))
  }
  new.exps <- list()
  new.smap <- list()
  for (e in e.use) {
    verboseMsg("Merging experiment: ", e)
    tmp.list <- lapply(X = SCMEs, FUN = experiment, e = e)
    if (length(x = tmp.list) == 1) {
      new.exps[[e]] <- tmp.list[[1]]
      if (any(assays %in% assayNames(x = new.exps[[e]]))) {
        a.use <- intersect(x = assays, y = assayNames(x = new.exps[[e]]))
        assays(x = new.exps[[e]]) <- assays(x = new.exps[[e]])[a.use]
      }
      if (any(reducedDims %in% reducedDimNames(x = new.exps[[e]]))) {
        r.use <- intersect(
          x = reducedDims,
          y = reducedDimNames(x = new.exps[[e]])
        )
        reducedDims(x = new.exps[[e]]) <- reducedDims(x = new.exps[[e]])[r.use]
      }
      if (any(altExps %in% altExpNames(x = new.exps[[e]]))) {
        alt.use <- intersect(x = altExps, y = altExpNames(x = new.exps[[e]]))
        altExps(x = new.exps[[e]]) <- altExps(x = new.exps[[e]])[altExps]
      }
      next
    }
    new.exps[[e]] <- merge(
      x = tmp.list[[1]],
      y = tmp.list[2:length(x = tmp.list)],
      assays = assays,
      reducedDims = reducedDims,
      altExps = altExps,
      altExps.assays = altExps.assays,
      label = NULL,
      add.idx = NULL,
      idx.collapse = "_",
      verbose = verbose
    )
    new.smap[[e]] <- SCMEs %>%
      lapply(FUN = function(x) mapToList(dfmap = sampleMap(x = x))[[e]]) %>%
      Reduce(f = rbind)
  }
  new.smap <- listToMap(listmap = new.smap)
  cdata <- cdata[rownames(x = cdata) %in% new.smap[['primary']], ]
  new.sce <- .merge_SCEs(
    SCEs = lapply(X = SCMEs, FUN = int_SCE),
    reducedDims = global.reducedDims,
    verbose = verbose
  )[, rownames(x = cdata)]
  new.SCME <- SingleCellMultiExperiment(
    experiments = new.exps,
    colData = cdata,
    sampleMap = new.smap,
    reducedDims = reducedDims(x = new.sce)
  )
  return(new.SCME)
}
