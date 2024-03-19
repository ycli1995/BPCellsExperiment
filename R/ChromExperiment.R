
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constructor ##################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The ChromExperiment class
#'
#' The \code{ChromExperiment} class is designed to represent single-cell
#' chromatin accessibility (such as ATAC-seq) data using \pkg{BPCells}. It
#' inherits from the \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' class and is used in the same manner. In addition, it supports storage of
#' barcoded and aligned fragments information via \code{\link{fragments}}, in
#' which the \code{\link[BPCells:IterableFragments-methods]{IterableFragments}}
#' objects are used as the underlying data.
#'
#' @param ... Arguments passed to the constructor function of
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} to fill the slots
#' of the base class.
#' @param sep Separators to use for strings encoding genomic coordinates. The
#' first element is used to separate the chromosome from the coordinates, the
#' second element is used to separate the start from end coordinate. Only used
#' if \code{ranges} is \code{NULL}
#' @param ranges A set of \code{\link[GenomicRanges]{GRanges}} corresponding to
#' the rows (peaks) of the input matrix
#' @param genome A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome used. Alternatively, the name of a UCSC genome
#' can be provided and the sequence information will be downloaded from UCSC.
#' @param annotations A \code{\link[GenomicRanges]{GRanges}} containing
#' annotations for the genome used
#' @param fragments Fragments data for the input matrix. Can be one of the
#' following:
#' \itemize{
#' \item Tabix-indexed fragment files (*.tsv.gz from 10x Genomics).
#' \item Directories created by \code{\link[BPCells]{write_fragments_dir}}.
#' \item HDF5 files created by \code{\link[BPCells]{write_fragments_hdf5}}.
#' \item One or a list of \code{IterableFragments} objects. Note that if the
#' \code{IterableFragments} is created by opening the 10x *.tsv.gz
#' (\code{\link[BPCells]{open_fragments_10x}}) directly, the cell names cannot
#' be accessed or modified by \code{\link[BPCells]{cellNames}}.
#' }
#'
#' @details
#' In this class, rows should represent genomic coordinates (e.g., peaks) while
#' columns represent samples generated from single cells.
#'
#' The extra arguments in the constructor (e.g., \code{fragments} represent the
#' main extensions implemented in the \code{ChromExperiment} class. Readers
#' are referred to the specific documentation pages for more details.
#'
#' @return
#' A \code{ChromExperiment} object
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @export
#' @docType class
ChromExperiment <- function(
    ...,
    sep = c("-", "-"),
    ranges = NULL,
    genome = NULL,
    annotations = NULL,
    fragments = NULL
) {
  SingleCellExperiment(...) %>%
    .sce_to_csce(
      sep = sep,
      ranges = ranges,
      genome = genome,
      annotations = annotations,
      fragments = fragments
    )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## fragments ###################################################################

#' Get or set fragments data
#'
#' Methods to get or set fragments data for a \code{\link{ChromExperiment}}
#' object, to hold information needed for working with fragment files.
#'
#' @param x A \code{\link{ChromExperiment}}.
#' @param ... Passed to \code{\link{openFragments}}
#'
#' @note Currently we only support storing fragments as a single
#' \code{\link[BPCells:IterableFragments-methods]{IterableFragments}} object.
#' If the input object inherits from \code{\link[Signac]{Fragment}}, it will be
#' converted into \code{IterableFragments} via \code{\link{openFragments}}.
#'
#' @returns
#' \itemize{
#' \item \code{fragments}: Retrieve the stored BPCells fragments.
#' \item \code{fragments<-}: \code{x} with updated \code{fragments}
#' }
#'
#' @seealso
#' \itemize{
#' \item \code{\link[BPCells:IterableFragments-methods]{IterableFragments}}
#' \item \code{\link[Signac]{Fragment}}
#' }
#'
#' @name fragments
NULL

#' @importFrom SingleCellExperiment int_metadata
#' @export
#' @rdname fragments
setMethod(
  f = "fragments",
  signature = "ChromExperiment",
  definition = function(x, ...) {
    return(int_metadata(x = x)[[.frag_key]])
  }
)

#' @param value \itemize{
#' \item A list of \code{IterableFragments} or \code{Fragment} objects.
#' \item A single \code{IterableFragments}  or \code{Fragment} object.
#' \item A string specifying the path to a tabix-indexed fragments file.
#' \item \code{NULL}. This will remove all \code{fragments} data.
#' }
#'
#' @importFrom SingleCellExperiment int_metadata<-
#' @export
#' @rdname fragments
setMethod(
  f = "fragments<-",
  signature = c("ChromExperiment", "IterableFragments"),
  definition = function(x, ..., value) {
    frag.cells <- cellNames(x = value)
    cells <- colnames(x = x)
    if (!any(frag.cells %in% cells)) {
      warning("No cell for the input IterableFragments is found in 'x'")
      return(x)
    }
    value <- select_cells(
      fragments = value,
      cell_selection = fastIntersect(x = frag.cells, y = cells)
    )
    cells.notfound <- setdiff(x = cells, y = cellNames(x = value))
    if (length(x = cells.notfound) > 0) {
      warning(
        length(x = cells.notfound), " cells were not found ",
        "in the input IterableFragments:\n  ",
        paste(head(x = cells.notfound), collapse = ", "), "...",
        immediate. = TRUE
      )
    }
    int_metadata(x = x)[[.frag_key]] <- value
    return(x)
  }
)

#' @importClassesFrom Signac Fragment
#' @export
#' @rdname fragments
setMethod(
  f = "fragments<-",
  signature = c("ChromExperiment", "Fragment"),
  definition = function(x, ..., value) {
    fragments(x = x) <- openFragments(x = value, ...)
    return(x)
  }
)

#' @export
#' @rdname fragments
setMethod(
  f = "fragments<-",
  signature = c("ChromExperiment", "character"),
  definition = function(x, ..., value) {
    frags <- list()
    for (i in seq_along(value)) {
      frags[[i]] <- openFragments(x = value[i], ...)
    }
    if (length(x = frags) == 0) {
      return(x)
    }
    if (length(x = frags) == 1) {
      fragments(x = x) <- frags[[1]]
      return(x)
    }
    fragments(x = x) <- Reduce(f = c, x = frags)
    return(x)
  }
)

#' @importClassesFrom S4Vectors list_OR_List
#' @export
#' @rdname fragments
setMethod(
  f = "fragments<-",
  signature = c("ChromExperiment", "list_OR_List"),
  definition = function(x, ..., value) {
    cells <- colnames(x = x)
    cells.loaded <- character()
    for (i in seq_along(value)) {
      check_inherits_for_func(
        value = value[[i]],
        classes = .support_frgas,
        func = "fragments<-",
        name = "value"
      )
      value[[i]] <- openFragments(x = value[[i]], ...)
      cells.toload <- setdiff(x = cellNames(x = value[[i]]), y = cells.loaded)
      value[[i]] <- select_cells(
        fragments = value[[i]],
        cell_selection = fastIntersect(x = cells, cells.toload)
      )
      cells.loaded <- c(cells.loaded, cellNames(x = value[[i]]))
    }
    cells.notfound <- setdiff(x = cells, y = cells.loaded)
    if (length(x = cells.notfound) > 0) {
      warning(
        length(x = cells.notfound), " cells not found in loaded fragments:\n  ",
        paste(head(cells.notfound), collapse = ", "), "...",
        immediate. = TRUE
      )
    }
    if (length(x = value) == 0) {
      return(x)
    }
    if (length(x = value) == 1) {
      fragments(x = x) <- value[[1]]
      return(x)
    }
    fragments(x = x) <- Reduce(f = c, x = value)
    return(x)
  }
)

#' @importFrom SingleCellExperiment int_metadata<-
#' @export
#' @rdname fragments
setMethod(
  f = "fragments<-",
  signature = c("ChromExperiment", "NULL"),
  definition = function(x, ..., value) {
    int_metadata(x = x)[[.frag_key]] <- NULL
    return(x)
  }
)

## annotations #################################################################

#' Get or set genomic annotation
#'
#' Methods to get or set a \code{\link[GenomicRanges:GRanges-class]{GRanges}}
#' specifying the genomic annotation in a \code{\link{ChromExperiment}}.
#'
#' @param x A \code{\link{ChromExperiment}} object.
#' @param ... Arguments passed to other methods.
#'
#' @seealso
#' \code{\link[Signac]{Annotation}}
#'
#' @name annotations
NULL

#' @importFrom SingleCellExperiment int_metadata
#' @export
#' @rdname annotations
setMethod(
  f = "annotations",
  signature = "ChromExperiment",
  definition = function(x, ...) {
    return(int_metadata(x = x)[[.annot_key]])
  }
)

#' @param value \itemize{
#' \item A \code{GRanges} object.
#' \item \code{NULL}. This will remove the existing \code{annotations}
#' }
#'
#' @importFrom SingleCellExperiment int_metadata<-
#' @importClassesFrom GenomicRanges GRanges
#' @export
#' @rdname annotations
setMethod(
  f = "annotations<-",
  signature = c("ChromExperiment", "GRanges"),
  definition = function(x, ..., value) {
    current.genome <- unique(x = genome(x = x))
    annot.genome <- unique(x = genome(x = value))
    if (!is.null(x = current.genome)) {
      if (!is.na(x = annot.genome) & (current.genome != annot.genome)) {
        stop("annotations genome does not match genome of the object")
      }
    }
    int_metadata(x = x)[[.annot_key]] <- value
    return(x)
  }
)

#' @importFrom SingleCellExperiment int_metadata<-
#' @export
#' @rdname annotations
setMethod(
  f = "annotations<-",
  signature = c("ChromExperiment", "NULL"),
  definition = function(x, ..., value) {
    int_metadata(x = x)[[.annot_key]] <- value
    return(x)
  }
)

## seqinfo #####################################################################

#' Access and modify sequence information for ChromExperiment
#'
#' Methods for accessing and modifying the \code{\link[GenomeInfoDb]{Seqinfo}}
#' object stored in a \code{\link{ChromExperiment}} object.
#'
#' @param x A \code{\link{ChromExperiment}} object.
#' @param value A \code{\link[GenomeInfoDb]{Seqinfo}} object or name of a UCSC
#' genome to store in the \code{\link{ChromExperiment}}
#'
#' @details
#' These methods are intended to be consistent with methods for
#' \code{\link[Signac]{ChromatinAssay-class}} in the \pkg{Signac} package.
#'
#' @seealso
#' \itemize{
#' \item \code{\link[GenomeInfoDb]{seqinfo}} in the \pkg{GenomeInfoDb} package.
#' \item \code{\link[Signac]{seqinfo-methods}}
#' }
#'
#' @name seqinfo-methods
NULL

#' @importFrom SingleCellExperiment int_metadata
#' @importFrom GenomeInfoDb seqinfo
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqinfo",
  signature = "ChromExperiment",
  definition = function(x) {
    return(int_metadata(x = x)[[.sinfo_key]])
  }
)

#' @importFrom SingleCellExperiment int_metadata<-
#' @importFrom GenomeInfoDb seqinfo<- Seqinfo
#' @importClassesFrom GenomeInfoDb Seqinfo
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqinfo<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    if (inherits(x = value, what = "Seqinfo")) {
      int_metadata(x = x)[[.sinfo_key]] <- value
    } else if (is(object = value, class2 = "character")) {
      int_metadata(x = x)[[.sinfo_key]] <- Seqinfo(genome = value)
    } else if (is.null(x = value)) {
      int_metadata(x = x)[[.sinfo_key]] <- NULL
    } else {
      stop(
        "Unknown object supplied. Choose a Seqinfo object ",
        "or the name of a UCSC genome"
      )
    }
    return(x)
  }
)

## seqlevels ###################################################################

#' @importFrom GenomeInfoDb seqlevels
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqlevels",
  signature = "ChromExperiment",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    }
    callGeneric()
  }
)

#' @importFrom GenomeInfoDb seqlevels<-
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqlevels<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    sinfo <- seqinfo(x = x)
    seqlevels(x = sinfo) <- value
    seqinfo(x = x) <- sinfo
    return(x)
  }
)

## seqnames ####################################################################

#' @importFrom GenomeInfoDb seqnames
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqnames",
  signature = "ChromExperiment",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    }
    callGeneric()
  }
)

#' @importFrom GenomeInfoDb seqnames<-
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqnames<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    sinfo <- seqinfo(x = x)
    seqnames(x = sinfo) <- value
    seqinfo(x = x) <- sinfo
    return(x)
  }
)

## seqlengths ##################################################################

#' @importFrom GenomeInfoDb seqlengths
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqlengths",
  signature = "ChromExperiment",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    }
    callGeneric()
  }
)

#' @importFrom GenomeInfoDb seqlengths<-
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "seqlengths<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    sinfo <- seqinfo(x = x)
    seqlengths(x = sinfo) <- value
    seqinfo(x = x) <- sinfo
    return(x)
  }
)

## genome ######################################################################

#' @importFrom GenomeInfoDb genome
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "genome",
  signature = "ChromExperiment",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    }
    callGeneric()
  }
)

#' @importFrom GenomeInfoDb genome<-
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "genome<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    seqinfo(x = x) <- value
    return(x)
  }
)

## isCircular ##################################################################

#' @importFrom GenomeInfoDb isCircular
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "isCircular",
  signature = "ChromExperiment",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    }
    callGeneric()
  }
)

#' @importFrom GenomeInfoDb isCircular<-
#' @export
#' @rdname seqinfo-methods
setMethod(
  f = "isCircular<-",
  signature = "ChromExperiment",
  definition = function(x, value) {
    sinfo <- seqinfo(x = x)
    isCircular(x = sinfo) <- value
    seqinfo(x = x) <- sinfo
    return(x)
  }
)

## Show ########################################################################

#' @param object A \code{ChromExperiment}
#'
#' @importFrom SingleCellExperiment altExpNames mainExpName reducedDimNames
#' @importFrom S4Vectors coolcat
#' @importFrom methods selectMethod show
#' @export
#' @rdname ChromExperiment
setMethod(
  f = "show",
  signature = "ChromExperiment",
  definition = function(object) {
    old_show <- selectMethod(f = "show", signature = "SummarizedExperiment")
    old_show(object)
    coolcat(
      fmt = "reducedDimNames(%d): %s\n",
      vals = reducedDimNames(x = object)
    )
    me <- mainExpName(x = object)
    if (is.null(x = me)) {
      me <- "NULL"
    }
    cat(sprintf(fmt = "mainExpName: %s\n", me))
    coolcat(fmt = "altExpNames(%d): %s\n", vals = altExpNames(x = object))
    cat("fragments: \n")
    if (!is.null(x = fragments(x = object))) {
      frag_paths <- getPath(object = fragments(x = object))
      for (i in seq_len(length.out = nrow(x = frag_paths))) {
        cat(" path:", frag_paths[i, 1], "\n")
        if (!is.na(x = frag_paths[i, 3])) {
          cat(" group:", frag_paths[i, 3], "\n")
        }
      }
    }
    .show_h5backed_path(object = object)
    return(invisible(x = NULL))
  }
)

## Subset ######################################################################

#' Subset ChromExperiment object
#'
#' Deal with the dimensions and subset of \code{\link{ChromExperiment}} object.
#' These methods are generally identical to those inherited operations, with
#' additional manipulation for \code{fragments(x)}.
#'
#' @param x A \code{\link{ChromExperiment}}.
#' @param ... Arguments passed to other methods.
#'
#' @name ChromExperiment-subset
NULL

#' @param value An object of a class specified in the S4 method signature.
#'
#' @examples
#' csce <- load_example_csce()
#' csce
#'
#' colnames(csce) <- paste0("Sorted_", colnames(csce))
#' csce
#' fragments(csce)
#'
#' @export
#' @rdname ChromExperiment-subset
setMethod(
  f = "dimnames<-",
  signature = c("ChromExperiment", "list"),
  definition = function(x, value) {
    old.colnames <- colnames(x = x)
    x <- callNextMethod()
    if (identical(x = old.colnames, y = value[[2L]])) {
      return(x)
    }
    if (is.null(x = fragments(x = x))) {
      return(x)
    }
    cellNames(x = fragments(x = x)) <- value[[2L]]
    return(x)
  }
)

#' @param i,j Indices specifying rows or columns to extract.
#' @param drop For matrices and arrays. If TRUE the result is coerced to the
#' lowest possible dimension.
#'
#' @examples
#' csce[1:10, 1:20]
#'
#' csce[, sample(colnames(csce), 20)]
#'
#' @export
#' @rdname ChromExperiment-subset
setMethod(
  f = "[",
  signature = c("ChromExperiment", "ANY", "ANY"),
  definition = function(x, i, j, ..., drop = TRUE) {
    old.frags <- fragments(x = x)
    x <- callNextMethod()
    if (ncol(x = x) == 0) {
      fragments(x = x) <- NULL
      return(x)
    }
    if (is.null(x = old.frags)) {
      return(x)
    }
    fragments(x = x) <- select_cells(
      fragments = old.frags,
      cell_selection = colnames(x = x)
    )
    return(x)
  }
)

## Validity ####################################################################

#' @importFrom SummarizedExperiment rowRanges
#' @importClassesFrom GenomicRanges GRanges
.valid_csce <- function(object) {
  msg <- NULL
  if (!inherits(x = rowRanges(x = object), what = "GRanges")) {
    msg <- msg %>%
      c(paste("'rowRanges' of", class(x = object), "must be GRanges"))
  }
  valid.frags <- is.null(x = fragments(x = object)) ||
    inherits(x = fragments(x = object), what = "IterableFragments")
  if (!valid.frags) {
    msg <- msg %>%
      c(paste("'fragments' of", class(x = object), "must be IterableFragments"))
  }
  if (length(x = msg) > 0) {
    return(msg)
  }
  return(TRUE)
}

#' @importFrom S4Vectors setValidity2
setValidity2(Class = "ChromExperiment", method = .valid_csce)

## Coerce ######################################################################

setAs(
  from = "SingleCellExperiment",
  to = "ChromExperiment",
  def = function(from) {
    ranges <- rowRanges(x = from)
    if (!inherits(x = ranges, what = "GenomicRanges")) {
      ranges <- NULL
    }
    fragments <- int_metadata(x = from)[[.frag_key]]
    genome <- int_metadata(x = from)[[.sinfo_key]]
    annotations <- int_metadata(x = from)[[.annot_key]]
    return(.sce_to_csce(
      sce = from,
      sep = c(":", "-"),
      ranges = ranges,
      genome = genome,
      annotations = annotations,
      fragments = fragments
    ))
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Format the row names of ChromExperiment
#'
#' Make sure the row names of \code{ChromExperiment} are formatted as
#' \code{"chr:start-end"}
#'
#' @importFrom SummarizedExperiment rowRanges
#'
#' @noRd
.format_csce_rownames <- function(csce) {
  rownames(x = csce) <- as.character(x = rowRanges(x = csce))
  return(csce)
}

#' Coerce SingleCellExperiment to ChromExperiment
#'
#' @importFrom SummarizedExperiment rowData rowData<- rowRanges<-
#' @importFrom GenomicRanges granges
#' @importFrom Signac StringToGRanges
#'
#' @noRd
.sce_to_csce <- function(
    sce,
    sep = c("-", "-"),
    ranges = NULL,
    genome = NULL,
    annotations = NULL,
    fragments = NULL,
    ...
) {
  # old <- S4Vectors:::disableValidity()
  # if (!isTRUE(x = old)) {
  #   S4Vectors:::disableValidity(disabled = TRUE)
  #   on.exit(S4Vectors:::disableValidity(disabled = old))
  # }
  #
  if (is.null(x = ranges)) {
    ranges <- ranges %||%
      StringToGRanges(regions = rownames(x = sce), sep = sep)
  }
  if (length(x = ranges) != nrow(x = sce)) {
    stop(
      "Length of 'ranges' does not match number of rows",
      " in SingleCellExperiment"
    )
  }
  if (!isDisjoint(x = ranges)) {
    warning(
      "Overlapping 'ranges' supplied. Ranges should be non-overlapping.",
      immediate. = TRUE, call. = FALSE
    )
  }
  row.data <- rowData(x = sce)
  rowRanges(x = sce) <- ranges
  rownames(x = sce) <- rownames(x = row.data)
  rowData(x = sce) <- row.data
  csce <- new(Class = "ChromExperiment", sce)
  csce <- .format_csce_rownames(csce = csce)
  fragments(x = csce) <- fragments
  annotations(x = csce) <- annotations
  genome(x = csce) <- genome
  return(csce)
}
