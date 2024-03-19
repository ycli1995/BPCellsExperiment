#' @importFrom methods as callGeneric callNextMethod getMethod is selectMethod
#' setGeneric setMethod slot slotNames
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 Generics ##################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Merge matrices ##############################################################

#' Merge Single-cell data matrices
#'
#' Single-cell data matrices are mostly in cell-by-feature formation. In this
#' function, matrices are merged through \code{\link{cbind}}. The union of rows
#' for input matrices are used.
#'
#' @param x A matrix
#' @param y One or more matrices of the same class or coercible to the same
#' class as \code{x}
#' @param ... Arguments passed to other methods
#'
#' @note
#' The column names along all the input matrices must be unique, or an error
#' will be raised.
#'
#' @seealso The manner of \code{mergeMatrices} is actually the same as
#' \pkg{Seurat}'s \code{\link[SeuratObject]{merge}}
#'
#' @return A single matrix of type \code{class(x)}
#'
#' @rdname mergeMatrices
#' @export mergeMatrices
mergeMatrices <- function(x, y, ...) {
  UseMethod(generic = "mergeMatrices", object = x)
}

## BPCells #####################################################################

#' Open fragments file
#'
#' @param x Can be one of the followings:
#' \itemize{
#' \item A \code{.tsv.gz} fragments file (must has \code{.tbi} index)
#' \item An HDF5 file created by \code{\link[BPCells]{write_fragments_hdf5}}.
#' \item A directory created by \code{\link[BPCells]{write_fragments_dir}}.
#' \item A \code{\link[Signac]{Fragment}} object from \pkg{Signac}.
#' }
#' @param group Control the link name of HDF5 group, where the fragments data
#' are stored. Default will be 'fragments'.
#' @param ... Arguments passed to other metheds.
#'
#' @return An \code{IterableFragments} object
#'
#' @seealso \code{\link[BPCells]{open_fragments_10x}} and
#' \code{\link[BPCells]{open_fragments_hdf5}}
#'
#' @rdname openFragments
#' @export openFragments
openFragments <- function(x, group = NULL, ...) {
  UseMethod(generic = "openFragments", object = x)
}

#' Fetch the file path information for an on-disk object
#'
#' Get the file path and the link (for HDF5-backed) name for an on-disk object.
#'
#' @param object A file-backed object.
#' @param ... Arguments passed to other metheds.
#'
#' @return A \code{\link{data.frame}} with three columns: path, type and group.
#'
#' @details
#' This function can be useful to check the file path of an object to avoid IO
#' conflicts.
#'
#' @rdname getPath
#' @export getPath
getPath <- function(object, ...) {
  UseMethod(generic = "getPath", object = object)
}

## HDF5 ########################################################################

#' Write a SingleCellExperiment to HDF5 file
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object
#' @param file Path of the HDF5 file
#' @param ... Arguments passed to other metheds.
#'
#' @rdname writeH5SCE
#' @export writeH5SCE
writeH5SCE <- function(object, file, ...) {
  UseMethod(generic = "writeH5SCE", object = object)
}

#' Write a SingleCellMultiExperiment to HDF5 file
#'
#' @param object A \code{\link{SingleCellMultiExperiment}} object
#' @param file Path of the HDF5 file
#' @param ... Arguments passed to other metheds.
#'
#' @rdname writeH5SCME
#' @export writeH5SCME
writeH5SCME <- function(object, file, ...) {
  UseMethod(generic = "writeH5SCME", object = object)
}

#' Write a SingleCellExperiment to H5AD file
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object
#' @param file Path of the H5AD file
#' @param ... Arguments passed to other metheds.
#'
#' @rdname exportH5AD
#' @export exportH5AD
exportH5AD <- function(object, file, ...) {
  UseMethod(generic = "exportH5AD", object = object)
}

## Diet object #################################################################

#' Slim down a SingleCellExperiment object
#'
#' Keep only certain aspects of the \code{SingleCellExperiment} object.
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object
#' @param ... Arguments passed to other metheds.
#'
#' @rdname dietSCE
#' @export dietSCE
dietSCE <- function(object, ...) {
  UseMethod(generic = "dietSCE", object = object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Generics ##################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Rename functions ############################################################

#' @export
#' @rdname set-dimnames
setGeneric(
  "setRownames",
  function(x, new.names, ...) standardGeneric("setRownames")
)

#' @export
#' @rdname set-dimnames
setGeneric(
  "setColnames",
  function(x, new.names, ...) standardGeneric("setColnames")
)

## SingleCellExperiment ########################################################

#' @export
#' @rdname fetchData
setGeneric(
  "fetchData",
  function(object, ...) standardGeneric("fetchData")
)

#' @export
#' @rdname fetchData
setGeneric(
  "variableFeatures",
  function(object, ...) standardGeneric("variableFeatures")
)

#' @export
#' @rdname fetchData
setGeneric(
  "variableFeatures<-",
  function(object, ..., value) standardGeneric("variableFeatures<-")
)

## ChromatinExperiment #########################################################

#' @export
#' @rdname fragments
setGeneric(
  "fragments",
  function(x, ...) standardGeneric("fragments")
)

#' @export
#' @rdname fragments
setGeneric(
  "fragments<-",
  function(x, ..., value) standardGeneric("fragments<-")
)

#' @export
#' @rdname annotations
setGeneric(
  "annotations",
  function(x, ...) standardGeneric("annotations")
)

#' @export
#' @rdname annotations
setGeneric(
  "annotations<-",
  function(x, ..., value) standardGeneric("annotations<-")
)

## MultiAssayExperiment ########################################################

#' @export
#' @rdname experiment
setGeneric(
  "experiment",
  function(x, e, ...) standardGeneric("experiment")
)

## SingleCellMultiExperiment ###################################################

#' @export
#' @rdname SingleCellMultiExperiment-internal
setGeneric(
  "int_SCE",
  function(x, ...) standardGeneric("int_SCE")
)

#' @export
#' @rdname SingleCellMultiExperiment-internal
setGeneric(
  "int_SCE<-",
  function(x, ..., value) standardGeneric("int_SCE<-")
)

#' @export
#' @rdname defaultExp
setGeneric(
  "defaultExp",
  function(x, ...) standardGeneric("defaultExp")
)

#' @export
#' @rdname defaultExp
setGeneric(
  "defaultExp<-",
  function(x, ..., value) standardGeneric("defaultExp<-")
)

