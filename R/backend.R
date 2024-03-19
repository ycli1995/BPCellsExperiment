
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## openFragments ###############################################################

#' @examples
#' frags_file <- system.file(
#'   "extdata", "pbmc_sorted_fragments.tsv.gz",
#'   package = "BPCellsExperiment"
#' )
#' tmp.file <- tempfile("pbmc_sorted_fragments", fileext = ".tsv.gz")
#' file.copy(frags_file, tmp.file)
#'
#' ## Open using file name
#' frags <- openFragments(tmp.file)
#'
#' @importFrom hdf5r is_hdf5
#' @importFrom hdf5r.Extra h5AbsLinkName
#' @importFrom dplyr filter select
#' @importFrom tibble add_row
#' @importFrom tools md5sum
#' @importFrom rlang .data
#' @export
#' @rdname openFragments
#' @method openFragments character
openFragments.character <- function(x, group = NULL, ...) {
  group <- group %||% "fragments"
  group <- h5AbsLinkName(name = group)
  path <- file_path_as_absolute(x = x)
  is.dir <- dir.exists(paths = path)
  is.tsv <- grepl(pattern = "tsv(\\.gz)?$", x = path)
  is.h5 <- is_hdf5(name = path)
  if (is.h5) {
    return(open_fragments_hdf5(path = path, group = group, ...))
  }
  if (is.dir) {
    return(open_fragments_dir(dir = path, ...))
  }
  if (!is.tsv) {
    stop("'x' must be an H5 file, a directory or an indexed .tsv")
  }
  hash <- md5sum(files = path) %>%
    paste(collapse = "_")
  if (hash %in% .BPCE.envs$frags_filemap$hash) {
    # The query .tsv file has been cached into an H5 file,
    # just open it.
    file.map <- .BPCE.envs$frags_filemap %>%
      filter(.data$hash == hash) %>%
      select(.data$bpcells_h5, .data$bpcells_group)
    return(open_fragments_hdf5(
      path = file.map$bpcells_h5,
      group = file.map$bpcells_group,
      ...
    ))
  }
  # Cache the query .tsv into an H5 file
  tmp.file <- basename(path = path) %>%
    gsub(pattern = "\\.tsv(\\.gz)?$", replacement = ".") %>%
    tempfile(tmpdir = dirname(path = path), fileext = ".h5") %>%
    normalizePath(mustWork = FALSE)
  frags <- open_fragments_10x(path = path, ...) %>%
    write_fragments_hdf5(path = tmp.file, group = group)
  .BPCE.envs$frags_filemap <- .BPCE.envs$frags_filemap %>%
    add_row(
      tsv = path,
      hash = hash,
      bpcells_h5 = file_path_as_absolute(x = tmp.file),
      bpcells_group = group
    )
  return(frags)
}

#' @examples
#' file.copy(paste0(frags_file, ".tbi"), paste0(tmp.file, ".tbi"))
#'
#' ## Open using Signac::Fragment
#' frags <- Signac::CreateFragmentObject(tmp.file, cells = cellNames(frags))
#' frags <- openFragments(frags)
#'
#' @importFrom rlang is_empty
#' @importClassesFrom Signac Fragment
#' @export
#' @rdname openFragments
#' @method openFragments Fragment
openFragments.Fragment <- function(x, group = NULL, ...) {
  file <- slot(object = x, name = "path")
  cells <- slot(object = x, name = "cells")
  frags <- openFragments(x = file, group = group, ...)
  if (!is_empty(x = cells)) {
    frags <- select_cells(fragments = frags, cell_selection = cells)
    cellNames(x = frags) <- names(x = cells)
  }
  return(frags)
}

#' @examples
#' ## Open an already existing Fragment
#' frags <- openFragments(frags)
#'
#' @export
#' @rdname openFragments
#' @method openFragments IterableFragments
openFragments.IterableFragments <- function(x, group = NULL, ...) {
  return(x)
}

#' @section Check H5 files for opened fragments:
#' When \code{openFragments} meets a \code{Fragments} object or a \code{.tsv.gz}
#' file, it will first check whether there is already an h5 file, which was
#' created by \code{\link[BPCells]{write_fragments_hdf5}} and cached the
#' fragments data to open. If so, it will directly open the h5 file. Otherwise,
#' it will calculate the \code{\link[tools]{md5sum}} of the \code{.tsv.gz}, then
#' cache it into a new h5 file. Use \code{showFragmentsFileMap()} to get the
#' file mapping table created by current R session.
#'
#' @export
#' @rdname openFragments
showFragmentsFileMap <- function() {
  .BPCE.envs$frags_filemap
}

#' @section Clear temporary fragments files:
#' When all needed fragments data in current R session have been stored in
#' proper destination paths, \code{clearFragmentsFileMap()} can be called to
#' delete all the temporary h5 files shown by \code{showFragmentsFileMap()}, so
#' that the disk storage can be saved.
#'
#' @export
#' @rdname openFragments
clearFragmentsFileMap <- function() {
  unlink(x = .BPCE.envs$frags_filemap$bpcells_h5, force = TRUE)
  .init_frags_filemap()
}

## getPath #####################################################################

#' @examples
#' ### For merged BPCells matrix
#' bpce <- load_example_sce()
#' getPath(counts(bpce))
#'
#' @export
#' @rdname getPath
#' @method getPath IterableFragments
getPath.IterableFragments <- function(object, ...) {
  return(.get_path_for_bpcells(object = object, ...))
}

#' @examples
#' ### For merged BPCells fragments
#' cbpce <- load_example_csce()
#' getPath(counts(cbpce))
#' getPath(fragments(cbpce))
#'
#' ### For merged BPCells object
#' cbpce2 <- load_example_csce("unsorted")
#' cbpce <- merge(cbpce, cbpce2)
#' getPath(counts(cbpce))
#' getPath(fragments(cbpce))
#'
#' @export
#' @rdname getPath
#' @method getPath IterableMatrix
getPath.IterableMatrix <- function(object, ...) {
  return(.get_path_for_bpcells(object = object, ...))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Get path for BPCells objects ################################################

#' Fetch all paths of backend files for BPCells object
#'
#' @importFrom hdf5r.Extra h5AbsLinkName
#' @noRd
.get_path_for_bpcells <- function(object, ...) {
  slots <- slotNames(x = object)
  use.slots <- c("matrix_list", "fragments_list", "path", "dir")
  empty.out <- data.frame(
    path = character(),
    type = character(),
    group = character()
  )
  out <- list(empty.out)
  if (!any(use.slots %in% slots)) {
    for (i in slots) {
      slot.obj <- slot(object = object, name = i)
      is.iter_frag <- inherits(x = slot.obj, what = "IterableFragments")
      is.iter_mat <- inherits(x = slot.obj,  what = "IterableMatrix")
      if (is.iter_mat || is.iter_frag) {
        out[[i]] <- .get_path_for_bpcells(object = slot.obj, ...)
      }
    }
    names(x = out) <- NULL
    return(do.call(what = rbind, args = out))
  }
  if ("path" %in% slots) {
    path <- slot(object = object, name = "path") %>%
      file_path_as_absolute()
    group <- NA
    type <- "tsv"
    if (inherits(x = object,  what = "10xMatrixH5")) {
      type <- "h5"
    }
    if ("group" %in% slots) {
      group <- h5AbsLinkName(name = slot(object = object, name = "group"))
      type <- "h5"
    }
    return(data.frame(path = path, type = type, group = group))
  }
  if ("dir" %in% slots) {
    path <- slot(object = object, name = "dir") %>%
      file_path_as_absolute()
    return(data.frame(path = path, type = "dir", group = NA))
  }
  iter_obj.list <- list()
  if ("fragments_list" %in% slots) {
    iter_obj.list <- c(
      iter_obj.list,
      slot(object = object, name = "fragments_list")
    )
  }
  if ("matrix_list" %in% slots) {
    iter_obj.list <- c(
      iter_obj.list,
      slot(object = object, name = "matrix_list")
    )
  }
  for (i in seq_along(iter_obj.list)) {
    out[[i]] <- .get_path_for_bpcells(object = iter_obj.list[[i]], ...)
  }
  names(x = out) <- NULL
  return(do.call(what = rbind, args = out, ...))
}

#' Check if a path has been occupied by BPCells object
#'
#' @param x The BPCells object to be check.
#' @param path Path to be check.
#' @param name Name of HDF5 link to be check.
#'
#' @noRd
#' @importFrom hdf5r.Extra h5AbsLinkName
.check_bpcells_file_occupy <- function(x, path, name = NULL) {
  path.use <- .get_path_for_bpcells(object = x)
  path <- normalizePath(path = path, mustWork = FALSE)
  name <- name %iff% h5AbsLinkName(name = name)
  path.use <- path.use[path.use$path %in% path, , drop = FALSE]
  if (any(name %in% path.use$group)) {
    stop(
      "\n  The destination has been occupied ",
      "by the iterable object to be written:",
      "\n  Path: ", path,
      "\n  Group: ", name
    )
  }
  if (nrow(x = path.use) > 0) {
    stop(
      "\n  The destination has been occupied ",
      "by the iterable object to be written:",
      "\n  Path: ", path
    )
  }
  return(invisible(x = NULL))
}

## Fragments file map ##########################################################

.init_frags_filemap <- function() {
  .BPCE.envs$frags_filemap <- data.frame(
    tsv = character(),
    hash = character(),
    bpcells_h5 = character(),
    bpcells_group = character()
  )
  return(invisible(x = NULL))
}
