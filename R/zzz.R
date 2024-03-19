
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constant #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Field keys ##################################################################
.int_cd <- "int_colData"

.int_emd <- "int_elementMetadata"

.int_md <- "int_metadata"

.colp_key <- "colPairs"

.rowp_key <- "rowPairs"

.alt_key <- "altExps"

.red_key <- "reducedDims"

.frag_key <- "fragments"

.annot_key <- "annotations"

.sinfo_key <- "seqinfo"

.path_key <- "filepath"

## BPCells version #############################################################

.bpcells_matrix_attr <- sapply(
  X = c("unpacked", "packed"),
  FUN = paste,
  c("uint", "float", "double"),
  sep = "-"
) %>%
  as.vector() %>%
  paste0("-matrix-v2")

## Supported classes ###########################################################

.support_assays <- c("IterableMatrix", "dgCMatrix")

.support_frgas <- c("IterableFragments", "Fragment")

.support_altExps <- "BPCExperiment"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global variables #############################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.BPCE.envs <- new.env(parent = emptyenv())

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package ######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @import BPCells
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom tools file_path_as_absolute
#' @importFrom rlang is_false is_true
#' @importFrom easy.utils fastIntersect verboseMsg
#' @importMethodsFrom Matrix t
#' @keywords internal
"_PACKAGE"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hooks ########################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  .init_frags_filemap()
  return(invisible(x = NULL))
}
