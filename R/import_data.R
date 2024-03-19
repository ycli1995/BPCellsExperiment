
#' Load in data from 10X
#'
#' Read the sparse data matrix provided by 10X genomics.
#'
#' @param path Directory containing the matrix.mtx, genes.tsv (or features.tsv),
#' and barcodes.tsv files provided by 10X.
#' @param gene.column Specify which column of genes.tsv or features.tsv to use
#' for gene names; default is 2
#' @param cell.column Specify which column of barcodes.tsv to use for cell
#' names; default is 1
#' @param use.BPCells Whether or not to use \pkg{BPCells} to hold the matrix.
#'
#' @importFrom utils read.delim read.table
#' @export
read10xMtx <- function(
    path,
    gene.column = 2,
    cell.column = 1,
    use.BPCells = FALSE
) {
  path <- path[1]
  has_dt <- requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("R.utils", quietly = TRUE)
  if (!dir.exists(paths = path)) {
    stop("Directory provided does not exist: ", path)
  }
  barcode.loc <- file.path(path, "barcodes.tsv")
  gene.loc <- file.path(path, "genes.tsv")
  features.loc <- file.path(path, "features.tsv.gz")
  matrix.loc <- file.path(path, "matrix.mtx")
  pre_ver_3 <- file.exists(gene.loc)
  if (!pre_ver_3) {
    barcode.loc <- paste0(barcode.loc, ".gz")
    matrix.loc <- paste0(matrix.loc, ".gz")
  }
  if (!file.exists(barcode.loc)) {
    stop("Barcode file missing: ", basename(path = barcode.loc))
  }
  if (!pre_ver_3 && !file.exists(features.loc)) {
    stop("Gene name or features file missing: ", basename(path = features.loc))
  }
  if (!file.exists(matrix.loc)) {
    stop("Expression matrix file missing: ", basename(path = matrix.loc))
  }
  ## Read cell barcodes
  if (has_dt) {
    cell.barcodes <- data.table::fread(
      input = barcode.loc,
      header = FALSE,
      data.table = FALSE
    )
  } else {
    cell.barcodes <- read.table(
      file = barcode.loc,
      header = FALSE,
      sep = "\t",
      row.names = NULL
    )
  }
  if (ncol(x = cell.barcodes) > 1) {
    rownames(x = cell.barcodes) <- cell.barcodes[, cell.column]
  } else {
    rownames(x = cell.barcodes) <- readLines(con = barcode.loc)
  }
  cell.barcodes <- data.frame(row.names = rownames(x = cell.barcodes))

  ## Read features
  features.file <- ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc)
  if (has_dt) {
    features <- data.table::fread(
      input = features.file,
      header = FALSE,
      data.table = FALSE
    )
  } else {
    features <- read.delim(
      file = features.file,
      header = FALSE,
      stringsAsFactors = FALSE
    )
  }
  if (anyNA(x = features[, gene.column])) {
    warning(
      "Replacing NA features with ID from the opposite column requested",
      call. = FALSE, immediate. = TRUE
    )
    na.features <- which(x = is.na(x = features[, gene.column]))
    repl.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
    features[na.features, gene.column] <- features[na.features, repl.column]
  }
  rownames(x = features) <- make.unique(names = features[, gene.column])
  if (ncol(x = features) >= 1) {
    colnames(x = features)[1] <- "ID"
  }
  if (ncol(x = features) >= 2) {
    colnames(x = features)[2] <- "Name"
  }
  if (ncol(x = features) >= 3) {
    colnames(x = features)[3] <- "Type"
    features <- features[, 1:3]
  }
  sce.func <- ifelse(
    test = use.BPCells,
    yes = .read_10x_mtx_bpcells,
    no = .read_10x_mtx_sce
  )
  sce <- sce.func(matrix.loc, features, colData = cell.barcodes)
  if (!"Type" %in% colnames(x = features)) {
    return(sce)
  }
  data_types <- factor(x = features$Type)
  all_types <- levels(x = data_types)
  if (length(x = all_types) == 1) {
    return(sce)
  }
  message(
    "10X data contains more than one type: ", paste(all_types, collapse = ", "),
    "\nReturn a list containing ", class(x = sce), " of each type."
  )
  expr_name <- "Gene Expression"
  if (expr_name %in% all_types) {
    all_types <- c(expr_name, all_types[-which(x = all_types == expr_name)])
  }
  sce <- lapply(X = sce, FUN = function(x) sce[data_types == x, , drop = FALSE])
  names(x = sce) <- all_types
  return(sce)
}

#' Read 10X hdf5 file
#'
#' Read count matrix from 10X CellRanger hdf5 file.
#'
#' @param path Path to HDF5 file]
#' @param use.names Label row names with feature names rather than ID numbers.
#' @param use.BPCells Whether or not to use \pkg{BPCells} to hold the matrix.
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom hdf5r.Extra h5Read h5List
#' @importFrom Matrix sparseMatrix
#' @export
read10xH5 <- function(path, use.names = TRUE, use.BPCells = FALSE) {
  path <- file_path_as_absolute(x = path)
  groups <- h5List(x = path)
  is.v3 <- "matrix" %in% groups
  read_features_func <- ifelse(
    test = is.v3,
    yes = .read_10x_h5_features_v3,
    no = .read_10x_h5_features_v2
  )
  output <- list()
  for (genome in groups) {
    features <- read_features_func(path, genome)
    barcodes <- h5Read(x = path, name = file.path(genome, "barcodes"))
    if (use.BPCells) {
      mat <- open_matrix_10x_hdf5(path = path)
    } else {
      mat <- sparseMatrix(
        i = h5Read(x = path, name = file.path(genome, "indices")),
        p = h5Read(x = path, name = file.path(genome, "indptr")),
        x = h5Read(x = path, name = file.path(genome, "data")),
        dims = h5Read(x = path, name = file.path(genome, "shape")),
        index1 = FALSE,
        repr = "C"
      )
    }
    if (use.names) {
      rownames(x = mat) <- make.unique(names = features$Name)
    } else {
      rownames(x = mat) <- make.unique(names = features$ID)
    }
    colnames(x = mat) <- barcodes
    if (use.BPCells) {
      mat <- write_matrix_dir(mat = mat, dir = tempfile("tenx_matrix_h5"))
      mat <- SingleCellExperiment(
        assays = list(counts = mat),
        rowData = DataFrame(features)
      )
    } else {
      mat <- SingleCellExperiment(
        assays = list(counts = mat),
        rowData = DataFrame(features)
      )
    }
    if ("Type" %in% colnames(x = features)) {
      types <- unique(x = features$Type)
      if (length(x = types) > 1) {
        message(
          "Genome ", genome, " has multiple modalities, ",
          "returning a list of SingleCellExperiment for this genome."
        )
        mat <- sapply(
          X = types,
          FUN = function(x) mat[features$Type == x, ],
          simplify = FALSE,
          USE.NAMES = TRUE
        )
        gc(verbose = FALSE)
      }
    }
    output[[genome]] <- mat
  }
  if (length(x = output) == 1) {
    return(output[[genome]])
  }
  return(output)
}

#' Read H5AD file
#'
#' Create a \code{\link[SingleCellExperiment]{SingleCellExperiment}} from a
#' H5AD file.
#'
#' @param path A H5AD file
#' @param name Name of the HDF5 group where the AnnData is stored.
#' @param layers Which \code{layers} to read. If \code{NULL}, read all layers.
#' @param obsm Which \code{obsm} to read. If \code{NULL}, read all \code{obsm}.
#' @param obsp,varp Which \code{obsp} (\code{varp}) to read. If \code{NULL},
#' read all \code{obsp} (\code{varp}).
#' @param uns Which slots of \code{uns} data to read. If \code{NULL}, read all
#' slots of \code{uns}.
#' @param raw Whether or not to read the \code{raw} of AnnData. By default is
#' \code{TRUE}.
#' @param use.BPCells Whether or not to use \pkg{BPCells} to hold the matrix.
#'
#' @importFrom S4Vectors metadata<-
#' @importFrom SummarizedExperiment assay<- rowData
#' @importFrom SingleCellExperiment altExp<- colPair<- LinearEmbeddingMatrix
#' reducedDim<- rowPair<-
#' @importFrom hdf5r.Extra h5Read h5List
#' @export
readH5AD <- function(
    path,
    name = "/",
    layers = NULL,
    obsm = NULL,
    obsp = NULL,
    varp = NULL,
    uns = NULL,
    raw = NULL,
    use.BPCells = FALSE
) {
  path <- file_path_as_absolute(x = path)
  all.slots <- h5List(x = path, name = name)
  obs <- h5Read(x = path, name = file.path(name, "obs"))
  sce <- .read_h5ad_raw(
    path = path,
    name = name,
    obs = obs,
    use.BPCells = use.BPCells
  )
  var <- rowData(x = sce)
  all.layers <- h5List(x = path, name = file.path(name, "layers"))
  layers <- layers %||% all.layers
  layers <- intersect(x = layers, y = all.layers)
  for (i in layers) {
    assay(x = sce, i = i) <- .read_h5ad_layer(
      path = path,
      name = file.path(name, "layers", i),
      obs = obs,
      var = var,
      use.BPCells = use.BPCells
    )
  }
  all.obsm <- h5List(x = path, name = file.path(name, "obsm"))
  obsm <- obsm %||% all.obsm
  obsm <- intersect(x = obsm, y = all.obsm)
  for (i in obsm) {
    obsm <- h5Read(x = path, name = file.path(name, "obsm", i))
    rownames(x = obsm) <- colnames(x = sce)
    colnames(x = obsm) <- paste0(i, "_", seq(ncol(x = obsm)))
    reducedDim(x = sce, type = i) <- LinearEmbeddingMatrix(
      sampleFactors = obsm,
      featureLoadings = matrix(data = 0, nrow = 0, ncol = ncol(x = obsm))
    )
  }
  all.obsp <- h5List(x = path, name = file.path(name, "obsp"))
  obsp <- obsp %||% all.obsp
  obsp <- intersect(x = obsp, y = all.obsp)
  for (i in obsp) {
    colPair(x = sce, type = i) <- h5Read(
      x = path,
      name = file.path(name, "obsm", i)
    )
  }
  all.varp <- h5List(x = path, name = file.path(name, "varp"))
  varp <- varp %||% all.varp
  varp <- intersect(x = varp, y = all.varp)
  for (i in varp) {
    rowPair(x = sce, type = i) <- h5Read(
      x = path,
      name = file.path(name, "varp", i)
    )
  }
  all.uns <- h5List(x = path, name = file.path(name, "uns"))
  uns <- uns %||% all.uns
  uns <- intersect(x = uns, y = all.uns)
  for (i in uns) {
    metadata(x = sce)[[i]] <- h5Read(x = path, name = file.path(name, "uns", i))
  }
  raw <- raw %||% TRUE
  if (!"raw" %in% all.slots) {
    return(sce)
  }
  raw.slots <- h5List(x = path, name = file.path(name, "raw"))
  if (!(length(x = raw.slots) > 0 && is_true(x = raw))) {
    return(sce)
  }
  altExp(x = sce, e = "raw") <- .read_h5ad_raw(
    path = path,
    name = file.path(name, "raw"),
    obs = obs,
    use.BPCells = use.BPCells
  )
  return(sce)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal #####################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom BPCells import_matrix_market
.read_10x_mtx_bpcells <- function(mtx_path, rowData, colData) {
  mat <- import_matrix_market(
    mtx_path = mtx_path,
    row_names = rownames(x = rowData),
    col_names = rownames(x = colData)
  )
  return(SingleCellExperiment(
    assays = list(counts = mat),
    rowData = DataFrame(rowData),
    colData = DataFrame(colData)
  ))
}

#' @importFrom Matrix readMM
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom Matrix CsparseMatrix
.read_10x_mtx_sce <- function(mat_path, rowData, colData) {
  mat <- readMM(file = mat_path)
  rownames(x = mat) <- rownames(x = rowData)
  colnames(x = mat) <- rownames(x = colData)
  mat <- as(object = mat, Class = "CsparseMatrix")
  gc(verbose = FALSE)
  return(SingleCellExperiment(
    assays = list(counts = mat),
    rowData = DataFrame(rowData),
    colData = DataFrame(colData)
  ))
}

#' @importFrom hdf5r.Extra h5Read
.read_10x_h5_features_v3 <- function(path, name) {
  features <- h5Read(x = path, name = file.path(name, "features"))
  features[["_all_tag_keys"]] <- NULL
  features <- as.data.frame(x = features, stringsAsFactors = FALSE)
  features <- features[, c("id", "name", "feature_type", "genome")]
  colnames(features)[1:3] <- c("ID", "Name", "Type")
  return(features)
}

#' @importFrom hdf5r.Extra h5Read
.read_10x_h5_features_v2 <- function(path, name) {
  features <- data.frame(
    ID = h5Read(x = path, name = file.path(name, "genes")),
    Name = h5Read(x = path, name = file.path(name, "gene_names")),
    stringsAsFactors = FALSE
  )
  return(features)
}

#' @importFrom HDF5Array HDF5Array HDF5ArraySeed
#' @importFrom S4Vectors new2
#' @importClassesFrom HDF5Array Dense_H5ADMatrixSeed
.read_h5ad_layer_dense <- function(path, name, obs, var) {
  dimnames <- list(rownames(x = var), rownames(x = obs))
  ans0 <- HDF5ArraySeed(filepath = path, name = name)
  if (length(x = dim(x = ans0)) == 2L) {
    return(new2("Dense_H5ADMatrixSeed", ans0, dimnames = dimnames))
  }
  warning(
    "HDF5 dataset '", name, "' in file '",
    path, "' does not have exactly 2 dimensions. ",
    "Using 'HDF5Array' to access this dataset.",
    immediate. = TRUE, call. = FALSE
  )
  mat <- HDF5Array(filepath = ans0)
  dimnames(x = mat) <- dimnames
  return(mat)
}

#' @importFrom HDF5Array HDF5Array
#' @importFrom hdf5r.Extra h5Read is.H5Group
.read_h5ad_layer <- function(path, name, obs, var, use.BPCells = FALSE) {
  dimnames <- list(rownames(x = var), rownames(x = obs))
  if (use.BPCells && !is.H5Group(file = path, name = name)) {
    warning(
      "'BPCells' doesn't support dense matrix: ", name,
      "\nUse 'H5ADMatrix' instead.",
      immediate. = TRUE, call. = FALSE
    )
    return(.read_h5ad_layer_dense(
      path = path,
      name = name,
      obs = obs,
      var = var
    ))
  }
  if (use.BPCells) {
    mat <- open_matrix_anndata_hdf5(path = path, group = name)
    dimnames(x = mat) <- dimnames
    return(write_matrix_dir(mat = mat, dir = tempfile("anndata_matrix_h5")))
  }
  mat <- h5Read(x = path, name = name)
  dimnames(x = mat) <- dimnames
  return(mat)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom hdf5r.Extra h5List h5Read
#' @importClassesFrom DelayedArray DelayedArray
.read_h5ad_raw <- function(path, name, obs, use.BPCells = FALSE) {
  var <- h5Read(x = path, name = file.path(name, "var"))
  mat <- .read_h5ad_layer(
    path = path,
    name = file.path(name, "X"),
    obs = obs,
    var = var,
    use.BPCells = use.BPCells
  )
  if (inherits(x = mat, what = "DelayedArray")) {
    use.BPCells <- FALSE
  }
  return(SingleCellExperiment(
    assays = list(X = mat),
    rowData = DataFrame(var),
    colData = DataFrame(obs)
  ))
}
