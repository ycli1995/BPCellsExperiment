
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Checking functions ###########################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom rlang inherits_any
check_inherits_for_vec <- function(x, class, name) {
  check.out <- vapply(
    X = x,
    FUN = inherits_any,
    FUN.VALUE = logical(length = 1L),
    class = class
  )
  if (all(check.out)) {
    return(invisible(x = NULL))
  }
  idx <- seq_along(along.with = x)[!check.out]
  stop(
    "The following element(s) in '", name, "'",
    " doesn't inherit from ", paste(class, collapse = ", "), ":\n  ",
    paste(idx, collapse = ", ")
  )
}

#' @importFrom rlang inherits_any
check_inherits_for_func <- function(value, classes, func, name) {
  check.class <- inherits_any(x = value, class = classes)
  if (check.class) {
    return(invisible(x = NULL))
  }
  stop(
    " Unsupported class: ", class(x = value), ". \n",
    " '", name, "' in '", func, "' should be the following classes: \n",
    paste(classes, collapse = ", "),
    call. = FALSE
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Printing functions ###########################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






