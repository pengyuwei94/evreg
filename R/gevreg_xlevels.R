gevreg_xlevels <- function(fo, data) {
  # Get all variable names used on RHSs of formulae
  getVars <- function(fo) {
    all.vars(update(fo, 0~.))
  }
  allVars <- unique(unlist(lapply(fo, getVars)))
  # Get rid of variables not in data, get classes, then get rid of non-factors
  data <- data[, allVars, drop=FALSE]
  classes <- sapply(data, class)
  wh <- classes %in% c("factor", "ordered", "character")
  data <- data[, wh, drop=FALSE]
  classes <- classes[wh]
  data[classes == "character"] <- lapply(data[classes == "character"],
                                         as.factor)
  # Get a single named list containing all levels
  xlevels <- lapply(data, levels)
  # Split it by formula
  res <- lapply(fo, getVars)
  return(lapply(res, function(X, wh) wh[names(wh) %in% X], wh=xlevels))
}
