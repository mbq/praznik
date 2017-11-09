#' @param X Attribute table, given as a data frame with only factor columns. \code{NA}s are not allowed.
#' @param Y Decision attribute; must be a factor. \code{NA}s are not allowed.
#' @param k Number of attributes to select. Must not exceed \code{ncol(X)}.
#' @return A list with two elements: \code{selection}, a vector of names of the selected features in the selection order (note that it may be shorter than \code{k}), and \code{scores}, a vector of feature scores.
