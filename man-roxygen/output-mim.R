#' @return A list with two elements: \code{selection}, a vector of indices of the selected features in the selection order, and \code{score}, a vector of corresponding feature scores.
#' Names of both vectors will correspond to the names of features in \code{X}.
#' Both vectors will be at most of a length \code{k}, as the selection will stop as soon as all the remaining features will have a score of zero.
#' This may happen during initial selection, in which case both vectors will be empty.


