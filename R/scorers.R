#' Calculate mutual information of all features
#'
#' Calculates mutual information between all attributes and the decision.
#' @template input
#' @return A numerical vector with mutual information scores, with names copied from \code{X}.
#' @examples
#' mi(iris[,-5],iris$Species)
#' @export
mi<-function(X,Y)
 .Call(C_mi,X,Y)
