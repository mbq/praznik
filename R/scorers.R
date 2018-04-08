# Feature scorers

#' Calculate mutual information of all features
#'
#' Calculates mutual information between each attribute and the decision, that is
#' \deqn{I(X,Y).}
#' @template input
#' @return A numerical vector with mutual information scores, with names copied from \code{X}.
#' @examples
#' miScores(iris[,-5],iris$Species)
#' @export
miScores<-function(X,Y,threads=0)
 .Call(C_mi,X,Y,as.integer(threads))

#' Calculate conditional mutual information of all features
#'
#' Calculates mutual information between each attributes and the decision, that is
#' \deqn{I(X,Y|Z).}
#' @template input
#' @param Z Condition; should be given as a factor, but other options are accepted, as for attributes. 
#' @return A numerical vector with conditional mutual information scores, with names copied from \code{X}.
#' @examples
#' cmiScores(iris[,-5],iris$Species,iris$Sepal.Length)
#' @export
cmiScores<-function(X,Y,Z,threads=0)
 .Call(C_cmi_jmi,X,Y,Z,791L,as.integer(threads))

#' Calculate joint mutual information of all features
#'
#' Calculated mutual information between each attribute joint with some other vector \code{Z} with the decision, that is
#' \deqn{I(X,Z;Y).}
#' This is the same as conditional mutual information between X and Y plus a constant that depends on Y and Z, that is
#' \deqn{I(X,Z;Y)=I(X;Y|Z)+I(Y;Z).}
#' @template input
#' @param Z Other vector; should be given as a factor, but other options are accepted, as for attributes. 
#' @return A numerical vector with joint mutual information scores, with names copied from \code{X}.
#' @examples
#' jmiScores(iris[,-5],iris$Species,iris$Sepal.Length)
#' @export
jmiScores<-function(X,Y,Z,threads=0)
 .Call(C_cmi_jmi,X,Y,Z,792L,as.integer(threads))

#' Calculate normalised joint mutual information of all features
#'
#' Calculated normalised mutual information between each attribute joint with some other vector \code{Z} with the decision, that is
#' \deqn{\frac{I(X,Z;Y)}{H(X,Y,Z)}.}
#' This is the same as in the criterion used by \code{\link{DISR}} and \code{\link{NJMIM}}.
#' @template input
#' @param Z Other vector; should be given as a factor, but other options are accepted, as for attributes. 
#' @return A numerical vector with the normalised joint mutual information scores, with names copied from \code{X}.
#' @examples
#' njmiScores(iris[,-5],iris$Species,iris$Sepal.Length)
#' @export
njmiScores<-function(X,Y,Z,threads=0)
 .Call(C_cmi_jmi,X,Y,Z,793L,as.integer(threads))

#' Calculate Gini impurity scores of all features
#'
#' Calculates Gini impurity between each attribute and the decision, that is
#' \deqn{G(X;Y)=\sum_{xy} \frac{p_{xy}^2}{p_x}-\sum_y p_y^2.}
#' @template input
#' @return A numerical vector with Gini impurity scores, with names copied from \code{X}.
#' @examples
#' impScores(iris[,-5],iris$Species)
#' @export
impScores<-function(X,Y,threads=0)
 .Call(C_im,X,Y,as.integer(threads))
