#' Calculate mutual information of all features
#'
#' Calculates mutual information between each attribute and the decision, that is
#; \deqn{I(X,Y).}
#' @template input
#' @return A numerical vector with mutual information scores, with names copied from \code{X}.
#' @examples
#' mi(iris[,-5],iris$Species)
#' @export
mi<-function(X,Y)
 .Call(C_mi,X,Y)


#' Calculate conditional mutual information of all features
#'
#' Calculates mutual information between each attributes and the decision, that is
#' \deqn{I(X,Y|Z).}
#' @template input
#' @param Z Condition; should be given as a factor, but other options are accepted, as for attributes. 
#' @return A numerical vector with conditional mutual information scores, with names copied from \code{X}.
#' @examples
#' mi(iris[,-5],iris$Species)
#' @export
cmi<-function(X,Y,Z)
 mi(X,factor(sprintf("%s%s",Y,Z)))-mi(X,Z) #TODO: Replace placeholder

#TODO: Remove this stuff?
#TODO: NJMI could be useful

#' Calculate joint mutual information of all features
#'
#' Calculated mutual information between each attribute joint with some other vector \code{Z} with the decision, that is
#' \deqn{I(XZ,Y).}
#' This is the same as conditional mutual information between X and Y plus a constant that depends on Y and Z, i.e.:
#' \deqn{I(XZ,Y)=I(X;Y|Z)+I(Y;Z).}
#' @template input
#' @param Z Other vector; should be given as a factor, but other options are accepted, as for attributes. 
#' @return A numerical vector with joint mutual information scores, with names copied from \code{X}.
#' @examples
#' mi(iris[,-5],iris$Species)
jmi<-function(X,Y,Z)
 mi(X,factor(sprintf("%s%s",Y,Z)))-mi(X,Z)-mi(data.frame(Y),Z) #TODO: Replace placeholder


