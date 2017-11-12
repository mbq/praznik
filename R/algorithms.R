#' Mutual information maximisation filter
#'
#' Calculates mutual information between all attributes and the decision, then returns top k.
#' @useDynLib praznik, .registration=TRUE
#' @template generic
#' @export
MIM<-function(X,Y,k=3){
 .Call(C_MIM,X,Y,as.integer(k))->ans
 names(ans)<-c("selection","scores")
 ans$selection<-colnames(X)[ans$selection]
 return(ans)
}

#' Minimal conditional mutual information maximisation filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\min_{W\in S} I(X;Y|W),}
#' where \eqn{S} is the set of already selected attributes.
#' The method stops either when no initial attribute with positive mutual information with \eqn{Y} can be found, or after \code{k} attributes are found.
#' @note CMIM is identical to the Informative Fragments (IF) method.
#' @references "Fast Binary Feature Selection using Conditional Mutual Information Maximisation" F. Fleuret, JMLR (2004)
#' @references "Object recognition with informative features and linear classification" M. Vidal-Naquet and S. Ullman, IEEE Conference on Computer Vision and Pattern Recognition (2003).
#' @template generic
#' @export
CMIM<-function(X,Y,k=3){
 .Call(C_CMIM,X,Y,as.integer(k))->ans
 names(ans)<-c("selection","scores")
 ans$selection<-colnames(X)[ans$selection]
 return(ans)
}

#' Minimum redundancy maximal relevancy filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=I(X;Y)-\frac{1}{|S|}\sum_{W\in S} I(X;W),}
#' where \eqn{S} is the set of already selected attributes.
#' The method stops either when no initial attribute with positive mutual information with \eqn{Y} can be found, or after \code{k} attributes are found.
#' @references "Feature Selection Based on Mutual Information: Criteria of Max-Dependency, Max-Relevance, and Min-Redundancy" H. Peng et al. IEEE Pattern Analysis and Machine Intelligence (PAMI) (2005)
#' @template generic
#' @export
MRMR<-function(X,Y,k=3){
 .Call(C_MRMR,X,Y,as.integer(k))->ans
 names(ans)<-c("selection","scores")
 ans$selection<-colnames(X)[ans$selection]
 return(ans)
}

#' Joint mutual information filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\sum_{W\in S} I(X,W;Y),}
#' where \eqn{S} is the set of already selected attributes.
#' The method stops either when no initial attribute with positive mutual information with \eqn{Y} can be found, or after \code{k} attributes are found.
#'
#' @note \code{\link{DISR}} is a normalised version of JMI; \code{\link{JMIM}} and \code{\link{NJMIM}} are modifications of JMI and DISR in which minimal joint information over already selected attributes is used instead of a sum.
#' @references "Data Visualization and Feature Selection: New Algorithms for Nongaussian Data H. Yang and J. Moody, NIPS (1999)
#' @template generic
#' @export
JMI<-function(X,Y,k=3){
 .Call(C_JMI,X,Y,as.integer(k))->ans
 names(ans)<-c("selection","scores")
 ans$selection<-colnames(X)[ans$selection]
 return(ans)
}

#' Double input symmetrical relevance filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\sum_{W\in S} \frac{I(X,W;Y)}{H(X,W,Y)},}
#' where \eqn{S} is the set of already selected attributes.
#' The method stops either when no initial attribute with positive mutual information with \eqn{Y} can be found, or after \code{k} attributes are found.
#'
#' @note DISR is a normalised version of \code{\link{JMI}}; \code{\link{JMIM}} and \code{\link{NJMIM}} are modifications of JMI and DISR in which minimal joint information over already selected attributes is used instead of a sum.
#' @references "On the Use of Variable Complementarity for Feature Selection in Cancer Classification" P. Meyer and G. Bontempi, (2006)
#' @template generic
#' @export
DISR<-function(X,Y,k=3){
 .Call(C_DISR,X,Y,as.integer(k))->ans
 names(ans)<-c("selection","scores")
 ans$selection<-colnames(X)[ans$selection]
 return(ans)
}

MI<-function(X,Y){
 .Call(C_MI,X,Y,0L)->ans
 names(ans)<-colnames(X)
 return(ans)
}

#' Minimal joint mutual information maximisation filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\min_{W\in S} I(X,W;Y),}
#' where \eqn{S} is the set of already selected attributes.
#' The method stops either when no initial attribute with positive mutual information with \eqn{Y} can be found, or after \code{k} attributes are found.
#'
#' @note \code{\link{NJMIM}} is a normalised version of JMIM; \code{\link{JMI}} and \code{\link{DISR}} are modifications of JMIM and NJMIM in which a sum of joint information over already selected attributes is used instead of a minimum.
#' @template generic
#' @references "Feature selection using Joint Mutual Information Maximisation" M. Bennasar, Y. Hicks and R. Setchi, (2015)
#' @export
JMIM<-function(X,Y,k=3){
 .Call(C_JMIM,X,Y,as.integer(k))->ans
 names(ans)<-c("selection","scores")
 ans$selection<-colnames(X)[ans$selection]
 return(ans)
}

#' Minimal normalised joint mutual information maximisation filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\min_{W\in S} \frac{I(X,W;Y)}{H(X,W,Y)},}
#' where \eqn{S} is the set of already selected attributes.
#' The method stops either when no initial attribute with positive mutual information with \eqn{Y} can be found, or after \code{k} attributes are found.
#'
#' @note NJMIM is a normalised version of \code{\link{JMIM}}; \code{\link{JMI}} and \code{\link{DISR}} are modifications of JMIM and NJMIM in which a sum of joint information over already selected attributes is used instead of a minimum.
#' @template generic
#' @references "Feature selection using Joint Mutual Information Maximisation" M. Bennasar, Y. Hicks and R. Setchi, (2015)
#' @export
NJMIM<-function(X,Y,k=3){
 .Call(C_NJMIM,X,Y,as.integer(k))->ans
 names(ans)<-c("selection","scores")
 ans$selection<-colnames(X)[ans$selection]
 return(ans)
}

