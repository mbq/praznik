# Feature selection algorithms

#' Mutual information maximisation filter
#'
#' Calculates mutual information between all attributes and the decision, then returns top k.
#' @template input
#' @template k
#' @template output-mim
#' @examples data(MadelonD)
#' MIM(MadelonD$X,MadelonD$Y,20)
#' @export
MIM<-function(X,Y,k=3,threads=0)
 .Call(C_MIM,X,Y,as.integer(k),as.integer(threads))

#' Minimal conditional mutual information maximisation filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\min_{W\in S} I(X;Y|W),}
#' where \eqn{S} is the set of already selected attributes.
#' @note CMIM is identical to the Informative Fragments (IF) method.
#' @references "Fast Binary Feature Selection using Conditional Mutual Information Maximisation" F. Fleuret, JMLR (2004)
#' @references "Object recognition with informative features and linear classification" M. Vidal-Naquet and S. Ullman, IEEE Conference on Computer Vision and Pattern Recognition (2003).
#' @template input
#' @template k
#' @template output-mim
#' @examples data(MadelonD)
#' CMIM(MadelonD$X,MadelonD$Y,20)
#' @export
CMIM<-function(X,Y,k=3,threads=0)
 .Call(C_CMIM,X,Y,as.integer(k),as.integer(threads))

#' Minimum redundancy maximal relevancy filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=I(X;Y)-\frac{1}{|S|}\sum_{W\in S} I(X;W),}
#' where \eqn{S} is the set of already selected attributes.
#' @references "Feature Selection Based on Mutual Information: Criteria of Max-Dependency, Max-Relevance, and Min-Redundancy" H. Peng et al. IEEE Pattern Analysis and Machine Intelligence (PAMI) (2005)
#' @template input
#' @template k
#' @template output
#' @examples data(MadelonD)
#' MRMR(MadelonD$X,MadelonD$Y,20)
#' @export
MRMR<-function(X,Y,k=3,threads=0)
 .Call(C_MRMR,X,Y,as.integer(k),as.integer(threads))

#' Joint mutual information filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\sum_{W\in S} I(X,W;Y),}
#' where \eqn{S} is the set of already selected attributes.
#' @note \code{\link{DISR}} is a normalised version of JMI; \code{\link{JMIM}} and \code{\link{NJMIM}} are modifications of JMI and DISR in which minimal joint information over already selected attributes is used instead of a sum.
#' @references "Data Visualization and Feature Selection: New Algorithms for Nongaussian Data H. Yang and J. Moody, NIPS (1999)
#' @template input
#' @template k
#' @template output
#' @examples data(MadelonD)
#' JMI(MadelonD$X,MadelonD$Y,20)
#' @export
JMI<-function(X,Y,k=3,threads=0)
 .Call(C_JMI,X,Y,as.integer(k),as.integer(threads))

#' Double input symmetrical relevance filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\sum_{W\in S} \frac{I(X,W;Y)}{H(X,W,Y)},}
#' where \eqn{S} is the set of already selected attributes.
#' @note DISR is a normalised version of \code{\link{JMI}}; \code{\link{JMIM}} and \code{\link{NJMIM}} are modifications of JMI and DISR in which minimal joint information over already selected attributes is used instead of a sum.
#' @references "On the Use of Variable Complementarity for Feature Selection in Cancer Classification" P. Meyer and G. Bontempi, (2006)
#' @template input
#' @template k
#' @template output
#' @examples data(MadelonD)
#' DISR(MadelonD$X,MadelonD$Y,20)
#' @export
DISR<-function(X,Y,k=3,threads=0)
 .Call(C_DISR,X,Y,as.integer(k),as.integer(threads))

#' Minimal joint mutual information maximisation filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\min_{W\in S} I(X,W;Y),}
#' where \eqn{S} is the set of already selected attributes.
#' @note \code{\link{NJMIM}} is a normalised version of JMIM; \code{\link{JMI}} and \code{\link{DISR}} are modifications of JMIM and NJMIM in which a sum of joint information over already selected attributes is used instead of a minimum.
#' @template input
#' @template k
#' @template output-mim
#' @examples data(MadelonD)
#' JMIM(MadelonD$X,MadelonD$Y,20)
#' @references "Feature selection using Joint Mutual Information Maximisation" M. Bennasar, Y. Hicks and R. Setchi, (2015)
#' @export
JMIM<-function(X,Y,k=3,threads=0)
 .Call(C_JMIM,X,Y,as.integer(k),as.integer(threads))

#' Minimal normalised joint mutual information maximisation filter
#'
#' The method starts with an attribute of a maximal mutual information with the decision \eqn{Y}.
#' Then, it greedily adds attribute \eqn{X} with a maximal value of the following criterion:
#' \deqn{J(X)=\min_{W\in S} \frac{I(X,W;Y)}{H(X,W,Y)},}
#' where \eqn{S} is the set of already selected attributes.
#' @note NJMIM is a normalised version of \code{\link{JMIM}}; \code{\link{JMI}} and \code{\link{DISR}} are modifications of JMIM and NJMIM in which a sum of joint information over already selected attributes is used instead of a minimum.
#' It stops returning features when the best score reaches 0.
#' @template input
#' @template k
#' @template output-mim
#' @examples data(MadelonD)
#' NJMIM(MadelonD$X,MadelonD$Y,20)
#' @references "Feature selection using Joint Mutual Information Maximisation" M. Bennasar, Y. Hicks and R. Setchi, (2015)
#' @export
NJMIM<-function(X,Y,k=3,threads=0)
 .Call(C_NJMIM,X,Y,as.integer(k),as.integer(threads))

