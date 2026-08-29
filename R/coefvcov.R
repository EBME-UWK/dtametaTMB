#' Extract model coefficients
#'
#' Returns the estimated model parameters.
#'
#' @param object A fitted model object.
#' @param ... Not currently used.
#'
#' @return
#' A named vector of parameter estimates.
#' 
#' @note For \code{ReitsmaSubgroup} models, \code{coef()} returns
#' parameter estimates from the cell-means parameterization rather
#' than treatment-contrast coefficients.
#'
#' @seealso
#' [vcov.dtametaTMB()]
#' @name coef.dtametaTMB


#' @rdname coef.dtametaTMB
#' @export
coef.DTAmodel <- function(object, ...) {
  as.matrix(object$sdreport2)[,"Estimate"]
}

#' @rdname coef.dtametaTMB
#' @export
coef.Reitsma <- function(object, ...){
  as.matrix(object$estimates)[,"Estimate"]
}

#' @rdname coef.dtametaTMB
#' @export
coef.ReitsmaSubgroup <- function(object, ...){
  as.matrix(object$estimates_mu)[,"Estimate"]
}

#' @rdname coef.dtametaTMB
#' @export
coef.ReitsmaLCA <- coef.DTAmodel

#' @rdname coef.dtametaTMB
#' @export
coef.ReitsmaSubgroupLCA <- coef.DTAmodel

#' Variance-covariance matrix
#'
#' Returns the variance-covariance matrix of the estimated model
#' parameters.
#'
#' @param object A fitted model object.
#' @param ... Not currently used.
#'
#' @return
#' A variance-covariance matrix corresponding to the parameters returned by 
#' \code{coef()}.
#'
#' @note For \code{ReitsmaSubgroup} models, the returned
#' variance-covariance matrix corresponds to the cell-means
#' parameterization returned by \code{coef()}.
#'
#' @seealso
#' [coef.dtametaTMB()]
#' @name vcov.dtametaTMB

#' @rdname vcov.dtametaTMB
#' @export
vcov.DTAmodel <- function(object, ...) {
  vc <- object$sdreport$cov
  colnames(vc) <- rownames(vc) <- rownames(object$sdreport2)
  vc
}

#' @rdname vcov.dtametaTMB
#' @export
vcov.Reitsma <- function(object, ...) {
  object$vcov
}

#' @rdname vcov.dtametaTMB
#' @export
vcov.ReitsmaSubgroup <- function(object, ...) {
  object$vcov_mu
}

#' @rdname vcov.dtametaTMB
#' @export
vcov.ReitsmaLCA         <- vcov.DTAmodel

#' @rdname vcov.dtametaTMB
#' @export
vcov.ReitsmaSubgroupLCA <- vcov.DTAmodel

