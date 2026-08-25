#' Extract model coefficients
#'
#' Returns the estimated model parameters
#'
#' @param object A fitted model object.
#' @param ... Not currently used.
#'
#' @return
#' A named vector of parameter estimates.
#'
#' @seealso
#' [vcov.dtametaTMB()]
#' @name coef.dtametaTMB


#' @rdname coef.dtametaTMB
#' @export
coef.RutterGatsonis <- function(object, ...) {
  as.matrix(object$sdreport2)[,"Estimate"]
}

#' @rdname coef.dtametaTMB
#' @export
coef.RutterGatsonisSubgroup <- function(object, ...) {
  as.matrix(object$sdreport2)[,"Estimate"]
}

#' @rdname coef.dtametaTMB
#' @export
coef.RutterGatsonisReg <- function(object, ...) {
  as.matrix(object$sdreport2)[,"Estimate"]
}

#' @rdname coef.dtametaTMB
#' @export
coef.HoyerAFT <- function(object, ...) {
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
coef.ReitsmaLCA <- function(object, ...) {
  as.matrix(object$sdreport2)[,"Estimate"]
}

#' @rdname coef.dtametaTMB
#' @export
coef.RutterGatsonisLCA <- function(object, ...) {
  as.matrix(object$sdreport2)[,"Estimate"]
}

#' @rdname coef.dtametaTMB
#' @export
coef.ReitsmaSubgroupLCA <- function(object, ...) {
  as.matrix(object$sdreport2)[,"Estimate"]
}

#' @rdname coef.dtametaTMB
#' @export
coef.RutterGatsonisSubgroupLCA <- function(object, ...) {
  as.matrix(object$sdreport2)[,"Estimate"]
}


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
#' @seealso
#' [coef.dtametaTMB()]
#' @name vcov.dtametaTMB


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
vcov.RutterGatsonis <- function(object, ...) {
  vc <- object$sdreport$cov
  colnames(vc) <- rownames(vc) <- rownames(object$sdreport2)
  vc
}

#' @rdname vcov.dtametaTMB
#' @export
vcov.RutterGatsonisSubgroup <- function(object, ...) {
  vc <- object$sdreport$cov
  colnames(vc) <- rownames(vc) <- rownames(object$sdreport2)
  vc
}

#' @rdname vcov.dtametaTMB
#' @export
vcov.RutterGatsonisReg <- function(object, ...) {
  vc <- object$sdreport$cov
  colnames(vc) <- rownames(vc) <- rownames(object$sdreport2)
  vc
}

#' @rdname vcov.dtametaTMB
#' @export
vcov.RutterGatsonisLCA <- function(object, ...) {
  vc <- object$sdreport$cov
  colnames(vc) <- rownames(vc) <- rownames(object$sdreport2)
  vc
}

#' @rdname vcov.dtametaTMB
#' @export
vcov.RutterGatsonisSubgroupLCA <- function(object, ...) {
  vc <- object$sdreport$cov
  colnames(vc) <- rownames(vc) <- rownames(object$sdreport2)
  vc
}

#' @rdname vcov.dtametaTMB
#' @export
vcov.ReitsmaLCA <- function(object, ...) {
  vc <- object$sdreport$cov
  colnames(vc) <- rownames(vc) <- rownames(object$sdreport2)
  vc
}

#' @rdname vcov.dtametaTMB
#' @export
vcov.ReitsmaSubgroupLCA <- function(object, ...) {
  vc <- object$sdreport$cov
  colnames(vc) <- rownames(vc) <- rownames(object$sdreport2)
  vc
}


#' @rdname vcov.dtametaTMB
#' @export
vcov.HoyerAFT <- function(object, ...) {
  vc <- object$sdreport$cov
  colnames(vc) <- rownames(vc) <- rownames(object$sdreport2)
  vc
}

