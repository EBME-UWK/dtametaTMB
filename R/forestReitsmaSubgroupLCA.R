#' Coupled Forest plot for latent class diagnostic test accuracy meta-analysis
#'
#' Study-specific sensitivities and specificities are empirical Bayes estimates
#' derived from the fitted latent class model. Confidence intervals are based on 
#' the conditional posterior variance of the study-specific random effects 
#' obtained from the TMB Laplace approximation.
#'
#' @param x Object of class \code{"ReitsmaSubgroupLCA"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param subgroup_label Column name for the subgroup. Defaults to \code{"Subgroup"}.
#' @param order Specifies the ordering of studies
#'   in the forest plot. Can be `"study"` (default), which orders
#'   observations by study identifier and then subgroup, or
#'   `"subgroup"`, which orders observations by subgroup and then
#'   study identifier.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest ReitsmaSubgroupLCA
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qnorm qlogis plogis
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.ReitsmaSubgroupLCA <- function(x, conflevel=0.95, subgroup_label="Subgroup", order=c("study","subgroup"), ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  order <- match.arg(order)
  ss <- getForestSensSpecLCA(x=x,conflevel=conflevel)

  XP <- ss$XP
  if(order == "study"){XP <- XP[order(XP$study, XP$subgroup), ]}
  if(order == "subgroup"){XP <- XP[order(XP$subgroup, XP$study), ]}
  dt <- XP[,c("study","subgroup","y11","y10","y01","y00","senslabel","speclabel")]
  dt$" "    <- " "
  dt$fsens  <- paste(rep(" ",18),collapse=" ")
  dt$a      <- " "
  dt$fspec  <- paste(rep(" ",18),collapse=" ")  
  colnames(dt) <- c("Study",subgroup_label,
                    "+/+","+/-","-/+","-/-",
                    ss$senslab,ss$speclab," ",ss$senslab," ",ss$speclab)
  
  getForestPlot(dt=dt,XP=XP,r=1,s=3)
}

