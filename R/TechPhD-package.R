#' TechPhD: Tests and Estimation of Covariance Change-Points for Hi-Dim Data
#'
#' An implementation of the procedures in Zhong et al. (2019) and Santo and Zhong (2020)
#' for testing the homogeneity of covariance matrices, and estimating
#' multiple change-points in high-dimensional (Hi-Dim) longitudinal/functional data with general
#' temporospatial dependence. The null hypothesis of the homogeneity test is that
#' all covariance matrices are equal at each time point. If the null hypothesis is rejected,
#' the procedure further identifies the locations of the change points.
#' Note: The package uses Open MP.  Mac OS X users may need to update clang compiler so that it supports Open MP.
#'
#' \tabular{ll}{
#'   Package: \tab TechPhD\cr
#'   Type: \tab package\cr
#'   Version: \tab 1.0.0\cr
#'   Date: \tab 2020-04-06\cr
#'   License: \tab GPL-2\cr
#'   }
#'
#' @section Functions:
#'
#' \itemize{
#'   \item test_covmat
#'   \item cpi_covmat
#' }
#'
#' @author \strong{Maintainer}: Ping-Shou Zhong \email{pszhong@uic.edu}
#'
#'   Authors:
#'     \itemize{
#'       \item Ping-Shou Zhong
#'       \item Shawn Santo
#'       \item Nurlan Abdukadyrov
#'       \item Bo Liu
#'     }
#'
#' @references
#' \emph{Zhong, Li, and Santo (2019). Homogeneity tests of covariance
#'   matrices with high-dimensional longitudinal data. Biometrika, 106, 619-634}
#'
#' \emph{Santo and Zhong (2020). Homogeneity tests of covariance and
#'    change-points identification for high-dimensional functional data. arXiv:2005.01895}
#'
#' @name TechPhD-package
#' @docType package
NULL
