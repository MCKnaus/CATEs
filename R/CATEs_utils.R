#' This function creates the nuisance parameters p(x), mu(x), and mu_d(x)
#' via cross-fitting using the \code{\link{glmnet}} package
#'
#' @param y Vector of outcome values
#' @param d Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix)
#' @param index List indicating indices for cross-fitting (e.g. obtained by \code{createFolds} of \code{\link{caret}} pkg)
#' @param args_p List of arguments passed to estimate propensity score model
#' @param args_y List of arguments passed to estimate outcome model
#' @param args_y1 List of arguments passed to estimate outcome model of treated
#' @param args_y0 List of arguments passed to estimate outcome model of non-treated
#' @import glmnet
#'
#' @return Returns n x 4 matrix containing the nuisance parameters
#'
#' @export

nuisance_cf_glmnet <- function(y,d,x,index,
                               args_p=list(),
                               args_y=list(),
                               args_y0=list(),
                               args_y1=list()) {


  np = matrix(NA,length(d),4)
  colnames(np) = c("p_hat","y_hat","y0_hat","y1_hat")

  for(i in 1:length(index)) {
    # P-score
    fit_p = do.call(cv.glmnet,c(list(x=x[-index[[i]],,drop=F],
                                     y=d[-index[[i]]],
                                     family="binomial"),
                                args_p))
    np[index[[i]],1] = predict(fit_p,x[index[[i]],,drop=F], s = "lambda.min", type = "response")

    # Outcome
    fit_y = do.call(cv.glmnet,c(list(x=x[-index[[i]],,drop=F],
                                     y=y[-index[[i]]]),
                                args_y))
    np[index[[i]],2] = predict(fit_y,x[index[[i]],,drop=F], s = "lambda.min")

    # Outcome of non-treated
    fit_y0 = do.call(cv.glmnet,c(list(x=x[-index[[i]],,drop=F][d[-index[[i]]] == 0,,drop=F],
                                      y=y[-index[[i]]][d[-index[[i]]] == 0]),
                                 args_y0))
    np[index[[i]],3] = predict(fit_y0,x[index[[i]],,drop=F], s = "lambda.min")

    # Outcome of non-treated
    fit_y1 = do.call(cv.glmnet,c(list(x=x[-index[[i]],,drop=F][d[-index[[i]]] == 1,,drop=F],
                                      y=y[-index[[i]]][d[-index[[i]]] == 1]),
                                 args_y1))
    np[index[[i]],4] = predict(fit_y1,x[index[[i]],,drop=F], s = "lambda.min")
  }
  return(np)
}



#' This function creates the nuisance parameters p(x), mu(x), and mu_d(x)
#' via cross-fitting using the \code{\link{grf}} package
#'
#' @param y Vector of outcome values
#' @param d Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix)
#' @param index List indicating indices for cross-fitting (e.g. obtained by \code{createFolds} of \code{\link{caret}} pkg)
#' @param args_p List of arguments passed to estimate propensity score model
#' @param args_y List of arguments passed to estimate outcome model
#' @param args_y1 List of arguments passed to estimate outcome model of treated
#' @param args_y0 List of arguments passed to estimate outcome model of non-treated
#' @import grf
#'
#' @return Returns n x 4 matrix containing the nuisance parameters
#'
#' @export


nuisance_cf_grf <- function(y,d,x,index,
                            args_p=list(),
                            args_y=list(),
                            args_y0=list(),
                            args_y1=list()) {

  np = matrix(NA,length(d),4)
  colnames(np) = c("p_hat","y_hat","y0_hat","y1_hat")

  for(i in 1:length(index)) {

    fit_p = do.call(regression_forest,c(list(X=x[-index[[i]],,drop=F],
                                             Y=d[-index[[i]]]),
                                        tune.parameters = TRUE,
                                        args_p))
    np[index[[i]],1] = predict(fit_p,x[index[[i]],,drop=F])$prediction

    fit_y = do.call(regression_forest,c(list(X=x[-index[[i]],,drop=F],
                                             Y=y[-index[[i]]]),
                                        tune.parameters = TRUE,
                                        args_y))
    np[index[[i]],2] = predict(fit_y,x[index[[i]],,drop=F])$prediction

    fit_y0 = do.call(regression_forest,c(list(X=x[-index[[i]],,drop=F][d[-index[[i]]] == 0,,drop=F],
                                              Y=y[-index[[i]]][d[-index[[i]]] == 0]),
                                         tune.parameters = TRUE,
                                         args_y0))
    np[index[[i]],3] = predict(fit_y0,x[index[[i]],,drop=F])$prediction

    fit_y1 = do.call(regression_forest,c(list(X=x[-index[[i]],,drop=F][d[-index[[i]]] == 1,,drop=F],
                                              Y=y[-index[[i]]][d[-index[[i]]] == 1]),
                                         tune.parameters = TRUE,
                                         args_y1))
    np[index[[i]],4] = predict(fit_y1,x[index[[i]],,drop=F])$prediction
    # Think about predicting also for the other treatment category in the other fold and take the average
  }

  return(np)
}


#' Implementation of MOM IPW using the \code{\link{glmnet}} package
#'
#' @param y Vector of outcome values
#' @param d Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix) for estimation
#' @param np Matrix of nuisance parameters obtained by by \code{nuisance_cf_glmnet} or \code{nuisance_cf_grf}
#' @param xnew Matrix of covariates (N x p matrix) for out-of-sample prediction
#' @param args_tau List of arguments passed to estimate IATEs
#' @import glmnet
#'
#' @return Returns vector containing the out-of-sample IATEs
#'
#' @export

mom_ipw_glmnet = function(y,d,x,np,xnew,args_tau=list()) {
  mo = y * (d-np[,"p_hat"]) / (np[,"p_hat"]*(1-np[,"p_hat"]))
  fit_tau = do.call(cv.glmnet,c(list(x=x,y=mo),args_tau))
  iate = predict(fit_tau,xnew, s = "lambda.min")
  return(iate)
}


#' Implementation of MOM DR using the \code{\link{glmnet}} package
#'
#' @param y Vector of outcome values
#' @param d Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix) for estimation
#' @param np Matrix of nuisance parameters obtained by by \code{nuisance_cf_glmnet} or \code{nuisance_cf_grf}
#' @param xnew Matrix of covariates (N x p matrix) for out-of-sample prediction
#' @param args_tau List of arguments passed to estimate IATEs
#' @import glmnet
#'
#' @return Returns vector containing the out-of-sample IATEs
#'
#' @export

mom_dr_glmnet = function(y,d,x,np,xnew,args_tau=list()) {
  mo = np[,"y1_hat"] - np[,"y0_hat"] + d * (y-np[,"y1_hat"]) / np[,"p_hat"] - (1-d) * (y-np[,"y0_hat"]) / (1-np[,"p_hat"])
  fit_tau = do.call(cv.glmnet,c(list(x=x,y=mo),args_tau))
  iate = predict(fit_tau,xnew, s = "lambda.min")
  return(iate)
}


#' Implementation of MCM using the \code{\link{glmnet}} package
#'
#' @param y Vector of outcome values
#' @param d Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix) for estimation
#' @param np Matrix of nuisance parameters obtained by by \code{nuisance_cf_glmnet} or \code{nuisance_cf_grf}
#' @param xnew Matrix of covariates (N x p matrix) for out-of-sample prediction
#' @param args_tau List of arguments passed to estimate IATEs
#' @import glmnet
#'
#' @return Returns vector containing the out-of-sample IATEs
#'
#' @export

mcm_glmnet = function(y,d,x,np,xnew,args_tau=list()) {
  mo = 2 * y * (2*d - 1)
  w =  (2*d - 1) * (d - np[,"p_hat"]) / (4 * np[,"p_hat"] * (1 - np[,"p_hat"]))
  fit_tau = do.call(cv.glmnet,c(list(x=x,y=mo,weights=w),args_tau))
  iate = predict(fit_tau,xnew, s = "lambda.min")
  return(iate)
}

#' Implementation of MCM with efficiency augmentation using the \code{\link{glmnet}} package
#'
#' @param y Vector of outcome values
#' @param d Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix) for estimation
#' @param np Matrix of nuisance parameters obtained by by \code{nuisance_cf_glmnet} or \code{nuisance_cf_grf}
#' @param xnew Matrix of covariates (N x p matrix) for out-of-sample prediction
#' @param args_tau List of arguments passed to estimate IATEs
#' @import glmnet
#'
#' @return Returns vector containing the out-of-sample IATEs
#'
#' @export

mcm_ea_glmnet = function(y,d,x,np,xnew,args_tau=list()) {
  mo = 2 * (y - np[,"y_hat"]) * (2*d - 1)
  w =  (2*d - 1) * (d - np[,"p_hat"]) / (4 * np[,"p_hat"] * (1 - np[,"p_hat"]))
  fit_tau = do.call(cv.glmnet,c(list(x=x,y=mo,weights=w),args_tau))
  iate = predict(fit_tau,xnew, s = "lambda.min")
  return(iate)
}


#' Implementation of R-learning using the \code{\link{glmnet}} package
#'
#' @param y Vector of outcome values
#' @param d Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix) for estimation
#' @param np Matrix of nuisance parameters obtained by by \code{nuisance_cf_glmnet} or \code{nuisance_cf_grf}
#' @param xnew Matrix of covariates (N x p matrix) for out-of-sample prediction
#' @param args_tau List of arguments passed to estimate IATEs
#' @import glmnet
#'
#' @return Returns vector containing the out-of-sample IATEs
#'
#' @export

rl_glmnet = function(y,d,x,np,xnew,args_tau=list()) {
  mo = (y - np[,"y_hat"]) / (d - np[,"p_hat"])
  w = (d - np[,"p_hat"])^2
  fit_tau = do.call(cv.glmnet,c(list(x=x,y=mo,weights=w),args_tau))
  iate = predict(fit_tau,xnew, s = "lambda.min")
  return(iate)
}


#' Implementation of MOM IPW using the \code{\link{grf}} package
#'
#' @param y Vector of outcome values
#' @param d Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix) for estimation
#' @param np Matrix of nuisance parameters obtained by by \code{nuisance_cf_glmnet} or \code{nuisance_cf_grf}
#' @param xnew Matrix of covariates (N x p matrix) for out-of-sample prediction
#' @param args_tau List of arguments passed to estimate IATEs
#' @import grf
#'
#' @return Returns vector containing the out-of-sample IATEs
#'
#' @export

mom_ipw_grf = function(y,d,x,np,xnew,args_tau=list()) {
  mo = y * (d-np[,"p_hat"]) / (np[,"p_hat"]*(1-np[,"p_hat"]))
  fit_tau = do.call(regression_forest,c(list(X=x,Y=mo),tune.parameters = TRUE,args_tau))
  iate = predict(fit_tau,xnew)$prediction
  return(iate)
}


#' Implementation of MOM DR using the \code{\link{grf}} package
#'
#' @param y Vector of outcome values
#' @param d Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix) for estimation
#' @param np Matrix of nuisance parameters obtained by by \code{nuisance_cf_glmnet} or \code{nuisance_cf_grf}
#' @param xnew Matrix of covariates (N x p matrix) for out-of-sample prediction
#' @param args_tau List of arguments passed to estimate IATEs
#' @import grf
#'
#' @return Returns vector containing the out-of-sample IATEs
#'
#' @export

mom_dr_grf = function(y,d,x,np,xnew,args_tau=list()) {
  mo = np[,"y1_hat"] - np[,"y0_hat"] + d * (y-np[,"y1_hat"]) / np[,"p_hat"] - (1-d) * (y-np[,"y0_hat"]) / (1-np[,"p_hat"])
  fit_tau = do.call(regression_forest,c(list(X=x,Y=mo),tune.parameters = TRUE,args_tau))
  iate = predict(fit_tau,xnew)$prediction
  return(iate)
}


#' This function implements the 50:50 cross-fitting
#'
#' @param est Vector of outcome values
#' @param d Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix)
#' @param index List indicating indices for cross-fitting (e.g. obtained by \code{createFolds} of \code{\link{caret}} pkg)
#' @param args_p List of arguments passed to estimate propensity score model
#' @param args_y List of arguments passed to estimate outcome model
#' @param args_y1 List of arguments passed to estimate outcome model of treated
#' @param args_y0 List of arguments passed to estimate outcome model of non-treated
#'
#' @return Returns n x 4 matrix containing the nuisance parameters
#' 
#' @export

cf_dml1 = function(est,y,d,x,np,xnew,index,args_tau=list()) {

  iate = matrix(0,length(d),1)

  for (i in 1:length(index)) {
    iate = iate + 1/length(index) *
                  do.call(est,list(y[index[[i]]],
                          d[index[[i]]],
                          x[index[[i]],,drop=F],
                          np[index[[i]],],
                          xnew,
                          args_tau=args_tau))
  }
  return(iate)
}

