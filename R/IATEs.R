#' This function produces the predicted values of CATE estimators on a validation set X
#'
#' @param y_t Vector of training outcome values
#' @param d_t Vector of training treament indicators
#' @param x_t Matrix of training covariates (N x p matrix)
#' @param x_v Matrix of validation covariates (N x p matrix)
#' @param tau_t (optional) Vector of training treatment effects
#' @param estimators Character vector of length e with names of estimators to fit - defaults to all supported
#' @import grf glmnet caret
#' @return Returns list with e vectors containing the IATEs of the validation sample
#' @export

IATEs <- function(y_t, d_t, x_t, x_v, tau_t = NULL,
    estimators = c("RF CMR", "RF MOM IPW", "RF MOM DR","CF","CF LC",
                   "LASSO CMR", "LASSO MOM IPW", "LASSO MOM DR",
                   "LASSO MCM", "LASSO MCM EA", "LASSO RL"),
    nFolds = 2
    ) {
  # return list with vectors - easy to coerce to matrix when necessary
  IATEs = list()
  if(!is.null(tau_t)) { # true TEs provided - simulation study
    # Infeasible RF
    rf_inf = regression_forest(x_t, tau_t)
    IATEs[['RF Inf']] = as.numeric(predict(rf_inf, x_v)$predictions)
    # Infeasible Lasso
    lasso_inf = cv.glmnet(x_t, tau_t)
    IATEs[['LASSO Inf']] = as.numeric(predict(lasso_inf, x_v))
  }

  # split for nuisance function estimation
  index = caret::createFolds(y_t, k = nFolds)

  ######################################################################
  ### RF Forest based methods
  ######################################################################

  if ("RF CMR" %in% estimators){
    # CMR Random Forest
    rf0 = regression_forest(x_t[d_t==0,], y_t[d_t==0])
    rf1 = regression_forest(x_t[d_t==1,], y_t[d_t==1])
    IATEs[['RF CMR']] = as.numeric(predict(rf1, x_v)$predictions - predict(rf0, x_v)$predictions)
  }
  if ("CF" %in% estimators){
    # Causal Forest
    cf = causal_forest(x_t, y_t, d_t, Y.hat = rep(0,n), W.hat = rep(0,n))
    IATEs[["CF"]] = as.numeric(predict(cf,x_v)$predictions)
  }

  ## modified outcome models  RF - fit nuisance fns first
  if (sum(c("RF MOM IPW", "RF MOM DR")  %in% estimators)){
        np = nuisance_cf_grf(y_t,d_t,x_t,index)
  }

  if ("RF MOM IPW" %in% estimators){
    IATEs[['RF MOM IPW']] = as.numeric(do.call(cf_dml1,
            list(mom_ipw_grf, y_t, d_t, x_t, np, x_v, index)
            ))
  }

  if ("RF MOM DR" %in% estimators){
    IATEs[["RF MOM DR"]] = as.numeric(do.call(cf_dml1,
            list(mom_dr_grf, y_t, d_t, x_t, np, x_v, index)
            ))
  }

  if ("CF LC" %in% estimators){
    # Causal Forest with local centering
    cf_lc = causal_forest(x_t, y_t, d_t, Y.hat = np[,2], W.hat = np[,1])
    IATEs[["CF LC"]] = as.numeric(predict(cf_lc,x_v)$predictions)
  }

  ######################################################################
  ### Lasso based methods
  ######################################################################
  if ("LASSO CMR" %in% estimators){
    # CMR Lasso
    lasso0 = cv.glmnet(x_t[d_t==0,], y_t[d_t==0])
    lasso1 = cv.glmnet(x_t[d_t==1,], y_t[d_t==1])
    IATEs[["LASSO CMR"]] = as.numeric(predict(lasso1, x_v) - predict(lasso0, x_v))
  }

  # need estimates of Lasso nuisance parameters for next 5 estimators
  if (sum(c("Lasso MOM IPW", "Lasso MOM DR", "Lasso MCM", "Lasso MCM EA", "Lasso RL")
        %in% estimators)){
    np = nuisance_cf_glmnet(y_t, d_t, x_t, index)
  }

  if ("LASSO MOM IPW" %in% estimators){
    IATEs[["LASSO MOM IPW"]] = as.numeric(do.call(cf_dml1,
                                list(mom_ipw_glmnet, y_t, d_t, x_t, np, x_v, index)
                              ))
  }

  if ("LASSO MOM DR" %in% estimators){
    IATEs[["LASSO MOM DR"]] = as.numeric(do.call(cf_dml1,
                                list(mom_dr_glmnet, y_t, d_t, x_t, np, x_v, index)
                              ))
  }

  if ("LASSO MCM" %in% estimators){
    IATEs[["LASSO MCM"]] = as.numeric(do.call(cf_dml1,
                                list(mcm_glmnet, y_t, d_t, x_t, np, x_v, index)
                              ))
  }

  if ("LASSO MCM EA" %in% estimators){
    IATEs[["LASSO MCM EA"]] = as.numeric(do.call(cf_dml1,
                                list(mcm_ea_glmnet, y_t, d_t, x_t, np, x_v, index)
                              ))
  }

  if ("LASSO RL" %in% estimators){
    IATEs[["LASSO RL"]] = as.numeric(do.call(cf_dml1,
                                list(rl_glmnet, y_t, d_t, x_t, np, x_v, index)
                              ))
  }

  # additional estimators go here
  return(IATEs)
}
