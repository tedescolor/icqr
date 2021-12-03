#################################
KMfunction = function(KMsurvfit){ # from survfit KM to function KM
  return(function(t){
    return(KMsurvfit$surv[ # return the $surv value
      which.max( # of the index of the max
        KMsurvfit$time[ # of the times
          KMsurvfit$time -rep(t,length(KMsurvfit$time)) <= 0  # among the times which are lower or equal then t
        ]
      )
    ])
  })
}

#given the data, returns the function data has to be minimized
createf = function(data,u,deltaName,timeName, covNames, IVNames, G){ return(
  function(beta){
    n = nrow(data)
    res = 0
    # OSS: weights can be computed only once:
    weigths = sapply(data[,timeName],FUN = G)
    # OSS: sapply(data$y,FUN = G) can return zero.
    #So we set them to NA and when we will sum (see for loop), we set na.rm = TRUE
    weigths[weigths==0] <- NA
    weigths = data[,deltaName] / weigths
    for(j in 1:n){
      #for optimization, we compute the indicator function for wi<=wj just once, and save it in the following wsj variable:
      wsj = as.numeric(apply( data.matrix(data[,IVNames]) <= data.matrix(data[rep(j, nrow(data)),IVNames]),1,all))
      # compute the statistic for the value wj
      res = res + ( (sum(
        as.numeric( data[,timeName]  <= exp(rowSums(matrix(beta, nrow = 1)[rep(1,n),] * data.matrix(data[,covNames])) )  ) *
          wsj * weigths, na.rm = TRUE) - u * sum(wsj) )/n )^2
    }
    return(res/n)
  })
}
#' icqr
#'
#'This function returns the coefficient vector beta 
#'in the model T = exp(z beta(u)) instrumental censored quantile regression. 
#'
#' @param data dataset as dataframe
#' @param attempts  attempts to estimate beta, to avoid local minima
#' @param method  estimation method, see optim function in stats package.
#' @param lower 1 dimension: (lower,upper) = interval where beta is searched (see optim function in stats package.),more dimensions: interval for uniform random starting point for beta, in each component
#' @param upper see "lower" parameter
#' @param deltaName variable name of the censor variable: 0 censored, 1 not censored. It is passed as list of string 
#' @param timeName variable name of the time variable.  It is passed as list of string 
#' @param covNames variable names of the covriate. It is passed as list of string 
#' @param IVNames variable names for the Instrumental variable, it has to include exogenuos variable name.  It is passed as list of string 
#' @param Nbootstrap number of bootstrap resampling. NA in case of not boostrap
#' @param bootstrapEpsilon when bootstrap, reaserch the minimum in (beta-epsilon, beta + epsilon), vector interval
#' @param verbose print during estimation
#' @return estimated beta vector parameter
#' @export
icqr = function(data, #dataset
                    u, #quantile
                    attempts = 10, #attempts to estimate beta, to avoid local minima
                    method = "Brent", # estimation method (see optim function).
                    lower = -5, upper = 5,
                    # 1 dimension: (lower,upper) = interval where beta is searched (see optim)
                    # more dimensions: interval for uniform random starting point for beta, in each component
                    deltaName = c("delta"),
                    timeName = c("y"),
                    covNames = c("x"),# covariate names
                    IVNames=c("w"), # Instrumental variables names
                    Nbootstrap = NA, # number of bootstrap resampling. By default not computed
                    bootstrapEpsilon = 1, #when bootstrap, reaserch the minimum in [beta-epsilon, beta + epsilon] (vector interval)
                    verbose = FALSE # print during estimation
){
  formulaSurv = formula(paste("Surv(",timeName[1],", 1 -", deltaName[1],") ~ 1", sep = ""))
  KMC = survfit(formulaSurv,data = data) #Kaplan Maier estimator for censor variable
  G = KMfunction(KMC) # Kaplan Maier estimator for censor variable as function
  # to avoid local minima, we minimize the function multiple times (attempts)
  # and return the value which obtains the minimum
  f = createf(data = data,u=u,deltaName=deltaName,timeName=timeName, covNames=covNames, IVNames=IVNames , G=G)
  betas = matrix(NA, nrow = attempts, ncol = length(covNames))
  for(i in 1:attempts){
    start =  runif(length(covNames),min = lower, max = upper)
    if(verbose){print(paste(i, " - starting point = ",toString(start)))}
    if(length(covNames) == 1){ #differentiate, as optim requires, 1 covarate or multiple
      betas[i,] = stats::optim(par = start,
                        fn  = f,
                        method = method,
                        lower = lower,
                        upper = upper
      )$par # we suppress the warning, passing the lower and upper parameters also when dim > 1
      
    }else{
      betas[i,] = stats::optim(par = start,
                        fn  = f,
                        method = method
      )$par
    }
  }
  indexMinBeta = which.min(apply(betas, 1, f))
  beta = betas[indexMinBeta,] # beta is the value which obtain the miniumm of f, to avoid local minima
  if(is.na(Nbootstrap)){
    if(verbose){betas = cbind(betas, apply(betas, 1, f));print(betas); print(paste("min beta at ", indexMinBeta, sep = ""));}
    return(beta)
  }else{
    boot_function <- function(original_vector, resample_vector) { # create bootstrap function that is the estimate function with no Nbootstrap and resampled data
      return( icqr(data = original_vector[resample_vector,],
                       attempts = attempts,
                       u = u,
                       method = method,
                       lower = beta - bootstrapEpsilon,
                       upper = beta + bootstrapEpsilon,
                       deltaName = deltaName,
                       timeName = timeName,
                       covNames = covNames,
                       IVNames = IVNames,
                       Nbootstrap = NA))
    }
    boot_results <- as.matrix( #transform to matrix
      broom::tidy( #use broom::tidy to obtain a dataframe as output.
        boot::boot(data, boot_function, R = Nbootstrap))) # bootstrap result obtained by boot: [original, bias, std.error] (see boot)
    return( boot_results )
  }
}




