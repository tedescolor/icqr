#' upperQuantileIdentification
#'
#'This function returns the a suggested upper bound for the quantile that can be identified 
#'depennding on the data.
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
#' @param us grid of quantiles that among which the suggested quantile is picked
#' @return estimated beta vector parameter
#' @export

upperQuantileIdentification = function(data, #dataset
                                       attempts = 1, #attempts to estimate beta, to avoid local minima
                                       method = "Brent", # estimation method (see optim function).
                                       lower = -5, upper = 5,
                                       # 1 dimension: (lower,upper) = interval where beta is searched (see optim)
                                       # more dimensions: interval for uniform random starting point for beta, in each component
                                       deltaName = c("delta"),
                                       timeName = c("y"),
                                       covNames = c("x"),# covariate names
                                       IVNames=c("w"), # Instrumental variables names
                                       us = seq(0.1,0.9, by= 0.1) #grid of quantile to consider 
){
  n = nrow(data) #size data
  censorIndex = data[,deltaName] == 0 #index of censor times
  namesResult = c("u", "#timesAboveSupportC", "%timesAboveSupportC", "suggestedLimit")
  result = matrix(NA,nrow = 0,ncol = length(namesResult) )
  for(u in us){
    print(u)
    beta = estimate(data = data,u = u,attempts = attempts,method = method,
                    lower = lower,upper =upper,
                    deltaName = deltaName,timeName = timeName,covNames = covNames,IVNames = IVNames) #estimate beta, qith paramenters given
    
    # check how many censor data times are less  that the estimated quantile value i.e. (censoreTime < exp(Z beta(u))) 
    numTimesAboveSupportC = sum(data[censorIndex, timeName] <  exp(rowSums(matrix(beta, nrow = 1)[rep(1,sum(censorIndex)),] * data.matrix(data[censorIndex,covNames])) ))
    
    result = rbind(result, c(u,numTimesAboveSupportC,numTimesAboveSupportC/n, " ") )  
    
    
  }
  result = as.data.frame(result)
  names(result) = namesResult
  
  #suggestedLimit is the biiger u for which timesAboveSupportC/n < 10%
  result[sum(result[,c("%timesAboveSupportC")] < 0.10), c("suggestedLimit")] = "*"
  return(result)
}
