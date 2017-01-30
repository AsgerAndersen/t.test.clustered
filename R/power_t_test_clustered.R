#' Power of the adjusted t test for cluster randomized designs
#' 
#' Calculates either the power of the cluster adjusted t test or the value 
#' of some other parameter needed to obtain a specified power of the test. 
#' The calculated power is only exact, if the design is balanced (same numbers
#' of clusters in each group, and same number of responses in each cluster). 
#' See XXXX for a description of the assumptions and construction of the 
#' cluster adjusted t test and how to calculate its power.
#' 
#' @param cluster_group_1 a vector that contains the sizes of the clusters in 
#' group 1. If group 1 for instance has one cluster with 10 measurements, and 
#' three clusters with 12 measurements each, then we have to set 
#' \code{clusters_group_1} to \code{c(10,12,12,12)}.
#' @param cluster_group_2 a vector that contains the sizes of the clusters in 
#' group 2.
#' @param delta the difference between the group means (the treatment effect).
#' @param sd the overall standard deviation of the responses, which is assumed 
#' to be equal in both groups.
#' @param rho the intracluster correlation coeffiecient, which is assumed to 
#' be equal in both groups.
#' @param power the power of the cluster adjusted t test.
#' @param alternative has to be set to one of two strings, either 
#' \code{"one.sided"} or \code{"two.sided"}, to say whether the alternative 
#' hypothesis is one-sided or two-sided.
#' @param sig_level the significance level of the t test. If it is not 
#' specified, it has default value \code{0.05}.
#' @param sample_size_calc has to be set to \code{TRUE} or \code{FALSE} to say 
#' whether the sample size needed to achieve a certain power should be 
#' calculated. This calculation assumes a balanced design. By default, 
#' \code{sample_size_calc} is \code{FALSE}.
#' @param clusters_pr_group the number of clusters in each group. It should 
#' only be specified, if \code{sample_size_calc} is \code{TRUE}, and 
#' \code{cluster_size} is not specified.
#' @param cluster_size the number of measurements in each cluster. It should 
#' only be specified, if \code{sample_size_calc} is \code{TRUE}, and 
#' \code{clusters_pr_group} is not specified.
#' 
#' @return A list - If \code{sample_size_calc} is set to \code{FALSE}, then one the 
#' parameters \code{delta, sd, rho, power} must be left unspecified. The 
#' unspecified parameter will then be calculated by the function. If 
#' \code{sample_size_calc} is set to \code{TRUE}, then all of the parameters 
#' \code{delta, sd, rho, power} must be specified, but \code{clusters_group_1} 
#' and \code{clusters_group_2} must not be specified. On top of that, either 
#' \code{clusters_pr_group} or \code{cluster_size} must be specified. The 
#' other one will then be calculated by the function. Whether or not 
#' \code{sample_size_calc} is \code{TRUE} or \code{FALSE}, the function returns 
#' a list of all of the arguments to the function and the calculated parameter 
#' (the argument left unspecified when calling the function)
#' 
#' @example examples/power.t.test.clustered.ex1.R
#' 
#' @seealso \code{\link{t.test.clustered.pval}} for p value of the cluster 
#' adjusted t test, \code{\link{t.test.clustered.stat}} for value of the 
#' cluster adjusted t statistic and \code{\link{simulate.power.t.test.clustered}}
#' for power simulation of the cluster adjusted t test.
#' 
#' @export 
power.t.test.clustered <- function (clusters_group_1=NULL, clusters_group_2=clusters_group_1, delta=NULL, sd = NULL, rho = NULL, 
                              power = NULL, alternative=c("one.sided","two.sided"), sig_level = 0.05,
                              sample_size_calc=FALSE, clusters_pr_group=NULL, cluster_size=NULL) {
  
  #------------------------------------------------------------------------------------------
  #checking if the arguments are valid
  if (sample_size_calc) {
    if (sum(sapply(list(clusters_pr_group,cluster_size), is.null)) != 1) {
      stop("to calculate sample-size, you need to specify the number of clusters per group 
           or the number of responses in each cluster (sample-size-calculation assumes a 
           balanced design)")
    }
  }
  else if (sum(sapply(list(delta, sd, rho, power, sig_level), is.null)) != 1) {
    stop("exactly one of 'delta', 'sd', 'power' or 'sig_level' must be NULL, 
         or sample_size_cal must be TRUE")
  }
  if (!is.null(clusters_group_1)) {
    if (length(clusters_group_1)<2) {
      stop("clusters_group_1 must be a vector with at least two arguments")
    }
  }
  if (!is.null(clusters_group_2)) {
    if (length(clusters_group_2)<2) {
      stop("clusters_group_2 must be a vector with at least two arguments")
    }
  }
  #------------------------------------------------------------------------------------------
  
  #------------------------------------------------------------------------------------------
  #Configuring according to whether it is one- or two-sided
  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 2 && !(is.null(delta))) {delta <- abs(delta)}
  #------------------------------------------------------------------------------------------
  
  #------------------------------------------------------------------------------------------
  #Finding the mean cluster-size og the mean number of clusters per group
  #If the design is balanced, this will just be the exact cluster-size and the exact number of clusters per group
  if (!sample_size_calc) {
    cluster_size <- mean(c(clusters_group_1,clusters_group_2))
    clusters_pr_group <- mean(c(length(clusters_group_1) , length(clusters_group_2)))
  }
  #------------------------------------------------------------------------------------------
  
  #------------------------------------------------------------------------------------------
  #The power of the test
  p.body <- quote({
    
    this_ncp <- ncp(total_obs = 2*(clusters_pr_group*cluster_size) ,cluster_size=cluster_size,
                      delta=delta,sd=sd,rho=rho) #Non-central parameter
    this_df <- 2*(clusters_pr_group-1) #Degrees of freedom
    this_q <- qt(sig_level/tside, df=this_df, lower.tail = FALSE) #Critical value under the null-hypothesis
  
    if (tside==1) {
    pt(q=this_q, df=this_df, ncp = this_ncp, lower.tail = FALSE) #Power, one-sided test
    }
    else if (tside==2) {
      (pt(q=this_q, df=this_df, ncp = this_ncp, lower.tail = FALSE)  +
        pt(q=-this_q, df=this_df, ncp = this_ncp, lower.tail = TRUE)) #Power, two-sided test
    }
  })
  #------------------------------------------------------------------------------------------
  
  #------------------------------------------------------------------------------------------
  #Evaluating the above expression for the test power according to which parameter should be determined 
  tol = .Machine$double.eps^0.25
  
  if (is.null(power)) {power <- eval(p.body)}
  
  else if (is.null(sd)) {
    sd <- uniroot(function(sd) eval(p.body) - power, delta * 
                    c(1e-07, 1e+07), tol = tol, extendInt = "downX")$root
  }
  
  else if (is.null(delta)) {
    delta <- uniroot(function(delta) eval(p.body) - power, 
                     sd * c(1e-07, 1e+07), tol = tol, extendInt = "upX")$root
  }
  
  else if (is.null(sig_level)) {
    sig_level <- uniroot(function(sig_level) eval(p.body) - 
                           power, c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root
  }
  
  else if (is.null(rho)) {
    rho <- uniroot(function(rho) eval(p.body) - 
                     power, c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root
  }
  
  else if (sample_size_calc && is.null(clusters_pr_group)) {
    clusters_pr_group <- uniroot(function(clusters_pr_group) eval(p.body) - 
                     power, c(2, 1e+10), tol = tol, extendInt = "upX")$root
  }
  
  else if (sample_size_calc && is.null(cluster_size)) {
    cluster_size <- uniroot(function(cluster_size) eval(p.body) - 
                     power, c(1, 1e+10), tol = tol, extendInt = "upX")$root
  }
  #------------------------------------------------------------------------------------------
  
  #------------------------------------------------------------------------------------------
  #Putting the arguments of the function and the calculated parameter in a structure, 
  #which is returned by the function
  METHOD <- paste(switch(tside, one.sample = "One-sided", two.sample = "Two-sided"), 
                  "t test power calculation")
  
  if (sample_size_calc) {
    structure(list(clusters_pr_group=clusters_pr_group, 
                   cluster_size = cluster_size,
                   delta = delta, sd = sd, sig_level = sig_level,
                   power = power, alternative = alternative, rho = rho,
                   method = METHOD), class = "power.htest")
  }
  else {
    structure(list(clusters_group_1=clusters_group_1, 
                 clusters_group_2 = clusters_group_2,
                 delta = delta, sd = sd, sig_level = sig_level,
                 power = power, alternative = alternative, rho = rho,
                 method = METHOD), class = "power.htest")
  }
  #------------------------------------------------------------------------------------------
}
#***********************************************************************************************

#***********************************************************************************************
#Calculates the non centrality parameter of the t-distribution of the t-statistic under
#alternative hypothesis
ncp <- function(total_obs,cluster_size,delta,sd,rho) {
  a <- 4*(1+(cluster_size-1)*rho)
  b <- sqrt(a/total_obs)
  d <- abs(delta)/sd
  return (d/b)
}
#***********************************************************************************************


