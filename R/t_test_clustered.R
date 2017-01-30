#' t statistic for cluster randomized designs
#' 
#' Calculates the value of a t statistic that has been adjusted for
#' a cluster randomized design. See XXXX for a description of the assumptions 
#' and construction of this t statistic.
#' 
#' @param data_experiment a dataframe with three columns: group, cluster and 
#' response. The response column holds the values of the measured responses. The 
#' group column assigns a group number (must be either 1 or 2) to each 
#' response. The cluster column assigns a cluster number to each response. 
#' Each different cluster must have a different number, even if they are in 
#' different groups. See the examples section for an example of a valid
#' \code{data_experiment}.
#' 
#' @return a number - the value of the cluster adjusted t statistic.
#' 
#' @example examples/validDataframe.R
#' @example examples/tTestClusteredStatEx1.R
#' 
#' @seealso \code{\link{t.test.clustered.pval}} for p value of the cluster 
#' adjusted t test, \code{\link{power.t.test.clustered}} for power 
#' calculation of the cluster adjusted t test and 
#' \code{\link{simulate.power.t.test.clustered}} for power simulation of 
#' the cluster adjusted t test.
#' 
#' @export

t.test.clustered.stat <- function(data_experiment) {
  group_1_mean <- mean(subset(data_experiment,group==1)$response)
  group_2_mean <- mean(subset(data_experiment,group==2)$response)
  this_SE_est <- SE_est(data_experiment)
  return ((group_2_mean-group_1_mean)/this_SE_est)
}

#' p value of the t test for cluster randomized designs
#' 
#' Calculates the p value of a t test that has been adjusted to work in a 
#' cluster randomized design. The t test uses the cluster adjusted t 
#' statistic, which has a t distribution with the total number of clusters 
#' minus 2 degrees of freedom. See XXXX for a description of the assumptions 
#' and construction of the cluster adjusted t statistic and test.
#' 
#' @inheritParams t.test.clustered.stat
#' @inheritParams power.t.test.clustered
#' 
#' @return a number - the p value of the t test for cluster randomized designs
#' 
#' @example examples/validDataframe.R
#' @example examples/tTestClusteredPvalEx1.R
#' 
#' @seealso \code{\link{t.test.clustered.stat}} for value of the cluster adjusted 
#' t statistic, \code{\link{power.t.test.clustered}} for power calculation of 
#' cluster adjusted t test and \code{\link{simulate.power.t.test.clustered}} 
#' for power simulation of the cluster adjusted t test.
#'  
#' @export  

t.test.clustered.pval <- function(data_experiment, alternative=c("one.sided","two.sided")) {
  
  num_clusters <- nlevels(data_experiment$cluster)
  t <- t.test.clustered.stat(data_experiment)
  if (alternative == "one.sided") {
    p <- pt(t,df=num_clusters-2,lower.tail = FALSE)
  }
  if (alternative == "two.sided") {
    t <- abs(t)
    p <- (pt(-t,df=num_clusters-2,lower.tail = TRUE) + 
            pt(t,df=num_clusters-2,lower.tail = FALSE))
  }
  return (p)
}

t.test.clustered <- function(data_experiment, alt, sig_level=0.05) {
  p <- t.test.clustered.pval(data_experiment,alt)
  if (p<sig_level) {return (0)}
  else {return (1)}
}
#**************************************************************************

#**************************************************************************
#Calculating the standard deviation of the clustered t-test

m_Ai <- function(data_experiment_group_i) {
  clus_sizes <- as.numeric(table(data_experiment_group_i$cluster))
  group_size <- length(data_experiment_group_i$response)
  return ((sum(clus_sizes^2))/group_size)
}

var_IF_group_i_est <- function(data_experiment_group_i,rho_est) {
  this_m_Ai <- m_Ai(data_experiment_group_i)
  return (1 + (this_m_Ai - 1)*rho_est)
}

pooled_sd <- function(data_experiment) {
  group_1_response <- subset(data_experiment,group==1)$response
  group_2_response <- subset(data_experiment,group==2)$response
  size_group_1 <- length(group_1_response)
  size_group_2 <- length(group_2_response)
  pooled_sd <- ((size_group_1-1)*sd(group_1_response)+(size_group_2-1)*sd(group_2_response))/(size_group_1+size_group_2-2)
  return (pooled_sd)
}

SE_est_fact <- function(data_experiment) {
  
  data_group_1 <- subset(data_experiment,group==1)
  data_group_2 <- subset(data_experiment,group==2)
  size_group_1 <- length(data_group_1$response)
  size_group_2 <- length(data_group_2$response)
  
  rho_est <- rho_est(data_experiment)
  
  var_IF_group_1 <- var_IF_group_i_est(data_group_1,rho_est)
  var_IF_group_2 <- var_IF_group_i_est(data_group_2,rho_est)
  
  a <- var_IF_group_1/size_group_1
  b <- var_IF_group_2/size_group_2
  
  return (sqrt(a+b))
  
}

SE_est <- function(data_experiment) {
  S_p <- pooled_sd(data_experiment)
  SE_fact <- SE_est_fact(data_experiment)
  return (S_p*SE_fact) 
}
#***************************************************************************

#***************************************************************************
#Using MSW (mean squared error within clusters) and MSC (mean squared error among clusters)
#to calculate an estimate of rho

#---------------------------------------------------------------------------
#Calculation of MSW
SS_within_from_one_clus <- function(cluster_response) {
  SS <- sum((cluster_response - mean(cluster_response))^2)
  return (SS)
}

SS_within <- function(data_experiment) {
  ss <- 0
  for (i in 1:nlevels(data_experiment$cluster)) {
    this_clus_response <- (subset(data_experiment,cluster==i))$response
    ss <- ss + SS_within_from_one_clus(this_clus_response)
  }
  return (ss)
}

MSW <- function(SS_within,num_obs_total,num_clus_total) {
  return (SS_within/(num_obs_total-num_clus_total))
}
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Calculation of MSC
SS_among_from_one_clus <- function(cluster_response,group_mean) {
  num_obs_cluster <- length(cluster_response)
  clus_mean <- mean(cluster_response)
  sqr_err <- (group_mean - clus_mean)^2
  return (num_obs_cluster*sqr_err)
}

SS_among <- function(data_experiment) {
  ss <- 0
  data_group_1 <- subset(data_experiment,group==1)
  mean_group_1 <- mean(data_group_1$response)
  num_clus_group_1 <- nlevels(droplevels(data_group_1$cluster))
  data_group_2 <- subset(data_experiment,group==2)
  mean_group_2 <- mean(data_group_2$response)
  first_clus_group_2 <- num_clus_group_1+1
  num_clus_total <- nlevels(data_experiment$cluster)
  for (i in 1:num_clus_group_1) {
    this_clus_response <- (subset(data_group_1,cluster==i))$response
    ss <- ss + SS_among_from_one_clus(this_clus_response,mean_group_1)
  }
  for (i in first_clus_group_2:num_clus_total) {
    this_clus_response <- (subset(data_group_2,cluster==i))$response
    ss <- ss + SS_among_from_one_clus(this_clus_response,mean_group_2)
  }
  return (ss)
}

MSC <- function(SS_among,num_clus_total) {
  return (SS_among/(num_clus_total-2))
}
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
#Estimate of rho
m_0 <- function(num_obs_total,m_A1,m_A2,num_clus_total) {
  a <- (num_obs_total - (m_A1 + m_A2))
  return (a/(num_clus_total-2))
}

rho_est <- function(data_experiment) {
  num_obs_total <- length(data_experiment$response)
  num_clus_total <- nlevels(data_experiment$cluster)
  
  SS_within <- SS_within(data_experiment)
  MSW <- MSW(SS_within,num_obs_total,num_clus_total)
  
  SS_among <- SS_among(data_experiment)
  MSC <- MSC(SS_among,num_clus_total)
  
  m_A1 <- m_Ai(subset(data_experiment,group==1))
  m_A2 <- m_Ai(subset(data_experiment,group==2))
  
  m_0 <- m_0(num_obs_total,m_A1,m_A2,num_clus_total)
  
  a <- (MSC - MSW)
  b <- (MSC + (m_0 -1)*MSW)
  
  return (a/b)
}
#--------------------------------------------------------------------------
#************************************************************************************************