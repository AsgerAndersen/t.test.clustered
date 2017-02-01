#' Simulated power of the t test for cluster randomized designs
#' 
#' The function simulates \code{num_experiments} instances of data according 
#' to the model explained in the vignette "Construction of the library 
#' functions" with parameters defined by the arguments of 
#' the function. For each one of the instances, the function uses the t test 
#' for cluster randomized designs to test the null hypothesis that the group 
#' means are equal. The function outputs the simulated power of the test, 
#' which is the proportion of test rejecting the null hypothesis out of the 
#' total number of test.
#' 
#' @param num_experiments is the number of experiments in the simulation.
#' @inheritParams power_t_test_clustered
#' 
#' @example examples/simulate.power.ex1.R
#' 
#' @return A number - simulated power of the t test for cluster randomized
#' designs 
#' 
#' @seealso \code{\link{t_test_clustered_pval}} for p value of the cluster 
#' adjusted t test, \code{\link{t_test_clustered_stat}} for value of the 
#' cluster adjusted t statistic and \code{\link{power_t_test_clustered}}
#' for power calculation of the cluster adjusted t test.
#' 
#' @export 
simulate_power_t_test_clustered <- function(num_experiments=1000, clusters_group_1, 
                                            clusters_group_2, delta, sd, rho, alternative, 
                                            sig_level=0.05) {
  accepted <- 0
  
  for (j in 1:num_experiments) {
    test_data <- data_generate_experiment(clusters_group_1,clusters_group_2,delta,sd,rho)
    accepted <- accepted + t_test_clustered(test_data,alternative,sig_level)
  }
  
  return (1 - (accepted/num_experiments))
  
}

#****************************************************************************************
#Simulating an experiment given the model described in the documentation

data_generate_group <- function(num_obs_clus, group_mean, sd_group, rho) {
  var_group <- sd_group^2
  var_among_clus <- rho*var_group #udregner variansen among clusters
  var_within_clus <- var_group - var_among_clus #udregner variansen within clusters
  clus_effect <- stats::rnorm(n=length(num_obs_clus),mean=0,sd=sqrt(var_among_clus)) #finder cluster-effekten for hvert cluster vha. among-clus-var
  data<-vector(mode='numeric',length = 0)
  for (i in 1:length(num_obs_clus)) {
    x <- ( (stats::rnorm(num_obs_clus[i], mean=group_mean, sd=sqrt(var_within_clus))) + clus_effect[i]) #Genererer normalfordelte m?linger med within-clus-var st?j og cluster-effekt
    data <- c(data,x)
  }
  return (data)
}

data_generate_experiment <- function(num_obs_clus_group_1, num_obs_clus_group_2=num_obs_clus_group_1, 
                                     delta, sd = 1, rho) {
  response_group_1 <- data_generate_group(num_obs_clus_group_1,group_mean=0,sd_group=sd,rho)
  response_group_2 <- data_generate_group(num_obs_clus_group_2,group_mean=delta,sd_group=sd,rho)
  response <- c(response_group_1,response_group_2)
  group <- as.factor(c(rep(1,sum(num_obs_clus_group_1)),rep(2,sum(num_obs_clus_group_2))))
  num_obs_clusters <- c(num_obs_clus_group_1,num_obs_clus_group_2)
  cluster <- vector(mode="numeric",length = 0)
  for (i in 1:length(num_obs_clusters)) {
    cluster <- c(cluster,rep(i,num_obs_clusters[i]))
  }
  cluster <- as.factor(cluster)
  data <- data.frame(group,cluster,response)
  return (data)
}
#***********************************************************************************************










