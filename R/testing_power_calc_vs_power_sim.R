test_calcpower_against_simpower <- function(num_test=10, num_experiments_pr_test=100, arg_clusters_group_1=NULL,
                                         arg_clusters_group_2=arg_clusters_group_1,arg_delta=NULL,
                                         arg_sd=NULL,arg_rho=NULL, arg_alt=NULL, arg_sig_level=0.05,
                                         balanced=FALSE,kind_of_balanced=FALSE, arg_cluster_size=NULL,arg_clusters_pr_group=NULL,
                                         bounds_cluster_size=c(10,100),bounds_cluster_pr_group=c(5,30),
                                         bound_rho=1) {
  
  results <- matrix(NA,nrow = num_test,ncol = 15)
  colnames(results) <- c("pow.cal", "pow.sim", "diff", "|diff|>3", "|diff|>8", "sdclussize1", "sd.clus.size.2", 
                         "mean.clus.size.1", "mean.clus.size.2", "clus.num.1", "clus.num.2", "delta", "sd" ,"rho", "alt")
  
  lower_clus_size <- bounds_cluster_size[1]
  upper_clus_size <- bounds_cluster_size[2]
  lower_clus_pr_group <- bounds_cluster_pr_group[1]
  upper_clus_pr_group <- bounds_cluster_pr_group[2]
  
  clusters <- list()
  browser()
  for (j in 1:num_test) {
    
    if (balanced) {
      if (is.null(arg_cluster_size)) {cluster_size <- sample(lower_clus_size:upper_clus_size,1)}
      else {cluster_size <- arg_cluster_size}
      if (is.null(arg_clusters_pr_group)) {clusters_pr_group <- sample(lower_clus_pr_group:upper_clus_pr_group,1)}
      else {clusters_pr_group <- arg_clusters_pr_group}
      clusters_group_1 <- c(rep(cluster_size,clusters_pr_group))
      clusters_group_2 <- c(rep(cluster_size,clusters_pr_group))
    }
    else if (kind_of_balanced) {
      clusters_pr_group <- sample(lower_clus_pr_group:upper_clus_pr_group,1)
      clusters_group_1 <- sample(lower_clus_size:upper_clus_size,clusters_pr_group,replace=TRUE)
      clusters_group_2 <- sample(lower_clus_size:upper_clus_size,clusters_pr_group,replace=TRUE)
    }
    else {
      num_clus_group_1 <- sample(lower_clus_pr_group:upper_clus_pr_group,1)
      num_clus_group_2 <- sample(lower_clus_pr_group:upper_clus_pr_group,1)
      clusters_group_1 <- sample(lower_clus_size:upper_clus_size,num_clus_group_1,replace=TRUE)
      clusters_group_2 <- sample(lower_clus_size:upper_clus_size,num_clus_group_2,replace=TRUE)
    }
    
    if (is.null(arg_sd)) {sd <- stats::runif(1,min=0,max=10)}
    else {sd <- arg_sd}
    if (is.null(arg_rho)) {rho <- stats::runif(1,min=0,max=bound_rho)}
    else {rho <- arg_rho}
    if (is.null(arg_alt)) {
      find_alt <- sample(1:2,1)
      if (find_alt==1) {
        alt <- "one.sided"
      }
      if (find_alt==2) {
        alt <- "two.sided"
      }
    }
    else {alt <- arg_alt}
    if (is.null(arg_delta)) {
      if (alt == "one.sided") {delta <- stats::runif(1,min=0,max=5)}
      if (alt == "two.sided") {delta <- stats::runif(1,min=-5,max=5)}
    }
    else {delta <- arg_delta}
    if (is.null(arg_sig_level)) {sig_level <- stats::runif(1,min=0.01,max=0.1)}
    else {sig_level <- arg_sig_level}
    
    simulated_power <- simulate_power_t_test_clustered(num_experiments = num_experiments_pr_test,
                                                       clusters_group_1,clusters_group_2, delta,
                                                       sd, rho, alt, sig_level)
    
    calculated_power <- power_t_test_clustered(clusters_group_1 = clusters_group_1, clusters_group_2 = clusters_group_2,
                            delta=delta,sd=sd,rho=rho,alternative=alt,sig_level=sig_level)$power
    
    results[j,1] <- 100*calculated_power
    results[j,2] <- 100*simulated_power
    diff <- 100*(calculated_power - simulated_power)
    results[j,3] <- diff
    if (abs(diff)>3) {
      results[j,4] <- 1
    }
    else {
      results[j,4] <- 0
    }
    if (abs(diff)>8) {
      results[j,5] <- 1
    }
    else {
      results[j,5] <- 0
    }
    results[j,6] <- stats::sd(clusters_group_1)
    results[j,7] <- stats::sd(clusters_group_2)
    results[j,8] <- mean(clusters_group_1)
    results[j,9] <- mean(clusters_group_2)
    results[j,10] <- length(clusters_group_1)
    results[j,11] <- length(clusters_group_2)
    results[j,12] <- delta
    results[j,13] <- sd
    results[j,14] <- rho
    results[j,15] <- find_alt
    
    these_clusters <- list(clusters_1 = clusters_group_1,clusters_2 = clusters_group_2)
    clusters[[j]] <- these_clusters
    
  }
  
  mean_diff <- mean(results[,3])
  sd_diff <- stats::sd(results[,3])
  
  structure(list(balanced = balanced, bounds_clusters_pr_group = bounds_cluster_pr_group, 
                 bounds_cluster_size = bounds_cluster_size, results=results,
                 mean_diff=mean_diff,sd_diff=sd_diff,clusters=clusters))
  
}

