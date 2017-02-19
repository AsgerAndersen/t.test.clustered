#****************************************************************************************
#Doing the testing

bounds1 <- list(c(25,30),c(50,70))
bounds2 <- list(c(50,70),c(25,30))
bounds3 <- list(c(5,10),c(50,70))
bounds4 <- list(c(50,70),c(5,10))
bounds5 <- list(c(5,10),c(20,40))
bounds6 <- list(c(20,40),c(5,10))
bounds7 <- list(c(5,10),c(5,10))
bounds8 <- list(c(5,10),c(10,100))
bounds9 <- list(c(10,100),c(5,10))
bounds10 <- list(c(30,60),c(30,60))

bounds <- list(bounds1,bounds2,bounds3,bounds4,bounds5,bounds6,bounds7,bounds8,bounds9,bounds10)

for (j in (1:2)) {
    if (j==1) {bool=TRUE}
    if (j==2) {bool=FALSE}
    for (i in (1:10)) {
      these_bounds = bounds[[i]]
      bounds_clusters_pr_group = these_bounds[[1]]
      bounds_cluster_size = these_bounds[[2]]
      name <- paste("testresults",j,sep="")
      name <- paste(name,i,sep="")
      testresults <- test_calcpower_against_simpower(num_test=20,num_experiments_pr_test = 1000, 
                                                     bounds_cluster_size = bounds_cluster_size,
                                                     bounds_cluster_pr_group = bounds_clusters_pr_group,
                                                     balanced = bool)
      assign(name,testresults)
      rm(testresults,name,these_bounds,bounds_cluster_size,bounds_clusters_pr_group)
      print( paste(Sys.time(), paste(j,i,sep=" "), sep=" "))
    }
    rm(bool)
}
#*********************************************************************************************
#*********************************************************************************************
#Gathering some descriptive statistics from the results

get_num_diff_3 <- function(results) {
  return (sum(results[,4]))
}

get_num_diff_8 <- function(results) {
  return (sum(results[,5]))
}

get_descriptive_stats <- function(testresultList) {
  stats_mat <- matrix(nrow=20,ncol=9)
  colnames(stats_mat) <- c("balanced ", "low-clus-pr-group ", "up-clus-pr-group ", "low-clus-size ", 
                           "up-clus-size ", "mean-diff ", "sd-diff ", "#(|diff|>3) ", "#(|diff|>8) ")
  for (i in (1:20)) {
    results <- testresultList[[i]]
    stats_mat[i,1] <- results$balanced
    stats_mat[i,2] <- results$bounds_clusters_pr_group[1]
    stats_mat[i,3] <- results$bounds_clusters_pr_group[2]
    stats_mat[i,4] <- results$bounds_cluster_size[1]
    stats_mat[i,5] <- results$bounds_cluster_size[2]
    stats_mat[i,6] <- results$mean_diff
    stats_mat[i,7] <- results$sd_diff
    stats_mat[i,8] <- get_num_diff_3(results$results)
    stats_mat[i,9] <- get_num_diff_8(results$results)
  }
  return (stats_mat)
}

allResults <- list(testresults11,testresults12,testresults13,testresults14,testresults15,testresults16,
                   testresults17, testresults18, testresults19, testresults110,testresults21,testresults22,
                   testresults23,testresults24,testresults25,testresults26,testresults27,testresults28,
                   testresults29,testresults210)

overviewTestResults <- get_descriptive_stats(allResults)
