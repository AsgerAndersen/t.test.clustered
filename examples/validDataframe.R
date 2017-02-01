#***************************************************************************
# a valid instance of a data_experiment.
# * this example fits with the underlying model - see XXXX - only in the 
#   uninteresting case, where the cluster correllation coefficient and the 
#   treatment effect are both 0
# * although this construction orders the rows of the data.frame according 
#   to groups and clusters, this ordering does not have to be present for the 
#   data.frame to be a valid data_experiment
clusters_group_1 <- c(60,50,55)
clusters_group_2 <- c(50,45,50,55)
clusters <- c(clusters_group_1,clusters_group_2)
group_1_size <- sum(clusters_group_1)
group_2_size <- sum(clusters_group_2)
total_size <- group_1_size+group_2_size
group <- as.factor(c(rep(1,group_1_size),rep(2,group_2_size))) 
cluster <- vector()
for (i in 1:length(clusters)) {
 cluster <- c(cluster,rep(i,clusters[i]))
}
cluster <- as.factor(cluster)
response <- rnorm(total_size)
test_data <- data.frame(group,cluster,response)
#***************************************************************************
