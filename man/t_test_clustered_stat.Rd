% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/t_test_clustered.R
\name{t_test_clustered_stat}
\alias{t_test_clustered_stat}
\title{t statistic for cluster randomized designs}
\usage{
t_test_clustered_stat(data_experiment)
}
\arguments{
\item{data_experiment}{a dataframe with three columns: group, cluster and 
response. The response column holds the values of the measured responses. The 
group column assigns a group number (must be either 1 or 2) to each 
response. The cluster column assigns a cluster number to each response. 
Each different cluster must have a different number, even if they are in 
different groups. See the examples section for an example of a valid
\code{data_experiment}.}
}
\value{
a number - the value of the cluster adjusted t statistic.
}
\description{
Calculates the value of a t statistic that has been adjusted for
a cluster randomized design. See the vignette "Construction of 
the library functions" for a description of the assumptions 
and construction of this t statistic.
}
\examples{
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
#***************************************************************************
# the value of the cluster adjusted t statistic on the data_experiment 
# constructed above
# * (since the intracluster correlation coefficient is 0 for this 
#   data_experiment, the cluster adjusted t statistic is actually 
#   unnecessary here...)
t_test_clustered_stat(data_experiment=test_data)
#***************************************************************************
}
\seealso{
\code{\link{t_test_clustered_pval}} for p value of the cluster 
adjusted t test, \code{\link{power_t_test_clustered}} for power 
calculation of the cluster adjusted t test and 
\code{\link{simulate_power_t_test_clustered}} for power simulation of 
the cluster adjusted t test.
}

