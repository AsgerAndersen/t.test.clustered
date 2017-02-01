#*************************************************************************
# power calculation in balanced design:
power_t_test_clustered(clusters_group_1=c(20,20,20), 
clusters_group_2=c(20,20,20), delta=2,sd=1.5,rho=0.5,alternative="one.sided")

# power calculation in unbalanced design:
power_t_test_clustered(clusters_group_1=c(22,20,24), 
clusters_group_2=c(19,21,20,23), delta=2,sd=1.5,rho=0.5,alternative="one.sided")

# calculation of the minimum treatment effect needed to 
# obtain a power of 0.8:
power_t_test_clustered(clusters_group_1=c(22,20,24), 
clusters_group_2=c(19,21,20,23), sd=1.5, rho=0.5, power=0.8, alternative="one.sided")

# calculation of the needed number of measurements in each cluster
# to obtain a power of 0.8:
power_t_test_clustered(delta=2, sd=1.5, rho=0.5, power=0.8, alternative="one.sided", 
sample_size_calc=TRUE, clusters_pr_group=5)

# calculation of the needed number of clusters in each group 
# to obtain a power of 0.8:
power_t_test_clustered(delta=2, sd=1.5, rho=0.5, power=0.8, alternative="one.sided", 
sample_size_calc=TRUE, cluster_size=10)
#*************************************************************************
