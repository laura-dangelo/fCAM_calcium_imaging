
## Example code for the data generated from sim_data_scenario1.R ##

library(Rcpp)
library(parallel)
library(RcppArmadillo)
sourceCpp('../fCAM_calcium_imaging.cpp')

load("data_scen1_seed1.Rdata") # this is an example for seed=1
# load the simulated data
# the data is a list named "out" with elements 
# y = observed fluorescence
# g = "group", visual stimulus. series of the same length as y with a categorical variable indicating the group
# cluster = initial guess for the spikes and cluster allocation (see sim_data_scenario1.R for a possible initialization)
# A_start = initial guess for the cluster parameters (related to cluster)

nrep = 4000
burnin = 1:2000
gammapar = 8 

run = calcium_gibbs(Nrep = nrep, 
                    y = out$y,
                    g = out$g,                      
                    cal = c(0,out$y),
                    clO = out$cluster, 
                    clD = 1:6,
                    A_start = out$A_start,
                    b_start = 0,
                    gamma_start = 0.5,
                    sigma2_start = 0.004,
                    tau2_start = 0.0003,
                    p_start = 0.001,
                    alpha_start = 1, 
                    beta_start = 1,
                    maxK_start = 7,
                    maxL_start = 20,
                    c0 = 0, varC0 = 0.1, 
                    hyp_A1 = gammapar, hyp_A2 = gammapar, 
                    hyp_b1 = 0, hyp_b2 = 1, 
                    hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                    hyp_tau21 = 1000, hyp_tau22 = 1, 
                    hyp_gamma1 = 1, hyp_gamma2 = 1,
                    hyp_p1 = 1, hyp_p2 = 999,
                    hyp_alpha1 = 3, hyp_alpha2 = 3,
                    hyp_beta1 = 3, hyp_beta2 = 3,
                    hyp_maxK1 = 4, hyp_maxK2 = 4, hyp_maxK3 = 3,
                    hyp_maxL1 = 4, hyp_maxL2 = 4, hyp_maxL3 = 3,
                    eps_alpha = 0.7, eps_beta = 0.7,
                    eps_gamma = 0.005, eps_A = 0.002,
                    eps_maxK = 7, eps_maxL = 10)

# just to free some memory
y = out$y
g = out$g
A = out$A
spp = which(out$s>0)
rm(out)
#----------------------------------------------------------------------------#
# check convergence of the chains

# spike probability in the "spike and slab" prior
plot(1:length(run$p), run$p, type = "l", ylab = "p")
lines(1:length(run$p), cumsum(run$p)/1:length(run$p), col =2)

# baseline
plot(1:length(run$b), run$b, type = "l", ylab = "b")
lines(1:length(run$b), cumsum(run$b)/1:length(run$p), col =2)

# variance on y_t 
plot(1:length(run$sigma2), run$sigma2, type = "l", ylab = "sigma2")
lines(1:length(run$sigma2), cumsum(run$sigma2)/1:length(run$sigma2), col =2)

# variance on c_t
plot(1:length(run$tau2), run$tau2, type = "l", ylab = "tau2")
lines(1:length(run$tau2), cumsum(run$tau2)/1:length(run$tau2), col =2)

# dirichlet parameters in the mixtures
plot(1:length(run$alpha), run$alpha, type = "l", ylab = "alpha")
lines(1:length(run$alpha), cumsum(run$alpha)/1:length(run$alpha), col =2)

plot(1:length(run$beta), run$beta, type = "l", ylab = "beta")
lines(1:length(run$beta), cumsum(run$beta)/1:length(run$beta), col =2)

# number of components of the mixtures
plot(1:length(run$maxL), run$maxL, type = "l", ylab = "maxL")
lines(1:length(run$maxL), cumsum(run$maxL)/1:length(run$maxL), col =2)

plot(1:length(run$maxK), run$maxK, type = "l", ylab = "maxK")
lines(1:length(run$maxK), cumsum(run$maxK)/1:length(run$maxL), col =2)


# just to free some memory
run$calcium = NULL
run$b = mean(run$b[-burnin,])
run$sigma2 = mean(run$sigma2[-burnin,])
run$tau2 = mean(run$tau2[-burnin,])

#----------------------------------------------------------------------------#
# extract the spikes
A_ext = t(sapply(1:ncol(run$clusterO[,-burnin]), function(x) run$A[,-burnin][run$clusterO[,x]+1,x]) )

est_spikes = colMeans(A_ext) 
est_spikes[which( apply(A_ext, 2, function(x) sum(x>0)) < (nrow(A_ext)/2) )] = 0 # I keep only the spikes which are present in at least half the iterations
times = which(est_spikes>0)
length(times) # number of detected spikes

# analysis of the clusters: compute Rand index to evaluate the accuracy of the clustering
library(mclust)
AA_cluster = apply(A_ext, 1, rank )
rm(A_ext)
rand_indexO = apply(AA_cluster, 2, function(x) adjustedRandIndex(x,rank(A)) )
summary(rand_indexO)
rand_indexD = apply(run$clusterD[,-burnin], 2, function(x) adjustedRandIndex(x, c(1,2,3,4,1,2) ) )
summary(rand_indexD)

false_negatives = sum(sapply(spp, function(x) !(x %in% times)))  # total number of false negatives
false_positives = sum(sapply(times, function(x) !(x %in% spp)))  # total number of false positives

(false_positives + false_negatives) / length(y) # misclassification error rate
