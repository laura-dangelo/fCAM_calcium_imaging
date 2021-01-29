library(Rcpp)
library(RcppDist)

library(ggplot2)
library(viridis)

sourceCpp('./fCAM_calcium_imaging.cpp')

 
data <- read.csv("cell_data.csv", header = FALSE)
y_real = c(data$V1)
rm(list = ("data"))
n = length(y_real)
load(file = "group.Rdata")




#---------------------------------------------------------#

library(FastLZeroSpikeInference)
fit <- FastLZeroSpikeInference::estimate_spikes(dat = y_real, gam = 0.35, lambda = 0.004)
clus = kmeans(y_real[fit$spikes], centers = 5)

nrep = 8500
start <- Sys.time()
run = calcium_gibbs(Nrep = nrep,
                    y = y_real,
                    g = g,
                    cal = c(0,y_real),
                    clO = clus$cluster,
                    clD = 1:4,
                    A_start = clus$centers,
                    b_start = 0,
                    gamma_start = 0.6,
                    sigma2_start = 0.001,
                    tau2_start = 0.00002,
                    p_start = 0.01,
                    alpha_start = 0.5,
                    beta_start = 0.5,
                    maxK_start = 5,
                    maxL_start = 10,
                    c0 = 0, varC0 = 0.1,
                    hyp_A1 = 5, hyp_A2 = 7,
                    hyp_b1 = 0, hyp_b2 = 1,
                    hyp_sigma21 = 1000, hyp_sigma22 = 1,
                    hyp_tau21 = 1000, hyp_tau22 = 1,
                    hyp_gamma1 = 1, hyp_gamma2 = 1,
                    hyp_p1 = 1, hyp_p2 = 999,
                    hyp_alpha1 = 3, hyp_alpha2 = 3,
                    hyp_beta1 = 3, hyp_beta2 = 6,
                    hyp_maxK1 = 2, hyp_maxK2 = 4, hyp_maxK3 = 3,
                    hyp_maxL1 = 2, hyp_maxL2 = 4, hyp_maxL3 = 3,
                    eps_alpha = 0.5, eps_beta = 0.7,
                    eps_gamma = 0.003,
                    eps_A = 0.02,
                    eps_maxK = 4, eps_maxL = 5)
end <- Sys.time()
end - start

burnin = 1:7000
run$calcium = run$calcium[,-burnin]
run$clusterO = run$clusterO[,-burnin]
run$clusterD = run$clusterD[,-burnin]
run$b = run$b[-burnin]
run$gamma = run$gamma[-burnin]
run$sigma2 = run$sigma2[-burnin]
run$tau2 = run$tau2[-burnin]
run$A = run$A[,-burnin]
run$p = run$p[-burnin]
run$alpha = run$alpha[-burnin]
run$beta = run$beta[-burnin]
run$maxK = run$maxK[-burnin]
run$maxL = run$maxL[-burnin]

run$calcium = NULL

#---------------------------------------------------------#
# check convergence

plot(1:length(run$p), run$p, type = "l", ylab = "p")
lines(1:length(run$p), cumsum(run$p)/1:length(run$p), col =2)

plot(1:length(run$sigma2), run$sigma2, type = "l", ylab = "sigma2")
lines(1:length(run$sigma2), cumsum(run$sigma2)/1:length(run$sigma2), col =2)

plot(1:length(run$tau), run$tau, type = "l", ylab = "tau2")
lines(1:length(run$tau), cumsum(run$tau)/1:length(run$tau), col =2)

plot(1:length(run$b), run$b, type = "l", ylab = "b")
lines(1:length(run$b), cumsum(run$b)/1:length(run$b), col =2)

plot(1:length(run$gamma), run$gamma, type = "l", ylab = "gamma")
lines(1:length(run$gamma), cumsum(run$gamma)/1:length(run$gamma), col =2)

plot(1:length(run$alpha), run$alpha, type = "l", ylab = "alpha")
lines(1:length(run$alpha), cumsum(run$alpha)/1:length(run$alpha), col =2)

plot(1:length(run$beta), run$beta, type = "l", ylab = "beta")
lines(1:length(run$beta), cumsum(run$beta)/1:length(run$beta), col =2)

plot(1:length(run$maxL), run$maxL, type = "l", ylab = "maxL")
lines(1:length(run$maxL), cumsum(run$maxL)/1:length(run$maxL), col =2)

plot(1:length(run$maxK), run$maxK, type = "l", ylab = "maxK")
lines(1:length(run$maxK), cumsum(run$maxK)/1:length(run$maxL), col =2)





AA_gMFM = t(sapply(1:ncol(run$clusterO), function(x) run$A[run$clusterO[,x]+1,x]) )

clusterO = run$clusterO
A_par = run$A
rm(run)

est_spikes = colMeans(AA_gMFM)
est_spikes[which( apply(t(clusterO), 2, function(x) mean(x != 0))<0.5)] = 0
times = which(est_spikes>0)
length(times)


#----------------# plot data + activity #----------------# 

df = data.frame(x = 1:n, y = y_real, g = g)

which(diff(g)!=0)

df_rect = data.frame(start = c(745, 16101, 39581, 54932, 70294, 80227, 97383),
                     end = c(15198, 30550, 54029, 69391, 79323, 96104, 113637),
                     Stimulus = as.factor(c("Static grating","Natural scene","Natural scene",
                                            "Static grating","Natural movie","Natural scene",
                                            "Static grating")) )


cols = c("#91ff00", "#00fffb","#ff3700")

df$x = df$x/30
df_rect$start = df_rect$start/30
df_rect$end = df_rect$end/30

df$AA = est_spikes

ggplot(data = df) +
  geom_rect(data = df_rect, inherit.aes = FALSE,
            aes(xmin = start, 
                xmax = end, 
                ymin = -Inf, 
                ymax = Inf, fill = Stimulus), alpha = 0.12 ) +
  scale_fill_manual(values = cols) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  theme(legend.position = "bottom",
        rect = element_rect(fill="transparent", colour=NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(name = "Time (seconds)") +
  scale_y_continuous(name = "Calcium level") +
  geom_line(aes(x = x, y = AA), col = "gold") 
