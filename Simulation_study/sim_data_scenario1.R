
#----------#         Scenario 1          #----------#
# J = 6
# number of distributional clusters = 4

## mixture ##
sim_data <- function(n, sigma2, tau2, time_spike, b, gamma, prob, par)
{
  c = rep(0,n)
  s = rep(0,n)
  A = rep(0,n)
  s[time_spike==1] = 1
  
  k <- sample(1:length(prob), sum(s), prob, replace = TRUE)
  A[time_spike==1] <- par[k]
  
  for(i in 2:n)
  {
    c[i] = gamma * c[i-1] + A[i] * s[i] + rnorm(1, 0, sqrt(tau2))
  }
  return(list("y" = b + c + rnorm(n, 0, sqrt(sigma2)), "c" = c, "s" = s, "A" = A, "k" = k))
}

times_spike <- function(n, ns, p1, p2)
{
  times = rbinom(n, 1, p1)
  for(i in 1:n)
  {
    if(times[i] == 1)
    {
      times[(i+1):(i+ns)] = rbinom(ns, 1, p2)
    }
  }
  times 
}

data_gen <- function(seedd)
{
  
  # general parameters
  sigma2 = 0.007
  tau2 = 0.0001
  n1 = n2 = n3 = n4 = n5 = n6 = 5000
  gamma = 0.6
  b = 0
  
  # spike probability (assigned uniformly over the interval)
  pp1 = 0.018
  pp2 = 0.010
  pp3 = 0.014
  pp4 = 0.012
  
  # number of consecutive frames after a positive spike I expect additional spikes
  m1 = 5
  m2 = 8
  m3 = 7
  m4 = 10
  
  # spike probability on the interval of length m after a spike
  pm1 = 0.17
  pm2 = 0.13
  pm3 = 0.13
  pm4 = 0.14
  
  
  prob1 = c(0.25, 0.25, 0.2, 0.15, 0.15)
  par1 = c(0.5, 1.2, 1.7, 2, 2.3)
  
  prob2 = rep(0.25, 4)
  par2 = c(0.9, 1.2, 1.7, 2.3)
  
  prob3 = c(0.4, 0.2, 0.4)
  par3 = c(0.5, 1.2, 1.7)
  
  prob4 = rep(1, 3)/3
  par4 = c(0.5, 0.9, 2)
  
  set.seed(seedd)
  spp = c( times_spike(n1, m1, pp1, pm1)[1:n1], # clD1
           times_spike(n2, m2, pp2, pm2)[1:n2], ## clD2
           times_spike(n3, m3, pp3, pm3)[1:n3], ### cld3
           times_spike(n4, m4, pp4, pm4)[1:n4], #### clD4
           times_spike(n5, m1, pp1, pm1)[1:n5], # clD1
           times_spike(n6, m2, pp2, pm2)[1:n6]) ## clD2
  
  
  ### n distr: 4
  set.seed(seedd)
  group1 <- sim_data(n = n1, sigma2 = sigma2, tau2 = tau2, 
                     time_spike = spp[1:n1],
                     gamma = gamma, b = b,
                     prob = prob1, par = par1)
  set.seed(seedd)
  group2 <- sim_data(n = n2, sigma2 = sigma2, tau2 = tau2, 
                     time_spike = spp[(n1+1):(n1 + n2)],
                     gamma = gamma, b = b,
                     prob = prob2, par = par2)
  set.seed(seedd)
  group3 <- sim_data(n = n3, sigma2 = sigma2, tau2 = tau2, 
                     time_spike = spp[(n1 + n2 +1):(n1 + n2 +n3)],
                     gamma = gamma, b = b,
                     prob = prob3, par = par3)
  set.seed(seedd)
  group4 <- sim_data(n = n4, sigma2 = sigma2, tau2 = tau2, 
                     time_spike = spp[(n1 + n2 + n3 +1):(n1 + n2 + n3 + n4)],
                     gamma = gamma, b = b,
                     prob = prob4, par = par4)
  set.seed(seedd)
  group5 <- sim_data(n = n5, sigma2 = sigma2, tau2 = tau2, 
                     time_spike = spp[(n1 + n2 + n3 + n4 +1):(n1 + n2 + n3 + n4 + n5)],
                     gamma = gamma, b = b,
                     prob = prob1, par = par1)
  set.seed(seedd)
  group6 <- sim_data(n = n6, sigma2 = sigma2, tau2 = tau2, 
                     time_spike = spp[(n1 + n2 + n3 + n4 + n5+1):(n1 + n2 + n3 + n4 + n5 + n6)],
                     gamma = gamma, b = b,
                     prob = prob2, par = par2)
  
  y = c(group1$y, group2$y, group3$y, group4$y, group5$y, group6$y)
  g = c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4), rep(5,n5), rep(6,n6))
  A = c(group1$A, group2$A, group3$A, group4$A, group5$A, group6$A)
  s = c(group1$s, group2$s, group3$s, group4$s, group5$s, group6$s)
  k = c(group1$k, group2$k, group3$k, group4$k, group5$k, group6$k)
  
  
  clus = kmeans(diff(y,1)[(diff(y,1))>0.5], centers = 8)
  A_start = rep(0,50)
  A_start[2:9] = c(clus$centers)
  cluster = numeric(length(y))
  cluster[-1][ (diff(y,1))>0.5] = clus$cluster
  
  out = list( "y" = y, "g" = g, "A" = A, "s" = s, "k" = k, "A_start" = A_start, "cluster" = cluster)
  
  filename = paste0("data_scen1_seed", seedd)
  save( out, file = paste0(filename, ".Rdata") )
  
  return(out)
}

sapply(1:50, function(x) data_gen(x))
























