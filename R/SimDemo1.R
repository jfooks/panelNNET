



library(mvtnorm)
library(panelNNET)
library(doParallel)
library(parallel)
library(mvtnorm)

rm(list=ls())
gc()
seed_x = 1
set.seed(seed_x) 

z_col <- 2
NumObs <- 400
time_periods <- 5
linear_const <- 1
probTreat = .50
N <- NumObs
t <- time_periods
ndt <- N/t
z_s <- z_col # columns_Z
id <- (1:N-1) %/% (N/ndt) + 1 #id_for_each_time_period
id.eff <- as.numeric(id) # distinct_effect

#mean
id_mean <- foreach(i = 1:ndt) %do% {
  rnorm(z_s, sd=1)
}

#cov_matrix
id_cov <- foreach(i = 1:ndt) %do% {
  A <- matrix(rnorm(z_s^2), z_s)
  t(A) %*% A
}

#dgp_from_distributions
z_1 <- foreach(i = 1:N, .combine = rbind) %do% {
  mvrnorm(1, id_mean[[id[i]]], id_cov[[id[i]]])
}

D <- rbinom(NumObs, 1, probTreat) #dummy vector
beta1 <- c(1, -1)

f_z <- z_1 %*% (beta1)
g_z <- dmvnorm(z_1, rep(0, z_s), diag(rep(1, z_s)))
g_z <- g_z*(sd(f_z)/sd(g_z))

e_1 <- rnorm(N, sd=sd(g_z))
y <- id.eff + f_z + g_z + e_1
id <- as.factor(id)

#############


X <- z_1
hidden_units <- c(8, 5, 2)
fe_var <- id
maxit = 100
lam = 0
time_var = NULL
param = NULL
parapen = rep(0, ncol(param))
parlist = NULL
verbose = TRUE
para_plot = FALSE
report_interval = 10
gravity = 1.01
convtol = 1e-3
bias_hlayers = TRUE
RMSprop = TRUE
start_LR = .01
activation = 'lrelu'
doscale = TRUE
treatment = D
batchsize = nrow(X)
maxstopcounter = 10
OLStrick = FALSE
initialization = 'enforce_normalization'
dropout_hidden = 1
dropout_input = 1
