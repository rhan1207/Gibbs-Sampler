#Implement Gibbs Sampler for autocorrelation sequential data
rm(list = ls())
library(MASS)
library(MCMCpack)

# theta contains all the intial values of parameters as well as 
# autocorrelation parameter rho
Gibbs_sampler <- function(data, theta, iter){
  
  X = data$X             # n * p matrix
  Y = data$Y             # n * 1 matrix
  
  beta_0 = theta$beta_0
  Sigma_0 = theta$Sigma_0
  Sigma_0_inv = solve(Sigma_0)
  sigma2_0 = theta$sigma2_0
  v0 = theta$v0
  
  
  sigma2 = sigma2_0
  
  #construct C_rho matrix
  rho = theta$rho
  n = dim(X)[1]
  
  #initialize C_rho
  C_rho = matrix(0, nrow = n, ncol = n)
  
  #fill in values of C_rho
  for (i in 1:n){
    for (j in 1:n){
      abs_diff = abs(i - j)
      C_rho[i,j] = rho^(abs_diff)
    }
  }
  #cache the value of product of beta_0 and Sigma_0_inv
  beta_Sigma_0 = Sigma_0_inv %*% beta_0
  
  # initialize empty tables to store beta and sigma2
  beta_tab = matrix(0, iter, dim(X)[2])
  sigma2_tab = matrix(0, iter, 1)
  
  for (t in 1:iter){
    print(t)
    #compute new Sigma and Sigma inverse
    Sigma = (sigma2 * diag(n)) %*% C_rho
    Sigma_inv = solve(Sigma)
    
    #compute Sigma_n and beta_n
    Sigma_n = solve(t(X)%*%Sigma_inv%*%X + Sigma_0_inv)
    beta_n = Sigma_n %*% (t(X)%*%Sigma_inv%*%Y + beta_Sigma_0)
    #sample and store beta  
    beta = mvrnorm(n = 1, beta_n, Sigma_n, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    beta_tab[t,] = beta
    
    temp = t(Y - X%*%beta) %*% solve(C_rho) %*% (Y - X%*%beta)
    #sample and store sigma^2
    sigma2 = rinvgamma(1, v0/2 + n/2, (v0*sigma2_0 + temp)/2)
    sigma2_tab[t,1] = sigma2
  }
  res = list(beta_tab, sigma2_tab)
  names(res) = c("beta", "sigma2")
  return(res)
}



# create random sample of X: 100 obs * 5 features
mu = c(0, 0, 0, 0, 0)
Sigma = diag(5)
n = 100
X = mvrnorm(n, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)


#fill in values of C_rho with rho = 0.3
rho = 0.3
C_rho = matrix(0, nrow = n, ncol = n)
for (i in 1:n){
  for (j in 1:n){
    abs_diff = abs(i - j)
    C_rho[i,j] = rho^(abs_diff)
  }
}

sigma2 = 1
Sigma_err = (sigma2 * diag(100)) %*% C_rho

err = mvrnorm(n=1, rep(0,100), Sigma_err, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
Y = matrix(0, nrow = n, ncol = 1)
Y[,1] = X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + err

# create data
data = list()
data$X = X
data$Y = Y

# create theta (initial values for Gibbs sampler)
theta = list()
beta = lm(Y ~ 0 + X)
theta$beta_0 = matrix(beta$coefficients, ncol = 1, nrow = 5)
theta$Sigma_0 = var(X)
theta$sigma2_0 = var(Y)[1,]
theta$v0 = 1  
theta$rho = 0.3
  
#test Gibbs Sampler
results = Gibbs_sampler(data, theta, iter = 10000)
  
