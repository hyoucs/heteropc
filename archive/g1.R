library(glasso)
set.seed(100)

# multi-variate Gaussian with 10 variables and 10000 sample points
# function generate positive definite matrix (generate eigenvalue)
Posdef <- function (n, ev = runif(n, 0, 10)){
Z <- matrix(ncol=n, rnorm(n^2))
decomp <- qr(Z)
Q <- qr.Q(decomp) 
R <- qr.R(decomp)
d <- diag(R)
ph <- d / abs(d)
O <- Q %*% diag(ph)
Z <- t(O) %*% diag(ev) %*% O
return(Z)
}

num <- 10 #number of genes
var_num <- 2 #number of observations
mu <- runif(num, min=0, max=10) # mean
Sigma <- Posdef(num) # covariance matrix is always positive semi-definite

#generate the normal distribution from mean and co-variance matrix
#n=no of observations of each of the p vectors - Z vector
data <- rmvnorm(n=var_num, mu, Sigma) 
#m x p = var_num x num
#m no of observations of p variate distribution

# Add constraints
#0 can mean no edge, 1 can mean edge exists, 2 means unknown?
J <- matrix(floor(runif(num*num,0,3)),ncol=num)
J[lower.tri(J)] = t(J)[lower.tri(J)]
diag(J) <- 0
print("Constraints matrix is:")
print(J)

theta <- matrix(10, ncol=num, nrow=num)
lasso_dim <- list()
count<-1

for (i in 1:num){
count <- 1
for (j in 1:num){
if (J[i,j] == 0) {
theta[i,j] = 0
} else if (J[i,j]==1){
#list of all dimensions to be considered for linear regression
theta[i,j] = 1
} else {
#list of all dimensions to be considered for lasso regression
if (count == 1){
lasso_dim <- append(lasso_dim,j)
count <- count + 1
} else {
lasso_dim[[i]][count]<-j
count <- count + 1
}
}
}
}

#Perform Lasso regression on all in Z - superimpose constraints on this
#What is the value of regularization parameter? Apply threshold from linear reg?
a <- glasso(Sigma, rho=0.09, approx=TRUE)
a_inv_cov <- a$wi

for (i in 1:num){
for (j in 1:num){
if ((theta[i,j]==0) || (theta[i,j]==1)){
a_inv_cov[i,j] = theta[i,j]
} 
}
}

# Generate a symmetric adjacency matrix (0s and 1s)
A <- round(matrix(runif(var_num * var_num), var_num, var_num))
A[lower.tri(A)] = t(A)[lower.tri(A)]
diag(A) <- 0
a
# TODO: pick some entries as the constraints