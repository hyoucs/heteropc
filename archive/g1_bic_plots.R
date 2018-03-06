library(glasso)
library(mvtnorm)
library(network)
library(Matrix)
library(matrixcalc)
library(matlib)
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
var_num <- 6 #number of observations
mu <- runif(num, min=0, max=10) # mean
Sigma <- Posdef(num) # covariance matrix is always positive semi-definite
print("Covariance matrix (gt) is:")
print(Sigma)

Sigma <- round(Sigma,10) #to make sure Sigma is symmetric
Sigma_inv <- inv(Sigma)
print("Inverse covariance matrix (gt) is:")
print(Sigma_inv)

threshold = 0.03
Sigma_inv[abs(Sigma_inv)<threshold] <- 0
print(Sigma_inv)

par(mfrow=c(2,2))
g1<-network(Sigma_inv)
plot(g1, main="Inv covariance matrix of data")

#generate the normal distribution from mean and co-variance matrix
#n=no of observations of each of the p vectors - Z vector
data <- rmvnorm(n=var_num, mu, Sigma) 
#m x p = var_num x num
#m no of observations of p variate distribution

# Add constraints
#0 can mean no edge, 1 can mean edge exists, 2 means unknown

#constraint_entries
c_entries <- matrix(floor(runif(num*num,0,2)),ncol=num)
c_entries[lower.tri(c_entries)] = t(c_entries)[lower.tri(c_entries)]

print("entries are:")
print(c_entries)

J <-matrix(0,ncol=num,nrow=num)
for(i in 1:num){
	for(j in 1:num){
		if(c_entries[i,j]==0){ #unknown
			J[i,j] <- 2
		}else{
			if(Sigma_inv[i,j]==0){
				J[i,j] <- 0
			}else{
				J[i,j] <- 1
			}
		}
	}
}

# J <- matrix(floor(runif(num*num,0,3)),ncol=num)
# J[lower.tri(J)] = t(J)[lower.tri(J)]
diag(J) <- 0

print("Constraints matrix is:")
print(J)

# J[J==2]<-0
# g0<-network(J)
# plot(g0, main="External constraints")

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

#bic to find rho
nr  <- 100
rho <- seq(0.1,1,length=nr)
bic <- rho
for(j in 1:nr){
a       <- glasso(Sigma,rho[j],approx=TRUE)
p_off_d <- sum(a$wi!=0 & col(Sigma)<row(Sigma))
bic[j]  <- -2*(a$loglik) + p_off_d*log(num*num)
}
best <- which.min(bic)
plot(rho,bic)
points(rho[best],bic[best],pch=15)

#Perform Lasso regression on all in Z - superimpose constraints on this
#What is the value of regularization parameter? Apply threshold from linear reg?
a <- glasso(Sigma, 0.4, approx=TRUE)
a_inv_cov <- a$wi
g2 <- network(a$wi)
plot(g2, main="Estimated inverse covariance matrix")

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
print(a)
# TODO: pick some entries as the constraints