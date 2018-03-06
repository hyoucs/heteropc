set.seed(100)
# x<-matrix(rnorm(50*20),ncol=20)
# s<- var(x)
# a<-glasso(s, rho=.01)
# aa<-glasso(s,rho=.02, w.init=a$w, wi.init=a$wi)

# example with structural zeros and no regularization,
# from Whittaker s Graphical models book  page xxx.

# s=c(10,1,5,4,10,2,6,10,3,10)
# S=matrix(0,nrow=4,ncol=4)
# S[row(S)>=col(S)]=s
# S=(S+t(S))
# diag(S)<-10
# zero<-matrix(c(1,3,2,4),ncol=2,byrow=TRUE)
# a<-glasso(S,0,zero=zero)


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

num <- 10000 #no of nodes?
var_num <-10
mu <- runif(var_num, min=0, max=10) # mean
Sigma <- Posdef(n=var_num) # covariance matrix

data <- rmvnorm(n=num, mu, Sigma)
#m x p = num x var_num
#m no of observations of p variate distribution

# Add constraints
#0 can mean no edge, 1 can mean edge exists, 2 means unknown?
J <- matrix(floor(runif(var_num*var_num,0,3)),ncol=var_num)
J[lower.tri(J)] = t(J)[lower.tri(J)]
diag(J) <- 0

theta <- matrix(0, ncol=var_num, nrow=num)

for i in 1:var_num{
k<-1
m<-1
for j in 1:var_num{
if (J[i,j]==0) {
theta[i,j]=0
} else if (J[i,j]==1){
#list of all dimensions to be considered for linear regression
defined_1[i,k] <- j
k<-k+1
} else {
#list of all dimensions to be considered for lasso regression
undefined[i,m] <- j
m<-m+1
}
}

for i in 1:var_num{

}

# Generate a symmetric adjacency matrix (0s and 1s)
A <- round(matrix(runif(var_num * var_num), var_num, var_num))
A[lower.tri(A)] = t(A)[lower.tri(A)]
diag(A) <- 0
a
# TODO: pick some entries as the constraints