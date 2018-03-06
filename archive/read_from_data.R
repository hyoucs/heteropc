library(glasso)
library(mvtnorm)
set.seed(101)

temp1 <- read.table("C:\\Users\\Nidhi1\\Dropbox\\glasso\\insilico_size10_1_knockdowns.tsv", sep="\t",  header=TRUE)
temp2 <- as.matrix(temp1)
obs <- nrow(temp2)
dimensions <- ncol(temp2)
mu <- matrix(dimensions, nrow=1, ncol=dimensions)
for (i in 1:dimensions){
mu[1,i]<-mean(temp2[,i])
}
Sigma<-cov(temp2)

data <- temp2
num <- dimensions
var_num <- obs

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
a <- glasso(Sigma, rho=0.03, approx=TRUE)
a_inv_cov <- a$wi

for (i in 1:num){
for (j in 1:num){
if ((theta[i,j]==0) || (theta[i,j]==1)){
a_inv_cov[i,j] = theta[i,j]
} 
}
}

# # Generate a symmetric adjacency matrix (0s and 1s)
# A <- round(matrix(runif(var_num * var_num), var_num, var_num))
# A[lower.tri(A)] = t(A)[lower.tri(A)]
# diag(A) <- 0
# a
# # TODO: pick some entries as the constraints