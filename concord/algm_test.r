library('R.matlab')

library(gconcordopt)
library(gconcord)
library(glasso)
library(space)

library(MASS)

p = 10
omega = diag(p)
omega[1,5] = omega[5,1] = .99
omega[2,6] = omega[6,2] = .99
# sigma = solve(omega)

# vectors = mvrnorm(200, mu=rep(0, p), Sigma=sigma)


# save(vectors,file='synthetic.Rdata')
load('synthetic.Rdata')
print(dim(vectors))

####### GROUND TRUTH #######
print(omega)

####### glasso #######
# cc <- cov(vectors)
# a <- glasso(cc,0.11,penalize.diagonal=FALSE)
# inv_cov <- a$wi
# print(inv_cov)
####### concord-ista #######
pMat <-matrix(0,nrow=dim(vectors)[2],ncol=dim(vectors)[2])
inv_cov <- gconcordopt::concordista(vectors, lam=0.5, pMat=pMat)
# print(inv_cov)
####### concord #######
# inv_cov <- concord(vectors,0.6)
# print(inv_cov)
####### space #######
# # Standardize
# vectors <- sweep(vectors, 2L, colMeans(vectors)) #col mean zero
# L2 <- function(x){return(sqrt(sum(x^2)))}
# col_norm <- apply(vectors, 2, L2)
# vectors <- sweep(vectors, 2L, col_norm, "/")	#col normalize

# alpha=0.1
# n = nrow(vectors)
# p = ncol(vectors)
# l1=1/sqrt(n)*qnorm(1-alpha/(2*p^2))
# iter=3
# a <- space.joint(vectors, lam1=l1*n*0.03, iter=iter)
# inv_cov <- a$ParCor
# print(inv_cov)

# Thresholding
tmp <- inv_cov
diag(tmp) <- 0
th <- abs(max(tmp))/100
inv_cov[abs(inv_cov) < th] <- 0
print(inv_cov)

nonzero <- function(x) sum(x != 0)
print(nonzero(inv_cov)-dim(vectors)[2])