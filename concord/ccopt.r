library('R.matlab')
library(gconcordopt)
# library(glasso)
# library(space)
library(MASS)
library(beepr)

mat_data <- readMat('~/Dropbox/glasso/data/HCP-V1/tfMRI-EMOTION.mat')
data <- mat_data$X
dim <- dim(data)
print(dim)

# Create 51 vetors, each with (83*82)/2 dims
vec_dim <- dim[1]*(dim[2]-1)/2
vectors <- matrix(0L,nrow=dim[3], ncol=vec_dim)

# Assign the vectors with data
counter <- 1
for(i in 1:(dim[1]-1)){
	for(j in (i+1):dim[2]){
		for(k in 1:dim[3]){
			vectors[k,counter]<-data[i,j,k]
		}
		counter <- counter+1
	}
}
print(dim(vectors))

# writeMat('~/Dropbox/glasso/concord_results/EMOTION/data.mat',M=vectors)

mat_data <- readMat('~/Dropbox/glasso/pMask.mat')
pMat <- mat_data$M

####### glasso #######
# cc <- cov(vectors)
# a <- glasso(cc,0.11,penalize.diagonal=FALSE)
# inv_cov <- a$wi
####### concord-ista #######
inv_cov <- gconcordopt::concordista(vectors, lam=0.23, pMat=pMat)
####### concord #######
# inv_cov <- concord(vectors,0.6)
####### space #######
# Standardize
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
#####################
tmp <- inv_cov
diag(tmp) <- 0
th <- abs(max(tmp))/100
inv_cov[abs(inv_cov) < th] <- 0

nonzero <- function(x) sum(x != 0)
print(nonzero(inv_cov)-vec_dim)
# print((nonzero(inv_cov)-vec_dim)/(vec_dim*vec_dim-vec_dim))
# save(inv_cov,file='~/Dropbox/glasso/concord_results/RELATIONAL/7.3e-3_.Rdata')
# writeMat('~/Dropbox/glasso/concord_results/EMOTION/0.15_.mat', M=inv_cov)
beep()