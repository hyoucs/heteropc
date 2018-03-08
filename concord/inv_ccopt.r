library('R.matlab')
library(gconcordopt)
library(MASS)

mat_data <- readMat('~/Dropbox/glasso/data/HCP-V1/tfMRI-GAMBLING.mat')
data <- mat_data$X
dim <- dim(data)
print(dim)

# Create 51 vetors, each with (83*82)/2 dims
vec_dim <- dim[1]*(dim[2]-1)/2
vectors <- matrix(0L,nrow=dim[3], ncol=vec_dim)

# Assign the vector with data
counter <- 1
for(i in 1:(dim[1]-1)){
	for(j in (i+1):dim[2]){
		for(k in 1:dim[3]){
			vectors[k,counter]<-data[i,j,k]
		}
		counter <- counter+1
	}
}

cc <- cov(vectors)

mat_data <- readMat('~/Dropbox/glasso/pMask.mat')
pMat <- mat_data$M
dim <- dim(pMat)
print(dim)

inv_cov<-gconcordopt::concordista(cc, lam=0.01, pMat=pMat)
nonzero <- function(x) sum(x != 0)
print(nonzero(inv_cov))
print((nonzero(inv_cov)-vec_dim)/(vec_dim*vec_dim-vec_dim))