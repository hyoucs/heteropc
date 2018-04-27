library('R.matlab')
library('MVN')
library('gconcord')
library(beepr)

mat_data <- readMat('~/Dropbox/glasso/data/HCP-V1/tfMRI-EMOTION.mat')
data <- mat_data$X
dim <- dim(data)
print(dim)

# Create 51 vectors, each with (83*82)/2 dims
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
print(dim(vectors))
# cc <- cov(vectors)
# print(dim(cc))
nonzero <- function(x) sum(x != 0)

inv_cov<-concord(vectors,0.6)
print(nonzero(inv_cov)-vec_dim)
writeMat('~/Dropbox/glasso/concord_results/EMOTION/0.42_.mat', M=inv_cov)
beep()
# print(nonzero(inv_cov))
# print((nonzero(inv_cov)-vec_dim)/(vec_dim*vec_dim-vec_dim))