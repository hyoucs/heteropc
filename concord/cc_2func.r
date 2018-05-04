library('R.matlab')
library(gconcordopt)

mat_data <- readMat('~/Dropbox/glasso/concord_results/GAMBLING/data.mat')
data1 <- mat_data$M
mat_data <- readMat('~/Dropbox/glasso/concord_results/EMOTION/data.mat')
data2 <- mat_data$M

D <- rbind(data1, data2)
print(dim(D))

nonzero <- function(x) sum(x != 0)

for(i in seq(0.311,0.315,by=0.001)){
	omega <- gconcordopt::concordista(D, lam=i)
	nz <- nonzero(omega)-dim(D)[2]
	print(nz)

	fname = paste0('result/',toString(i),'_',toString(nz),'.mat')
	writeMat(fname, M=omega)
}