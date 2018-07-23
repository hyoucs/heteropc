library('R.matlab')

library(gconcordopt)

mat_data <- readMat('lw/data_lw1.mat')
data1 <- mat_data$M
mat_data <- readMat('lw/data_lw2.mat')
data2 <- mat_data$M

nonzero <- function(x) sum(x != 0)

# for(i in seq(0.5,0.4,by=-0.01)){
# 	omega1 <- gconcordopt::concordista(data1, lam=i)
# 	nz <- nonzero(omega1)-dim(data1)[2]
# 	# print(nz)

# 	fname = paste0('sw1/',toString(i),'_',toString(nz),'.mat')
# 	writeMat(fname, M=omega1)
# }

# for(i in seq(0.5,0.4,by=-0.01)){
# 	omega2 <- gconcordopt::concordista(data2, lam=i)
# 	nz <- nonzero(omega2)-dim(data2)[2]
# 	# print(nz)

# 	fname = paste0('sw2/',toString(i),'_',toString(nz),'.mat')
# 	writeMat(fname, M=omega2)
# }

D <- rbind(data1, data2)
# writeMat('outlier/data_sw.mat', M=D)

for(i in seq(0.9,0.2,by=-0.1)){
	omega <- gconcordopt::concordista(D, lam=i)
	nz <- nonzero(omega)-dim(D)[2]
	# print(nz)

	fname = paste0('c/',toString(i),'_',toString(nz),'.mat')
	writeMat(fname, M=omega)
}
