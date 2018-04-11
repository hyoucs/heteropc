load('~/Dropbox/glasso/concord_results/MOTOR/lam3e-3_nz216.Rdata')
a <- inv_cov
num_diag <- dim(a)[1]
idx <- which(a!=0,arr.ind=T)

num_nz <- dim(idx)[1]
print(num_nz-num_diag)
# num_cover <- 0.0

for(i in 1:num_nz){
	x<-idx[i,1]
	y<-idx[i,2]
	a[x,y] <- 1
}

load('~/Dropbox/glasso/concord_results/LANGUAGE/lam5.2e-3_nz220.Rdata')
b <- inv_cov
idx2 <- which(b!=0,arr.ind=T)
num_nz2 <- dim(idx2)[1]
print(num_nz2-num_diag)

for(i in 1:num_nz2){
	x<-idx2[i,1]
	y<-idx2[i,2]
	b[x,y] <- 1
}

idx3 <- which((a-b)!=0.0,arr.ind=T)
num_nz3 <- dim(idx3)[1]
# print(num_nz3)
print(num_nz+num_nz2-num_diag*2-num_nz3)