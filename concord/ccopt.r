library('R.matlab')
library(gconcordopt)
# library(glasso)
# library(space)
library(MASS)
library(beepr)

load_real_data <-function(dir1,dir2=''){
	# prepare f_data & s_data
	if(dir2 != ''){
		mat_data <- readMat(dir1) #'~/Dropbox/glasso/data/HCP-V1/tfMRI-EMOTION.mat'
		f_data <- mat_data$X

		mat_data <- readMat(dir2) #'~/Dropbox/glasso/data/HCP-V1/Diffusion-q.mat'
		s_data <- mat_data$X
	}else{
		mat_data <- readMat(dir1) #'~/Dropbox/glasso/data/Bassette/data_matrices.mat'
		f_mat <- mat_data$Fs 
		s_mat <- mat_data$Ss 

		dim <- dim(f_mat[[1]][[1]])

		f_data <- array(0,c(dim[1],dim[2],length(f_mat)))
		s_data <- array(0,c(dim[1],dim[2],length(f_mat)))

		for(i in 1:length(f_mat)){
			f_data[,,i] = f_mat[[i]][[1]]
			s_data[,,i] = s_mat[[i]][[1]]
		}
	}
	dim <- dim(f_data)
	print(dim)
	# build data vectors
	# Create 51 vetors, each with (83*82)/2 dims
	p <- dim[1]*(dim[2]-1)/2
	f <- matrix(0L,nrow=dim[3], ncol=p) # function vector
	s <- matrix(0L,nrow=dim[3], ncol=p)

	# Assign the vectors with data
	counter <- 1
	for(i in 1:(dim[1]-1)){
		for(j in (i+1):dim[2]){
			for(k in 1:dim[3]){
				f[k,counter]<-f_data[i,j,k]
				s[k,counter]<-s_data[i,j,k]
			}
			counter <- counter+1
		}
	}
	returnlist <- list('f'=f, 's'=s)
	returnlist
}

load_syn_data <- function(dir_='synthetic.mat'){
	mat_data <- readMat(dir_)
	f <- mat_data$F 
	s <- mat_data$S
	returnlist <- list('f'=f, 's'=s)
	returnlist
}

nonzero <- function(x) sum(x != 0)

glasso_omega <- function(vectors,root='result/'){
	cc <- cov(vectors)
	a <- glasso(cc,0.11,penalize.diagonal=FALSE)
	inv_cov <- a$wi

	nz <- (nonzero(inv_cov)-p)/2
	print(nz)

	fname = paste0(root,'glasso_',toString(nz),'.mat')
	writeMat(fname, M=inv_cov)
	# save(inv_cov,file='***.Rdata')
}

space_par_cor <- function(vectors,root='result/'){
	# Standardize
	vectors <- sweep(vectors, 2L, colMeans(vectors)) #col mean zero
	L2 <- function(x){return(sqrt(sum(x^2)))}
	col_norm <- apply(vectors, 2, L2)
	vectors <- sweep(vectors, 2L, col_norm, "/")	#col normalize

	alpha=0.1
	n <- nrow(vectors)
	l1 <- 1/sqrt(n)*qnorm(1-alpha/(2*p^2))
	iter <- 3
	a <- space.joint(vectors, lam1=l1*n*0.03, iter=iter)
	par_cor <- a$ParCor

	nz <- (nonzero(par_cor)-p)/2
	print(nz)

	fname = paste0(root,'space_',toString(nz),'.mat')
	writeMat(fname, M=par_cor)
}

concord_omega <- function(vectors,lam,root='result/'){
	inv_cov <- concord(vectors,lam)

	nz <- (nonzero(inv_cov)-p)/2
	print(nz)

	fname = paste0(root,toString(lam),'_concord_',toString(nz),'.mat')
	writeMat(fname, M=inv_cov)
}

concord_ista_crossval <- function(vectors,s){
	vectors <- sign(vectors)*(abs(vectors)-s)
	vectors <- vectors-s
	fold <- 10.0
	n <- nrow(vectors)
	t <- floor(n/fold)
	print(t)
	m <- 0
	for(k in seq(0, 9, by=1)){
		d1<-head(vectors,k*t)
		d2<-tail(vectors,n-(k+1)*t)
		train <- rbind(d1,d2)

		inv_cov <-gconcordopt::concordista(train,lam=0.29)
		nz <- (nonzero(inv_cov)-p)/2
		print(nz)
		m <- m+nz

		fname = paste0('syn/sf/',toString(k),'.mat')
		writeMat(fname, S=inv_cov)
	}
	print(m/10.0)
}

concord_ista <- function(vectors, thresh=FALSE){
	for(i in seq(0.259,0.259,by=-0.0002)){
		inv_cov <- gconcordopt::concordista(vectors, lam=i, pMat=pMat)
		# inv_cov <- gconcordopt::concordista(vectors, lam=i)
		
		if(thresh){
			tmp <- inv_cov
			diag(tmp) <- 0
			th <- abs(max(tmp))/100
			inv_cov[abs(inv_cov) < th] <- 0
		}

		nz <- (nonzero(inv_cov)-p)/2
		print(nz)

		fname = paste0('result/',toString(i),'_',toString(nz),'.mat')
		writeMat(fname, M=inv_cov)
	}

}

main <- function(dir1, dir2='', pMat='', normalize_s=FALSE){
	data <- load_real_data(dir1,dir2)
	f <- data$f 
	s <- data$s

	if(normalize_s){
		print(min(s))
		print(max(s))
		s <- s-min(s)
		s <- s/max(s)
	}

	vectors <- cbind(f,s)
	# vectors <- s
	p <- ncol(vectors)
	print(dim(vectors))

	browser()

	if(pMat!=''){
		# mat_data <- readMat('~/Dropbox/glasso/pMask.mat')
		mat_data <- readMat(pMat)
		pMat <- mat_data$M
	}else{
		pMat <- matrix(1L,p,p)
	}
	beep()
}

main('~/Dropbox/glasso/data/Bassette/data_matrices.mat','','pMask/pMask_ss.mat')