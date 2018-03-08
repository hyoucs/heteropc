wthresh <- function(X, type='s', L){
    sign(X) * pmax(abs(X)-L, 0)
}
tr <- function(X){
    sum(diag(X))
}

concordonkar <- function(dat, lambda, eps=1e-8, c=0.5, ...){

    p <- ncol(dat)
    n <- nrow(dat)
    S <- (t(dat)%*%dat)/(n-1)
    
    lambdaMat = 0.5*lambda*(matrix(1,p,p) - diag(p))
    X = diag(diag(S)) + 0.5*lambda*diag(p)
    W = S %*% X
    Df = -diag(1/diag(X)) + 0.5*(W + t(W))
    subg = (X!=0)*(Df + lambdaMat*sign(X)) + (X==0)*wthresh(Df, 's', lambdaMat)

    h = -sum(log(diag(X))) + 0.5 * tr(X%*%(S%*%X))
    epss = 1e-8; c = 0.5; itr = 0; loop = 1; fistory=NULL; tau_n = 1; Dcollect = NULL

    while (loop){

        itr = itr + 1
        tau = tau_n
        Xn = wthresh ( X - tau*Df , 's' , tau*lambdaMat )
        while ( min ( diag ( Xn ) ) < 1e-8 ){
            tau = tau * c
            Xn = wthresh ( X - tau * Df , 's' , tau * lambdaMat )
        }

        Qn = h + tr( (Xn-X)%*%Df ) + (1/(2*tau))*norm(Xn-X, "F")^2
        hn = -sum(log(diag(Xn))) + 0.5*tr(Xn%*%(S%*%Xn))

        while ( hn > Qn ) {
            tau = tau * c
            Xn = wthresh ( X - tau * Df , 's' , tau * lambdaMat ) 
            Qn =  h +  tr ( (Xn-X) %*% Df ) + (1/(2*tau))*norm(Xn-X,"F")^2
            hn =  - sum( log( diag ( Xn ) ) ) + 0.5 * tr( Xn %*% ( S %*% Xn ) )    
        }

        Wn = S %*% Xn;
        Dfn = - ( diag ( 1 / diag ( Xn ) ) ) + 0.5 * ( Wn + t(Wn) );
        tau_n = min( tr ( t( Xn - X )%*%( Xn - X ) ) / tr( t( Xn - X )%*%( Dfn - Df ) ) );

        subg = (Xn!=0)*(Dfn + lambdaMat*sign(Xn)) + (Xn==0)*wthresh(Dfn, 's', lambdaMat)
        loop = ( norm(subg, "F") / norm(Xn, "F") > epss );
        X = Xn; h = hn; Df = Dfn;
        
        f = h + sum(lambdaMat*abs(Xn));
        info = data.frame(itr=itr, f=f, tau=tau, subg_Xn=log10(norm(subg,"F")/norm(Xn,"F")))
        fistory = rbind(fistory, info)

    }

    Xn

}
    
