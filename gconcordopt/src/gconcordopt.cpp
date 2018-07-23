#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;

double sgn(double val) {
  return (double(0) < val) - (val < double(0));
}

double sthresh(double x, double t ){
  return sgn(x) * max(abs(x)-t, 0.0);
}

// SEXP sthreshmat(Eigen::Map<Eigen::MatrixXd> & x,
// 		double tau,
// 		Eigen::Map<Eigen::MatrixXd> & t){
void sthreshmat(Eigen::MatrixXd & x,
		double tau,
		Eigen::MatrixXd & t){

  Eigen::MatrixXd tmp1(x.cols(), x.cols());
  Eigen::MatrixXd tmp2(x.cols(), x.cols());
  
  tmp1 = x.array().unaryExpr(ptr_fun(sgn));
  tmp2 = (x.cwiseAbs() - tau*t).cwiseMax(0.0);

  x = tmp1.cwiseProduct(tmp2);

  return;
}

// another simple example: outer product of a vector, 
// returning a matrix
//
// [[Rcpp::export]]
SEXP cc_ista(const Eigen::Map<Eigen::MatrixXd> & D, double lam, 
       const Eigen::Map<Eigen::MatrixXd> & pMat,
	     const Eigen::MappedSparseMatrix<double> & X0_,
	     int penalize_diag = 0, int BBstep = 1,
	     double tol = 1e-5, int maxit = 100, int info = 0, 
	     int trkeigs = 0, int DisS = 0) {
	     // const Eigen::Map<Eigen::MatrixXd> & X0_,


  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;

  string convcond = "";
  NumericVector fhist(maxit+1);
  NumericVector tauehist(maxit+1); // effective tau (decreased by some power of c)
  NumericVector tauohist(maxit+1); // original tau computed by BB step
  NumericVector subgXnhist(maxit+1);
  NumericVector maxdiffhist(maxit+1);
  NumericVector eigprodhist(maxit+1);
  IntegerVector diagitrhist(maxit);
  IntegerVector backitrhist(maxit);

  int p = D.cols(); // problem size

  // standardize
  MatrixXd S;

  if (DisS) {
    S = D;
  } else {
    // standardize
    S = D.rowwise() - D.colwise().mean();                 // mean zero
    S = (S.transpose() * S) / (S.rows() - 1);
  }

  MatrixXd LambdaMat(p, p);

  for (int j=0; j<p; j++){
    for(int i=0; i<p; i++){
      if(pMat(i,j)==1){
        LambdaMat(i,j) = lam;
        // if(i==j){
        //   LambdaMat(i,j) = lam/10;
        // }
      }else{
        LambdaMat(i,j) = lam*100; //*100
      }
    }
  }
  // LambdaMat.setConstant(lam);
  if (penalize_diag == 0) {
    LambdaMat.diagonal().setZero().eval();
  }


  // // initial guess is a diagonal matrix
  // SparseMatrix<double> X(p, p);      // Omega
  // X = MatrixXd (S.diagonal().asDiagonal()).sparseView();
  // for (int i=0; i < p; i++) { 
  //   X.coeffRef(i,i) += 0.5 * lam;
  // }
  // // cout << X << endl;

  
  // SparseMatrix<double> X = MatrixXd::Identity(p,p).sparseView();
  SparseMatrix<double> X(X0_);
  SparseMatrix<double> Xn(p, p);     // next omega
  SparseMatrix<double> Step(p, p);   // Omega_n - Omega

  MatrixXd W = S * X;
  MatrixXd Wn(p, p);

  MatrixXd G(p, p);
  MatrixXd Gn(p, p);
  MatrixXd subg(p, p); // computesubg here
  MatrixXd tmp(p, p);
  
  double h = - X.diagonal().array().log().sum() + 0.5*(X*W).trace();
  double hn = 0; 
  double Qn = 0;
  double f = 0;
  double subgnorm, Xnnorm;
  double maxdiff;

  double tau;
  double taun = 1.0;
  double c = 0.5;
  int itr = 0;
  int loop = 1;
  int diagitr = 0;
  int backitr = 0;

  G = 0.5 * (W + W.transpose());
  G += - MatrixXd((MatrixXd((1/X.diagonal().array()))).asDiagonal());

  while (loop != 0){

    tau = taun;
    // cout << "starting loop: " << tau << endl;
    
    diagitr = 0;
    backitr = 0;

    while ( 1 ) {

      if (diagitr != 0 || backitr != 0) { tau = tau * c; } // decrease tau only if needed

      tmp = MatrixXd(X) - tau*G;
      sthreshmat(tmp, tau, LambdaMat);
      Xn = tmp.sparseView();

      if (Xn.diagonal().minCoeff() < 1e-8 && diagitr < 10) { 
	// cout << "diagonal scaling" << endl;
	diagitr += 1;
	continue;
      }

      Step = Xn - X;
      Wn = S * Xn;
      Qn = h + Step.cwiseProduct(G).sum() + (1/(2*tau))*pow(Step.norm(),2);
      hn = - Xn.diagonal().array().log().sum() + 0.5*(Xn*Wn).trace();

      if (hn > Qn) { 
	// cout << "backtracking" << endl;
	backitr += 1;
      } else {
	// cout << "terminating" << endl;
	break;
      }

    }

    if (info) {
      tauehist(itr) = tau;  // effective tau
      tauohist(itr) = taun; // original tau
    }

    Gn = 0.5 * (Wn + Wn.transpose());
    Gn += - MatrixXd((MatrixXd((1/Xn.diagonal().array()))).asDiagonal());

    if ( BBstep == 0 ) {
      taun = 1;
    } else if ( BBstep == 1 ) {
      taun = ( Step * Step ).eval().diagonal().array().sum() / (Step * ( Gn - G )).trace();
    }

    // compute subg
    tmp = MatrixXd(Xn).array().unaryExpr(ptr_fun(sgn));   // sign term
    tmp = Gn + tmp.cwiseProduct(LambdaMat);               // first term is in "tmp"
    subg = Gn;                                            // second term is in "subg"
    sthreshmat(subg, 1.0, LambdaMat);
    subg = (MatrixXd(Xn).array() != 0).select(tmp, subg); // select terms

    subgnorm = subg.norm();
    Xnnorm = Xn.norm();

    maxdiff = 0;
    for (int k=0; k<Step.outerSize(); ++k) {
      for (SparseMatrix<double>::InnerIterator it(Step,k); it; ++it) {
	maxdiff = max(abs(it.value()), maxdiff);
      }
    }

    X = Xn; h = hn; G = Gn;
    f = h + Xn.cwiseAbs().cwiseProduct(LambdaMat).sum();
    // cout<<"1st term: "<<- Xn.diagonal().array().log().sum()<<"; 2nd term: "<<0.5*(Xn*Wn).trace()<<endl;
    // cout<<"Likelihood: "<<h<<endl;
    // cout<<"L1 penality: "<<f-h<<endl;

    if (info) {
      fhist(itr) = f;    
      subgXnhist(itr) = subgnorm/Xnnorm;
      maxdiffhist(itr) = maxdiff;
      diagitrhist(itr) = diagitr;
      backitrhist(itr) = backitr;
      if (trkeigs) {
	eigprodhist(itr) = double((MatrixXd(Xn).eigenvalues()).prod().real());
      } else {
	eigprodhist(itr) = -1.0;
      }
    }
    itr += 1;

    // loop = int((itr < maxit) && (maxdiff > tol) && (subgnorm/Xnnorm > tol));
    // loop = int((itr < maxit) && (maxdiff > tol));
    loop = int((itr < maxit) && (subgnorm/Xnnorm > tol));

    if (info) {
      if (subgnorm/Xnnorm < tol) { convcond += "subg_2/Xn_2"; }
      if (maxdiff < tol)         { convcond += "maxdiff"; }
      if (itr > maxit)           { convcond += "maxitr"; }
    }

  }
  cout<<"Likelihood "<<h<<endl;
  cout<<"L1 penality "<<f-h<<endl;
  
  if (info) {
    return List::create(Named("omega") = wrap(Xn),
			Named("itr") = itr,
			Named("convcond") = wrap(convcond),
			Named("info") = 
			DataFrame::create(Named("f") = fhist[Range(0,itr-1)],
					  Named("taue") = tauehist[Range(0,itr-1)],
					  Named("tauo") = tauohist[Range(0,itr-1)],
					  Named("subgXn") = subgXnhist[Range(0,itr-1)],
					  Named("maxdiff") = maxdiffhist[Range(0,itr-1)],
					  Named("eigprod") = eigprodhist[Range(0,itr-1)],
					  Named("diagitr") = diagitrhist[Range(0,itr-1)],
					  Named("backitr") = backitrhist[Range(0,itr-1)]));
  } else {
    return List::create(Named("omega") = wrap(Xn),
			Named("itr") = itr);
  }
}

// another simple example: outer product of a vector, 
// returning a matrix
//
// [[Rcpp::export]]
SEXP cc_fista(const Eigen::Map<Eigen::MatrixXd> & D, double lam,
        const Eigen::Map<Eigen::MatrixXd> & pMat,
	      const Eigen::MappedSparseMatrix<double> & X0_,
	      int penalize_diag = 0, int steptype = 1,
	      double tol = 1e-5, int maxit = 100, int info = 0, 
	      int trkeigs = 0, int DisS = 0) {
  
  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;

  string convcond = "";
  NumericVector fhist(maxit+1);
  NumericVector tauehist(maxit+1); // effective tau (decreased by some power of c)
  NumericVector tauohist(maxit+1); // original tau computed by BB step
  NumericVector subgXnhist(maxit+1);
  NumericVector maxdiffhist(maxit+1);
  NumericVector eigprodhist(maxit+1);
  IntegerVector diagitrhist(maxit);
  IntegerVector backitrhist(maxit);

  int p = D.cols(); // problem size

  MatrixXd S;

  if (DisS) {
    S = D;
  } else {
    // standardize
    S = D.rowwise() - D.colwise().mean();                 // mean zero
    S = (S.transpose() * S) / (S.rows() - 1);
  }

  MatrixXd LambdaMat(p, p);
  for (int j=0; j<p; j++){
    for(int i=0; i<p; i++){
      if(pMat(i,j)==1){
        LambdaMat(i,j) = lam;
      }else{
        LambdaMat(i,j) = lam*100;
      }
    }
  }
  // LambdaMat.setConstant(lam);
  if (penalize_diag == 0) {
    LambdaMat.diagonal().setZero().eval();
  }

  // // initial guess is a diagonal matrix
  // SparseMatrix<double> X(p, p);      // Omega
  // X = MatrixXd (S.diagonal().asDiagonal()).sparseView();
  // for (int i=0; i < p; i++) { 
  //   X.coeffRef(i,i) += 0.5 * lam;
  // }
  // // cout << X << endl;

  // SparseMatrix<double> X =  MatrixXd::Identity(p, p).sparseView();
  // SparseMatrix<double> Theta =  MatrixXd::Identity(p, p).sparseView();
  SparseMatrix<double> X(X0_);
  SparseMatrix<double> Theta(X0_);
  SparseMatrix<double> Xn(p, p);     // next omega
  SparseMatrix<double> Step(p, p);   // Omega_n - Omega

  MatrixXd W = S * X;
  MatrixXd Wn(p, p);
  MatrixXd WTh = S * Theta;

  MatrixXd G(p, p);
  MatrixXd Gn(p, p);
  MatrixXd subg(p, p); // computesubg here
  MatrixXd tmp(p, p);
  
  double h = - X.diagonal().array().log().sum() + 0.5*(X*W).trace();
  double hn = 0; 
  double hTh = 0;
  double Qn = 0;
  double f = 0;
  double subgnorm, Xnnorm;
  double maxdiff;

  double tau;
  double taun = 1.0;
  double alpha = 1.0;
  double alphan;
  double c = 0.9;
  int itr = 0;
  int loop = 1;
  int diagitr = 0;
  int backitr = 0;

  G = 0.5 * (WTh + WTh.transpose());
  G += - MatrixXd((MatrixXd((1/Theta.diagonal().array()))).asDiagonal());

  while (loop != 0){

    tau = taun;
    // cout << "starting loop: " << tau << endl;
    
    diagitr = 0;
    backitr = 0;

    while ( 1 ) {

      if (diagitr != 0 || backitr != 0) { tau = tau * c; } // decrease tau only if needed

      tmp = MatrixXd(Theta) - tau*G;
      sthreshmat(tmp, tau, LambdaMat);
      Xn = tmp.sparseView();

      if (Xn.diagonal().minCoeff() < 1e-8 && diagitr < 50) { 
	// cout << "diagonal scaling" << endl;
	diagitr += 1;
	continue;
      }

      Step = Xn - Theta;
      Wn = S * Xn;
      hTh = - Theta.diagonal().array().log().sum() + 0.5*(Theta*WTh).trace();
      Qn = hTh + Step.cwiseProduct(G).sum() + (1/(2*tau))*pow(Step.norm(),2);
      hn = - Xn.diagonal().array().log().sum() + 0.5*(Xn*Wn).trace();

      if (hn > Qn) { 
	// cout << "backtracking" << endl;
	backitr += 1;
      } else {
	// cout << "terminating" << endl;
	break;
      }

    }

    if (info) {
      tauehist(itr) = tau;  // effective tau
      tauohist(itr) = taun; // original tau
    }
    
    alphan = (1 + sqrt(1 + 4*pow(alpha,2)))/2;

    Theta = Xn + ((alpha - 1)/alphan) * (Xn - X);
    
    WTh = S * Theta;
    Gn = 0.5 * (WTh + WTh.transpose());
    Gn += - MatrixXd((MatrixXd((1/Theta.diagonal().array()))).asDiagonal());

    if ( steptype == 0 ) {
      taun = 1;
    } else if ( steptype == 1 ) {
      taun = tau;
    } else if ( steptype == 2 ) {
      taun = ( Step * Step ).eval().diagonal().array().sum() / (Step * ( Gn - G )).trace();
      if (taun < 0.0) { taun = tau; }
    }

    // compute subg
    tmp = MatrixXd(Xn).array().unaryExpr(ptr_fun(sgn));   // sign term
    tmp = Gn + tmp.cwiseProduct(LambdaMat);               // first term is in "tmp"
    subg = Gn;                                            // second term is in "subg"
    sthreshmat(subg, 1.0, LambdaMat);
    subg = (MatrixXd(Xn).array() != 0).select(tmp, subg); // select terms

    subgnorm = subg.norm();
    Xnnorm = Xn.norm();

    maxdiff = 0;
    for (int k=0; k<Step.outerSize(); ++k) {
      for (SparseMatrix<double>::InnerIterator it(Step,k); it; ++it) {
	maxdiff = max(abs(it.value()), maxdiff);
      }
    }

    alpha = alphan;
    X = Xn; h = hn; G = Gn;
    f = h + Xn.cwiseAbs().cwiseProduct(LambdaMat).sum();

    if (info) {
      fhist(itr) = f;    
      subgXnhist(itr) = subgnorm/Xnnorm;
      maxdiffhist(itr) = maxdiff;
      diagitrhist(itr) = diagitr;
      backitrhist(itr) = backitr;
      if (trkeigs) {
	eigprodhist(itr) = double((MatrixXd(Xn).eigenvalues()).prod().real());
      } else {
	eigprodhist(itr) = -1.0;
      }
    }
    itr += 1;

    // loop = int((itr < maxit) && (maxdiff > tol) && (subgnorm/Xnnorm > tol));
    // loop = int((itr < maxit) && (maxdiff > tol));
    loop = int((itr < maxit) && (subgnorm/Xnnorm > tol));

    if (info) {
      if (subgnorm/Xnnorm < tol) { convcond += "subg_2/Xn_2"; }
      if (maxdiff < tol)         { convcond += "maxdiff"; }
      if (itr > maxit)           { convcond += "maxitr"; }
    }

  }

  if (info) {
    return List::create(Named("omega") = wrap(Xn),
			Named("itr") = itr,
			Named("convcond") = wrap(convcond),
			Named("info") = 
			DataFrame::create(Named("f") = fhist[Range(0,itr-1)],
					  Named("taue") = tauehist[Range(0,itr-1)],
					  Named("tauo") = tauohist[Range(0,itr-1)],
					  Named("subgXn") = subgXnhist[Range(0,itr-1)],
					  Named("maxdiff") = maxdiffhist[Range(0,itr-1)],
					  Named("eigprod") = eigprodhist[Range(0,itr-1)],
					  Named("diagitr") = diagitrhist[Range(0,itr-1)],
					  Named("backitr") = backitrhist[Range(0,itr-1)]));
  } else {
    return List::create(Named("omega") = wrap(Xn),
			Named("itr") = itr);
  }
}

// another simple example: outer product of a vector, 
// returning a matrix
//
// [[Rcpp::export]]
SEXP cc_elastic_ista(const Eigen::Map<Eigen::MatrixXd> & D, double lam,
         const Eigen::Map<Eigen::MatrixXd> & pMat,
		     int penalize_diag = 0, int BBstep = 1,
		     double tol = 1e-5, int maxit = 100, int info = 0, 
		     int trkeigs = 0, int DisS = 0, double en = 1.0) {

  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;

  string convcond = "";
  NumericVector fhist(maxit+1);
  NumericVector tauehist(maxit+1); // effective tau (decreased by some power of c)
  NumericVector tauohist(maxit+1); // original tau computed by BB step
  NumericVector subgXnhist(maxit+1);
  NumericVector maxdiffhist(maxit+1);
  NumericVector eigprodhist(maxit+1);
  IntegerVector diagitrhist(maxit);
  IntegerVector backitrhist(maxit);

  int p = D.cols(); // problem size

  // standardize
  MatrixXd S;

  if (DisS) {
    S = D;
  } else {
    // standardize
    S = D.rowwise() - D.colwise().mean();                 // mean zero
    S = (S.transpose() * S) / (S.rows() - 1);
  }

  MatrixXd LambdaMat(p, p);
  for (int j=0; j<p; j++){
    for(int i=0; i<p; i++){
      if(pMat(i,j)==1){
        LambdaMat(i,j) = lam;
      }else{
        LambdaMat(i,j) = lam*100;
      }
    }
  }
  // LambdaMat.setConstant(lam);
  if (penalize_diag == 0) {
    LambdaMat.diagonal().setZero().eval();
  }

  // // initial guess is a diagonal matrix
  // SparseMatrix<double> X(p, p);      // Omega
  // X = MatrixXd (S.diagonal().asDiagonal()).sparseView();
  // for (int i=0; i < p; i++) { 
  //   X.coeffRef(i,i) += 0.5 * lam;
  // }
  // // cout << X << endl;

  SparseMatrix<double> X =  MatrixXd::Identity(p, p).sparseView();
  SparseMatrix<double> Xn(p, p);     // next omega
  SparseMatrix<double> Step(p, p);   // Omega_n - Omega

  MatrixXd W = S * X;
  MatrixXd Wn(p, p);

  MatrixXd G(p, p);
  MatrixXd Gn(p, p);
  MatrixXd subg(p, p); // computesubg here
  MatrixXd tmp(p, p);
  
  double h = (- X.diagonal().array().log().sum()) + (0.5*(X*W).trace()) + (en * pow(X.norm(),2));
  double hn = 0; 
  double Qn = 0;
  double f = 0;
  double subgnorm, Xnnorm;
  double maxdiff;

  double tau;
  double taun = 1.0;
  double c = 0.5;
  int itr = 0;
  int loop = 1;
  int diagitr = 0;
  int backitr = 0;

  G = (0.5 * (W + W.transpose())) + (en * 2.0 * MatrixXd(X));
  G += - MatrixXd((MatrixXd((1/X.diagonal().array()))).asDiagonal());

  while (loop != 0){

    tau = taun;
    // cout << "starting loop: " << tau << endl;
    
    diagitr = 0;
    backitr = 0;

    while ( 1 ) {

      if (diagitr != 0 || backitr != 0) { tau = tau * c; } // decrease tau only if needed

      tmp = MatrixXd(X) - tau*G;
      sthreshmat(tmp, tau, LambdaMat);
      Xn = tmp.sparseView();

      if (Xn.diagonal().minCoeff() < 1e-8 && diagitr < 10) { 
	// cout << "diagonal scaling" << endl;
	diagitr += 1;
	continue;
      }

      Step = Xn - X;
      Wn = S * Xn;
      Qn = h + Step.cwiseProduct(G).sum() + (1/(2*tau))*pow(Step.norm(),2);
      // hn = - Xn.diagonal().array().log().sum() + 0.5*(Xn*Wn).trace() + (1.0 * (pow(Xn.norm(),2)/pow(itr + 1.0,4)));
      hn = - Xn.diagonal().array().log().sum() + 0.5*(Xn*Wn).trace() + (en * pow(Xn.norm(),2));

      if (hn > Qn) { 
	// cout << "backtracking" << endl;
	backitr += 1;
      } else {
	// cout << "terminating" << endl;
	break;
      }

    }

    if (info) {
      tauehist(itr) = tau;  // effective tau
      tauohist(itr) = taun; // original tau
    }

    // Gn = 0.5 * (Wn + Wn.transpose()) + (1.0 * 2 * MatrixXd(Xn)/pow(itr + 1.0,4));
    Gn = 0.5 * (Wn + Wn.transpose()) + (en * 2 * MatrixXd(Xn));
    Gn += - MatrixXd((MatrixXd((1/Xn.diagonal().array()))).asDiagonal());

    if ( BBstep == 0 ) {
      taun = 1;
    } else if ( BBstep == 1 ) {
      taun = ( Step * Step ).eval().diagonal().array().sum() / (Step * ( Gn - G )).trace();
    }

    // compute subg
    tmp = MatrixXd(Xn).array().unaryExpr(ptr_fun(sgn));   // sign term
    tmp = Gn + tmp.cwiseProduct(LambdaMat);               // first term is in "tmp"
    subg = Gn;                                            // second term is in "subg"
    sthreshmat(subg, 1.0, LambdaMat);
    subg = (MatrixXd(Xn).array() != 0).select(tmp, subg); // select terms

    subgnorm = subg.norm();
    Xnnorm = Xn.norm();

    maxdiff = 0;
    for (int k=0; k<Step.outerSize(); ++k) {
      for (SparseMatrix<double>::InnerIterator it(Step,k); it; ++it) {
	maxdiff = max(abs(it.value()), maxdiff);
      }
    }

    X = Xn; h = hn; G = Gn;
    f = h + Xn.cwiseAbs().cwiseProduct(LambdaMat).sum();

    if (info) {
      fhist(itr) = f;    
      subgXnhist(itr) = subgnorm/Xnnorm;
      maxdiffhist(itr) = maxdiff;
      diagitrhist(itr) = diagitr;
      backitrhist(itr) = backitr;
      if (trkeigs) {
	eigprodhist(itr) = double((MatrixXd(Xn).eigenvalues()).prod().real());
      } else {
	eigprodhist(itr) = -1.0;
      }
    }
    itr += 1;

    // loop = int((itr < maxit) && (maxdiff > tol) && (subgnorm/Xnnorm > tol));
    // loop = int((itr < maxit) && (maxdiff > tol));
    loop = int((itr < maxit) && (subgnorm/Xnnorm > tol));

    if (info) {
      if (subgnorm/Xnnorm < tol) { convcond += "subg_2/Xn_2"; }
      if (maxdiff < tol)         { convcond += "maxdiff"; }
      if (itr > maxit)           { convcond += "maxitr"; }
    }

  }

  cout<<"Likelihood: "<<h<<endl;
  cout<<"L1 penality: "<<f-h<<endl;

  if (info) {
    return List::create(Named("omega") = wrap(Xn),
			Named("itr") = itr,
			Named("convcond") = wrap(convcond),
			Named("info") = 
			DataFrame::create(Named("f") = fhist[Range(0,itr-1)],
					  Named("taue") = tauehist[Range(0,itr-1)],
					  Named("tauo") = tauohist[Range(0,itr-1)],
					  Named("subgXn") = subgXnhist[Range(0,itr-1)],
					  Named("maxdiff") = maxdiffhist[Range(0,itr-1)],
					  Named("eigprod") = eigprodhist[Range(0,itr-1)],
					  Named("diagitr") = diagitrhist[Range(0,itr-1)],
					  Named("backitr") = backitrhist[Range(0,itr-1)]));
  } else {
    return List::create(Named("omega") = wrap(Xn),
			Named("itr") = itr);
  }
}
