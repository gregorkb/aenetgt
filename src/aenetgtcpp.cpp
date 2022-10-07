# include <RcppArmadillo.h>
// [[ Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector EYgibbs(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z,
                        NumericVector se, NumericVector sp, 
                        int na, int GI) {
  int g, i, np, l, j, Zj, cj, ybar, tid, id, t;
  float pi1, pi2, pistar, sej, spj, u;
  NumericVector WW(N);
  
  for(g=0;g<GI;g++){
  for(i=0;i<N;i++){
    pi1=p(i);
    pi2=1-p(i);
    np=Y(i,1);
    for(l=0;l<np;l++){
      j=Y(i,(l+2));
      Zj=Z(j-1,0);
      cj=Z(j-1,1);
      tid=Z(j-1,2);
      sej=se(tid-1);
      spj=sp(tid-1);
      ybar=0;
      Y(i,0)=0;
      for(t=0;t<cj;t++){
        id=Z(j-1,(t+3));
        ybar=ybar+Y(id-1,0);
      }
      pi1=pi1*(sej*Zj + (1-sej)*(1-Zj));
      if(ybar > 0){
        pi2=pi2*(sej*Zj + (1-sej)*(1-Zj));
      }else{
        pi2=pi2*((1-spj)*Zj + spj*(1-Zj));
      }
    }
    pistar=(pi1/(pi1+pi2));
    u = R::runif(0,1);
//    u=rand() / double(RAND_MAX);
    if(u<pistar){
      Y(i,0)=1;
    }else{Y(i,0)=0;}
    WW(i)=WW(i)+Y(i,0);
  }}  

  return WW;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix EYiYjgibbs_slow(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z,
                        NumericVector se, NumericVector sp, 
                        int na, int GI) {
  int g, i, np, l, j, Zj, cj, ybar, tid, id, t;
  float pi1, pi2, pistar, sej, spj, u;
  NumericMatrix EYiYj(N,N);
  
  for(g=0;g<GI;g++){
  for(i=0;i<N;i++){
    pi1=p(i);
    pi2=1-p(i);
    np=Y(i,1);
    for(l=0;l<np;l++){
      j=Y(i,(l+2));
      Zj=Z(j-1,0);
      cj=Z(j-1,1);
      tid=Z(j-1,2);
      sej=se(tid-1);
      spj=sp(tid-1);
      ybar=0;
      Y(i,0)=0;
      for(t=0;t<cj;t++){
        id=Z(j-1,(t+3));
        ybar=ybar+Y(id-1,0);
      }
      pi1=pi1*(sej*Zj + (1-sej)*(1-Zj));
      if(ybar > 0){
        pi2=pi2*(sej*Zj + (1-sej)*(1-Zj));
      }else{
        pi2=pi2*((1-spj)*Zj + spj*(1-Zj));
      }
    }
    pistar=(pi1/(pi1+pi2));
    u = R::runif(0,1);
    //u=rand() / double(RAND_MAX);
    if(u<pistar){
      Y(i,0)=1;
    }else{Y(i,0)=0;}
    // WW(i)=WW(i)+Y(i,0);
    for(j=std::max(0,i - Z.nrow()*Z.nrow());j<=i;j++){
	    // If individuals are ordered according to the groups they are in,  
	    // then we only need to compute EYiYj for |i-j|< cj, where cj is the
	    // size of the group. This speeds things MASSIVELY.
	EYiYj(i,j) = EYiYj(i,j) + Y(i,0)*Y(j,0);
    }}}  

  return EYiYj;
}



// [[Rcpp::export]]
Rcpp::NumericMatrix CovYiYjgibbs(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z, NumericMatrix W,
                                 NumericVector se, NumericVector sp, NumericVector EY,
                                 int na, int GI) {
  int g, i, np, l, j, Zj, cj, ybar, tid, id, t, k;
  float pi1, pi2, pistar, sej, spj, u;
  NumericMatrix EYiYj(N,N);
  NumericMatrix CovYiYj(N,N);
  
  for(g=0;g<GI;g++){
    for(i=0;i<N;i++){
      pi1=p(i); // the value P_{\alpha,\beta}(Y_i=1|X_i)
      pi2=1-p(i);// the value P_{\alpha,\beta}(Y_i=0|X_i)
      np=Y(i,1); //number of assays in which individual i participated
      k=0;
      for(l=0;l<np;l++){
        j=Y(i,(l+2)); // number of the lth assay in which individual i participated
        Zj=Z(j-1,0); // outcome of the lth assay in which individual i participated
        cj=Z(j-1,1); // number of individuals participating in the lth assay in which individual i participated
        tid=Z(j-1,2); // an index to the sens or spec of the lth assay in which individual i participated
        sej=se(tid-1); // sensitivity of the lth assay in which individual i participated
        spj=sp(tid-1); // specificity of the same.
        ybar=0;
        Y(i,0)=0;
        for(t=0;t<cj;t++){  
          id=Z(j-1,(t+3)); // number of the tth individual in the lth assay in which individual i participated
          ybar=ybar+Y(id-1,0); // ybar > 0 gives the true disease status of the group of individuals in the lth assay of individual i
        }
        pi1=pi1*(sej*Zj + (1-sej)*(1-Zj)); // product of independent assay probabilities, given Y_i = 1 (group is positive), times P_{\alpha,\beta}(Y_i=1|X_i)
        if(ybar > 0){// since Y(i,0)=0, ybar > 0 only if OTHER individuals in the group are positive.
          pi2=pi2*(sej*Zj + (1-sej)*(1-Zj)); // product of independent assay probabilities, given Y_i = 0 but group is positive, times P_{\alpha,\beta}(Y_i=0|X_i)
        }else{
          pi2=pi2*((1-spj)*Zj + spj*(1-Zj));// product of independent assay probabilities, given Y_i = 0 and group is negative, times P_{\alpha,\beta}(Y_i=0|X_i)
        }
      }
      pistar=(pi1/(pi1+pi2));// probability that Y_i = 1 given all other disease statuses and outcomes of all assays in which individual i participated
      u = R::runif(0,1);
      //u=rand() / double(RAND_MAX);
      if(u<pistar){
        Y(i,0)=1;
      }else{Y(i,0)=0;}
      k = W(i,0);
      for(j=0;j<k;j++){
        t = W(i,1+j)-1;
        EYiYj(i,t) = EYiYj(i,t) + Y(i,0)*Y(t,0);
      }}}
  
  // Now compute the covariance matrix
  for(i=0;i<N;i++){
    k = W(i,0);
    for(j=0;j<k;j++){
      t = W(i,1+j)-1;
      CovYiYj(i,t) = EYiYj(i,t)/GI - EY(i)*EY(t) ;
    }
  }
  
  return CovYiYj;
}

//' Compute the elastic net estimator for logistic regression
//' 
//' @param Yr Response vector of 1s and 0s
//' @param Xr A design matrix with the first column a column of 1s
//' @param lambda The tuning parameter governing the strength of the elastic net penalty
//' @param gammar A vector of length \code{ncol(X) - 1} giving the weights applied to each covariate in the elastic net penalization
//' @param theta Value controlling the relative strength of the ridge and lasso penalties; 1 gives lasso.
//' @param tol Convergence tolerance
//' @return a list with the estimated coefficients, etc.
//' 
//' @examples
//' # generate some data
//' n <- 5000
//' p <- 40
//' b <- c(0,3,0,1,-2,0,rep(0,p-5)) # first is intercept
//' X <- cbind(rep(1,n),scale(matrix(rnorm(n*p),nrow=n),TRUE,TRUE))
//' eta <- X %*% b
//' Y <- rbinom(n,1,1/(1 + exp(-eta)))
//' 
//' # compute elastic net estimator  
//' logistic_enet(Y, X, lambda = 30, gammar = rep(1,p), theta = 0.5, tol = 0.0001)$b
// [[Rcpp::export()]]
Rcpp::List logistic_enet(Rcpp::NumericVector Yr, 
                         Rcpp::NumericMatrix Xr,
                         float lambda,
                         Rcpp::NumericVector gammar,
                         float theta,
                         float tol){
                           
		  
		  int n = Xr.nrow(), p = Xr.ncol()-1, i, k;
		  float uj, vj, wj, sj;
		  bool conv00, conv0;
		  
		  arma::mat X(Xr.begin(), n, p+1, false); 
		  arma::colvec Y(Yr.begin(),Yr.size(), false);
		  arma::colvec gamma(gammar.begin(),gammar.size(),false);
		  
		  arma::colvec b = arma::zeros(p+1);
		  arma::colvec b0(p+1);
		  arma::colvec b00(p+1);
		  arma::colvec diff(p+1);
		  
		  arma::colvec eta = arma::zeros(n);
		  arma::colvec pr(n);
		  arma::colvec w(n);
		  arma::colvec z(n);
		  arma::colvec zj(n);
		  
		  i = 0;
		  conv00 = false;
		  while( (i < 500) & (conv00 == false)){
		    
		    b00 = b;
		    
		    pr = 1 / (1 + exp(-eta));
		    w = pr % (1 - pr);
		    z = eta + (Y - pr) / w;
		    
		    k = 0;
		    conv0 = false;
		    while( (k < 500) & (conv0 == false)){
		      
		      b0 = b; 
		      
		      b(0) = sum( w % ( z - (eta - X.col(0) * b(0) ))) / sum(w);
		      
		      eta = eta + X.col(0) * ( b(0) - b0(0) );
		      
		      for(int j=1; j < (p+1) ; j++){
		        
		        zj = z - (eta - X.col(j) * b(j));
		        uj = sum( w % zj % X.col(j) );
		        vj = theta * lambda * arma::as_scalar(gamma(j-1));
		        wj = sum( w % pow(X.col(j),2)) + lambda * (1 - theta); 
		        
		        // soft-threshold
		        if( (uj > 0) & (vj < std::abs(uj)) ){ 
		          
		          sj = uj - vj;
		          
		        } else if( (uj < 0) & (vj < std::abs(uj))){
		          
		          sj = uj + vj;
		          
		        } else {
		          
		          sj = 0;
		          
		        }
		        
		        b(j) = sj / wj;
		        eta = eta + X.col(j) * ( b(j) - b0(j) );
		        
		      }
		      
		      diff = abs(b - b0);
		      conv0 = diff.max() < tol;
		      k++;
		      
		    }
		    
		    diff = abs(b - b00);
		    conv00 = diff.max() < tol;
		    i++;
		  
		  }
		  
		  if(b.has_nan()) {
		    
		    Rcpp::Rcout << "warning: failure to converge due to complete or quasi-complete separation" << std::endl;
		    
		  }
		  
		  return Rcpp::List::create(Named("b") = b,
                              Named("lambda") = lambda,
                              Named("theta") = theta,
                              Named("gamma") = gamma,
                              Named("eta") = eta);
		}

// [[Rcpp::export()]]
Rcpp::List llj_array(	Rcpp::IntegerVector Zjr, 
			Rcpp::IntegerVector Zjc,
			Rcpp::IntegerVector Yji,
			Rcpp::IntegerVector whichjretest,
			Rcpp::NumericVector pxji,
			Rcpp::NumericVector Se,
			Rcpp::NumericVector Sp,
			int B){

arma::ivec Zjrow(Zjr.begin(),Zjr.size(),false);
arma::ivec Zjcol(Zjc.begin(),Zjc.size(),false);
arma::ivec Yj_retested(Yji.begin(),Yji.size(),false);			
arma::ivec retestint(whichjretest.begin(), whichjretest.size(),false);
arma::uvec retest = arma::conv_to<arma::uvec>::from(retestint);	
retest = retest - 1;
arma::colvec px(pxji.begin(),pxji.size(),false);
arma::colvec Sens(Se.begin(),Se.size(),false);
arma::colvec Spec(Sp.begin(),Sp.size(),false);				


int cj = Zjr.size(),
	M = whichjretest.size();

arma::uword cjsq = cj*cj;

float pZjrZjcYji_tYj, pZjrow_tYj, pZjcol_tYj, pYji_tYj, p1, p0, llj, Lj;

arma::colvec U(cjsq);
arma::vec tYj = arma::zeros(cjsq);			
arma::vec tYj_retested(M);
				
arma::umat Array(cj,cj);
arma::uvec tZjrow(cj);
arma::uvec tZjcol(cj);

Lj = 0;

for(int b=0;b<B;b++)								
{
	
	U.randu();
	
	for(arma::uword i=0;i < cjsq ;i++)
	{ 	
		if(U(i) <= px(i)) 
		{
			tYj(i) = 1;
			
		} else {
			
			tYj(i) = 0 ; 
		}
	
	}
	
	for(int row=0;row<cj;row++)
		for(int col = 0 ; col < cj ; col++)
			{
				
				Array(row,col) = tYj( col*cj + row);
				
			}
								
	pZjrow_tYj = 1;
	pZjcol_tYj = 1;
	
	for(int k=0;k<cj;k++)
	{
		
		tZjrow(k) = max(Array.row(k));
		tZjcol(k) = max(Array.col(k));
		
		p0 = (Sens(0)*Zjrow(k) + (1-Sens(0))*(1 - Zjrow(k))) * tZjrow(k) + ((1-Spec(0))*Zjrow(k) + Spec(0)*(1 - Zjrow(k))) * (1-tZjrow(k));
		p1 = (Sens(0)*Zjcol(k) + (1-Sens(0))*(1 - Zjcol(k))) * tZjcol(k) + ((1-Spec(0))*Zjcol(k) + Spec(0)*(1 - Zjcol(k))) * (1-tZjcol(k));
		
		pZjrow_tYj = pZjrow_tYj * p0 ;
		pZjcol_tYj = pZjcol_tYj * p1	 ;					
		
	}
		
	tYj_retested = tYj.rows(retest);
	
	pYji_tYj = 1;

	for(arma::uword i=0;i<tYj_retested.n_elem ;i++)
	{
		
		p1 = Sens(1)*Yj_retested(i) + (1-Sens(1))*(1-Yj_retested(i));
		p0 = (1-Spec(1))*Yj_retested(i)+ Spec(1)*(1-Yj_retested(i));
		
		pYji_tYj = pYji_tYj *  ( p1 *  tYj_retested(i) + p0 * ( 1 -  tYj_retested(i)) ) ;
				
	}
	
	
		pZjrZjcYji_tYj = pZjrow_tYj * pZjcol_tYj * pYji_tYj	;						
	
		Lj = Lj + pZjrZjcYji_tYj ;
						
}	

llj = log(Lj/B);
				
return Rcpp::List::create(Named("retest") = retest,
					Named("Yj_retested") = Yj_retested,
					Named("tYj_retested") = tYj_retested,
					Named("tYj") = tYj,
					Named("Se") = Sens,
					Named("Sp") = Spec,
					Named("Array") = Array,
					Named("tZjrow") = tZjrow,
					Named("tZjcol") = tZjcol,
					Named("pYji_tYj") = pYji_tYj,
					Named("pZjrZjcYji_tYj") = pZjrZjcYji_tYj,
					Named("Lj") = Lj,
					Named("llj") = llj,
					Named("U") = U
					);
			}
			

			
//' Generate all possible sequences of 0s and 1s of a given length
//' @param a the length of the sequences.
//' @return a matrix containing in its rows the sequences of 0s and 1s.
// [[Rcpp::export]]
arma::mat all_binary_sequences(int a){
  
  int n_rows = (int)(pow(2,a) + .5);
  int i,n,m;
  arma::mat M = arma::zeros(n_rows,a);
  
  // count from 0 to 2^a - 1 in binary
  for( n = 0; n < n_rows ; n++ ){
    
    m = n;
    i = 0;
    while( m > 0){
      
      M(n,i) = m % 2;    
      m = m / 2;  
      i++;
      
    }  
    
  }
  
  return M;
  
}

//' Computes conditional expectations of individual disease statuses for individual, master pool, or Dorfman testing
//'   
//' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}.
//' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}.
//' @param eta the value of the linear predictor
//' @param Se A vector of testing sensitivities of length \code{max(Z[,3])}.
//' @param Sp A vector of testing specificities of length \code{max(Z[,3])}.
//' @return The vector of conditional expectations.
//' 
//' This function computes the conditional expectations of each individual disease status given the observed assay data.
//' 
//' @examples
//' # generate individual covariate values and disease statuses
//' N <- 100
//' data <- model1(N)
//' X <- data$X
//' Y.true <- data$Yi
//' Se <- c(.95,.92) # set master pool and individual assay sensitivity
//' Sp <- c(.97,.98) # set master pool and individual assay specificity
//' cj <- 4 # set size of master pools
//' # subject individuals to Dorfman testing
//' assay.data <- dorfman.assay.gen(Y.true,Se,Sp,cj)
//' Z <- assay.data$Z
//' Y <- assay.data$Y
//' b <- data$b
//' eta <- X %*% b
//' EY <- EYexact(Z,Y,eta,Se,Sp)
// [[Rcpp::export]]
arma::colvec EYexact(IntegerMatrix Z, IntegerMatrix Y, NumericVector eta, NumericVector Se, NumericVector Sp){
  
  int n = Y.nrow(), max_c = Z.ncol() - 3;
  int assay_in_which_involved;
  int group_assay;
  int i, j, k, l, cj, zj;
  int which_SeSp;
  int n_left,n_group;
  
  double Sej, Spj, prd, max_Y;
  double prob_whole_group_negative;
  
  // arma versions of inputs
  arma::imat aY(Y.begin(),n,4,false);
  arma::imat aZ(Z.begin(),Z.nrow(),Z.ncol(),false);
  arma::colvec aeta(eta.begin(),n,false);
  
  // initialize some armadillo vectors
  arma::ivec group_long(max_c);
  arma::ivec U(max_c);
  arma::uvec individual_assays(max_c);
  arma::uvec group(max_c);
  arma::uvec to_drop(1);
  arma::colvec Sej_i(max_c);
  arma::colvec Spj_i(max_c);
  arma::colvec p_group(max_c);
  arma::colvec Aj(max_c);
  arma::colvec Bj(max_c);
  arma::colvec EY(n);
  arma::colvec aa(max_c);
  arma::colvec bb(max_c);
  arma::colvec cc(max_c);
  arma::colvec aa_bb_cc(max_c);
  arma::colvec Y_no_i(max_c-1);
  arma::umat cj_by_cj(max_c-1,max_c);
  arma::uvec k_uvec(1);
  
  // get the fitted probs
  arma::colvec p_all = 1 / ( 1 + exp( - eta )) ;
  
  // make lists of individuals tested in a single assay and in two assays
  arma::uvec involved_in_only_one_assay = arma::find( aY.col(1) == 1 );
  arma::uvec involved_in_two_assays = arma::find( aY.col(1) == 2 );
  
  // generate matrices with rows giving all sequences of 0s and 1s of lengths cj
  arma::field<arma::mat> YY(max_c - 1);
  for(i = 0; i < max_c - 1; i++){
    
    YY(i) = all_binary_sequences(i + 2);
    
  }
  
  //----------------------------------------------------------
  // handle individual, masterpool, and negative dorfman pools
  //----------------------------------------------------------
  
  n_left = involved_in_only_one_assay.size();
  while(n_left > 0){
    
    // in which assay was he/she involved?
    assay_in_which_involved = aY(involved_in_only_one_assay[0],2) - 1;
    
    // all involved in this assay
    group_long = aZ(assay_in_which_involved,arma::span(3,3 + max_c - 1)).t() - 1;
    group.set_size(max_c);
    cj = 0;
    for(i = 0; i < max_c ; i++){
      
      if(group_long(i) == -100) break;
      
      group(i) = group_long(i);
      cj++;
      
    }
    
    // set size of some arma vectors
    group.set_size(cj);
    p_group.set_size(cj);
    Aj.set_size(cj);
    Bj.set_size(cj);
    
    // get probs for the group
    p_group = p_all(group);
    
    // get Se and Sp for this assay
    which_SeSp = aZ(assay_in_which_involved,2) - 1;
    Sej = Se[which_SeSp];
    Spj = Sp[which_SeSp];
    
    // prepare to compute the group probabilities
    prob_whole_group_negative = prod(1 - p_group);
    zj = aZ(assay_in_which_involved,0);
    Aj = arma::ones(cj,1) * (Sej*zj + (1-Sej)*(1-zj));
    Bj = ((1-Spj)*zj + Spj*(1-zj) - Aj ) * prob_whole_group_negative / (1 - p_group) + Aj;
    
    // compute group probabilities
    EY(group) = Aj % p_group /( Aj % p_group + Bj % (1 - p_group) );
    
    // remove the "processed" individuals from the waiting list
    n_group = group.size();
    
    for(i = 0; i < n_group; i++){
      
      to_drop = find(involved_in_only_one_assay == group(i));
      involved_in_only_one_assay.shed_row(to_drop[0]);
      
    }
    
    n_left = involved_in_only_one_assay.size();
    
  }
  
  //--------------------------------------------
  // handle positive pools under dorfman testing
  //--------------------------------------------
  
  n_left = involved_in_two_assays.size();
  while(n_left > 0){
    
    // in which assays was he/she involved? ASSUMING THE GROUP ASSAY COMES FIRST!!! 
    group_assay = aY(involved_in_two_assays[0],2) - 1;
    
    // all involved in this group : ASSUMING THE GROUP ASSAY COMES FIRST!!! 
    // Which it should, unless the data are encoded in the Z matrix in a different way.
    group_long = aZ(group_assay,arma::span(3,3 + max_c - 1)).t() - 1;
    group.set_size(max_c);
    cj = 0;
    for(i = 0; i < max_c ; i++){
      
      if(group_long(i) == -100) break;
      
      group(i) = group_long(i);
      cj++;
      
    }
    
    // reset size of some arma vectors/matrices
    group.set_size(cj);
    p_group.set_size(cj);
    Aj.set_size(cj);
    Bj.set_size(cj);
    individual_assays.set_size(cj);
    Sej_i.set_size(cj);
    Spj_i.set_size(cj);
    U.set_size(cj);
    aa.set_size(cj);
    bb.set_size(cj);
    cc.set_size(cj);
    aa_bb_cc.set_size(cj);
    cj_by_cj.set_size(cj-1,cj);
    U.set_size(cj);
    
    // get probs for the group
    p_group = p_all(group);
    
    // get Se and Sp for the group assay:
    which_SeSp = aZ(group_assay,2) - 1;
    Sej = Se[which_SeSp];
    Spj = Sp[which_SeSp];
    
    // get individual assays and Se and Sp for individual assays
    for( i = 0; i < cj ; i ++){
      
      individual_assays(i) = aY(group(i),3) - 1;
      which_SeSp = aZ(individual_assays(i),2) - 1;
      Sej_i(i) = Se[which_SeSp];
      Spj_i(i) = Sp[which_SeSp];
      U(i) = aZ(individual_assays(i),0);
      
    }
    
    Aj = arma::zeros(cj);
    Bj = arma::zeros(cj);
    
    for(k = 0; k < YY(cj-2).n_rows ; k++){
      
      k_uvec(0) = k;
      
      int yi;
      
      for(i = 0; i < cj ; i++){
        
        yi = YY(cj-2)(k,i);
        
        if( yi == 1){
          
          aa(i) = (Sej_i(i) * U(i) + ( 1 - Sej_i(i)) * (1-U(i)));
          bb(i) = 1;
          cc(i) = p_group(i);
          
        } else {
          
          aa(i) = 1;
          bb(i) = ((1-Spj_i(i))*U(i) + Spj_i(i)*(1-U(i)));
          cc(i) = 1 - p_group(i);
        }
        
        l = 0;
        j = 0;
        while(l < cj - 1){
          
          cj_by_cj(l,i) = j;
          
          if(j != i){
            l++;
          }
          
          j++;
          
        }
        
      }
      
      aa_bb_cc = aa % bb % cc;
      
      for(i = 0; i < cj; i++){
        
        prd = prod(aa_bb_cc(cj_by_cj.col(i)));
        Y_no_i = YY(cj-2)(k_uvec,cj_by_cj.col(i)).t();
        max_Y = max(Y_no_i);
        
        Aj(i) = Aj(i) + prd;
        Bj(i) = Bj(i) + (Sej*max_Y + (1-Spj)*(1 - max_Y)) * prd;
        
      }
      
    }
    
    Aj = Sej * (Sej_i % U + ( 1 - Sej_i) % (1 - U) ) % Aj;
    Bj = ( (1-Spj_i) % U + Spj_i % (1-U) ) % Bj;
    
    EY(group) = Aj % p_group /( Aj % p_group + Bj % (1 - p_group) );
    
    // remove the "processed" individuals from the waiting list
    n_group = group.size();
    
    for(i = 0; i < n_group; i++){
      
      to_drop = find(involved_in_two_assays == group(i));
      involved_in_two_assays.shed_row(to_drop[0]);
      
    }
    
    n_left = involved_in_two_assays.size();
    
  }
  
  return EY;
  
}


// //' Computes conditional expectations of individual disease statuses for individual, master pool, or Dorfman testing
// //'   
// //' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}.
// //' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}.
// //' @param X Design matrix with first column a column of 1s.
// //' @param b Parameter values at which to compute the conditional expectations.
// //' @param Se A vector of testing sensitivities of length \code{max(Z[,3])}.
// //' @param Sp A vector of testing specificities of length \code{max(Z[,3])}.
// //' @return The vector of conditional expectations.
// //' 
// //' This function computes the conditional expectations of each individual disease status given the observed assay data.
// //' 
// //' @examples
// //' # generate individual covariate values and disease statuses
// //' N <- 100
// //' data <- model1(N)
// //' X <- data$X
// //' Y.true <- data$Yi
// //' Se <- c(.95,.92) # set master pool and individual assay sensitivity
// //' Sp <- c(.97,.98) # set master pool and individual assay specificity
// //' cj <- 4 # set size of master pools
// //' # subject individuals to Dorfman testing
// //' assay.data <- dorfman.assay.gen(Y.true,Se,Sp,cj)
// //' Z <- assay.data$Z
// //' Y <- assay.data$Y
// //' b <- data$b
// //'   
// //' EY <- EYexact(Z,Y,X,b,Se,Sp)
// // [[Rcpp::export]]
// arma::colvec EYexact(IntegerMatrix Z, IntegerMatrix Y, NumericMatrix X, NumericVector b, NumericVector Se, NumericVector Sp){
//   
//   int n = X.nrow(), p = X.ncol(), max_c = Z.ncol() - 3;
//   int assay_in_which_involved;
//   int group_assay;
//   int i, j, k, l, cj, zj;
//   int which_SeSp;
//   int n_left,n_group;
//   
//   double Sej, Spj, prd, max_Y;
//   double prob_whole_group_negative;
//   
//   // arma versions of inputs
//   arma::imat aY(Y.begin(),n,4,false);
//   arma::imat aZ(Z.begin(),Z.nrow(),Z.ncol(),false);
//   
//   // initialize some armadillo vectors
//   arma::ivec group_long(max_c);
//   arma::ivec U(max_c);
//   arma::uvec individual_assays(max_c);
//   arma::uvec group(max_c);
//   arma::uvec to_drop(1);
//   arma::colvec Sej_i(max_c);
//   arma::colvec Spj_i(max_c);
//   arma::colvec p_group(max_c);
//   arma::colvec Aj(max_c);
//   arma::colvec Bj(max_c);
//   arma::colvec EY(n);
//   arma::colvec aa(max_c);
//   arma::colvec bb(max_c);
//   arma::colvec cc(max_c);
//   arma::colvec aa_bb_cc(max_c);
//   arma::colvec Y_no_i(max_c-1);
//   arma::umat cj_by_cj(max_c-1,max_c);
//   arma::uvec k_uvec(1);
//   
//   // get the fitted probs
//   arma::mat aX(X.begin(),n,p,false);
//   arma::colvec ab(b.begin(),p,false);
//   arma::colvec p_all = 1 / ( 1 + exp( - aX * ab )) ;
//   
//   // make lists of individuals tested in a single assay and in two assays
//   arma::uvec involved_in_only_one_assay = arma::find( aY.col(1) == 1 );
//   arma::uvec involved_in_two_assays = arma::find( aY.col(1) == 2 );
//   
//   // generate matrices with rows giving all sequences of 0s and 1s of lengths cj
//   arma::field<arma::mat> YY(max_c - 1);
//   for(i = 0; i < max_c - 1; i++){
//     
//     YY(i) = all_binary_sequences(i + 2);
//     
//   }
//   
//   //----------------------------------------------------------
//   // handle individual, masterpool, and negative dorfman pools
//   //----------------------------------------------------------
//   
//   n_left = involved_in_only_one_assay.size();
//   while(n_left > 0){
//     
//     // in which assay was he/she involved?
//     assay_in_which_involved = aY(involved_in_only_one_assay[0],2) - 1;
//     
//     // all involved in this assay
//     group_long = aZ(assay_in_which_involved,arma::span(3,3 + max_c - 1)).t() - 1;
//     group.set_size(max_c);
//     cj = 0;
//     for(i = 0; i < max_c ; i++){
//       
//       if(group_long(i) == -100) break;
//       
//       group(i) = group_long(i);
//       cj++;
//       
//     }
//     
//     // set size of some arma vectors
//     group.set_size(cj);
//     p_group.set_size(cj);
//     Aj.set_size(cj);
//     Bj.set_size(cj);
//     
//     // get probs for the group
//     p_group = p_all(group);
//     
//     // get Se and Sp for this assay
//     which_SeSp = aZ(assay_in_which_involved,2) - 1;
//     Sej = Se[which_SeSp];
//     Spj = Sp[which_SeSp];
//     
//     // prepare to compute the group probabilities
//     prob_whole_group_negative = prod(1 - p_group);
//     zj = aZ(assay_in_which_involved,0);
//     Aj = arma::ones(cj,1) * (Sej*zj + (1-Sej)*(1-zj));
//     Bj = ((1-Spj)*zj + Spj*(1-zj) - Aj ) * prob_whole_group_negative / (1 - p_group) + Aj;
//     
//     // compute group probabilities
//     EY(group) = Aj % p_group /( Aj % p_group + Bj % (1 - p_group) );
//     
//     // remove the "processed" individuals from the waiting list
//     n_group = group.size();
//     
//     for(i = 0; i < n_group; i++){
//       
//       to_drop = find(involved_in_only_one_assay == group(i));
//       involved_in_only_one_assay.shed_row(to_drop[0]);
//       
//     }
//     
//     n_left = involved_in_only_one_assay.size();
//     
//   }
//   
//   //--------------------------------------------
//   // handle positive pools under dorfman testing
//   //--------------------------------------------
//   
//   n_left = involved_in_two_assays.size();
//   while(n_left > 0){
//     
//     // in which assays was he/she involved? ASSUMING THE GROUP ASSAY COMES FIRST!!! 
//     group_assay = aY(involved_in_two_assays[0],2) - 1;
//     
//     // all involved in this group : ASSUMING THE GROUP ASSAY COMES FIRST!!! 
//     // Which it should, unless the data are encoded in the Z matrix in a different way.
//     group_long = aZ(group_assay,arma::span(3,3 + max_c - 1)).t() - 1;
//     group.set_size(max_c);
//     cj = 0;
//     for(i = 0; i < max_c ; i++){
//       
//       if(group_long(i) == -100) break;
//       
//       group(i) = group_long(i);
//       cj++;
//       
//     }
//     
//     // reset size of some arma vectors/matrices
//     group.set_size(cj);
//     p_group.set_size(cj);
//     Aj.set_size(cj);
//     Bj.set_size(cj);
//     individual_assays.set_size(cj);
//     Sej_i.set_size(cj);
//     Spj_i.set_size(cj);
//     U.set_size(cj);
//     aa.set_size(cj);
//     bb.set_size(cj);
//     cc.set_size(cj);
//     aa_bb_cc.set_size(cj);
//     cj_by_cj.set_size(cj-1,cj);
//     U.set_size(cj);
//     
//     // get probs for the group
//     p_group = p_all(group);
//     
//     // get Se and Sp for the group assay:
//     which_SeSp = aZ(group_assay,2) - 1;
//     Sej = Se[which_SeSp];
//     Spj = Sp[which_SeSp];
//     
//     // get individual assays and Se and Sp for individual assays
//     for( i = 0; i < cj ; i ++){
//       
//       individual_assays(i) = aY(group(i),3) - 1;
//       which_SeSp = aZ(individual_assays(i),2) - 1;
//       Sej_i(i) = Se[which_SeSp];
//       Spj_i(i) = Sp[which_SeSp];
//       U(i) = aZ(individual_assays(i),0);
//       
//     }
//     
//     Aj = arma::zeros(cj);
//     Bj = arma::zeros(cj);
//     
//     for(k = 0; k < YY(cj-2).n_rows ; k++){
//       
//       k_uvec(0) = k;
//       
//       int yi;
//       
//       for(i = 0; i < cj ; i++){
//         
//         yi = YY(cj-2)(k,i);
//         
//         if( yi == 1){
//           
//           aa(i) = (Sej_i(i) * U(i) + ( 1 - Sej_i(i)) * (1-U(i)));
//           bb(i) = 1;
//           cc(i) = p_group(i);
//           
//         } else {
//           
//           aa(i) = 1;
//           bb(i) = ((1-Spj_i(i))*U(i) + Spj_i(i)*(1-U(i)));
//           cc(i) = 1 - p_group(i);
//         }
//         
//         l = 0;
//         j = 0;
//         while(l < cj - 1){
//           
//           cj_by_cj(l,i) = j;
//           
//           if(j != i){
//             l++;
//           }
//           
//           j++;
//           
//         }
//         
//       }
//       
//       aa_bb_cc = aa % bb % cc;
//       
//       for(i = 0; i < cj; i++){
//         
//         prd = prod(aa_bb_cc(cj_by_cj.col(i)));
//         Y_no_i = YY(cj-2)(k_uvec,cj_by_cj.col(i)).t();
//         max_Y = max(Y_no_i);
//         
//         Aj(i) = Aj(i) + prd;
//         Bj(i) = Bj(i) + (Sej*max_Y + (1-Spj)*(1 - max_Y)) * prd;
//         
//       }
//       
//     }
//     
//     Aj = Sej * (Sej_i % U + ( 1 - Sej_i) % (1 - U) ) % Aj;
//     Bj = ( (1-Spj_i) % U + Spj_i % (1-U) ) % Bj;
//     
//     EY(group) = Aj % p_group /( Aj % p_group + Bj % (1 - p_group) );
//     
//     // remove the "processed" individuals from the waiting list
//     n_group = group.size();
//     
//     for(i = 0; i < n_group; i++){
//       
//       to_drop = find(involved_in_two_assays == group(i));
//       involved_in_two_assays.shed_row(to_drop[0]);
//       
//     }
//     
//     n_left = involved_in_two_assays.size();
//     
//   }
//   
//   return EY;
//   
// }
