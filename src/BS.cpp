/*

	ARAB Marouane,
	LUQUET Jérémy,
	RAIMBAULT Victor


*/

#ifndef DEBUG
#define DEBUG 1
#endif

#include "assert.h"


#include "BS.h"
#include <cmath>


BS::BS(double r, const PnlVect *sigma, double rho){

	/*
		mean_ = (r-sigma²)/2
	*/
	mean_=pnl_vect_copy(sigma);
	pnl_vect_mult_vect_term(mean_,sigma);
	pnl_vect_minus(mean_);
	pnl_vect_plus_double(mean_,r);
	pnl_vect_div_double(mean_,2.0);

	/*
		Matrice de cholesky
	*/
	#if DEBUG
		assert(sigma->size >0);
		
	#endif
	cholesky_ = pnl_mat_create_from_double(sigma->size,sigma->size,rho);
	pnl_mat_set_diag(cholesky_,1.0,0);
	pnl_mat_chol(cholesky_);	
	
	/*
		Sigma
	*/
	sigma_ = pnl_vect_copy(sigma);
	
};

BS::~BS(){
	pnl_vect_free(&mean_);
	pnl_vect_free(&sigma_);
	pnl_mat_free(&cholesky_);
};


void BS::asset (PnlMat *path, double T, int N, PnlRng *rng,
		PnlVect *R, PnlVect *W){

	int nt = path->n;	
	int nd = path->m;
	double dt=T/(N*1.0);

	#if DEBUG
		assert(sigma_->size == nd);
		assert(W->size == nd);
		assert(path->n == N+1);
	#endif
	int i,j;
	for(i=1;i< nt;i++){
		pnl_vect_rng_normal(R,nd,rng);
		for(j=0;j<nd;j++){
			pnl_mat_get_row(W,cholesky_,j);
			MLET(path,j,i) = MGET(path,j,i-1)*exp(
				GET(mean_,j)+
				GET(sigma_,j)*
				pnl_vect_scalar_prod(R,W));
		}
	}	
};


