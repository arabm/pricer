/*

	ARAB Marouane,
	LUQUET Jérémy,
	RAIMBAULT Victor.

	Teide

*/

#include "BS.h"

#include <iostream>


#include "time.h"

using namespace std;

int main(){
	
	int nd = 3;
	int nt = 4;

	double arr[] = {0.2,0.5,0.3};
	PnlVect *S = pnl_vect_create_from_ptr(nd,arr);
	double r =0.2;
	double dt = 4;
	double rho = 0.5;
	double T = 2.0;
	BS bs(r,S,rho);


	cout << "Cholesky" << endl;
	pnl_mat_print(bs.cholesky_);

	cout << "mean_" << endl;
	pnl_vect_print(bs.mean_);


	PnlVect *R = pnl_vect_create(nd);
	PnlVect *W = pnl_vect_create(nd);
	PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
	pnl_rng_sseed(rng, time(NULL));

	PnlMat *path = pnl_mat_create_from_double(nd,nt+1,100.0);	
	bs.asset(path,T,nt,rng,R,W);	
	
	cout << "Path" << endl;
	pnl_mat_print(path);
	
	pnl_vect_free(&R);
	pnl_vect_free(&W);
	pnl_rng_free(&rng);
	pnl_mat_free(&path);	
	return 0;
}
