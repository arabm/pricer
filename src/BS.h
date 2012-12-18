/*
	ARAB Marouane, 
	LUQUET Jérémy,
	RAIMBAULT Victor.

	BS.h

	Modification:
		

*/

#include "pnl_vector.h"
#include "pnl_matrix.h"
#include "pnl_random.h"

class BS {

public:
	PnlVect *mean_;
	PnlMat *cholesky_;
	PnlVect *sigma_;	

public:
	BS(double r, const PnlVect *sigma,
		double rho);

	~BS();


	void asset (PnlMat *path, double T, int N, PnlRng *rng,
		PnlVect *R, PnlVect *W);

};

