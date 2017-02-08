/*
 * @file matrixes.h
 *
 * @author Anna Subbotina
 */

#ifndef _SLB_MATRIXES
#define _SLB_MATRIXES

// porosity
#define POR 0.238

extern "C" {
	// LU decomoposition of a general matrix
	void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

	// generate inverse of a matrix given its LU decomposition
	void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

static void ll_inverse(double* A, int N)
{
	int *IPIV = new int[N+1];
	int LWORK = N*N;
	double *WORK = new double[LWORK];
	int INFO;

	dgetrf_(&N,&N,A,&N,IPIV,&INFO);
	dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

	delete IPIV;
	delete WORK;
}

// coefficients for the problem of the single phase isothermal flow transport
// Al(x) dp(t,x)/dt - d(K(x) dp(t,x)/dx)/dx + D(x) p(t,x) = 0

static void get_Al(vector<mat> &Al_arr)
{
	int n = Al_arr[0].n_rows;
	int l = Al_arr.size();
	mat Al(n, n, fill::eye);

	for (int i = 0; i < l; i++)
		Al_arr[i] = Al;
}

static void get_K(vector<mat> &K_arr)
{
	int n = K_arr[0].n_rows;
	int l = K_arr.size();

	mat K(n, n, fill::zeros);
	mat K0(n, n, fill::zeros);

	K(0, 0) = 700.7; K(1, 0) = -90.5; K(2, 0) = -51.8; K(3, 0) = -20.1;
	K(4, 0) = -16.7; K(5, 0) = -36.7; K(6, 0) = -12.6; K(7, 0) = -11.1;

	K(0, 1) = -90.5; K(1, 1) = 410.5; K(2, 1) = 10.8; K(3, 1) = -6.0;
	K(4, 1) = 7.8; K(5, 1) = 4.9; K(6, 1) = 13.6; K(7, 1) = 6.6;

	K(0, 2) = -51.8; K(1, 2) = 10.8; K(2, 2) = 162.8; K(3, 2) = 24.6;
	K(4, 2) = 12.9; K(5, 2) = 16.8; K(6, 2) = 16.6; K(7, 2) = 20.7;

	K(0, 3) = -20.1; K(1, 3) = -6.0; K(2, 3) = 24.6; K(3, 3) = 110.0;
	K(4, 3) = 7.2; K(5, 3) = 14.9; K(6, 3) = 10.6; K(7, 3) = 18.4;

	K(0, 4) = -16.7; K(1, 4) = 7.8; K(2, 4) = 12.9; K(3, 4) = 7.2;
	K(4, 4) = 86.0; K(5, 4) = 10.6; K(6, 4) = 10.8; K(7, 4) = 13.1;

	K(0, 5) = -36.7; K(1, 5) = 4.9; K(2, 5) = 16.8; K(3, 5) = 14.9;
	K(4, 5) = 10.6; K(5, 5) = 92.0; K(6, 5) = 13.6; K(7, 5) = 18.6;

	K(0, 6) = -12.6; K(1, 6) = 13.6; K(2, 6) = 16.6; K(3, 6) = 10.6;
	K(4, 6) = 10.8; K(5, 6) = 13.6; K(6, 6) = 85.8; K(7, 6) = 23.8;

	K(0, 7) = -11.1; K(1, 7) = 6.6; K(2, 7) = 20.7; K(3, 7) = 18.4;
	K(4, 7) = 13.1; K(5, 7) = 18.6; K(6, 7) = 23.8; K(7, 7) = 68.8;

	K0(0, 0) = 700.7; K0(1, 0) = -90.5; K0(2, 0) = -51.8; K0(3, 0) = -20.1;
	K0(4, 0) = -16.7; K0(5, 0) = -36.7; K0(6, 0) = 0.00; K0(7, 0) = 0.00;

	K0(0, 1) = -90.5; K0(1, 1) = 410.5; K0(2, 1) = 10.8; K0(3, 1) = -6.0;
	K0(4, 1) = 7.8; K0(5, 1) = 4.9; K0(6, 1) = 0.00; K0(7, 1) = 0.00;

	K0(0, 2) = -51.8; K0(1, 2) = 10.8; K0(2, 2) = 162.8; K0(3, 2) = 24.6;
	K0(4, 2) = 12.9; K0(5, 2) = 16.8; K0(6, 2) = 0.00; K0(7, 2) = 0.00;

	K0(0, 3) = -20.1; K0(1, 3) = -6.0; K0(2, 3) = 24.6; K0(3, 3) = 110.1;
	K0(4, 3) = 7.2; K0(5, 3) = 14.9; K0(6, 3) = 0.00; K0(7, 3) = 0.00;

	K0(0, 4) = -16.7; K0(1, 4) = 7.8; K0(2, 4) = 12.9; K0(3, 4) = 7.2;
	K0(4, 4) = 86.0; K0(5, 4) = 10.6; K0(6, 4) = 0.00; K0(7, 4) = 0.00;

	K0(0, 5) = -36.7; K0(1, 5) = 4.9; K0(2, 5) = 16.8; K0(3, 5) = 14.9;
	K0(4, 5) = 10.6; K0(5, 5) = 92.0; K0(6, 5) = 0.00; K0(7, 5) = 0.00;

	K0(0, 6) = 0.00; K0(1, 6) = 0.00; K0(2, 6) = 0.00; K0(3, 6) = 0.00;
	K0(4, 6) = 0.00; K0(5, 6) = 0.00; K0(6, 6) = 85.8; K0(7, 6) = 0.00;

	K0(0, 7) = 0.00; K0(1, 7) = 0.00; K0(2, 7) = 0.00; K0(3, 7) = 0.00;
	K0(4, 7) = 0.00; K0(5, 7) = 0.00; K0(6, 7) = 0.00; K0(7, 7) = 68.8;


	for (int i = 0; i < l; i++)
		K_arr[i] = K * POR;

/*
	for (int i = 0; i < l* 0.2; i++)
		K_arr[i] = K * POR;
	for (int i = l * 0.2; i < l * 0.8; i++)
		K_arr[i] = K0 * POR;
	for (int i = l * 0.8; i < l; i++)
		K_arr[i] = K * POR;
*/
}

static void get_D(vector<mat> &D_arr)
{
	int n = D_arr[0].n_rows;
	int l = D_arr.size();
	mat D(n, n, fill::zeros);
	mat D0(n, n, fill::zeros);

	D(0, 0) = 58.0; D(1, 0) = -13.5; D(2, 0) = -6.0; D(3, 0) = -5.1;
	D(4, 0) = -5.4; D(5, 0) = -7.2; D(6, 0) = -9.7; D(7, 0) = -11.1;

	D(0, 1) = -13.5; D(1, 1) = 46.0; D(2, 1) = -4.0; D(3, 1) = -3.8;
	D(4, 1) = -4.1; D(5, 1) = -5.3; D(6, 1) = -7.0; D(7, 1) = -8.3;

	D(0, 2) = -6.0; D(1, 2) = -4.0; D(2, 2) = 24.3; D(3, 2) = -1.8;
	D(4, 2) = -1.9; D(5, 2) = -2.8; D(6, 2) = -3.5; D(7, 2) = -4.3;

	D(0, 3) = -5.1; D(1, 3) = -3.8; D(2, 3) = -1.8; D(3, 3) = 21.9;
	D(4, 3) = -1.7; D(5, 3) = -2.5; D(6, 3) = -3.1; D(7, 3) = -3.9;

	D(0, 4) = -5.4; D(1, 4) = -4.1; D(2, 4) = -1.9; D(3, 4) = -1.7;
	D(4, 4) = 22.9; D(5, 4) = -2.4; D(6, 4) = -3.3; D(7, 4) = -4.1;

	D(0, 5) = -7.2; D(1, 5) = -5.3; D(2, 5) = -2.8; D(3, 5) = -2.5;
	D(4, 5) = -2.4; D(5, 5) = 30.3; D(6, 5) = -4.5; D(7, 5) = -5.6;

	D(0, 6) = -9.7; D(1, 6) = -7.0; D(2, 6) = -3.5; D(3, 6) = -3.1;
	D(4, 6) = -3.3; D(5, 6) = -4.5; D(6, 6) = 38.2; D(7, 6) = -7.1;

	D(0, 7) = -11.1; D(1, 7) = -8.3; D(2, 7) = -4.3; D(3, 7) = -3.9;
	D(4, 7) = -4.1; D(5, 7) = -5.6; D(6, 7) = -7.1; D(7, 7) = 44.4;

	D0(0, 0) = 37.2; D0(1, 0) = -13.5; D0(2, 0) = -6.0; D0(3, 0) = -5.1;
	D0(4, 0) = -5.4; D0(5, 0) = -7.2; D0(6, 0) = 0.00; D0(7, 0) = 0.00;

	D0(0, 1) = -13.5; D0(1, 1) = 30.7; D0(2, 1) = -4.0; D0(3, 1) = -3.8;
	D0(4, 1) = -4.1; D0(5, 1) = -5.3; D0(6, 1) = 0.00; D0(7, 1) = 0.00;

	D0(0, 2) = -6.0; D0(1, 2) = -4.0; D0(2, 2) = 16.5; D0(3, 2) = -1.8;
	D0(4, 2) = -1.9; D0(5, 2) = -2.8; D0(6, 2) = 0.00; D0(7, 2) = 0.00;

	D0(0, 3) = -5.1; D0(1, 3) = -3.8; D0(2, 3) = -1.8; D0(3, 3) = 14.9;
	D0(4, 3) = -1.7; D0(5, 3) = -2.5; D0(6, 3) = 0.00; D0(7, 3) = 0.00;

	D0(0, 4) = -5.4; D0(1, 4) = -4.1; D0(2, 4) = -1.9; D0(3, 4) = -1.7;
	D0(4, 4) = 15.5; D0(5, 4) = -2.4; D0(6, 4) = 0.00; D0(7, 4) = 0.00;

	D0(0, 5) = -7.2; D0(1, 5) = -5.3; D0(2, 5) = -2.8; D0(3, 5) = -2.5;
	D0(4, 5) = -2.4; D0(5, 5) = 20.2; D0(6, 5) = 0.00; D0(7, 5) = 0.00;


	D0(0, 6) = -0.00; D0(1, 6) = -0.00; D0(2, 6) = -0.00; D0(3, 6) = -0.00;
	D0(4, 6) = -0.00; D0(5, 6) = -0.00; D0(6, 6) = 0.00; D0(7, 6) = -0.00;

	D0(0, 7) = -0.00; D0(1, 7) = -0.00; D0(2, 7) = -0.00; D0(3, 7) = -0.00;
	D0(4, 7) = -0.00; D0(5, 7) = -0.00; D0(6, 7) = -0.00; D0(7, 7) = 0.00;



	// * on the exchange rate (may be from 0.1 to 10)
	for (int i = 0; i < l; i++)
		D_arr[i] = D * POR * 0.1;


/*
	for (int i = 0; i < l * 0.2; i++)
		D_arr[i] = D * POR * 0.1;
	for (int i = l * 0.2; i < l * 0.8; i++)
		D_arr[i] = D0 * POR * 0.1;
	for (int i = l * 0.8; i < l; i++)
		D_arr[i] = D * POR * 0.1;
*/
}


//edge conditions
//	1 - fixed point (free edge)
//	0 - dp/dx = 0 (closed edge)
static void get_left_edge(vec &left)
{

//	left(0) = 1;
//	left(1) = 1;
	left(2) = 1;
	left(3) = 1;
//	left(4) = 1;
	left(5) = 1;
//	left(6) = 1;
//	left(7) = 1;

}

static void get_right_edge(vec &right)
{

//	right(0) = 1;
//	right(1) = 1;
//	right(2) = 1;
//	right(3) = 1;
	right(4) = 1;
//	right(5) = 1;
	right(6) = 1;
	right(7) = 1;

}

#endif // _SLB_MATRIXES
