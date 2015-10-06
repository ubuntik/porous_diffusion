/*
 * @file solvers.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <assert.h>
#include <iostream>
#include <vector>
#include <math.h>

#include "solvers.h"
#include "lalgebra.h"
#include "matrixes.h"

progonka::progonka(uint n_size)
{
	n = n_size;
	assert(n >= 0);

	A_ini = new matrix(n, 1.0);
	A = new matrix(n);
	B = new matrix(n);
	C = new matrix(n);
	matrix K_ini(n);
	matrix D_ini(n);
	get_K(K_ini);
	get_D(D_ini);

#if 0
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			K_ini(i, j) = (i == j) ? K_DIAG : K_ELEM;
			D_ini(i, j) = (i == j) ? 1 : 0;

			/* for start conditions != 0
			 * It must decrease exponentially
			 */
//			D_ini(i, j) = (i == j) ? K_DIAG : K_ELEM;
//			K_ini(i, j) = 0;
		}
	}
#endif

	std::cout << "Initian matrix (A, K, D):" << std::endl;
	(*A_ini).print();
	std::cout << "******************" << std::endl;
	K_ini.print();
	std::cout << "******************" << std::endl;
	D_ini.print();

	*A = K_ini * (-1.0 / h / h);
	*B = ((*A_ini) * (1.0 / t)) + K_ini * (2.0 / h / h) + D_ini;
	*C = K_ini * (-1.0 / h / h);
};

progonka::~progonka()
{
	delete A_ini;
	delete A;
	delete B;
	delete C;
};


void progonka::edge_conditions(const std::vector<vector>& u, std::vector<vector>& u1)
{
	u1[0] = u[0];
	u1[L - 1] = u[L - 1];
};


void progonka::calculate(const std::vector<vector>& u, std::vector<vector>& u1)
{
	assert(A_ini);
	assert(A);
	assert(B);
	assert(C);

	/* right part */
	vector Fi(n, 1);

	/* coefficients */
	std::vector<matrix> Ps(L);
	std::vector<vector> Qs(L);
	matrix G(n);

	/* due to the specific realization of std containers */
	for (int i = 0; i < L; i++) {
		Ps[i].init(n);
		Qs[i].init(n);
	}

	/* left edge */
	uint *u_left = (uint *)calloc(sizeof(uint), n);
/*
	u_left[1] = 1;
	u_left[3] = 1;
	u_left[5] = 1;
	u_left[7] = 1;
*/
	for (int i = 0; i < n; i++) {
		Ps[0](i, i) = u_left[i] ? P_LEFT : 0;
		Qs[0](i) = u_left[i] ? 0 : 1;
	}

	/* edge conditions */
	//edge_conditions(u, u1);

	/* There */
	for (int i = 1; i < L; i++) {
		Fi = ((*A_ini) * (1.0 / t)) * u[i];
		G = ((*A * Ps[i - 1]) + *B).inverse();

		Ps[i] = (G * (*C)) * (-1);
		Qs[i] = G * (Fi - (*A * Qs[i - 1]));
	}

	/* and Back Again */
	for (int i = L - 2; i >= 0; i--) {
		u1[i] = (Ps[i] * u1[i + 1]) + Qs[i];
	}

	/* right edge */
	uint *u_right = (uint *)calloc(sizeof(uint), n);
/*
	u_right[0] = 1;
	u_right[2] = 1;
	u_right[4] = 1;
	u_right[6] = 1;
*/

	//Qs[L - 1] = u1[L - 2] - Qs[L - 1];

	double *x = (double *)calloc(sizeof(double), n);
	double *a = (double *)calloc(sizeof(double), n * n);
	double *b = (double *)calloc(sizeof(double), n);
	memcpy(a, Ps[L - 1].get_ptr(), sizeof(double) * n * n);
	memcpy(b, Qs[L - 1].get_ptr(), sizeof(double) * n);

	//solve_eq(a, b, x, n);

	for (int i = 0; i < n; i++) {
		u1[L - 1](i) = u_right[i] ? x[i] : u1[L - 1](i);
	}

	Ps.clear();
	Qs.clear();
};

corner::corner(uint n_size)
{
	n = n_size;
	assert(n >= 0);

	matrix K_ini(n);
	get_K(K_ini);

	K = new matrix(n);
	*K = K_ini * (1.0 / ETA);
};

corner::~corner()
{
	delete K;
};

void corner::edge_conditions(const std::vector<vector>& C, std::vector<vector>& C1)
{
	C1[0] = C[0];
	C1[L] = C[L];
};

void corner::calculate(	std::vector<vector>& v_med,
			const std::vector<vector>& C,
			std::vector<vector>& C1, double dt)
{
	uint n = power(P, GAMMA - 1);
	vector v(n);
	vector crt(n);
	vector crt_1(n);
	/* edge conditions */
	edge_conditions(C, C1);

	for (int i = 1; i < L; i++) {
		v = (*K) * v_med[i];
		crt = C[i];
		crt_1 = C[i + 1];
		C1[i] = C[i];
		C1[i] = C1[i] - ((v.dif_abs()).mult(crt_1 - C[i]) +
			(v.add_abs()).mult(crt - C[i - 1])) * (dt / 2.0 / h);
	}
};



