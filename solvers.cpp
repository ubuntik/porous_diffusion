/*
 * @file solvers.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <assert.h>
#include <iostream>
#include <vector>

#include "lalgebra.h"
#include "solvers.h"

progonka::progonka(uint n_size)
{
	n = n_size;
	assert(n >= 0);

	A_ini = new matrix(n, 1.0 / t);
	A = new matrix(n);
	B = new matrix(n);
	C = new matrix(n);
	matrix K_ini(n);
	matrix D_ini(n);

	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			K_ini(i, j) = (i == j) ? K_DIAG : K_ELEM;
			D_ini(i, j) = 0;
		}
	}
	*A = K_ini * (-1.0 / h / h);
	*B = (*A_ini) + K_ini * (2.0 / h / h) + D_ini;
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
//	for (std::vector<vector>::iterator it = Ps.begin(); it != Ps.end(); ++it) {
	for (int i = 0; i < L; i++) {
		Ps[i].init(n);
		Qs[i].init(n);
	}

	/* progonka coeffs start cond */
	Fi = (*A_ini) * u[0];
	Ps[0] = ((*B).inverse() * (*C)) * (-1);
	Qs[0] = (*B).inverse() * Fi;

	/* edge conditions */
	edge_conditions(u, u1);

	/* There */
	for (int i = 1; i < L - 1; i++) {
		Fi = (*A_ini) * u[i];
		G = ((*A * Ps[i - 1]) + *B).inverse();

		Ps[i] = (G * (*C)) * (-1);
		Qs[i] = G * (Fi - (*A * Qs[i - 1]));
	}

	u1[L - 2] = Qs[L - 2];

	/* and Back Again */
	for (int i = L - 2; i > 0; i--) {
		u1[i] = (Ps[i] * u1[i + 1]) + Qs[i];
	}
}


