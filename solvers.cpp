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

#define ETA 1
#define h 1.0

progonka::progonka(uint n_size, uint length)
{
	n = n_size;
	l = length;
	assert(n >= 0);
};

void progonka::calculate(std::vector<vector>& u1,
		std::vector<matrix> &A,
		std::vector<matrix> &B,
		std::vector<matrix> &C,
		std::vector<vector> &F)
{
	/* coefficients */
	std::vector<matrix> Ps(l);
	std::vector<vector> Qs(l);
	matrix G(n);

	/* due to the specific realization of std containers */
	for (int i = 0; i < l; i++) {
		Ps[i].init(n);
		Qs[i].init(n);
	}

	// Ps[0] and Qs[0] are not used
	Ps[1] = B[0].inverse() * C[0];
	Qs[1] = B[0].inverse() * F[0] * (-1);

	/* There */
	for (int i = 2; i < l; i++) {
		G = (B[i - 1] - (A[i - 1] * Ps[i - 1])).inverse();

		Ps[i] = G * C[i - 1];
		Qs[i] = G * ((A[i - 1] * Qs[i - 1]) - F[i - 1]);
	}

	G = (B[l - 1] - (A[l - 1] * Ps[l - 1])).inverse();

//	for unit_tests check.cpp without edge conditions
//	u1[l - 1] = G * ((A[l - 1] * Qs[l - 1]) - F[l - 1]);

	u1[l] = G * ((A[l - 1] * Qs[l - 1]) - F[l - 1]);

	/* and Back Again */

//	for unit_tests check.cpp without edge conditions
/*
	for (int i = l - 2; i >= 0; i--) {
		u1[i] = (Ps[i + 1] * u1[i + 1]) + Qs[i + 1];
*/
	for (int i = l - 1; i > 0; i--) {
		u1[i] = (Ps[i] * u1[i + 1]) + Qs[i];
	}

	Ps.clear();
	Qs.clear();
};

corner::corner(uint n_size, uint length)
{
	n = n_size;
	l = length;
	assert(n >= 0);

	matrix K_ini(n);
	get_K(K_ini);

	K = new matrix(n);
	*K = K_ini * (1.0 / ETA);

	D = new matrix(n);
	get_D(*D);
};

corner::~corner()
{
	delete K;
	delete D;
};

void corner::left_edge(std::vector<vector>& C1)
{
	uint n = C1[0].size();
	uint *u_left = (uint *)calloc(sizeof(uint), n);
	get_left_edge(u_left);

	for (int i = 0; i < n; i++)
		C1[0](i) = (u_left[i] == 1) ? 0 : 1;

	free(u_left);
};

void corner::right_edge(std::vector<vector>& C1)
{
	/* for all classes of pores: dC/dx = 0 */
	uint n = C1[0].size();
	for (int i = 0; i < n; i++)
		C1[l](i) = C1[l - 1](i);
};

void corner::calculate(	const std::vector<vector>& v_med,
			const std::vector<vector>& p,
			const std::vector<vector>& C,
			std::vector<vector>& C1, double dt)
{
	vector v(n);
	vector crt(n);
	vector crt_1(n);
	vector q(n);

	left_edge(C1);

	for (int x = 1; x < l; x++) {
		/* split on physical processes */
		/* 1. calculate dC/dt + U dC/dx = 0 */
		v = (*K) * v_med[x];
		crt = C[x];
		crt_1 = C[x + 1];
		C1[x] = C[x];
		C1[x] = C1[x] - ((v.dif_abs()).mult(crt_1 - C[x]) +
			(v.add_abs()).mult(crt - C[x - 1])) * (dt / 2.0 / h);

		/* 2. calculate additional mass transfer dC/dt + U dC/dx = q */
		/* q = -C(i|Pi >= Pj, i != j / j) D(i,j) (Pi - Pj) */

		for (int i = 0; i < n; i++) {
			q(i) = 0;
			for (int j = 0; j < n; j++) {
				if (i == j)
					continue;
				q(i) += ((p[x](i) >= p[x](j)) ? C1[x](i) : C1[x](j)) *
					((*D)(i,j)) * (p[x](i) - p[x](j));
			}

			C1[x](i) += q(i) * dt;
			if (C1[x](i) < 0)
				C1[x](i) = 0;
		}

	}

	right_edge(C1);
};



