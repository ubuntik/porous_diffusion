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

#undef L
#define L 8

progonka::progonka(uint n_size)
{
	n = n_size;
	assert(n >= 0);

/*	A_ini = new matrix(n, 1.0);
	A = new matrix(n);
	B = new matrix(n);
	C = new matrix(n);
	matrix D_ini(n);
	matrix K_ini(n);
	get_K(K_ini);
	get_D(D_ini);
*/
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
/*
	std::cout << "Initian matrix (A, K, D):" << std::endl;
	(*A_ini).print();
	std::cout << "******************" << std::endl;
	K_ini.print();
	std::cout << "******************" << std::endl;
	D_ini.print();

	*A = K_ini * (-1.0 / h / h);
	*B = ((*A_ini) * (1.0 / t)) + K_ini * (2.0 / h / h) + D_ini;
	*C = K_ini * (-1.0 / h / h);
*/

	A = new matrix(n, 1);
	B = new matrix(n, 3);
	C = new matrix(n, 2);
	std::cout << "******************" << std::endl;
	A->print();
	std::cout << "******************" << std::endl;
	B->print();
	std::cout << "******************" << std::endl;
	C->print();
	std::cout << "******************" << std::endl;

};

progonka::~progonka()
{
	delete A_ini;
	delete A;
	delete B;
	delete C;
};

void progonka::calculate(const std::vector<vector>& u, std::vector<vector>& u1)
{
//	assert(A_ini);
	assert(A);
	assert(B);
	assert(C);

	/* coefficients */
	std::vector<matrix> Ps(L);
	std::vector<vector> Qs(L);
	matrix G(n);

	std::vector<vector> F(L);

	/* due to the specific realization of std containers */
	for (int i = 0; i < L; i++) {
		Ps[i].init(n);
		Qs[i].init(n);
	}

	F[0].init(n, 1);
	F[1].init(n, 1);
	F[2].init(n, 1);
	F[3].init(n, -3);
	F[4].init(n, -1);
	F[5].init(n, -1);
	F[6].init(n, -1);
	F[7].init(n, 1);

	/* edge conditions */
//	left_edge(Ps, Qs);

	Ps[1] = (*B).inverse() * (*C);
	Qs[1] = (*B).inverse() * F[0] * (-1);

	/* There */
	for (int i = 2; i < L; i++) {
		G = (*B - (*A * Ps[i - 1])).inverse();

		Ps[i] = G * (*C);
		Qs[i] = G * ((*A * Qs[i - 1]) - F[i - 1]);
	}

	/* edge conditions */
//	right_edge(Ps, Qs);

	G = ((*B) * (-1) - (*A * Ps[L - 1])).inverse();
	u1[L - 1] = G * ((*A * Qs[L - 1]) - F[L - 1]);

	/* and Back Again */
	for (int i = L - 2; i >= 0; i--) {
		u1[i] = (Ps[i + 1] * u1[i + 1]) + Qs[i + 1];
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
		C1[L](i) = C1[L - 1](i);
};

void corner::calculate(	const std::vector<vector>& v_med,
			const std::vector<vector>& p,
			const std::vector<vector>& C,
			std::vector<vector>& C1, double dt)
{
	uint n = power(P, GAMMA - 1);
	vector v(n);
	vector crt(n);
	vector crt_1(n);
	vector q(n);

	left_edge(C1);

	for (int x = 1; x < L; x++) {
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



