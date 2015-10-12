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
	matrix D_ini(n);
	matrix K_ini(n);
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


void progonka::left_edge(std::vector<matrix>& Ps, std::vector<vector>& Qs)
{
	uint n = Ps[0].size();
	uint *u_left = (uint *)calloc(sizeof(uint), n);

	u_left[1] = 1;
	u_left[3] = 1;
	u_left[5] = 1;
	u_left[7] = 1;

	for (int i = 0; i < n; i++) {
		Ps[0](i, i) = u_left[i] ? 1 : 0;
		Qs[0](i) = u_left[i] ? 0 : P_LEFT;
	}

	free(u_left);
};

void progonka::right_edge(std::vector<matrix>& Ps, std::vector<vector>& Qs)
{
	uint n = Ps[0].size();
	uint *u_right = (uint *)calloc(sizeof(uint), n);

	u_right[0] = 1;
	u_right[2] = 1;
	u_right[4] = 1;
	u_right[6] = 1;

	/* Here we should to solve u = Ps u + Qs for u in L - 1 (which is equal to L - 2)
	 * For some edges we have constant valuse, for others we should to solve equation
	 * So, reorganize equations with already known valuse (constant edge condition)
	 * f.i:
	 * (Ps11 - 1) u1 + Ps12 u2 + Ps13 u3 = -Qs1
	 * Ps21 u1 + (Ps22 - 1) u2 + Ps23 u3 = -Qs2
	 * Ps31 u1 + Ps32 u2 + (Ps33 - 1) u3 = -Qs3
	 * Imagine for the second class of pores we have constant edge condition.
	 * Thus, exclude the second equation from the system and use known u2 = P_RIGHT
	 * (Ps11 - 1) u1 + Ps13 u3 = -Qs1 - Ps12 * P_RIGHT
	 * Ps31 u1 + (Ps33 - 1) u3 = -Qs3 - Ps32 * P_RIGHT
	 * Than, we use lapack to solve this system of linear equations
	 */

	matrix Ai(n); Ai = Ps[L - 2];
	vector Bi(n); Bi = Qs[L - 2] * (-1);

	for (int j = 0; j < n; j++) {
		if (u_right[j] == 1) {
			Ai(j, j) -= 1;
			continue;
		} // else (u_right[j] == 0)
		for (int i = 0; i < n; i++)
			Bi(i) -= P_RIGHT * Ai(j, i);
		Bi(j) = 0;
	}

	uint m = 0;
	for (int i = 0; i < n; i++)
		m += u_right[i];

	matrix A(m);
	vector B(m);
	int cnt_i = 0;
	for (int i = 0; i < n; i++) {
		if (u_right[i] == 0)
			continue;
		int cnt_j = 0;
		for(int j = 0; j < n; j++) {
			if (u_right[j] == 0)
				continue;
			A(cnt_i, cnt_j) = Ai(i, j);
			cnt_j++;
		}
		B(cnt_i) = Bi(i);
		cnt_i++;
	}

	double *x = (double *)calloc(sizeof(double), m);
	solve_eq(A.get_ptr(), B.get_ptr(), x, m);

	cnt_i = 0;
	for (int i = 0; i < n; i++) {
		if (u_right[i] == 1) {
			Qs[L - 1](i) = x[cnt_i];
			cnt_i++;
		} else
			Qs[L - 1](i) = P_RIGHT;
	}

	free(x);
	free(u_right);
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

	/* edge conditions */
	left_edge(Ps, Qs);

	/* There */
	for (int i = 1; i < L - 1; i++) {
		Fi = ((*A_ini) * (1.0 / t)) * u[i];
		G = ((*A * Ps[i - 1]) + *B).inverse();

		Ps[i] = (G * (*C)) * (-1);
		Qs[i] = G * (Fi - (*A * Qs[i - 1]));
	}

	/* edge conditions */
	right_edge(Ps, Qs);

	u1[L - 1] = Qs[L - 1];

	/* and Back Again */
	for (int i = L - 2; i >= 0; i--) {
		u1[i] = (Ps[i] * u1[i + 1]) + Qs[i];
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



