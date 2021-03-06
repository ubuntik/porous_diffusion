/*
 * @file solvers.cpp
 * @author Anna Subbotina
 *
 */

#include <assert.h>
#include <iostream>
#include <math.h>

#include "solvers.h"
#include "matrixes.h"

progonka::progonka(uint n_size, uint length)
{
	n = n_size;
	l = length;
	assert(n >= 0);
};

void progonka::calculate(vector<vec>& u1,
		vector<mat> &A,
		vector<mat> &B,
		vector<mat> &C,
		vector<vec> &F)
{
	/* coefficients */
	vector<mat> Ps (l, mat(n, n, fill::zeros));
	vector<vec> Qs (l, vec(n, fill::zeros));
	mat G(n, n, fill::zeros);

	// Ps[0] and Qs[0] are not used
	mat inver(n, n, fill::zeros);
	try {
		inver = inv(B[0]);
	} catch (...) {
		inver = pinv(B[0]);
	}
	Ps[1] = inver * C[0];
	Qs[1] = -inver * F[0];

	/* There */
	for (int i = 2; i < l; i++) {
		G = (B[i - 1] - (A[i - 1] * Ps[i - 1])).i();

		Ps[i] = G * C[i - 1];
		Qs[i] = G * (A[i - 1] * Qs[i - 1] - F[i - 1]);
	}

	G = (B[l - 1] - (A[l - 1] * Ps[l - 1])).i();

//	for unit_tests check.cpp without edge conditions
//	u1[l - 1] = G * ((A[l - 1] * Qs[l - 1]) - F[l - 1]);

	u1[l] = G * ((A[l - 1] * Qs[l - 1]) - F[l - 1]);

	/* and Back Again */

//	for unit_tests check.cpp without edge conditions
/*
	for (int i = l - 2; i >= 0; i--)
		u1[i] = (Ps[i + 1] * u1[i + 1]) + Qs[i + 1];
*/
	for (int i = l - 1; i > 0; i--)
		u1[i] = (Ps[i] * u1[i + 1]) + Qs[i];

};

secondord::secondord(uint n_size, uint length,
			vector<mat> *Ki, vector<mat> *Di)
{
	n = n_size;
	l = length;
	assert(n >= 0);
	K = Ki;
	D = Di;
};

secondord::~secondord()
{};

void secondord::left_edge(vector<vec>& C1)
{
	uint n = C1[0].size();
	vec left(n, fill::zeros);
	get_left_edge(left);

	for (int i = 0; i < n; i++) {
		C1[1](i) = left(i) ? 0 : 1;
		C1[0](i) = left(i) ? 0 : 1;
	}
};

void secondord::right_edge(vector<vec>& C1)
{
	/* for all classes of pores: dC/dx = 0 */
	uint n = C1[0].size();
	for (int i = 0; i < n; i++) {
		C1[l - 1](i) = C1[l - 2](i);
		C1[l](i) = C1[l - 1](i);
	}
};

void secondord::calculate(const vector<vec>& v,
			const vector<vec>& p,
			const vector<vec>& C, vector<vec>& C1,
			double dt, double h)
{
	double cur = dt / h;
	vec c1(n, fill::zeros);
	vec c2(n, fill::zeros);

	vec c_i(n, fill::zeros);
	vec c_im1(n, fill::zeros);
	vec c_ip1(n, fill::zeros);
	vec c_im2(n, fill::zeros);
	vec c_ip2(n, fill::zeros);

	vec Dm(n, fill::zeros);
	vec Dp(n, fill::zeros);
	vec Dmm(n, fill::zeros);
	vec Dpp(n, fill::zeros);

	vec Qp(n, fill::zeros);
	vec Qm(n, fill::zeros);

	vec q(n, fill::zeros);

	for (int x = 3; x < l - 3; x++) {
		/* split on physical processes */
		/* 1. calculate dC/dt + U dC/dx = 0 */
		c1 = C[x] - ((C[x] - C[x - 1]) % v[x]) * cur;
		c2 = C[x + 1] - ((C[x + 1] - C[x]) % v[x + 1]) * cur;
		c_i = (c1 + C[x] - ((c2 - c1) % v[x]) * (dt / h)) * 0.5;

		c1 = C[x - 1] - ((C[x - 1] - C[x - 2]) % v[x - 1]) * cur;
		c2 = C[x] - ((C[x] - C[x - 1]) % v[x]) * (dt / h);
		c_im1 = (c1 + C[x - 1] - ((c2 - c1) % v[x - 1]) * cur) * 0.5;

		c1 = C[x + 1] - ((C[x + 1] - C[x]) % v[x + 1]) * cur;
		c2 = C[x + 2] - ((C[x + 2] - C[x + 1]) % v[x + 2]) * cur;
		c_ip1 = (c1 + C[x + 1] - ((c2 - c1) % v[x + 1]) * cur) * 0.5;

		c1 = C[x - 2] - ((C[x - 2] - C[x - 3]) % v[x - 2]) * cur;
		c2 = C[x - 1] - ((C[x - 1] - C[x - 2]) % v[x - 1]) * cur;
		c_im2 = (c1 + C[x - 2] - ((c2 - c1) % v[x - 2]) * cur) * 0.5;

		c1 = C[x + 2] - ((C[x + 2] - C[x + 1]) % v[x + 2]) * cur;
		c2 = C[x + 3] - ((C[x + 3] - C[x + 2]) % v[x + 3]) * cur;
		c_ip2 = (c1 + C[x + 2] - ((c2 - c1) % v[x + 2]) * cur) * 0.5;

		Dm = c_i - c_im1;
		Dp = c_ip1 - c_i;
		Dmm = c_im1 - c_im2;
		Dpp = c_ip2 - c_ip1;

		for (int i = 0; i < n; i++) {
			Qp(i) = ((Dm(i) * Dp(i) <= 0) || (Dp(i) * Dpp(i) <= 0)) ? Dp(i) : 0;
			Qm(i) = ((Dm(i) * Dp(i) <= 0) || (Dm(i) * Dmm(i) <= 0)) ? Dm(i) : 0;
		}

		C1[x] = c_i + (Qp - Qm) * 0.1;

		/* 2. calculate additional mass transfer dC/dt + U dC/dx = q */
		/* q = -C(i|Pi >= Pj, i != j / j) D(i,j) (Pi - Pj) */

		for (int i = 0; i < n; i++) {
			q(i) = 0;
			for (int j = 0; j < n; j++) {
				if (i == j)
					continue;
				q(i) += ((p[x](i) >= p[x](j)) ? C1[x](i) : C1[x](j)) *
					((*D)[x](i,j)) * (p[x](i) - p[x](j));
			}

			C1[x](i) += q(i) * dt;
			if (C1[x](i) < 0)
				C1[x](i) = 0;
		}

	}
};

