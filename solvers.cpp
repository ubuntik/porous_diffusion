/*
 * @file solvers.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <assert.h>
#include <iostream>
#include <math.h>

#include "solvers.h"
//#include "matrixes.h"

#define ETA 1
#define h 1.0

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
	mat inver = inv(B[0]);
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
	u1[l - 1] = G * ((A[l - 1] * Qs[l - 1]) - F[l - 1]);

//	u1[l] = G * ((A[l - 1] * Qs[l - 1]) - F[l - 1]);

	/* and Back Again */

//	for unit_tests check.cpp without edge conditions

	for (int i = l - 2; i >= 0; i--)
		u1[i] = (Ps[i + 1] * u1[i + 1]) + Qs[i + 1];
/*
	for (int i = l - 1; i > 0; i--)
		u1[i] = (Ps[i] * u1[i + 1]) + Qs[i];
*/
	Ps.clear();
	Qs.clear();
};

#if 0
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

	perm = new vector(n);
	for (int i = 0; i < n; i++)
		(*perm)(i) = (*K)(i, i);
};

corner::~corner()
{
	delete K;
	delete D;
};

void corner::left_edge(vector<vec>& C1)
{
	uint n = C1[0].size();
	vector left(n);
	get_left_edge(left);

	for (int i = 0; i < n; i++)
		C1[0](i) = left(i) ? 0 : 1;
};

void corner::right_edge(vector<vec>& C1)
{
	/* for all classes of pores: dC/dx = 0 */
	uint n = C1[0].size();
	for (int i = 0; i < n; i++)
		C1[l](i) = C1[l - 1](i);
};

void corner::calculate(	const vector<vec>& v_med,
			const vector<vec>& p,
			const vector<vec>& C,
			vector<vec>& C1, double dt)
{
	vector v(n);
	vector crt(n);
	vector crt_1(n);
	vector q(n);

	for (int x = 1; x < l; x++) {
		/* split on physical processes */
		/* 1. calculate dC/dt + U dC/dx = 0 */
		v = (*perm).mult(v_med[x]);

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
};

secondord::secondord(uint n_size, uint length)
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

	perm = new vector(n);
	for (int i = 0; i < n; i++)
		(*perm)(i) = (*K)(i, i);
};

secondord::~secondord()
{
	delete K;
	delete D;
};

void secondord::left_edge(vector<vec>& C1)
{
	uint n = C1[0].size();
	vector left(n);
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

void secondord::calculate(const vector<vec>& v_med,
			const vector<vec>& p,
			const vector<vec>& C,
			vector<vec>& C1, double dt)
{
	vector v(n);
	vector c1(n);
	vector c2(n);

	vector c_i(n);
	vector c_im1(n);
	vector c_ip1(n);
	vector c_im2(n);
	vector c_ip2(n);

	vector Dm(n);
	vector Dp(n);
	vector Dmm(n);
	vector Dpp(n);

	vector Qp(n);
	vector Qm(n);

	vector q(n);

	vector crt(n);

	for (int x = 3; x < l - 3; x++) {
		/* split on physical processes */
		/* 1. calculate dC/dt + U dC/dx = 0 */
		v = (*perm).mult(v_med[x]);
		// v = (*K) * v_med[x];

		crt = C[x];
		c1 = C[x];
		c1 = c1 - (crt - C[x - 1]).mult(v) * (dt / h);
		crt = C[x + 1];
		c2 = C[x + 1];
		c2 = c2 - (crt - C[x]).mult(v) * (dt / h);
		c_i = (c1 + C[x] - (c2 - c1).mult(v) * (dt / h)) * 0.5;

		crt = C[x - 1];
		c1 = C[x - 1];
		c1 = c1 - (crt - C[x - 2]).mult(v) * (dt / h);
		crt = C[x];
		c2 = C[x];
		c2 = c2 - (crt - C[x - 1]).mult(v) * (dt / h);
		c_im1 = (c1 + C[x - 1] - (c2 - c1).mult(v) * (dt / h)) * 0.5;

		crt = C[x + 1];
		c1 = C[x + 1];
		c1 = c1 - (crt - C[x]).mult(v) * (dt / h);
		crt = C[x + 2];
		c2 = C[x + 2];
		c2 = c2 - (crt - C[x + 1]).mult(v) * (dt / h);
		c_ip1 = (c1 + C[x + 1] - (c2 - c1).mult(v) * (dt / h)) * 0.5;

		crt = C[x - 2];
		c1 = C[x - 2];
		c1 = c1 - (crt - C[x - 3]).mult(v) * (dt / h);
		crt = C[x - 1];
		c2 = C[x - 1];
		c2 = c2 - (crt - C[x - 2]).mult(v) * (dt / h);
		c_im2 = (c1 + C[x - 2] - (c2 - c1).mult(v) * (dt / h)) * 0.5;

		crt = C[x + 2];
		c1 = C[x + 2];
		c1 = c1 - (crt - C[x + 1]).mult(v) * (dt / h);
		crt = C[x + 3];
		c2 = C[x + 3];
		c2 = c2 - (crt - C[x + 2]).mult(v) * (dt / h);
		c_ip2 = (c1 + C[x + 2] - (c2 - c1).mult(v) * (dt / h)) * 0.5;

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
					((*D)(i,j)) * (p[x](i) - p[x](j));
			}

			C1[x](i) += q(i) * dt;
			if (C1[x](i) < 0)
				C1[x](i) = 0;
		}

	}
};
#endif

