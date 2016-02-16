/*
 * @file main.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <iostream>
#include <vector>

#include <math.h>
#include <assert.h>

#include "lalgebra.h"
#include "solvers.h"
#include "vtk.h"
#include "matrixes.h"

#define PRNT 100
#define TIME 10000
#define t 1
#define L 100
#define h 1.0
#define P_LEFT 1.0
#define P_RIGHT 0.0

void start_cond(std::vector<vector>& u)
{
	uint n = u[0].size();
	for (int i = 0; i < n; i++) {
		u.at(0)(i) = P_LEFT;
		u.at(L - 1)(i) = P_RIGHT;
	}
};

void direct_problem(uint n, std::vector<vector>& u, std::vector<vector>& u1)
{
	char buf[256];
	double time = (double)TIME / t;
	uint l = L - 2; // exept edge points
	matrix Al(n, 1.0/t);
	matrix K(n);
	get_K(K);
	matrix D(n);
	get_D(D);

	start_cond(u);
	sprintf(buf, "res/data_000000.vtk");
	write_to_vtk1(u, buf, n, L);

	vector left(n);
	vector right(n);

//	left(0) = 1;
	left(1) = 1;
//	left(2) = 1;
	left(3) = 1;
//	left(4) = 1;
	left(5) = 1;
//	left(6) = 1;
	left(7) = 1;

	right(0) = 1;
	right(1) = 1;
	right(2) = 1;
	right(3) = 1;
	right(4) = 1;
	right(5) = 1;
	right(6) = 1;
	right(7) = 1;

	std::vector<matrix> A(l);
        std::vector<matrix> B(l);
	std::vector<matrix> C(l);
	std::vector<vector> F(l);

	for (int i = 0; i < l; i++) {
		A[i].init(n);
		A[i] = K * (1.0 / h / h);
		B[i].init(n);
		B[i] = Al + K * (2.0 / h / h) + D;
		C[i].init(n);
		C[i] = K * (1.0 / h / h);
		F[i].init(n);
	}

	matrix E(n);
	vector Fl(n);
	for (int i = 0; i < n; i++) {
		E(i, i) = left(i);
		if (left(i) == 1)
			continue;
		vector l(n);
		l(i) = P_LEFT;
		Fl = Fl + A[0] * l;
	}
	matrix Bl(n);
	Bl = A[0] * E;

	vector Fr(n);
	for (int i = 0; i < n; i++) {
		E(i, i) = right(i);
		if (right(i) == 1)
			continue;
		vector r(n);
		r(i) = P_RIGHT;
		Fr = Fr + C[l - 1] * r;
	}
	matrix Br(n);
	Br = C[l - 1] * E;

	B[0] = B[0] - Bl;
	B[l - 1] = B[l - 1] - Br;

	progonka method(n, l);

	for (int i = 1; i < time; i++) {
		for (int i = 0; i < l; i++)
			F[i] = Al * u[i + 1] * (-1);    // u[i + 1] -> start from 1-st index, not 0
		F[0] = F[0] - Fl;
		F[l - 1] = F[l - 1] - Fr;

		method.calculate(u1, A, B, C, F);
		for (int i = 0; i < n; i++) {
			u1[0](i) = left(i) ? u1[1](i) : P_LEFT;
			u1[L - 1](i) = right(i) ? u1[L - 2](i) : P_RIGHT;
		}

		if (i % PRNT == 0) {
			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk1(u1, buf, n, L);
		}

		u = u1;
	}
}

int main(int argc, char **argv)
{
	// the problem size
	uint n = 8;
	std::vector<vector> u(L);
	std::vector<vector> u1(L);

	/* due to the specific realization of std containers */
	for (int i = 0; i < L; i++) {
		u[i].init(n);
		u1[i].init(n);
	}

	std::cout << "The dimension of prodlem: " << n << std::endl;
	std::cout << "The length = " << L << ", the time = " << TIME << std::endl;
	std::cout << "Space step = " << h << ", time step = " << t << std::endl;

	direct_problem(n, u, u1);

	std::cout << "DONE" << std::endl;
	return 0;
}

