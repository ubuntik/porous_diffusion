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

#define PRNT 1
#define TIME 100
#define t 1
#define L 100
#define h 1.0
#define P_LEFT 1.0
#define P_RIGHT 0.0

void start_cond_1(std::vector<vector>& u)
{
	uint n = u[0].size();
	for (int j = 0; j < n; j++) {
		u.at(0)(j) = P_LEFT;
		for (int i = 1; i < L - 1; i++)
			u.at(i)(j) = 1;
		u.at(L - 1)(j) = P_RIGHT;
	}
};

void simple_ODE(uint n, std::vector<vector>& u, std::vector<vector>& u1)
{
	char buf[256];
	double time = (double)TIME / t;
	uint l = L - 2; // exept edge points
	matrix Al(n, 1.0/t);
	matrix K(n);
	matrix D(n, 1);

	/* This structure keep all data (not using in general) */
	std::vector<std::vector<vector> > u_all(time, std::vector<vector>(L));

	/* due to the specific realization of std containers */
	for (int i = 0; i < time; i++)
		for (int j = 0; j < L; j++)
			u_all[i][j].init(n);

	start_cond_1(u);

	// A(n) u(n - 1) - B(n) u(n) + C(n) u(n + 1) = F(n)
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
	vector left(n, P_LEFT);
	vector right(n, P_RIGHT);

	progonka method(n, l);

	for (int i = 0; i < time; i++) {
		for (int i = 0; i < l; i++)
			F[i] = Al * u[i + 1] * (-1);	// u[i + 1] -> start from 1-st index, not 0
		F[0] = F[0] + K * left * (-1.0 / h / h);
		F[l - 1] = F[l - 1] + K * right * (-1.0 / h / h);

		method.calculate(u1, A, B, C, F);
		u1[0] = left;
		u1[L - 1] = right;

		if (i % PRNT == 0) {
			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk1(u, buf, n, L);
		}

		u = u1;
		u_all[i] = u1;
	}

	/* This printing is u(TIME) throughout Length exept edges
	 * (not using in general)
	 */
	for (int i = 0; i < L; i++) {
		sprintf(buf, "res/time_%06d.vtk", i);
		write_to_vtk3(u_all, buf, n, i, time);
	}
}

int main(int argc, char **argv)
{
	// the problem size
	uint n = 1;
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

	simple_ODE(n, u, u1);

	std::cout << "DONE" << std::endl;
	return 0;
}

