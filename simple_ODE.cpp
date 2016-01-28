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
#include "config.h"

#define PRNT 1
#undef TIME
#define TIME 100
#undef t
#define t 1
#undef L
#define L 100

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
	/* This structure keep all data (not using in general) */
	std::vector<std::vector<vector> > u_all((double)TIME / t, std::vector<vector>(L));

	/* due to the specific realization of std containers */
	for (int i = 0; i < ((double)TIME / t); i++)
		for (int j = 0; j < L; j++)
			u_all[i][j].init(n);

	start_cond_1(u);

	progonka method(n);

	for (int i = 0; i < ((double)TIME / t); i++) {
		if (i % PRNT == 0) {
			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk1(u, buf, n);
		}

		method.calculate(u, u1);

		u = u1;
		u_all[i] = u1;
	}

	/* This printing is u(TIME) throughout Length exept edges
	 * (not using in general)
	 */
	for (int i = 1; i < L - 1; i++) {
		sprintf(buf, "res/time_%06d.vtk", i);
		write_to_vtk3(u_all, buf, n, i);
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

