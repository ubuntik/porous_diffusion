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

#define PRNT 100

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
	start_cond(u);

	progonka method(n);

	for (int i = 1; i < ((double)TIME / t); i++) {
		if (i % PRNT == 0) {
			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk1(u, buf, n);
		}

		method.calculate(u, u1);
		u = u1;
	}
}

int main(int argc, char **argv)
{
	// the problem size
	uint n = power(P, GAMMA - 1);
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

