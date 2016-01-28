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

std::vector<ptype> acc(L);

ptype erf(ptype x)
{
	ptype ret = 0;
	double a1 = 0.278393;
	double a2 = 0.230389;
	double a3 = 0.000972;
	double a4 = 0.078108;

	ret = (1.0 - (1.0 / (1.0 + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x) /
			(1.0 + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x) /
			(1.0 + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x) /
			(1.0 + a1*x + a2*x*x + a3*x*x*x + a4*x*x*x*x)));
	return ret;
}

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

void compare_accurate(uint n, std::vector<vector>& u, std::vector<vector>& u1)
{
	char buf[256];
	start_cond_1(u);

	progonka method(n);

	for (int i = 1; i < ((double)TIME / t); i++) {
		if (i % PRNT == 0) {
			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk2(u, buf, n, i*t);
		}

		method.calculate(u, u1);
		u = u1;
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

	compare_accurate(n, u, u1);

	std::cout << "DONE" << std::endl;
	return 0;
}

