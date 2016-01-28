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

#undef L
#define L 8

void start_cond(std::vector<vector>& u)
{
	uint n = u[0].size();
	for (int i = 0; i < n; i++) {
		u.at(0)(i) = P_LEFT;
		u.at(L - 1)(i) = P_RIGHT;
	}
};

void check_progonka(uint n, std::vector<vector>& u, std::vector<vector>& u1)
{
	start_cond(u);
	progonka method(n);

	std::cout << "was" << std::endl;
	for (int i = 0; i < L; i++)
		u[i].print();

	method.calculate(u, u1);

	std::cout << "now" << std::endl;
	for (int i = 0; i < L; i++)
		u1[i].print();

	std::cout << "should be:" << std::endl;
	std::cout << "( 1 1 )\n( 2 2 )\n( 3 3 )\n( 4 4 )" << std::endl;
	std::cout << "( 3 3 )\n( 2 2 )\n( 1 1 )\n( 0 0 )" << std::endl;
}

int main(int argc, char **argv)
{
	// the problem size
	uint n = 2;
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

	check_progonka(n, u, u1);

	std::cout << "DONE" << std::endl;
	return 0;
}

