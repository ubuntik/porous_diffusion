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

#include "../lalgebra.h"
#include "../solvers.h"
#include "../vtk.h"

// edge conditions
#define P_LEFT 1
#define P_RIGHT 0

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

	std::cout << "was" << std::endl;
	for (int i = 0; i < L; i++)
		u[i].print();

	// A(n) u(n - 1) - B(n) u(n) + C(n) u(n + 1) = F(n)
	std::vector<matrix> A(L);
	std::vector<matrix> B(L);
	std::vector<matrix> C(L);
	std::vector<vector> F(L);

	for (int i = 0; i < L; i++) {
		A[i].init(n, 1);
		B[i].init(n, 3);
		C[i].init(n, 2);
	}

	F[0].init(n, 1);
	F[1].init(n, 1);
	F[2].init(n, 1);
	F[3].init(n, -3);
	F[4].init(n, -1);
	F[5].init(n, -1);
	F[6].init(n, -1);
	F[7].init(n, 1);

	progonka method(n, L, A, B, C, F);

	method.calculate(u1);

	std::cout << "now" << std::endl;
	for (int i = 0; i < L; i++)
		u1[i].print();

	std::cout << "should be:" << std::endl;
	std::cout << "( 1 1 )\n( 2 2 )\n( 3 3 )\n( 4 4 )" << std::endl;
	std::cout << "( 3 3 )\n( 2 2 )\n( 1 1 )\n( 0 0 )" << std::endl;

	A.clear();
	B.clear();
	C.clear();
	F.clear();
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
	std::cout << "The length = " << L << std::endl;

	check_progonka(n, u, u1);

	std::cout << "DONE" << std::endl;

	u.clear();
	u1.clear();
	return 0;
}

