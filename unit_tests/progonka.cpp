/*
 * @file main.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <iostream>
#include <math.h>
#include <assert.h>

#include "../solvers.h"

// edge conditions
#define P_LEFT 1
#define P_RIGHT 0

#define L 8
#define n 2

void start_cond(vector<vec>& u)
{
	for (int i = 0; i < n; i++) {
		u.at(0)(i) = P_LEFT;
		u.at(L - 1)(i) = P_RIGHT;
	}
};

void check_progonka(vector<vec>& u, vector<vec>& u1)
{
	start_cond(u);

	cout << "was" << endl;
	for (int i = 0; i < L; i++)
		cout << u[i];

	// A(n) u(n - 1) - B(n) u(n) + C(n) u(n + 1) = F(n)
	vector<mat> A(L, mat(n, n, fill::eye));
	vector<mat> B(L, mat(n, n, fill::zeros));
	vector<mat> C(L, mat(n, n, fill::zeros));
	vector<vec> F(L, vec(n, fill::zeros));


	for (int i = 0; i < L; i++) {
		mat e(n, n, fill::eye);
		B[i] = 3 * e;
		C[i] = 2 * e;
	}

	F[0].fill(1);
	F[1].fill(1);
	F[2].fill(1);
	F[3].fill(-3);
	F[4].fill(-1);
	F[5].fill(-1);
	F[6].fill(-1);
	F[7].fill(1);

	progonka method(n, L);

	method.calculate(u1, A, B, C, F);

	cout << "now" << endl;
	for (int i = 0; i < L; i++)
		cout << u1[i];

	cout << "should be:" << endl;
	cout << "( 1 1 )\n( 2 2 )\n( 3 3 )\n( 4 4 )" << endl;
	cout << "( 3 3 )\n( 2 2 )\n( 1 1 )\n( 0 0 )" << endl;

	A.clear();
	B.clear();
	C.clear();
	F.clear();
}

int main(int argc, char **argv)
{
	// the problem size
	vector<vec> u (L, vec(n, fill::zeros));
	vector<vec> u1 (L, vec(n, fill::zeros));

	std::cout << "The dimension of prodlem: " << n << std::endl;
	std::cout << "The length = " << L << std::endl;

	check_progonka(u, u1);

	std::cout << "DONE" << std::endl;

	u.clear();
	u1.clear();
	return 0;
}

