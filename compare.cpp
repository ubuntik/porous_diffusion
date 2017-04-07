/*
 * @file compare.cpp
 * @brief Comparison with accurate solution for D = 0
 * @author Anna Subbotina
 *
 */

#include <iostream>

#include <math.h>
#include <assert.h>

#include "solvers.h"
#include "vtk.h"

#define PRNT 1
#define TIME 100
#define t 1
#define L 100
#define h 1.0
#define P_LEFT 0.0
#define P_RIGHT 1.0

vector<double> acc(L);

double erf(double x)
{
	double ret = 0;
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

void start_cond_1(vector<vec>& u)
{
	uint n = u[0].size();
	for (int j = 0; j < n; j++) {
		u.at(0)(j) = P_LEFT;
		for (int i = 1; i < L - 1; i++)
			u.at(i)(j) = 1;
		u.at(L - 1)(j) = P_RIGHT;
	}
};

void compare_accurate(uint n, vector<vec>& u, vector<vec>& u1)
{
	char buf[256];
	double time = (double)TIME / t;
	uint l = L - 2; // exept edge points
	mat e(n, n, fill::eye);
	mat Al = e * (1.0/t);
	mat K(n, n, fill::eye);
	mat D(n, n, fill::zeros);

	start_cond_1(u);

	vector<mat> A(l, mat(n, n, fill::zeros));
	vector<mat> B(l, mat(n, n, fill::zeros));
	vector<mat> C(l, mat(n, n, fill::zeros));
	vector<vec> F(l, vec(n, fill::zeros));

	for (int i = 0; i < l; i++) {
		A[i] = K * (1.0 / h / h);
		B[i] = Al + K * (2.0 / h / h) + D;
		C[i] = K * (1.0 / h / h);
	}
	vec ve(n, fill::ones);
	vec left = ve * P_LEFT;
	vec right = ve * P_RIGHT;

	progonka method(n, l);

	for (int i = 1; i < time; i++) {
		for (int i = 0; i < l; i++)
			F[i] = Al * u[i + 1] * (-1);    // u[i + 1] -> start from 1-st index, not 0
		F[0] = F[0] + K * left * (-1.0 / h / h);
		F[l - 1] = F[l - 1] + K * right * (-1.0 / h / h);

		method.calculate(u1, A, B, C, F);
		u1[0] = left;
		u1[L - 1] = right;

		if (i % PRNT == 0) {
			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk2(u, buf, n, L, i*t);
		}

		u = u1;
	}
}

int main(int argc, char **argv)
{
	// the problem size
	uint n = 1;
	vector<vec> u(L, vec(n, fill::zeros));
	vector<vec> u1(L, vec(n, fill::zeros));

	cout << "The dimension of prodlem: " << n << endl;
	cout << "The length = " << L << ", the time = " << TIME << endl;
	cout << "Space step = " << h << ", time step = " << t << endl;

	compare_accurate(n, u, u1);

	cout << "DONE" << endl;
	return 0;
}

