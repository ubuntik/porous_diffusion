/*
 * @file main.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <iostream>

#include <math.h>
#include <assert.h>

#include "solvers.h"
#include "vtk.h"
#include "matrixes.h"

#define PRNT 1
#define TIME 100
#define t 1
#define L 100
#define h 1.0
#define P_LEFT 1.0
#define P_RIGHT 0.0

void start_cond(vector<vec>& u)
{
	uint n = u[0].size();
	for (int i = 0; i < n; i++) {
		u.at(0)(i) = P_LEFT;
		u.at(L - 1)(i) = P_RIGHT;
	}
};

void direct_problem(uint n, vector<vec>& u, vector<vec>& u1)
{
	char buf[256];
	double time = (double)TIME / t;
	uint l = L - 2; // exept edge points
	mat e(n, n, fill::eye);
	mat Al = e * (1.0/t);
	mat K(n, n, fill::zeros);
	get_K(K);
	mat D(n, n, fill::zeros);
	get_D(D);

	start_cond(u);
	sprintf(buf, "res/data_000000.vtk");
	write_to_vtk1(u, buf, n, L);

	vec left(n, fill::zeros);
	get_left_edge(left);
	vec right(n, fill::zeros);
	get_right_edge(right);

	vector<mat> A(l, mat(n, n, fill::zeros));
        vector<mat> B(l, mat(n, n, fill::zeros));
	vector<mat> C(l, mat(n, n, fill::zeros));
	vector<vec> F(l, vec(n, fill::zeros));

	for (int i = 0; i < l; i++) {
		A[i] = K * (1.0 / h / h);
		B[i] = Al + K * (2.0 / h / h) + D;
		C[i] = K * (1.0 / h / h);
	}

	mat E(n, n, fill::zeros);
	vec Fl(n, fill::zeros);
	for (int i = 0; i < n; i++) {
		E(i, i) = left(i);
		if (left(i) == 1)
			continue;
		vec l(n, fill::zeros);
		l(i) = P_LEFT;
		Fl = Fl + A[0] * l;
	}
	mat Bl(n, n, fill::zeros);
	Bl = A[0] * E;

	vec Fr(n, fill::zeros);
	for (int i = 0; i < n; i++) {
		E(i, i) = right(i);
		if (right(i) == 1)
			continue;
		vec r(n, fill::zeros);
		r(i) = P_RIGHT;
		Fr = Fr + C[l - 1] * r;
	}
	mat Br(n, n, fill::zeros);
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
	vector<vec> u(L, vec(n, fill::zeros));
	vector<vec> u1(L, vec(n, fill::zeros));

	cout << "The dimension of prodlem: " << n << endl;
	cout << "The length = " << L << ", the time = " << TIME << endl;
	cout << "Space step = " << h << ", time step = " << t << endl;

	direct_problem(n, u, u1);

	cout << "DONE" << endl;
	return 0;
}

