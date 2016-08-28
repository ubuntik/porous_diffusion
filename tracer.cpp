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

#include "solvers.h"
#include "vtk.h"
#include "matrixes.h"

#define TR_START 10

#define PRNT 1
#define TIME 100
#define t 0.5
#define L 100
#define h 1.0
#define P_LEFT 1.0
#define P_RIGHT 0.0
// coefficient for tracer calculation (0 <= THETA <= 1)
#define THETA 0.5
// curant
#define CURANT 0.001

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void start_cond(vector<vec>& u)
{
	uint n = u[0].size();
	for (int i = 0; i < n; i++) {
		u.at(0)(i) = P_LEFT;
		u.at(L - 1)(i) = P_RIGHT;
	}
};

double maximum(vector<vec>& v_med)
{
	double ret = 0;
	for (int i = 0; i <= L; i++)
		for (int j = 0; j < v_med[i].size(); j++)
			ret = MAX(fabs(v_med[i](j)), ret);
	return ret;
}

void grad(uint n, vector<vec>& p,
		  vector<vec>& p1,
		  vector<vec>& v_med)
{
	vec v(n, fill::zeros);
	vec v1(n, fill::zeros);

	for (int i = 1; i < L; i++) {
		v = p[i - 1] - p[i];
		v = v * (1.0 / h);
		v1 = p1[i - 1] - p1[i];
		v1 = v1 * (1.0 / h);
		v_med[i] = v1 * THETA + v * (1.0 - THETA);
	}
	v_med[0] = v_med[1];
	v_med[L] = v_med[L - 1];
}


void start_cond_con(vector<vec>& C)
{
	uint n = C[0].size();
	vec left(n, fill::zeros);
	get_left_edge(left);
	/* free edge -> close C=0, fixed edge -> pressure C=1 */
	for (int i = 0; i < n; i++)
		C.at(TR_START)(i) = 1; //left(i) ? 0 : 1;
};

void tracer_problem(uint n, vector<vec>& p, vector<vec>& p1)
{
	char buf[256];
	double time = (double)TIME / t;
	uint l = L - 2; // exept edge points

	start_cond(p);

	double hh[2] = {h, 40.0}; /* Шаг сетки */
	int N[2] = {L, n}; /* Число точек расчетной области по осям */

	vector<vec> v_med(L + 1, vec(n, fill::zeros));
	vector<vec> c(L + 1, vec(n, fill::zeros));
	vector<vec> c1(L + 1, vec(n, fill::zeros));
	vector<vec> c_slice(L + 1, vec(n, fill::zeros));

	vector<mat> Al(L, mat(n, n, fill::zeros));
	get_Al(Al);
	vector<mat> K(L, mat(n, n, fill::zeros));
	get_K(K);
	vector<mat> D(L, mat(n, n, fill::zeros));
	get_D(D);

	start_cond_con(c);

	vec left(n, fill::zeros);
	get_left_edge(left);
	vec right(n, fill::zeros);
	get_right_edge(right);

	vector<mat> A(l, mat(n, n, fill::zeros));
	vector<mat> B(l, mat(n, n, fill::zeros));
	vector<mat> C(l, mat(n, n, fill::zeros));
	vector<vec> F(l, vec(n, fill::zeros));

	for (int i = 0; i < l; i++) {
		mat Kp(n, n, fill::zeros);
		Kp = (2 * K[i + 1] * K[i + 2]) * ((K[i + 1] + K[i + 2]).i());
		mat Km(n, n, fill::zeros);
		Km = (2 * K[i] * K[i + 1]) * ((K[i] + K[i + 1]).i());
		A[i] = Km * (1.0 / h / h);
		B[i] = Al[i + 1] * (1.0 / t) +
			Km * (1.0 / h / h) +
			Kp * (1.0 / h / h) + D[i + 1];
		C[i] = Kp * (1.0 / h / h);
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

	progonka gone(n, l);
	secondord turn(n, L, &K, &D);

	for (int i = 0; i < time; i++) {
		for (int i = 0; i < l; i++)
			// u[i + 1] -> start from 1-st index, not 0
			F[i] = Al[i + 1] * p[i + 1] * (-1 / t);
		F[0] = F[0] - Fl;
		F[l - 1] = F[l - 1] - Fr;

		gone.calculate(p1, A, B, C, F);
		for (int i = 0; i < n; i++) {
			p1[0](i) = left(i) ? p1[1](i) : P_LEFT;
			p1[L - 1](i) = right(i) ? p1[L - 2](i) : P_RIGHT;
		}

		grad(n, p, p1, v_med);
		// prepare to calculate next step by corner
		double c_time = (double)CURANT * h / maximum(v_med);
		uint steps = (double)t / c_time + 1;

		for (int i = 0; i < n; i++)
			c1[0](i) = c1[1](i); //left(i) ? 0 : 1;

		for (int j = 0; j < steps; j++) {
			turn.calculate(v_med, p1, c, c1, c_time);
		}

		for (int i = 0; i < n; i++)
			c1[L](i) = c1[L - 1](i);

		if (i % PRNT == 0) {
			N[0] = L;

			sprintf(buf, "res/pres_%06d.vtk", i);
			write_to_vtk2d(p, buf, "pressure", N, hh);

			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk1(p, buf, n, L);

			N[0] = L + 1;

			sprintf(buf, "res/conc_%06d.vtk", i);
			write_to_vtk2d(c, buf, "concentration", N, hh);
		}

		c = c1;
		p = p1;
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

	tracer_problem(n, u, u1);

	cout << "DONE" << endl;
	return 0;
}

