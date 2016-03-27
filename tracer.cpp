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
#include "matrixes.h"

#define TR_START 10

#define PRNT 100
#define TIME 20000
#define t 1.0
#define L 500
#define h 1.0
#define P_LEFT 1.0
#define P_RIGHT 0.0
// coefficient for tracer calculation (0 <= THETA <= 1)
#define THETA 0.5
// curant
#define CURANT 0.001

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void start_cond(std::vector<vector>& u)
{
	uint n = u[0].size();
	for (int i = 0; i < n; i++) {
		u.at(0)(i) = P_LEFT;
		u.at(L - 1)(i) = P_RIGHT;
	}
};

ptype maximum(std::vector<vector>& v_med)
{
	ptype ret = 0;
	for (int i = 0; i <= L; i++)
		for (int j = 0; j < v_med[i].size(); j++)
			ret = MAX(fabs(v_med[i](j)), ret);
	return ret;
}

void grad(uint n, std::vector<vector>& p,
		  std::vector<vector>& p1,
		  std::vector<vector>& v_med)
{
	vector v(n);
	vector v1(n);

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


void start_cond_con(std::vector<vector>& C)
{
	uint n = C[0].size();
	vector left(n);
	get_left_edge(left);
	/* free edge -> close C=0, fixed edge -> pressure C=1 */
	for (int i = 0; i < n; i++)
		C.at(TR_START)(i) = 1; //left(i) ? 0 : 1;
};

void tracer_problem(uint n, std::vector<vector>& p, std::vector<vector>& p1)
{
	char buf[256];
	double time = (double)TIME / t;
	uint l = L - 2; // exept edge points
	matrix Al(n, 1.0/t);
	matrix K(n);
	get_K(K);
	matrix D(n);
	get_D(D);

	start_cond(p);

	double hh[2] = {h, 40.0}; /* Шаг сетки */
	int N[2] = {L, n}; /* Число точек расчетной области по осям */

	std::vector<vector> v_med(L + 1);
	std::vector<vector> c(L + 1);
	std::vector<vector> c1(L + 1);
	std::vector<vector> c_slice(L + 1);

	/* due to the specific realization of std containers */
	for (int i = 0; i <= L; i++) {
		v_med[i].init(n);
		c[i].init(n);
		c1[i].init(n);
		c_slice[i].init(n);
	}

	start_cond_con(c);

	vector left(n);
	get_left_edge(left);
	vector right(n);
	get_right_edge(right);

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

	matrix E(n);
	vector Fl(n);
	for (int i = 0; i < n; i++) {
		E(i, i) = left(i);
		if (left(i) == 1)
			continue;
		vector l(n);
		l(i) = P_LEFT;
		Fl = Fl + A[0] * l;
	}
	matrix Bl(n);
	Bl = A[0] * E;

	vector Fr(n);
	for (int i = 0; i < n; i++) {
		E(i, i) = right(i);
		if (right(i) == 1)
			continue;
		vector r(n);
		r(i) = P_RIGHT;
		Fr = Fr + C[l - 1] * r;
	}
	matrix Br(n);
	Br = C[l - 1] * E;

	B[0] = B[0] - Bl;
	B[l - 1] = B[l - 1] - Br;

	progonka gone(n, l);
	corner turn(n, L);

	for (int i = 0; i < time; i++) {
		for (int i = 0; i < l; i++)
			F[i] = Al * p[i + 1] * (-1);    // u[i + 1] -> start from 1-st index, not 0
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

	tracer_problem(n, u, u1);

	std::cout << "DONE" << std::endl;
	return 0;
}

