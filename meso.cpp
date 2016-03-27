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

#define PRNT 1
#define TIME 20//15000
#define STABLE 1000
#define t 1.0
#define L 200
#define h 1.0
#define P_LEFT 1.0
#define P_RIGHT 0.0
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
		  std::vector<vector>& v_med)
{
	for (int i = 1; i < L; i++) {
		v_med[i] = p[i - 1] - p[i];
		v_med[i] = v_med[i] * (1.0 / h);
	}
	v_med[0] = v_med[1];
	v_med[L] = v_med[L - 1];
}


void start_cond_con(std::vector<vector>& C, uint dif)
{
	uint n = C[0].size();
	vector left(n);
	get_left_edge(left);
	/* free edge -> close C=0, fixed edge -> pressure C=1 */
	for (int i = 0; i < n; i++)
		C.at(dif)(i) = 1;
};

void concentration(uint n, std::vector<vector>& p1, uint dif)
{
	char buf[256];
	double time = (double)TIME / t;
	corner turn(n, L);

	double hh[2] = {h, 10.0}; /* Шаг сетки */
	int N[2] = {L, n}; /* Число точек расчетной области по осям */

	std::vector<vector> v_med(L + 1);
	std::vector<vector> c(L + 1);
	std::vector<vector> c1(L + 1);

	/* due to the specific realization of std containers */
	for (int i = 0; i <= L; i++) {
		v_med[i].init(n);
		c[i].init(n);
		c1[i].init(n);
	}

	start_cond_con(c, dif);

	grad(n, p1, v_med);
	// prepare to calculate next step by corner
	double c_time = (double)CURANT * h / maximum(v_med);
	uint steps = (double)t / c_time + 1;

	for (int j = 0; j < time; j++) {
		for (int i = 0; i < n; i++)
			c1[0](i) = 0; //c1[1](i);

		turn.calculate(v_med, p1, c, c1, c_time);

		for (int i = 0; i < n; i++)
			c1[L](i) = c1[L - 1](i);

		N[0] = L + 1;

		if (j % PRNT == 0) {
			sprintf(buf, "res/%03d_%06d.vtk",dif, j);
			write_to_vtk2d(c, buf, "concentration", N, hh);
		}
		c = c1;
	}
};

void pressure(uint n, std::vector<vector>& p, std::vector<vector>& p1)
{
	char buf[256];
	double stable = (double)STABLE / t;
	uint l = L - 2; // exept edge points
	matrix Al(n, 1.0/t);
	matrix K(n);
	get_K(K);
	matrix D(n);
	get_D(D);

	start_cond(p);

	double hh[2] = {h, 10.0}; /* Шаг сетки */
	int N[2] = {L, n}; /* Число точек расчетной области по осям */

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

	// wait for pressure becomes stable
	for (int i = 0; i < stable; i++) {
		for (int j = 0; j < l; j++)
			F[j] = Al * p[j + 1] * (-1);    // u[i + 1] -> start from 1-st index, not 0
		F[0] = F[0] - Fl;
		F[l - 1] = F[l - 1] - Fr;

		gone.calculate(p1, A, B, C, F);
		for (int j = 0; j < n; j++) {
			p1[0](j) = left(j) ? p1[1](j) : P_LEFT;
			p1[L - 1](j) = right(j) ? p1[L - 2](j) : P_RIGHT;
		}
		p = p1;
	}

	N[0] = L;
	sprintf(buf, "res/pressure.vtk");
	write_to_vtk1(p, buf, n, L);

	for (uint i = 1; i < 200; i += 2)
		concentration(n, p, i);

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

	pressure(n, u, u1);

	std::cout << "DONE" << std::endl;
	return 0;
}

