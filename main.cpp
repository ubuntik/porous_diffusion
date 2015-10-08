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

#define PRNT 20

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
			u.at(i)(j) = START;
		u.at(L - 1)(j) = P_RIGHT;
	}
};

void start_cond(std::vector<vector>& u)
{
	uint n = u[0].size();
	for (int i = 0; i < n; i++) {
		u.at(0)(i) = P_LEFT;
		u.at(L - 1)(i) = P_RIGHT;
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

void simple_ODE(uint n, std::vector<vector>& u, std::vector<vector>& u1)
{
	char buf[256];
	/* This structure keep all data (not using in general) */
	std::vector<std::vector<vector> > u_all((double)TIME / t, std::vector<vector>(L));

	/* due to the specific realization of std containers */
	for (int i = 0; i < ((double)TIME / t); i++)
		for (int j = 0; j < L; j++)
			u_all[i][j].init(n);

	start_cond_1(u);

	progonka method(n);

	for (int i = 0; i < ((double)TIME / t); i++) {
		if (i % PRNT == 0) {
			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk1(u, buf, n);
		}

		method.calculate(u, u1);

		u = u1;
		u_all[i] = u1;
	}

	/* This printing is u(TIME) throughout Length exept edges
	 * (not using in general)
	 */
	for (int i = 1; i < L - 1; i++) {
		sprintf(buf, "res/time_%06d.vtk", i);
		write_to_vtk3(u_all, buf, n, i);
	}
}

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
	for (int i = 0; i < n; i++) {
		C.at(0)(i) = 1;
		C.at(L)(i) = 0;
	}
};

void tracer_problem(uint n, std::vector<vector>& p, std::vector<vector>& p1)
{
	char buf[256];
	start_cond(p);
	double hh[2] = {h, 10}; /* Шаг сетки */
	double o[2] = {0.0, 0.0}; /* Положение в пространстве */
	int N[2] = {L, n}; /* Число точек расчетной области по осям */

	std::vector<vector> v_med(L + 1);
	std::vector<vector> C(L + 1);
	std::vector<vector> C1(L + 1);
	std::vector<vector> C_slice(L + 1);

	/* due to the specific realization of std containers */
	for (int i = 0; i <= L; i++) {
		v_med[i].init(n);
		C[i].init(n);
		C1[i].init(n);
		C_slice[i].init(n);
	}

	start_cond_con(C);

	progonka gone(n);
	corner turn(n);

	for (int i = 0; i < ((double)TIME / t); i++) {
		gone.calculate(p, p1);
		grad(n, p, p1, v_med);
		// prepare to calculate next step by corner
		double c_time = (double)CURANT * h / maximum(v_med);
		uint steps = (double)t / c_time + 1;

		//printf(">>> max = %g, time = %g, steps = %d\n", maximum(v_med), c_time, steps);

		for (int j = 0; j < steps; j++) {
			turn.calculate(v_med, C, C1, c_time);
		}

		if (i % PRNT == 0) {

			N[0] = L;

			sprintf(buf, "res/pres_%06d.vtk", i);
			write_to_vtk2d(p, buf, "pressure", N, o, hh);

			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk1(p, buf, n);

			N[0] = L + 1;

			sprintf(buf, "res/velo_%06d.vtk", i);
			write_to_vtk2d(v_med, buf, "velocity", N, o, hh);

			sprintf(buf, "res/conc_%06d.vtk", i);
			write_to_vtk2d(C, buf, "concentration", N, o, hh);

			for (int k = 0; k < L + 1; k++) {
				C_slice[k] = C[k].mult(v_med[k]) * (1.0 / v_med[k].sum());
			}
			sprintf(buf, "res/scon_%06d.vtk", i);
			write_to_vtk1(C_slice, buf, n);
		}

		C = C1;
		p = p1;
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

//	direct_problem(n, u, u1);
//	simple_ODE(n, u, u1);
//	compare_accurate(n, u, u1);
	tracer_problem(n, u, u1);

	std::cout << "DONE" << std::endl;
	return 0;
}

