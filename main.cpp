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

#include "config.h"

#define START 1
#define PRNT 1

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

void write_to_vtk1(std::vector<vector>& u, const char *path, uint n)
{
	FILE *f = fopen(path, "w");
	if (f == NULL) {
		std::cerr << "Cannot open file " << path << std::endl;
		exit(-1);
	}
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk1\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %d 1 1\n", L);
	fprintf(f, "SPACING %g 0.0 0.0\n", (double)h);
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %d\n", L);

	for (int j = 0; j < n; j++) {
		fprintf(f, "SCALARS %s%d float 1\n", "p", j);
		fprintf(f, "LOOKUP_TABLE %s%d_table\n", "p", j);
		for (int i = 0; i < L; i++) {
			fprintf(f, "%f\n", u[i](j));
		}
	}
	fclose(f);
}

void write_to_vtk2(std::vector<vector>& u, const char *path, uint n, double crt_time)
{
	FILE *f = fopen(path, "w");
	if (f == NULL) {
		std::cerr << "Cannot open file " << path << std::endl;
		exit(-1);
	}
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk2\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %d 1 1\n", L);
	fprintf(f, "SPACING %g 0.0 0.0\n", (double)h);
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %d\n", L);

	for (int j = 0; j < n; j++) {
		fprintf(f, "SCALARS %s%d float 1\n", "p", j);
		fprintf(f, "LOOKUP_TABLE %s%d_table\n", "p", j);
		for (int i = 0; i < L; i++) {
			fprintf(f, "%f\n", u[i](j));
		}
	}
	fprintf(f, "SCALARS accurate float 1\n");
	fprintf(f, "LOOKUP_TABLE accurate_table\n");
	for (int i = 0; i < L; i++) {
		ptype point = erf((ptype)i / 2.0 / sqrt(crt_time));
		fprintf(f, "%f\n", point);
	}
	fclose(f);
}

void write_to_vtk3(std::vector< std::vector<vector> >& u, const char *path, uint n, uint l)
{
	FILE *f = fopen(path, "w");
	if (f == NULL) {
		std::cerr << "Cannot open file " << path << std::endl;
		exit(-1);
	}
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk3\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %g 1 1\n", ((double)TIME / t));
	fprintf(f, "SPACING %g 0.0 0.0\n", (double)t);
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %g\n", ((double)TIME / t));

	for (int j = 0; j < n; j++) {
		fprintf(f, "SCALARS %s%d float 1\n", "p", j);
		fprintf(f, "LOOKUP_TABLE %s%d_table\n", "p", j);
		for (int i = 0; i < ((double)TIME / t); i++) {
			fprintf(f, "%f\n", u[i][l](j));
		}
	}
	fclose(f);
}

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
	start_cond_1(u);

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

	direct_problem(n, u, u1);
//	simple_ODE(n, u, u1);
//	compare_accurate(n, u, u1);

	std::cout << "DONE" << std::endl;
	return 0;
}

