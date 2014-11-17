/*
 * @file main.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <iostream>
#include <vector>

#include <math.h>

#include "lalgebra.h"
#include "solvers.h"

#include "config.h"

#define PRNT 10

ptype power(ptype base, ptype exponent)
{
	ptype half_pow = 0.0, int_part = 0.0;
	ptype fr_part = modf(exponent, &int_part);

	if ((fr_part != 0) && (fr_part != 0.5)) {
		return pow(base, exponent);
	} else if (fr_part == 0.5) {
		return (sqrt(base) * power(base, int_part));
	}

	if (exponent == 0)
		return (ptype)1;
	else if (exponent < 0)
		return (ptype)1 / power(base, -exponent);
	else if (fmod(exponent, 2) == 0) {
		half_pow = power(base, exponent / (ptype)2);
		return half_pow * half_pow;
	} else
		return base * power(base, exponent - 1);
};

void start_cond(std::vector<vector>& u)
{
	uint n = u[0].size();
	for (int i = 0; i < n; i++) {
		u.at(0)(i) = P_LEFT;
		u.at(L - 1)(i) = P_RIGHT;
	}
};

void write_to_vtk2(std::vector<vector>& u, const char *path, uint n)
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
	fclose(f);
}


int main(int argc, char **argv)
{
	// the problem size
	uint n = power(P, GAMMA - 1);
	char buf[256];

	progonka method(n);

	std::vector<vector> u(L);
	std::vector<vector> u1(L);

	/* due to the specific realization of std containers */
	for (int i = 0; i < L; i++) {
		u[i].init(n);
		u1[i].init(n);
	}

	start_cond(u);

	for (int i = 0; i < ((double)TIME / t); i++) {
		if (i % PRNT == 0) {
			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk2(u, buf, n);
		}

		method.calculate(u, u1);
		u = u1;
	}
	std::cout << "DONE" << std::endl;
	return 0;
}

