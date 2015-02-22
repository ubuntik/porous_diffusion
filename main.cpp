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

void start_cond_1(std::vector<vector>& u)
{
	for (int i = 0; i < u.size(); i++)
		for (int j = 0; j < u.at(i).size(); j++)
			u.at(i)(j) = START;
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

void write_to_vtk2(std::vector< std::vector<vector> >& u, const char *path, uint n, uint l)
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
	fprintf(f, "DIMENSIONS %d 1 1\n", TIME);
	fprintf(f, "SPACING %g 0.0 0.0\n", (double)t);
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %d\n", TIME);

	for (int j = 0; j < n; j++) {
		fprintf(f, "SCALARS %s%d float 1\n", "p", j);
		fprintf(f, "LOOKUP_TABLE %s%d_table\n", "p", j);
		for (int i = 0; i < TIME; i++) {
			fprintf(f, "%f\n", u[i][l](j));
		}
	}
	fclose(f);
}

int main(int argc, char **argv)
{
	// the problem size
	uint n = power(P, GAMMA - 1);
	char buf[256];

	std::cout << "The dimension of prodlem: " << n << std::endl;
	std::cout << "The length = " << L << ", the time = " << TIME << std::endl;
	std::cout << "Space step = " << h << ", time step = " << t << std::endl;

	progonka method(n);

	std::vector<vector> u(L);
	std::vector<vector> u1(L);

	/* due to the specific realization of std containers */
	for (int i = 0; i < L; i++) {
		u[i].init(n);
		u1[i].init(n);
	}


	/* This structure keep all data (not using in general) */
/*
	std::vector<std::vector<vector> > u_all(TIME, std::vector<vector>(L));

	for (int i = 0; i < TIME; i++)
		for (int j = 0; j < L; j++)
			u_all[i][j].init(n);
*/

//	start_cond_1(u);
	start_cond(u);

	for (int i = 0; i < ((double)TIME / t); i++) {
		if (i % PRNT == 0) {
			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk1(u, buf, n);
		}

		method.calculate(u, u1);
		u = u1;
//		u_all[i] = u1;
	}

	/* This printing is u(TIME) throughout Length exept edges
	 * (not using in general)
	 */
/*
	for (int i = 1; i < L - 1; i++) {
		sprintf(buf, "res/time_%06d.vtk", i);
		write_to_vtk2(u_all, buf, n, i);
	}
*/

	std::cout << "DONE" << std::endl;
	return 0;
}

