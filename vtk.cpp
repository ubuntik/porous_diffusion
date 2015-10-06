/*
 * @file main.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <iostream>
#include <math.h>

#include "vtk.h"

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
	fprintf(f, "DIMENSIONS %g 1 1\n", (double)L);
	fprintf(f, "SPACING %g 0.0 0.0\n", (double)h);
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %g\n", (double)L);

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
	fprintf(f, "DIMENSIONS %g 1 1\n", (double)L);
	fprintf(f, "SPACING %g 0.0 0.0\n", (double)h);
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %g\n", (double)L);

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

void write_to_vtk2d(std::vector<vector>& u, const char *path, const char *data,
			const int N[2], const double o[2], const double hh[2])
{
	FILE *f;
	f = fopen(path, "w");
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk2d\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %d %d 1\n", N[0], N[1]);
	fprintf(f, "SPACING %f %f 0.0\n", hh[0], hh[1]);
	fprintf(f, "ORIGIN %f %f 0.0\n", o[0], o[1]);
	fprintf(f, "POINT_DATA %d\n", N[0] * N[1]);

	fprintf(f, "SCALARS %s float 1\n", data);
	fprintf(f, "LOOKUP_TABLE %s_table\n", data);
	for (int i2 = 0; i2 < N[1]; i2++) {
		for (int i1 = 0; i1 < N[0]; i1++) {
			fprintf(f, "%f\n", u[i1](i2));
		}
	}
        fclose(f);
}

