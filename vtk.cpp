/*
 * @file vtk.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <iostream>
#include <math.h>

#include "vtk.h"

using namespace arma;

void write_plain_data(vector<double>& u, const char *path, uint l)
{
	FILE *f = fopen(path, "w+");
	if (f == NULL) {
		cerr << "Cannot open file " << path << endl;
		exit(-1);
	}
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk1\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %g 1 1\n", (double)l);
	fprintf(f, "SPACING 1.0 0.0 0.0\n");
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %g\n", (double)l);
	fprintf(f, "SCALARS %s float 1\n", "wc");
	fprintf(f, "LOOKUP_TABLE %s_table\n", "wc");

	for (int i = 0; i < l; i++) {
		fprintf(f, "%0.6f\n", u[i]);
	}
	fclose(f);
}

void write_to_vtk1(vector<vec>& u, const char *path, uint n, uint l)
{
	FILE *f = fopen(path, "w+");
	if (f == NULL) {
		cerr << "Cannot open file " << path << endl;
		exit(-1);
	}
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk1\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %g 1 1\n", (double)l);
	fprintf(f, "SPACING 1.0 0.0 0.0\n");
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %g\n", (double)l);

	for (int j = 0; j < n; j++) {
		fprintf(f, "SCALARS %s%d float 1\n", "p", j);
		fprintf(f, "LOOKUP_TABLE %s%d_table\n", "p", j);
		for (int i = 0; i < l; i++) {
			fprintf(f, "%g\n", u[i](j));
		}
	}
	fclose(f);
}

void write_to_vtk2(vector<vec>& u, const char *path, uint n, uint l, double crt_time)
{
	FILE *f = fopen(path, "w+");
	if (f == NULL) {
		cerr << "Cannot open file " << path << endl;
		exit(-1);
	}
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk2\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %g 1 1\n", (double)l);
	fprintf(f, "SPACING 1.0 0.0 0.0\n");
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %g\n", (double)l);

	for (int j = 0; j < n; j++) {
		fprintf(f, "SCALARS %s%d float 1\n", "p", j);
		fprintf(f, "LOOKUP_TABLE %s%d_table\n", "p", j);
		for (int i = 0; i < l; i++) {
			fprintf(f, "%g\n", u[i](j));
		}
	}
	fprintf(f, "SCALARS accurate float 1\n");
	fprintf(f, "LOOKUP_TABLE accurate_table\n");
	for (int i = 0; i < l; i++) {
		double point = erf((double)i / 2.0 / sqrt(crt_time));
		fprintf(f, "%g\n", point);
	}
	fclose(f);
}

void write_to_vtk3(vector< vector<vec> >& u, const char *path, uint n, uint l, uint t)
{
	FILE *f = fopen(path, "w+");
	if (f == NULL) {
		cerr << "Cannot open file " << path << endl;
		exit(-1);
	}
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk3\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %g 1 1\n", (double)t);
	fprintf(f, "SPACING 1.0 0.0 0.0\n");
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %g\n", (double)t);

	for (int j = 0; j < n; j++) {
		fprintf(f, "SCALARS %s%d float 1\n", "p", j);
		fprintf(f, "LOOKUP_TABLE %s%d_table\n", "p", j);
		for (int i = 0; i < t; i++) {
			fprintf(f, "%g\n", u[i][l](j));
		}
	}
	fclose(f);
}

void write_to_vtk2d(vector<vec>& u, const char *path, const char *data,
			const int N[2], const double hh[2])
{
	FILE *f = fopen(path, "w+");
	if (f == NULL) {
		cerr << "Cannot open file " << path << endl;
		exit(-1);
	}
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk2d\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %d %d 1\n", N[0], N[1]);
	fprintf(f, "SPACING %g %g 0.0\n", hh[0], hh[1]);
	fprintf(f, "ORIGIN 0.0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %d\n", N[0] * N[1]);

	fprintf(f, "SCALARS %s float 1\n", data);
	fprintf(f, "LOOKUP_TABLE %s_table\n", data);
	for (int i2 = 0; i2 < N[1]; i2++) {
		for (int i1 = 0; i1 < N[0]; i1++) {
			fprintf(f, "%0.6f\n", u[i1](i2));
		}
	}
        fclose(f);
}

