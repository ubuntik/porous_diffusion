#include <iostream>
#include <vector>

#include "lalgebra.h"
#include "solvers.h"

#include "config.h"

void start_cond(std::vector<vector>& u)
{
	uint n = u[0].size();
	for (int i = 0; i < n; i++) {
		u[0](i) = P_LEFT;
		u[N - 1](i) = P_RIGHT;
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
	fprintf(f, "DIMENSIONS %d 1 1\n", N);
	fprintf(f, "SPACING %g 0.0 0.0\n", (double)h);
	fprintf(f, "ORIGIN 0 0.0 0.0\n");
	fprintf(f, "POINT_DATA %d\n", N);

	for (int j = 0; j < n; j++) {
		fprintf(f, "SCALARS %s%d float 1\n", "p", j);
		fprintf(f, "LOOKUP_TABLE %s%d_table\n", "p", j);
		for (int i = 0; i < N; i++) {
			fprintf(f, "%f\n", u[i](j));
		}
	}
	fclose(f);
}


int main(int argc, char **argv)
{
	// the problem size
	uint n = GAMMA;
	char buf[256];

	progonka method(n);

	std::vector<vector> u(N, vector(n));
	std::vector<vector> u1(N, vector(n));

	start_cond(u);

//	for (int i = 0; i < TIME; i++) {
	for (int i = 0; i < 2; i++) {
		method.calculate(u, u1);
		if (i % 10 == 0) {
			sprintf(buf, "data_%06d.vtk", i);
			write_to_vtk2(u1, buf, n);
		}
		u = u1;
	}
	return 0;
}

