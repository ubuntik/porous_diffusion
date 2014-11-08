#include <iostream>
#include <vector>

#include "lalgebra.h"
#include "solvers.h"

#include "config.h"

#define PRNT 2

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
	uint n = GAMMA;
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

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < L; j++) {
			std::cout << u[j](i) << " ";
		}
		std::cout << std::endl;
	}

	for (int i = 0; i < TIME; i++) {
		method.calculate(u, u1);
		if (i % PRNT == 0) {
			sprintf(buf, "res/data_%06d.vtk", i);
			write_to_vtk2(u1, buf, n);
		}
		u = u1;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < L; j++) {
				std::cout << u[j](i) << " ";
			}
			std::cout << std::endl;

	}

	}
	return 0;
}

