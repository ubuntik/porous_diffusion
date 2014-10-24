#include <iostream>
#include <cblas.h>

#include "config.h"

#define MTX_SZ (sizeof(double) * GAMMA * GAMMA)

static matrix A_ini;
static matrix A;
static matrix B;
static matrix C;

static matrix E;
static p_vector e;

void init_matrix()
{
	matrix K_ini;
	matrix D_ini;

	bzero(A, MTX_SZ);
	bzero(A_ini, MTX_SZ);
	bzero(B, MTX_SZ);
	bzero(C, MTX_SZ);

	for (int i = 0; i < GAMMA; i++) {
		for (int j = 0; j < GAMMA; j++) {
			A_ini[i][j] = (i == j) ? 1 : 0;
			K_ini[i][j] = (i == j) ? 0.1 : 0.001;
			D_ini[i][j] = (i == j) ? 1 : P;
			E[i][j] = (i == j) ? 1 : 0;
		}
		e[i] = 1;
	}

	for (int i = 0; i < GAMMA; i++) {
		for (int j = 0; j < GAMMA; j++) {
			A_ini[i][j] = A_ini[i][j] / t;
			A[i][j] = - K_ini[i][j] / h / h;
			B[i][j] = A_ini[i][j] + 2 * K_ini[i][j] / h / h + D_ini[i][j];
			C[i][j] = - K_ini[i][j] / h / h;
		}
	}
}

void progonka(const p_vector *u, p_vector *u1)
{
	/* right part */
	p_vector Fi;

	/* coefficients */
	matrix *Ps = (matrix *)calloc(L, MTX_SZ);
	p_vector *Qs = (p_vector *)calloc(L, sizeof(p_vector));

	/* edge conditions */
	memcpy(&u1[0], &u[0], sizeof(double) * GAMMA);
	memcpy(&u1[L - 1], &u[L - 1], sizeof(double) * GAMMA);

	init_matrix();

	/* left edge */
	memcpy(Ps[0], &E, MTX_SZ);

	/* There */
	for (int i = 0; i < L; i++) {
//		Fi = (matrix_A_ini * u[i]);
//		matrix G = (A * Ps[i - 1] + B)(^ -1);
//		Ps[i] = - G * C;
//		Qs[i] = G * (Fi - A * Qs[i - 1]);
	}

	/* and Back Again */
	for (int i = L - 2; i >= 0; i--) {
//		u1[i] = Ps[i] * u1[i + 1] + Qs[i];
	}

	free(Ps);
	free(Qs);
}

int main(int argc, char **argv)
{
	p_vector *u = NULL;
	p_vector *u1 = NULL;
	p_vector *temp = NULL;

	u = (p_vector *)calloc(L, sizeof(double) * GAMMA);
	u1 = (p_vector *)calloc(L, sizeof(double) * GAMMA);
	if (u == NULL || u1 == NULL) {
		fprintf(stderr, "Cannot alloch enough memory\n");
		exit(-1);
	}

	start_cond(u);
	print_vec(u);
	std::cout << "************************************" << std::endl;
	print_vec(u1);

	std::cout << "************************************" << std::endl;
	std::cout << "************************************" << std::endl;

	for (int i = 0; i < TIME; i++) {
		progonka(u, u1);
		temp = u;
		u = u1;
		u1 = temp;
	}
	print_vec(u);
	std::cout << "************************************" << std::endl;
	print_vec(u1);
	return 0;
}

