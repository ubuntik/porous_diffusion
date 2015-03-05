/*
 * @file solvers.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <assert.h>
#include <iostream>
#include <vector>
#include <math.h>

#include "lalgebra.h"
#include "solvers.h"

int p_norm(ptype num)
{
	int cnt = 0;
	while((num / P) > 1) {
		num /= P;
		cnt++;
	}
	return cnt;
};

progonka::progonka(uint n_size)
{
	n = n_size;
	assert(n >= 0);

	A_ini = new matrix(n, 1.0);
	A = new matrix(n);
	B = new matrix(n);
	C = new matrix(n);
	matrix K_ini(n);
	matrix D_ini(n);

	K_ini(0, 0) = 700.72; K_ini(1, 0) = -90.791; K_ini(2, 0) = -51.8508; K_ini(3, 0) = -20.9356;
	K_ini(4, 0) = -16.6301; K_ini(5, 0) = -36.7135; K_ini(6, 0) = -12.6038; K_ini(7, 0) = -11.0965;

	K_ini(0, 1) = -90.4545; K_ini(1, 1) = 410.508; K_ini(2, 1) = 10.7879; K_ini(3, 1) = -6.02379;
	K_ini(4, 1) = 7.86314; K_ini(5, 1) = 4.83955; K_ini(6, 1) = 13.3205; K_ini(7, 1) = 6.58391;

	K_ini(0, 2) = -51.7445; K_ini(1, 2) = 10.9115; K_ini(2, 2) = 162.8; K_ini(3, 2) = 24.564;
	K_ini(4, 2) = 12.9702; K_ini(5, 2) = 16.7895; K_ini(6, 2) = 16.4923; K_ini(7, 2) = 20.6626;

	K_ini(0, 3) = -20.952; K_ini(1, 3) = -5.76781; K_ini(2, 3) = 24.6143; K_ini(3, 3) = 110.027;
	K_ini(4, 3) = 7.23522; K_ini(5, 3) = 14.8657; K_ini(6, 3) = 10.6427; K_ini(7, 3) = 18.3762;

	K_ini(0, 4) = -16.8055; K_ini(1, 4) = 7.7814; K_ini(2, 4) = 13.0678; K_ini(3, 4) = 7.29436;
	K_ini(4, 4) = 85.9633; K_ini(5, 4) = 10.5721; K_ini(6, 4) = 10.8133; K_ini(7, 4) = 13.1033;

	K_ini(0, 5) = -36.7031; K_ini(1, 5) = 4.94953; K_ini(2, 5) = 16.7641; K_ini(3, 5) = 14.8475;
	K_ini(4, 5) = 10.5869; K_ini(5, 5) = 91.9372; K_ini(6, 5) = 13.9441; K_ini(7, 5) = 18.5697;

	K_ini(0, 6) = -12.6047; K_ini(1, 6) = 13.253; K_ini(2, 6) = 16.6121; K_ini(3, 6) = 10.6599;
	K_ini(4, 6) = 10.821; K_ini(5, 6) = 13.9652; K_ini(6, 6) = 85.7859; K_ini(7, 6) = 23.7649;

	K_ini(0, 7) = -11.1; K_ini(1, 7) = 6.64; K_ini(2, 7) = 20.7; K_ini(3, 7) = 18.4;
	K_ini(4, 7) = 13.1; K_ini(5, 7) = 18.6; K_ini(6, 7) = 23.8; K_ini(7, 7) = 68.8;

	D_ini(0, 0) = 0.00579; D_ini(1, 0) = -0.00135; D_ini(2, 0) = -0.00060; D_ini(3, 0) = -0.00051;
	D_ini(4, 0) = -0.00054; D_ini(5, 0) = -0.00072; D_ini(6, 0) = -0.00097; D_ini(7, 0) = -0.00111;

	D_ini(0, 1) = -0.00135; D_ini(1, 1) = 0.00462; D_ini(2, 1) = -0.00040; D_ini(3, 1) = -0.00038;
	D_ini(4, 1) = -0.00041; D_ini(5, 1) = -0.00053; D_ini(6, 1) = -0.00070; D_ini(7, 1) = -0.00083;

	D_ini(0, 2) = -0.00060; D_ini(1, 2) = -0.00040; D_ini(2, 2) = 0.00242; D_ini(3, 2) = -0.00018;
	D_ini(4, 2) = -0.00019; D_ini(5, 2) = -0.00028; D_ini(6, 2) = -0.00035; D_ini(7, 2) = -0.00043;

	D_ini(0, 3) = -0.00051; D_ini(1, 3) = -0.00038; D_ini(2, 3) = -0.00018; D_ini(3, 3) = 0.00219;
	D_ini(4, 3) = -0.00017; D_ini(5, 3) = -0.00025; D_ini(6, 3) = -0.00031; D_ini(7, 3) = -0.00039;

	D_ini(0, 4) = -0.00054; D_ini(1, 4) = -0.00041; D_ini(2, 4) = -0.00019; D_ini(3, 4) = -0.00017;
	D_ini(4, 4) = 0.00228; D_ini(5, 4) = -0.00024; D_ini(6, 4) = -0.00033; D_ini(7, 4) = -0.00041;

	D_ini(0, 5) = -0.00072; D_ini(1, 5) = -0.00053; D_ini(2, 5) = -0.00028; D_ini(3, 5) = -0.00025;
	D_ini(4, 5) = -0.00024; D_ini(5, 5) = 0.00304; D_ini(6, 5) = -0.00045; D_ini(7, 5) = -0.00056;

	D_ini(0, 6) = -0.00097; D_ini(1, 6) = -0.00070; D_ini(2, 6) = -0.00035; D_ini(3, 6) = -0.00031;
	D_ini(4, 6) = -0.00033; D_ini(5, 6) = -0.00045; D_ini(6, 6) = 0.00383; D_ini(7, 6) = -0.00071;

	D_ini(0, 7) = -0.00112; D_ini(1, 7) = -0.00083; D_ini(2, 7) = -0.00043; D_ini(3, 7) = -0.00039;
	D_ini(4, 7) = -0.00041; D_ini(5, 7) = -0.00056; D_ini(6, 7) = -0.00071; D_ini(7, 7) = 0.00445;

#if 0
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			K_ini(i, j) = (i == j) ? K_DIAG : K_ELEM;
			D_ini(i, j) = (i == j) ? 1 : 0;

			/* for start conditions != 0
			 * It must decrease exponentially
			 */
//			D_ini(i, j) = (i == j) ? K_DIAG : K_ELEM;
//			K_ini(i, j) = 0;
		}
	}
#endif

	std::cout << "Initian matrix (A, K, D):" << std::endl;
	(*A_ini).print();
	std::cout << "******************" << std::endl;
	K_ini.print();
	std::cout << "******************" << std::endl;
	D_ini.print();

	*A = K_ini * (-1.0 / h / h);
	*B = ((*A_ini) * (1.0 / t)) + K_ini * (2.0 / h / h) + D_ini;
	*C = K_ini * (-1.0 / h / h);
};

progonka::~progonka()
{
	delete A_ini;
	delete A;
	delete B;
	delete C;
};


void progonka::edge_conditions(const std::vector<vector>& u, std::vector<vector>& u1)
{
	u1[0] = u[0];
	u1[L - 1] = u[L - 1];
};


void progonka::calculate(const std::vector<vector>& u, std::vector<vector>& u1)
{
	assert(A_ini);
	assert(A);
	assert(B);
	assert(C);

	/* right part */
	vector Fi(n, 1);

	/* coefficients */
	std::vector<matrix> Ps(L);
	std::vector<vector> Qs(L);
	matrix G(n);

	/* due to the specific realization of std containers */
//	for (std::vector<vector>::iterator it = Ps.begin(); it != Ps.end(); ++it) {
	for (int i = 0; i < L; i++) {
		Ps[i].init(n);
		Qs[i].init(n);
	}

	/* progonka coeffs start cond */

/*
 * Something wrong, this gets left edge falling to null
 * But this is done according to prof. Lobanov book
 *
	Fi = ((*A_ini) * (1.0 / t)) * u[0];
	Ps[0] = ((*B).inverse() * (*C)) * (-1);
	Qs[0] = (*B).inverse() * Fi;

*/

	Ps[0].init(n);
	Qs[0].init(n, P_LEFT);

	/* edge conditions */
	edge_conditions(u, u1);

	/* There */
	for (int i = 1; i < L - 1; i++) {
		Fi = ((*A_ini) * (1.0 / t)) * u[i];
		G = ((*A * Ps[i - 1]) + *B).inverse();

		Ps[i] = (G * (*C)) * (-1);
		Qs[i] = G * (Fi - (*A * Qs[i - 1]));
	}

// I don't remember why this is here
//	u1[L - 2] = Qs[L - 2];


	/* and Back Again */
	for (int i = L - 2; i >= 0; i--) {
		u1[i] = (Ps[i] * u1[i + 1]) + Qs[i];
	}

	Ps.clear();
	Qs.clear();
}

