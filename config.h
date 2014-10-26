/*
 * DECLARATION
 */

// Degree of hierarchy of porous media
#define P 3
// Number of porosity levels
#define GAMMA 4
// Length of considerable sample (1 meter)
#define LENGTH 10
// Step by space
#define h 1
// Number of steps by space
#define L (LENGTH / h)
// Considerable time
#define TIME 1000
// Step by time
#define t 1

// types
typedef unsigned int uint;
typedef double ptype;

typedef double p_vector[GAMMA];
typedef double matrix[GAMMA][GAMMA];

void print_vec(p_vector *u)
{
	for (int j = 0; j < GAMMA; j++) {
		for (int i = 0; i < L; i++) {
			std::cout << u[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void start_cond(p_vector *u)
{
	// left edge
	for (int j = 0; j < GAMMA; j++)
		u[0][j] = 0;

	for (int i = 1; i < (L - 1); i++) {
		for (int j = 0; j < GAMMA; j++)
			u[i][j] = 0;
	}

	// right edge
	for (int j = 0; j < GAMMA; j++)
		u[L - 1][j] = 1;
}

