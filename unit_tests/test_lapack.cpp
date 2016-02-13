#include <iostream>
#include <vector>
#include <stdlib.h>
#include "../lalgebra.h"

using namespace std;

int main()
{
	int n = 4;
        double *x = (double *)calloc(sizeof(double), n);
        double *a = (double *)calloc(sizeof(double), n * n);
        double *b = (double *)calloc(sizeof(double), n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			a[i * n + j] = random() % 9;
		b[i] = random() % 9;
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			std::cout << a[i * n + j] << ", ";
		std::cout << std::endl;
	}
	std::cout << std::endl;

	for (int i = 0; i < n; i++)
		std::cout << b[i] << ", ";
	std::cout << std::endl;

	solve_eq(a, b, x, n);

	std::cout << "solution is: [";
	for (int i = 0; i < n; i++)
		std::cout << x[i] << ", ";
	std::cout << "]" << std::endl;

	return 0;

}
