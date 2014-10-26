/*
 * DECLARATION
 */

#include <iostream>

#include <cblas.h>

#include "lalgebra.h"

matrix::matrix(uint n_size)
{
	n = n_size;
	size = n * n * sizeof(ptype);
	array = (ptype *)malloc(size);
	bzero(array, size);
	// TODO: check memory
};

matrix::matrix(uint n_size, ptype a)
{
	n = n_size;
	size = n * n *sizeof(ptype);
	array = (ptype *)malloc(size);
	// TODO: check memory

	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			array[j * n + i] = (i == j) ? a : 0;
};

matrix::~matrix()
{
	if (array)
		free(array);
};

void matrix::print()
{
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			std::cout << array[j * n + i] << " ";
		}
		std::cout << std::endl;
	}
};

matrix matrix::operator*(const matrix &b)
{
	matrix c(n);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < n; k++) {
				c.array[j * n + i] += array[j * n + k] * b(i, k);
			}
		}
	}
	return c;
};

vector matrix::operator*(const vector &b)
{
	vector c(n);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			c.array[j] += array[j * n + i] * b(i);
		}
	}
	return c;
};

matrix matrix::operator*(const ptype b)
{
	matrix c(n);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			c.array[j * n + i] = b * array[j * n + i];
		}
	}
	return c;
};

matrix matrix::operator+(const matrix &b)
{
	matrix c(n);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			c.array[j * n + i] = array[j * n + i] + b.array[j * n + i];
		}
	}
	return c;
};

matrix matrix::operator-(const matrix &b)
{
	matrix c(n);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			c.array[j * n + i] = array[j * n + i] - b.array[j * n + i];
		}
	}
	return c;
};

matrix& matrix::operator=(const matrix& b)
{
	n = b.n;
	size = b.size;
	memcpy(array, b.array, size);
	return *this;
};

void __gen_add_matrix(matrix& m, const matrix& a,
			const int j, const uint n)
{
	for (int i = 1; i < n; i++) {
		int j2 = 0;
		for (int j1 = 0; j1 < n; j1++) {
			if (j == j1)
				continue;
			m(i - 1, j2) = a(i, j1);
			j2++;
		}
	}
}

ptype matrix::det()
{
	ptype det = 0;
	matrix m(n - 1);

	if (n == 1)
		return array[0];
	if (n == 2)
		return array[0] * array[3] - array[1] * array[2];

	for (int j = 0; j < n; j++) {
		__gen_add_matrix(m, *this, j, n);
		det += pow(-1.0, j + 2.0) * array[j * n] * m.det();
	}
	return(det);
};

void __gen_adj_matrix(matrix& c, const matrix& a,
		const int i, const int j, const uint n)
{
	int i1 = 0;
	for (int ii = 0; ii < n; ii++) {
		if (ii == i)
			continue;
		int j1 = 0;
		for (int jj = 0; jj < n; jj++) {
			if (jj == j)
				continue;
			c(j1, i1) = a(jj, ii);
			j1++;
		}
		i1++;
	}
};

matrix matrix::coFactor()
{
	matrix b(n);
	matrix c(n - 1);
	ptype det;

	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			__gen_adj_matrix(c, *this, i, j, n);
			det = c.det();
			b(j, i) = pow(-1.0, i + j + 2.0) * det;
		}
	}
	return b;
}

matrix matrix::transpose()
{
	matrix c(n);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			c(i, j) = array[i * n + j];
		}
	}
	return c;
};

matrix matrix::inverse()
{
	ptype idet = 1.0 / det();
	std::cout << "det = " << det() << "; 1/det = " << idet << std::endl;
	return coFactor().transpose() * idet;
};

vector::vector(uint n_size)
{
	n = n_size;
	size = n * sizeof(ptype);
	array = (ptype *)malloc(size);
	bzero(array, size);
	// TODO: check memory
};

vector::vector(uint n_size, ptype a)
{
	n = n_size;
	size = n *sizeof(ptype);
	array = (ptype *)malloc(size);
	// TODO: check memory

	for (int i = 0; i < n; i++)
		array[i] = a;
};

vector::~vector()
{
	if (array)
		free(array);
};

void vector::print()
{
	std::cout << "( ";
	for (int i = 0; i < n; i++)
		std::cout << array[i] << " ";
	std::cout << ")" << std::endl;
};

vector vector::operator+(const vector &b)
{
	vector c(n);
	for (int i = 0; i < n; i++) {
		c.array[i] = array[i] + b.array[i];
	}
	return c;
};

vector vector::operator-(const vector &b)
{
	vector c(n);
	for (int i = 0; i < n; i++) {
		c.array[i] = array[i] - b.array[i];
	}
	return c;
};

vector& vector::operator=(const vector& b)
{
	n = b.n;
	size = b.size;
	memcpy(array, b.array, size);
	return *this;
};

vector vector::operator*(const ptype b)
{
	vector c(n);
	for (int i = 0; i < n; i++) {
		c.array[i] = array[i] * b;
	}
	return c;
};

#if 0
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
#endif
