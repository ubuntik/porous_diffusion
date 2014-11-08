/*
 * DECLARATION
 */

#include <assert.h>
#include <iostream>
#include <math.h>

#include "lalgebra.h"

matrix::matrix()
{
	n = 0;
	size_in_bytes = 0;
	array = NULL;
};

matrix::matrix(uint n_size)
{
	n = n_size;
	assert(n >= 0);
	size_in_bytes = n * n * sizeof(ptype);
	array = (ptype *)malloc(size_in_bytes);
	bzero(array, size_in_bytes);
	assert(array != NULL);
};

matrix::matrix(uint n_size, ptype a)
{
	n = n_size;
	assert(n >= 0);
	size_in_bytes = n * n *sizeof(ptype);
	array = (ptype *)malloc(size_in_bytes);
	assert(array != NULL);

	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			array[j * n + i] = (i == j) ? a : 0.0;
};

matrix::~matrix()
{
//	if (array)
//		free(array);
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
	assert(n == b.n);
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
	assert(n == b.n);
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
	assert(n == b.n);
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
	assert(n == b.n);
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
	assert(n == b.n);
	n = b.n;
	size_in_bytes = b.size_in_bytes;
	memcpy(array, b.array, size_in_bytes);
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
	return coFactor().transpose() * idet;
};

vector::vector()
{
	n = 0;
	size_in_bytes = 0;
	array = NULL;
};

vector::vector(uint n_size)
{
	n = n_size;
	assert(n >= 0);
	size_in_bytes = n * sizeof(ptype);
	array = (ptype *)malloc(size_in_bytes);
	bzero(array, size_in_bytes);
	assert(array != NULL);
};

vector::vector(uint n_size, ptype a)
{
	n = n_size;
	assert(n >= 0);
	size_in_bytes = n * sizeof(ptype);
	array = (ptype *)malloc(size_in_bytes);
	assert(array != NULL);

	for (int i = 0; i < n; i++)
		array[i] = a;
};

vector::~vector()
{
//	if (array)
//		free(array);
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
	assert(n == b.n);
	vector c(n);
	for (int i = 0; i < n; i++) {
		c.array[i] = array[i] + b.array[i];
	}
	return c;
};

vector vector::operator-(const vector &b)
{
	assert(n == b.n);
	vector c(n);
	for (int i = 0; i < n; i++) {
		c.array[i] = array[i] - b.array[i];
	}
	return c;
};

vector& vector::operator=(const vector& b)
{
	assert(n == b.n);
	n = b.n;
	size_in_bytes = b.size_in_bytes;
	memcpy(array, b.array, size_in_bytes);
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

