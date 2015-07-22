/*
 * @file lalgebra.cpp
 *
 * @author Anna Subbotina
 *
 */

#include <assert.h>
#include <iostream>
#include <cstdio>
#include <math.h>

#include "lalgebra.h"

extern "C" {
	// LU decomoposition of a general matrix
	void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

	// generate inverse of a matrix given its LU decomposition
	void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

ptype power(ptype base, ptype exponent)
{
	ptype half_pow = 0.0, int_part = 0.0;
	ptype fr_part = modf(exponent, &int_part);

	if ((fr_part != 0) && (fr_part != 0.5)) {
		return pow(base, exponent);
	} else if (fr_part == 0.5) {
		return (sqrt(base) * power(base, int_part));
	}

	if (exponent == 0)
		return (ptype)1;
	else if (exponent < 0)
		return (ptype)1 / power(base, -exponent);
	else if (fmod(exponent, 2) == 0) {
		half_pow = power(base, exponent / (ptype)2);
		return half_pow * half_pow;
	} else
		return base * power(base, exponent - 1);
};

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

void matrix::init(uint n_size)
{
	n = n_size;
	assert(n >= 0);
	size_in_bytes = n * n * sizeof(ptype);
	if (array)
		free(array);
	array = (ptype *)malloc(size_in_bytes);
	bzero(array, size_in_bytes);
	assert(array != NULL);
};

void matrix::init(uint n_size, ptype a)
{
	n = n_size;
	assert(n >= 0);
	size_in_bytes = n * n *sizeof(ptype);
	if (array)
		free(array);
	array = (ptype *)malloc(size_in_bytes);
	assert(array != NULL);

	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			array[j * n + i] = (i == j) ? a : 0.0;
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
	matrix *b = new matrix(n);
	matrix *c = new matrix(n - 1);
	ptype det;

	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			__gen_adj_matrix(*c, *this, i, j, n);
			det = (*c).det();
			(*b)(j, i) = pow(-1.0, i + j + 2.0) * det;
		}
	}
	delete c;
	return *b;
}

matrix matrix::transpose()
{
	if (n == 1) {
		return *this;
	}
	matrix *c = new matrix(n);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			(*c)(i, j) = array[i * n + j];
		}
	}
	return *c;
};

static void ll_inverse(double* A, int N)
{
	int *IPIV = new int[N+1];
	int LWORK = N*N;
	double *WORK = new double[LWORK];
	int INFO;

	dgetrf_(&N,&N,A,&N,IPIV,&INFO);
	dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

	delete IPIV;
	delete WORK;
}

matrix matrix::inverse()
{
	if (n == 1) {
		matrix *c = new matrix(1, 1.0/array[0]);
		return *c;
	}

	matrix *c = new matrix(n);
	memcpy((*c).array, array, size_in_bytes);
	ll_inverse((*c).array, (*c).n);

	return *c;
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


void vector::init(uint n_size)
{
	n = n_size;
	assert(n >= 0);
	size_in_bytes = n * sizeof(ptype);
	if (array)
		free(array);
	array = (ptype *)malloc(size_in_bytes);
	bzero(array, size_in_bytes);
	assert(array != NULL);
};

void vector::init(uint n_size, ptype a)
{
	n = n_size;
	assert(n >= 0);
	size_in_bytes = n * sizeof(ptype);
	if (array)
		free(array);
	array = (ptype *)malloc(size_in_bytes);
	assert(array != NULL);

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

vector vector::add_abs()
{
	vector c(n);
	for (int i = 0; i < n; i++) {
		c.array[i] = (array[i] >= 0) ? 2 * array[i] : 0;
	}
	return c;
};

vector vector::dif_abs()
{
	vector c(n);
	for (int i = 0; i < n; i++) {
		c.array[i] = (array[i] >= 0) ? 0 : 2 * array[i];
	}
	return c;
};

vector vector::mult(const vector& b)
{
	assert(n == b.n);
	vector c(n);
	for (int i = 0; i < n; i++) {
		c.array[i] = array[i] * b.array[i];
	}
	return c;
};

