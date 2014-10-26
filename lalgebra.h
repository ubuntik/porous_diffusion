/*
 * DECLARATION
 */

#include <iostream>

// types
typedef unsigned int uint;
typedef double ptype;

class matrix;

class vector
{
friend class matrix;
public:
	vector(uint n_size);
	vector(uint n_size, ptype a);
	~vector();
	ptype& operator()(int i) const { return array[i]; };
	vector operator+(const vector &b);
	vector operator-(const vector &b);
	vector operator*(const ptype b);
	vector& operator=(const vector& b);
	void print();
private:
	ptype *array;
	uint n;
	uint size;
};

class matrix
{
friend class vector;
public:
	matrix(uint n_size);
	matrix(uint n_size, ptype a);
	~matrix();
	ptype& operator()(int i, int j) const { return array[j * n + i]; };
	matrix operator*(const matrix &b);
	vector operator*(const vector &b);
	matrix operator*(ptype b);
	matrix operator+(const matrix &b);
	matrix operator-(const matrix &b);
	matrix& operator=(const matrix& b);
	ptype det();
	matrix coFactor();
	matrix transpose();
	matrix inverse();
	void print();

private:
	ptype *array;
	uint n;
	uint size;
};
