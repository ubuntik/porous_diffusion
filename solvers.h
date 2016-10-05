/*
 * @file solvers.h
 *
 * @brief header of tompson's method of numerical solution
 *
 * @author Anna Subbotina
 *
 */

#ifndef _SLB_SOLVERS_LIB
#define _SLB_SOLVERS_LIB

#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

typedef unsigned int uint;

class progonka
{
public:
	progonka(uint n_size, uint length);
	~progonka() {};
	void calculate(vector<vec>& u1,
		vector<mat> &a, vector<mat> &b,
		vector<mat> &c, vector<vec> &f);
private:
	uint n;
	uint l;
};

class secondord
{
public:
	secondord(uint n_size, uint length, vector<mat> *Ki, vector<mat> *Di);
	~secondord();
	void calculate( const vector<vec>& v,
			const vector<vec>& p,
			const vector<vec>& C, vector<vec>& C1,
			double dt, double h);
private:
	void left_edge(vector<vec>& C1);
	void right_edge(vector<vec>& C1);
	vector<mat> *K;
	vector<mat> *D;
	uint n;
	uint l;
};
#endif // _SLB_SOLVERS_LIB
