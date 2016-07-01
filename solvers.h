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
	secondord(uint n_size, uint length);
	~secondord();
	void calculate( const vector<vec>& v_med,
			const vector<vec>& p,
			const vector<vec>& C,
			vector<vec>& C1, double dt);
private:
	void left_edge(vector<vec>& C1);
	void right_edge(vector<vec>& C1);
	uint n;
	uint l;
	mat *K;
	mat *D;
	vec *perm;
};
#endif // _SLB_SOLVERS_LIB
