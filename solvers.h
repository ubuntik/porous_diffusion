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

#include "lalgebra.h"

class progonka
{
public:
	progonka(uint n_size, uint length);
	~progonka() {};
	void calculate(std::vector<vector>& u1,
		std::vector<matrix> &a, std::vector<matrix> &b,
		std::vector<matrix> &c, std::vector<vector> &f);
private:
	uint n;
	uint l;
};

class corner
{
public:
	corner(uint n_size, uint length);
	~corner();
	void calculate( const std::vector<vector>& v_med,
			const std::vector<vector>& p,
			const std::vector<vector>& C,
			std::vector<vector>& C1, double dt);
private:
	void left_edge(std::vector<vector>& C1);
	void right_edge(std::vector<vector>& C1);
	uint n;
	uint l;
	matrix *K;
	matrix *D;
	vector *perm;
};

class secondord
{
public:
	secondord(uint n_size, uint length);
	~secondord();
	void calculate( const std::vector<vector>& v_med,
			const std::vector<vector>& p,
			const std::vector<vector>& C,
			std::vector<vector>& C1, double dt);
private:
	void left_edge(std::vector<vector>& C1);
	void right_edge(std::vector<vector>& C1);
	uint n;
	uint l;
	matrix *K;
	matrix *D;
	vector *perm;
};

#endif // _SLB_SOLVERS_LIB
