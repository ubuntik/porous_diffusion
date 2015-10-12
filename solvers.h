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
#include "config.h"

class progonka
{
public:
	progonka(uint n_size);
	~progonka();
	void calculate(const std::vector<vector>& u, std::vector<vector>& u1);
private:
	void edge_conditions(const std::vector<vector>& u, std::vector<vector>& u1,
				std::vector<matrix>& Ps, std::vector<vector>& Qs);
	uint n;
	matrix *A_ini;
	matrix *A;
	matrix *B;
	matrix *C;
};

class corner
{
public:
	corner(uint n_size);
	~corner();
	void calculate( std::vector<vector>& v_med,
			const std::vector<vector>& C,
			std::vector<vector>& C1, double dt);
private:
	void edge_conditions(const std::vector<vector>& C, std::vector<vector>& C1);
	uint n;
	matrix *K;
};


#endif // _SLB_SOLVERS_LIB
