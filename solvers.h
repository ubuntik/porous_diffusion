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
	void edge_conditions(const std::vector<vector>& u, std::vector<vector>& u1);
	uint n;
	matrix *A_ini;
	matrix *A;
	matrix *B;
	matrix *C;
};

#endif // _SLB_SOLVERS_LIB
