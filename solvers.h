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
	progonka(uint n_size, uint length,
		const std::vector<matrix> &a, const std::vector<matrix> &b,
		const std::vector<matrix> &c, const std::vector<vector> &f);
	~progonka() {};
	void calculate(std::vector<vector>& u1);
private:
	uint n;
	uint l;
	std::vector<matrix> A;
	std::vector<matrix> B;
	std::vector<matrix> C;
	std::vector<vector> F;
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
};

static void get_left_edge(uint* u_left)
{

//	u_left[0] = 1;
//	u_left[1] = 1;
//	u_left[2] = 1;
//	u_left[3] = 1;
//	u_left[4] = 1;
//	u_left[5] = 1;
//	u_left[6] = 1;
//	u_left[7] = 1;

}

static void get_right_edge(uint* u_right)
{

//	u_right[0] = 1;
//	u_right[1] = 1;
//	u_right[2] = 1;
//	u_right[3] = 1;
//	u_right[4] = 1;
//	u_right[5] = 1;
//	u_right[6] = 1;
//	u_right[7] = 1;

}

#endif // _SLB_SOLVERS_LIB
