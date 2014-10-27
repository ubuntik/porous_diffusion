/*
 * DECLARATION
 */
#ifndef _SLB_CONFIG_LIB
#define _SLB_CONFIG_LIB

// Degree of hierarchy of porous media
#define P 3
// Number of porosity levels
#define GAMMA 4
// Length of considerable sample (1 meter)
#define LENGTH 10
// Step by space
#define h 1
// Number of steps by space
#define L (LENGTH / h)
// Considerable time
#define TIME 1000
// Step by time
#define t 1
// matrix size N x N
#define N 4
// edge conditions
#define P_LEFT 2
#define P_RIGHT 1


// types
typedef unsigned int uint;
typedef double ptype;

#endif // _SLB_CONFIG_LIB
