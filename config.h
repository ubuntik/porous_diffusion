/*
 * @file config.h
 *
 * @brief file with problem main parameters
 *
 * @author Anna Subbotina
 *
 */

#ifndef _SLB_CONFIG_LIB
#define _SLB_CONFIG_LIB

// Degree of hierarchy of porous media
#define P 2
// Number of porosity levels
#define GAMMA 4
// Length of considerable sample (1 meter)
#define LENGTH 250
// Step by space
#define h 1.0
// Number of steps by space
//#define L (LENGTH / h)
#define L 8
// Considerable time
#define TIME 8000
// Step by time
#define t 0.1
// coefficient for tracer calculation (0 <= THETA <= 1)
#define THETA 0.5
// curant
#define CURANT 0.0001
// v = - (K / eta)(dP/dx)
#define ETA 1

// edge conditions
#define P_LEFT 1
#define P_RIGHT 0

#define K_DIAG 1
#define K_ELEM 0

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// types
typedef unsigned int uint;
typedef double ptype;

#endif // _SLB_CONFIG_LIB
