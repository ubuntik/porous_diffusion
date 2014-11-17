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
#define P 3
// Number of porosity levels
#define GAMMA 2
// Length of considerable sample (1 meter)
#define LENGTH 100
// Step by space
#define h 1
// Number of steps by space
#define L (LENGTH / h)
// Considerable time
#define TIME 2000
// Step by time
#define t 0.1
// edge conditions
#define P_LEFT 0
#define P_RIGHT 1

#define K_DIAG 1
#define K_ELEM 0.1

// types
typedef unsigned int uint;
typedef double ptype;

#endif // _SLB_CONFIG_LIB
