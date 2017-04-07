/*
 * @file vtk.h
 * @brief Functions for dumping data to files for paraview
 * @author Anna Subbotina
 *
 */

#ifndef _SLB_VTKPRINT_LIB
#define _SLB_VTKPRINT_LIB

#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

typedef unsigned int uint;

void write_plain_data(vector<double>& u, const char *path, uint l);

void write_to_vtk1(vector<vec>& u, const char *path, uint n, uint l);

void write_to_vtk2(vector<vec>& u, const char *path, uint n, uint l, double crt_time);

void write_to_vtk3(vector< vector<vec> >& u, const char *path, uint n, uint l, uint t);

void write_to_vtk2d(vector<vec>& u, const char *path, const char *data,
			const int N[2], const double hh[2]);

#endif // _SLB_VTKPRINT_LIB
