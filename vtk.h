/*
 * @file vtk.h
 *
 * @author Anna Subbotina
 *
 */

#ifndef _SLB_VTKPRINT_LIB
#define _SLB_VTKPRINT_LIB

#include <vector>
#include "lalgebra.h"

#define START 1

void write_to_vtk1(std::vector<vector>& u, const char *path, uint n);

void write_to_vtk2(std::vector<vector>& u, const char *path, uint n, double crt_time);

void write_to_vtk3(std::vector< std::vector<vector> >& u, const char *path, uint n, uint l);

void write_to_vtk2d(std::vector<vector>& u, const char *path, const char *data,
			const int N[2], const double o[2], const double hh[2]);

#endif // _SLB_VTKPRINT_LIB
