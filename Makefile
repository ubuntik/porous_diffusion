CC = g++
HOME = /Users/admin/Documents/Schlum/
CFLAGS = -g -llapack -lblas
#BLAS = ${HOME}/OpenBLAS/
#CFLAGS = -I${BLAS} -L${BLAS} -lopenblas

all: main

main:
	$(CC) lalgebra.cpp solvers.cpp vtk.cpp main.cpp $(CFLAGS) -o main

clean:
	rm -rf main main.dSYM res/*
