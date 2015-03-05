CC = g++
HOME = /Users/admin/Documents/Schlum/
CFLAGS = -g -llapack
#BLAS = ${HOME}/OpenBLAS/
#CFLAGS = -I${BLAS} -L${BLAS} -lopenblas

all: main

main:
	$(CC) $(CFLAGS) lalgebra.cpp solvers.cpp main.cpp -o main

clean:
	rm -rf main main.dSYM res/*
