CC = g++
HOME = /Users/admin/Documents/Schlum/
CFLAGS = -g -llapack

all: main

main:
	$(CC) lalgebra.cpp solvers.cpp vtk.cpp main.cpp $(CFLAGS) -o main

clean:
	rm -rf main main.dSYM res; mkdir res
