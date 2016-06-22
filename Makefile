CC = g++
HOME = /Users/admin/Documents/Schlum/
CFLAGS = -g -larmadillo -llapack -lblas

all: tracer

simple_ODE:
	$(CC) solvers.cpp vtk.cpp simple_ODE.cpp $(CFLAGS) -o simple_ODE

compare:
	$(CC) solvers.cpp vtk.cpp compare.cpp $(CFLAGS) -o compare

direct:
	$(CC) solvers.cpp vtk.cpp direct.cpp $(CFLAGS) -o direct

tracer:
	$(CC) solvers.cpp vtk.cpp tracer.cpp $(CFLAGS) -o tracer

meso:
	$(CC) solvers.cpp vtk.cpp meso.cpp $(CFLAGS) -o meso

clean:
	rm -rf simple_ODE compare check direct tracer meso *.dSYM res; mkdir res
