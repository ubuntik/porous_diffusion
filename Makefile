CC = g++
CFLAGS = -g -larmadillo -llapack -lblas

all: tracer direct compare simple_ODE

simple_ODE:
	$(CC) solvers.cpp vtk.cpp simple_ODE.cpp $(CFLAGS) -o simple_ODE

compare:
	$(CC) solvers.cpp vtk.cpp compare.cpp $(CFLAGS) -o compare

direct:
	$(CC) solvers.cpp vtk.cpp direct.cpp $(CFLAGS) -o direct

tracer:
	$(CC) solvers.cpp vtk.cpp tracer.cpp $(CFLAGS) -o tracer

clean:
	rm -rf simple_ODE compare check direct tracer meso *.dSYM res; mkdir res
