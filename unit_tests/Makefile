all: progonka lapack lalgebra
flags = -larmadillo -llapack -lblas -g

progonka:
	g++ -o progonka ../solvers.cpp ../lalgebra.cpp progonka.cpp ${flags}

lapack:
	g++ -o lapack ../solvers.cpp ../lalgebra.cpp test_lapack.cpp ${flags}

lalgebra:
	g++ -o lalgebra ../lalgebra.cpp lalgebra_tests.cpp ${flags}

clean:
	rm -rf lalgebra lapack progonka *.dSYM
