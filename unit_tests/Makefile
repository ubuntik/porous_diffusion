all: progonka lapack lalgebra

progonka:
	g++ -o progonka ../solvers.cpp ../lalgebra.cpp check.cpp -llapack

lapack:
	g++ -o lapack ../solvers.cpp ../lalgebra.cpp test_lapack.cpp -llapack

lalgebra:
	g++ -o lalgebra ../lalgebra.cpp lalgebra_tests.cpp -llapack

clean:
	rm -rf lalgebra lapack progonka
