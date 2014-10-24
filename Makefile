CC = gcc
CP = g++
HOME = /Users/admin/Documents/Schlum/
BLAS = ${HOME}/OpenBLAS/
CFLAGS = -I${BLAS} -L${BLAS} -lopenblas

all: main

main:
	$(CP) $(CFLAGS) main.cpp -o main

clean:
	rm -rf main
