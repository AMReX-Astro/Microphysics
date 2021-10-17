test: eos17.f90 eos_c.cpp
	g++ -o eos_c.o -c eos_c.cpp
	gfortran -o test eos_c.o eos17.f90 -lstdc++

run:
	./test

