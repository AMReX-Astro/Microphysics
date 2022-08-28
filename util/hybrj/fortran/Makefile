fsources = $(wildcard *.f)
f90sources = $(wildcard *.f90)

objects = $(fsources:.f=.o) $(f90sources:.f90=.o)

%.o : %.f
	gfortran -c $<

%.o : %.f90
	gfortran -c $<

test: $(objects)
	gfortran -o test $(objects)
