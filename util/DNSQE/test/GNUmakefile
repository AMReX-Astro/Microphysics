include ../GPackage.mak

ALL: test

vpath %.f90 . ..
vpath %.f . ..

f90sources += test.f90

objects := $(fsources:.f=.o)
objects += $(f90sources:.f90=.o)

%.o: %.f
	gfortran -c $<

%.o: %.f90
	gfortran -c $<

test: $(objects)
	gfortran -o test $(objects)



