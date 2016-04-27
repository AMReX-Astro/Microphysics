FC := gfortran
FFLAGS = -c

f90sources = stiff_ode.f90 bl_types.f90
f90sources += $(notdir $(wildcard ../LAPACK/*.f90))

fsources += $(notdir $(wildcard ../LAPACK/*.f))

vpath %.f90 . ..
vpath %.f . .. ../LAPACK


ALL: test.ex

odir := _build

OBJECTS := $(addprefix $(odir)/, $(f90sources:.f90=.o))
OBJECTS += $(addprefix $(odir)/, $(fsources:.f=.o))

print-%: ; @echo $* is $($*)

$(odir)/stiff_ode.o: $(odir)/bl_types.o

$(odir)/%.o: %.f90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	$(FC) $(FFLAGS) -o $@ $<

$(odir)/%.o: %.f
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	$(FC) $(FFLAGS) -o $@ $<


test.ex: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

