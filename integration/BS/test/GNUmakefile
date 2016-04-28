FC := gfortran
FFLAGS = -c

LA_DIRS := $(CASTRO_DIR)/Util/VODE $(CASTRO_DIR)/Util/BLAS

f90sources = stiff_ode.f90 bl_types.f90
f90sources += $(notdir $(wildcard ../LAPACK/*.f90))

fsources += $(notdir $(foreach dir, $(LA_DIRS), $(wildcard $(dir)/*.f)))

vpath %.f90 . ..
vpath %.f . .. $(LA_DIRS)


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

