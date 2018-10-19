F90sources += actual_eos.F90
F90sources += weaklib_type.F90

F90sources += PhysicalConstantsModule.f90
F90sources += UnitsModule.f90

# In weaklib/Distributions/Library
F90sources += wlKindModule.f90
F90sources += wlGridModule.f90
F90sources += wlThermoStateModule.f90
F90sources += wlDependentVariablesModule.f90
F90sources += wlInterpolationModule.f90
F90sources += wlIOModuleHDF.f90

# In weaklib/Distributions/EOSSource
F90sources += wlEquationOfStateTableModule.f90
F90sources += wlEOSIOModuleHDF.f90

# In weaklib/Distributions/OpacitySource
F90sources += wlOpacityFieldsModule.f90
F90sources += wlOpacityTableModule.f90
F90sources += wlOpacityTableIOModuleHDF.f90
