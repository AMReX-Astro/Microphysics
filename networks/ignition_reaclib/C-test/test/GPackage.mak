# we need a dummy GPackage.mak here just to have this directory be
# part of the search path.  But since all the files here are a 'main',
# they don't get added to the general objects variable -- they are
# brought in one-by-one for their respective executable's linking

f90sources += probin.f90
f90sources += constants_cgs.f90
