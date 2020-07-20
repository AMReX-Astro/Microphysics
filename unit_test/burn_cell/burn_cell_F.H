#ifndef BURN_CELL_F_H_
#define BURN_CELL_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

void init_unit_test(const int* name, const int* namlen);

void burn_cell();

#ifdef __cplusplus
}
#endif

#endif
