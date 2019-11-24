#ifndef BURN_CELL_F_H_
#define BURN_CELL_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

void burn_cell(const int* name, const int* namlen);

#ifdef __cplusplus
}
#endif

#endif
