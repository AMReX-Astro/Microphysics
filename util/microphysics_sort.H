#ifndef MICROPHYSICS_SORT_H
#define MICROPHYSICS_SORT_H
template <int l, int m>
AMREX_INLINE
void sort_Array1D(amrex::Array1D<int, l, m>& array_1d){

  // a simple bubble sort, for sorting Array1D
  // assumes index start from int l, and end at int m

  int scratch;

  for (int i = l; i <= m; ++i){
    for (int j = l; j <= m-i; ++j){
      if (array_1d(j) > array_1d(j+1)){
	scratch = array_1d(j);
	array_1d(j) = array_1d(j+1);
	array_1d(j+1) = scratch;
      }
    }
  }
}

#endif
