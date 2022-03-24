#ifndef _network_H_
#define _network_H_

#include <string>
#include <network_properties.H>
#include <react_data.H>

void network_init();

AMREX_INLINE int network_spec_index(const std::string name) {
  int idx = -1;
  for (int n = 0; n < NumSpec; n++) {
    if (name == spec_names_cxx[n]) {
      idx = n;
      break;
    }
  }
  return idx;
}

#if NAUX_NET > 0
AMREX_INLINE int network_aux_index(const std::string name) {
  int idx = -1;
  for (int n = 0; n < NumAux; n++) {
    if (name == aux_names_cxx[n]) {
      idx = n;
      break;
    }
  }
  return idx;
}
#endif

#endif
