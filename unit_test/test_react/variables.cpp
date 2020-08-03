#include <variables.H>
#include <network.H>

plot_t init_variables() {

  plot_t p;

  p.irho = p.next_index(1);
  p.itemp = p.next_index(1);

  p.ispec = p.next_index(NumSpec);
  p.ispec_old = p.next_index(NumSpec);

#if NAUX_NET > 0
  p.iaux = p.next_index(NumAux);
  p.iaux_old = p.next_index(NumAux);
#endif

  p.irodot = p.next_index(NumSpec);
  p.irho_hnuc = p.next_index(1);

  return p;

}
