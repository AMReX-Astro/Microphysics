#include <set>
#include <tuple>

#include <variables.H>
#include <network.H>
#include <rhs_type.H>

struct reactants_t {
  int species_A;
  int species_B;
  int exponent_A;
  int exponent_B;

  bool operator<(const reactants_t &rhs) const {
    return std::tie(species_A, species_B, exponent_A, exponent_B) <
      std::tie(rhs.species_A, rhs.species_B, rhs.exponent_A, rhs.exponent_B);
  }
};

plot_t init_variables(amrex::Vector<std::string>& names) {

  plot_t p;

  p.irho = p.next_index(1);
  names.push_back("density");
  p.itemp = p.next_index(1);
  names.push_back("temperature");
  p.ispec = p.next_index(NumSpec);

  for (int n = 0; n < NumSpec; n++) {
    names.push_back("X_" + spec_names_cxx[n]);
  }

  // keep track of what reactants we've added to avoid duplication
  std::set<reactants_t> seen;

  for (int rate = 1; rate <= Rates::NumRates; ++rate) {
    RHS::rhs_t data = RHS::rhs_data(rate);
    if (data.screen_forward_reaction == 0 && data.screen_reverse_reaction == 0) {
      continue;
    }
    reactants_t key{data.species_A, data.species_B, data.exponent_A, data.exponent_B};
    const auto it = seen.find(key);
    if (it != seen.end()) {
      continue;
    }
    if (data.exponent_A == 1 && data.exponent_B == 1 && data.exponent_C == 0) {
      p.iscn(rate).value = p.next_index(1);
      names.push_back("scn_" + short_spec_names_cxx[data.species_A-1] +
                      "_" + short_spec_names_cxx[data.species_B-1]);
      seen.insert(key);
    }
    if (data.exponent_A == 2 && data.exponent_B == 0 && data.exponent_C == 0) {
      p.iscn(rate).value = p.next_index(1);
      names.push_back("scn_" + short_spec_names_cxx[data.species_A-1] +
                      "_" + short_spec_names_cxx[data.species_A-1]);
      seen.insert(key);
    }
    if (data.exponent_A == 3 && data.exponent_B == 0 && data.exponent_C == 0) {
      p.iscn(rate).value = p.next_index(1);
      names.push_back("scn_triple-" + short_spec_names_cxx[data.species_A-1]);
      p.iscn(rate).aux_value = p.next_index(1);
      names.push_back("scn_triple-" + short_spec_names_cxx[data.species_A-1] + "_aux");
      seen.insert(key);
    }
  }

  for (int rate = 1; rate <= Rates::NumRates; ++rate) {
    if (p.iscn(rate).value != -1) {
      p.iscn(rate).dt = p.next_index(1);
      names.push_back(names[p.iscn(rate).value] + "_dt");
    }
    if (p.iscn(rate).aux_value != -1) {
      p.iscn(rate).aux_dt = p.next_index(1);
      names.push_back(names[p.iscn(rate).aux_value] + "_dt");
    }
  }

  return p;
}
