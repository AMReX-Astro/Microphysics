#ifndef _actual_eos_H_
#define _actual_eos_H_
#include <string>
#include <eos.H>

const std::string eos_name = "gamma_law_general";

void actual_eos_cxx_init();

void actual_eos_cxx(const eos_input_t input, eos_t& state);

void actual_eos_cxx_finalize();

#endif
