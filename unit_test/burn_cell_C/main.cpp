#include <iostream>
#include <cstring>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
using namespace amrex;

#include "burn_cell_F.H"

int main(int argc, char *argv[]) {

  std::cout << "starting the single zone burn..." << std::endl;

  ParmParse ppa("amr");

  std::string probin_file = "probin";

  ppa.query("probin_file", probin_file);

  const int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  burn_cell(probin_file_name.dataPtr(), &probin_file_length);

}
