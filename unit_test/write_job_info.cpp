#include <iostream>
#include <unistd.h>
#include <iomanip>
#include <string>
#include <cstring>
#include <AMReX_Utility.H>
#include <AMReX_buildInfo.H>
#include <AMReX_GpuDevice.H>
#include <extern_parameters.H>
#include <network.H>
#include <unit_test.H>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

void write_job_info(const std::string& dir) {

  std::ofstream jobInfoFile;
  std::string FullPathJobInfoFile = dir;
  FullPathJobInfoFile += "/job_info";
  jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

  std::string PrettyLine = std::string(78, '=') + "\n";
  std::string OtherLine = std::string(78, '-') + "\n";
  std::string SkipSpace = std::string(8, ' ');

  // job information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Microphysics Job Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "number of MPI processes: " << amrex::ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
  jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

  jobInfoFile << "\n\n";

  // plotfile information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Plotfile Information\n";
  jobInfoFile << PrettyLine;

  const std::time_t now = time(nullptr);
  char buf[64];
  if (strftime(buf, sizeof buf, "%c\n", std::localtime(&now))) {
      jobInfoFile << "output date / time: " << buf << "\n";
  }

  char currentDir[FILENAME_MAX];
  if (getcwd(currentDir, FILENAME_MAX)) {
    jobInfoFile << "output dir:         " << currentDir << "\n";
  }

  jobInfoFile << "\n\n";

#ifdef AMREX_USE_GPU
  // This output assumes for simplicity that every rank uses the
  // same type of GPU.

  jobInfoFile << PrettyLine;
  jobInfoFile << "GPU Information:       " << "\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "GPU model name: " << amrex::Gpu::Device::deviceName() << "\n";
  jobInfoFile << "Number of GPUs used: " << amrex::Gpu::Device::numDevicesUsed() << "\n";

  jobInfoFile << "\n\n";
#endif

  // build information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Build Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "build date:    " << amrex::buildInfoGetBuildDate() << "\n";
  jobInfoFile << "build machine: " << amrex::buildInfoGetBuildMachine() << "\n";
  jobInfoFile << "build dir:     " << amrex::buildInfoGetBuildDir() << "\n";
  jobInfoFile << "AMReX dir:     " << amrex::buildInfoGetAMReXDir() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "COMP:          " << amrex::buildInfoGetComp() << "\n";
  jobInfoFile << "COMP version:  " << amrex::buildInfoGetCompVersion() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "C++ compiler:  " << amrex::buildInfoGetCXXName() << "\n";
  jobInfoFile << "C++ flags:     " << amrex::buildInfoGetCXXFlags() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "Fortran comp:  " << amrex::buildInfoGetFName() << "\n";
  jobInfoFile << "Fortran flags: " << amrex::buildInfoGetFFlags() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "Link flags:    " << amrex::buildInfoGetLinkFlags() << "\n";
  jobInfoFile << "Libraries:     " << amrex::buildInfoGetLibraries() << "\n";

  jobInfoFile << "\n";

  for (int n = 1; n <= amrex::buildInfoGetNumModules(); n++) {
    jobInfoFile << amrex::buildInfoGetModuleName(n) << ": " << amrex::buildInfoGetModuleVal(n) << "\n";
  }

  jobInfoFile << "\n";

  const char* githash1 = amrex::buildInfoGetGitHash(1);
  const char* githash2 = amrex::buildInfoGetGitHash(2);
  if (strlen(githash1) > 0) {
    jobInfoFile << "Microphysics git describe: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    jobInfoFile << "AMReX        git describe: " << githash2 << "\n";
  }

  jobInfoFile << "\n\n";


  // species info
  int mlen = 20;

  jobInfoFile << PrettyLine;
  jobInfoFile << " Species Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile <<
    std::setw(6) << "index" << SkipSpace <<
    std::setw(mlen+1) << "name" << SkipSpace <<
    std::setw(7) << "A" << SkipSpace <<
    std::setw(7) << "Z" << "\n";
  jobInfoFile << OtherLine;

  for (int i = 0; i < NumSpec; i++)
    {
      jobInfoFile <<
        std::setw(6) << i << SkipSpace <<
        std::setw(mlen+1) << std::setfill(' ') << short_spec_names_cxx[i] << SkipSpace <<
        std::setw(7) << aion[i] << SkipSpace <<
        std::setw(7) << zion[i] << "\n";
    }
  jobInfoFile << "\n\n";


  // runtime parameters
  jobInfoFile << PrettyLine;
  jobInfoFile << " Inputs File Parameters\n";
  jobInfoFile << PrettyLine;

#include <extern_job_info_tests.H>

  jobInfoFile.close();


}
