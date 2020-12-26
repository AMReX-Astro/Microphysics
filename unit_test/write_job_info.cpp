#include <iostream>
#include <unistd.h>
#include <iomanip>
#include <string>
#include <cstring>
#include <AMReX_Utility.H>
#include <AMReX_buildInfo.H>

#include <extern_parameters.H>
#include <network.H>
#include <unit_test.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace amrex;

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

  jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
  jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

  jobInfoFile << "\n\n";

  // plotfile information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Plotfile Information\n";
  jobInfoFile << PrettyLine;

  time_t now = time(0);

  // Convert now to tm struct for local timezone
  tm* localtm = localtime(&now);
  jobInfoFile   << "output date / time: " << asctime(localtm);

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

  jobInfoFile << "GPU model name: " << Gpu::Device::deviceName() << "\n";
  jobInfoFile << "Number of GPUs used: " << Gpu::Device::numDevicesUsed() << "\n";

  jobInfoFile << "\n\n";
#endif

  // build information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Build Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
  jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
  jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
  jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
  jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

  jobInfoFile << "\n";
  
  jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
  jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
  jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
  jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

  jobInfoFile << "\n";

  for (int n = 1; n <= buildInfoGetNumModules(); n++) {
    jobInfoFile << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
  }

  jobInfoFile << "\n";

  const char* githash1 = buildInfoGetGitHash(1);
  const char* githash2 = buildInfoGetGitHash(2);
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

  jobInfoFile.close();

  // now the external parameters
  const int jobinfo_file_length = FullPathJobInfoFile.length();
  Vector<int> jobinfo_file_name(jobinfo_file_length);

  for (int i = 0; i < jobinfo_file_length; i++) {
    jobinfo_file_name[i] = FullPathJobInfoFile[i];
  }

  runtime_pretty_print(jobinfo_file_name.dataPtr(), &jobinfo_file_length);

}
