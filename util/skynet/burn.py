#!/usr/bin/env python

from SkyNet import *
import numpy as np
import shutil
import multiprocessing

num_done = multiprocessing.Value('i', 0)
output_dir = 'isochoric_tests/'

def run_skynet(args):
    global num_done
    # unpack arguments
    rho, total_size = args

    # set up network and options
    nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot
      + "/data/webnucleo_nuc_v2.0.xml")

    opts = NetworkOptions()
    opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
    opts.MassDeviationThreshold = 1.0E-10
    opts.IsSelfHeating = True
    opts.EnableScreening = True
    opts.DisableStdoutOutput = False

    screen = SkyNetScreening(nuclib)
    helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")

    strongReactionLibrary = REACLIBReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType.Strong, True, LeptonMode.TreatAllAsDecayExceptLabelEC,
	  "Strong reactions", nuclib, opts, True)
    symmetricFission = REACLIBReactionLibrary(SkyNetRoot
      + "/data/netsu_panov_symmetric_0neut", ReactionType.Strong, False,
      LeptonMode.TreatAllAsDecayExceptLabelEC,
      "Symmetric neutron induced fission with 0 neutrons emitted", nuclib, opts,
      False)
    spontaneousFission = REACLIBReactionLibrary(SkyNetRoot +
      "/data/netsu_sfis_Roberts2010rates", ReactionType.Strong, False,
      LeptonMode.TreatAllAsDecayExceptLabelEC, "Spontaneous fission", nuclib, opts,
      False)

    # use only REACLIB weak rates
    weakReactionLibrary = REACLIBReactionLibrary(SkyNetRoot + "/data/reaclib",
	    ReactionType.Weak, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
	    "Weak reactions", nuclib, opts, True)
    reactionLibraries = [strongReactionLibrary, symmetricFission,
            spontaneousFission, weakReactionLibrary]

    net = ReactionNetwork(nuclib, reactionLibraries, helm, screen, opts)
    nuclib = net.GetNuclideLibrary()

    # set up initial conditions
    Y0 = np.zeros(nuclib.NumNuclides())

    #URCA
    Y0[nuclib.NuclideIdsVsNames()["c12"]]  = 0.4998 / 12.0
    Y0[nuclib.NuclideIdsVsNames()["o16"]]  = 0.4998 / 16.0
    #Y0[nuclib.NuclideIdsVsNames()["ne23"]] = 4e-4 / 23.0 # above threshold
    Y0[nuclib.NuclideIdsVsNames()["na23"]] = 4e-4 / 23.0  # below threshold

    #XRB 
    #Y0[nuclib.NuclideIdsVsNames()["he4"]]  = 1.0 / 4.0

    #SubCh
    #Y0[nuclib.NuclideIdsVsNames()["c12"]]  = 0.08 / 12.0
    #Y0[nuclib.NuclideIdsVsNames()["o16"]]  = 0.08 / 16.0
    #Y0[nuclib.NuclideIdsVsNames()["n14"]] =  0.01 / 14.0
    #Y0[nuclib.NuclideIdsVsNames()["he4"]] =  0.83 / 4.0

    # T0 = 0.55 URCA, 0.2 XRB, 0.3 SubCh
    T0 = 0.55    # initial temperature in GK
    Ye = 0.5     # initial Ye
    #s = 10.0    # initial entropy in k_B / baryon
    #tau = 7.1   # expansion timescale in ms

    time = np.linspace(0.0, 1e9) 
    density_vs_time = ConstantFunction(rho) #isochoric

    # run through network, will save as "prefix".h5 and .log files
    output = net.EvolveSelfHeatingWithInitialTemperature(Y0, 0.0, 1.0E1,
		T0, density_vs_time, "urca_isoc_belowthresh")

    # save some additional information, not necessary as all will be
    # contained in h5 file
    YvsA = np.array(output.FinalYVsA())
    A = np.arange(len(YvsA))

    np.savetxt("urca_isoc_belowthresh.txt", np.array([A, YvsA]).transpose(),
	"%6i  %30.20E")

    return


if __name__ == '__main__':

  num_cores = multiprocessing.cpu_count()
  print "Running with %i worker threads" % num_cores
  pool = multiprocessing.Pool(num_cores)
  num_done.value = 0

  # use if you have trajectory files to read
  #files = os.listdir(input_dir)

  # otherwise, give the initial density
  #densities = [1e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9] # multiple
  densities = [1.0e9] # URCA
  #densities = [2e6] #XRB Aprox13
  #densities = [1e6] #SubCh
  
  # number of jobs to run
  size = len(densities)
  args = [(r, size) for r in densities]
  
  # run skynet in parallel
  pool.map_async(run_skynet, args)

  # done submitting jobs
  pool.close()
  pool.join()

# monitor progress with
# $ cd work_dir
# $ watch "bash -c 'tail -n 3 *.log'"


