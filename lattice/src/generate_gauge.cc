// generate_gauge.cc -- pure SU(3) Wilson HMC with CLI beta and seed, for INDEPENDENT parallel
// streams (Eternal Dawn lattice). Adapted verbatim from Grid's tests/hmc/Test_hmc_WilsonGauge.cc
// (proven to build and run on the GH200); the ONLY changes are: beta and the RNG seeds are read
// from the command line (--beta, --seed) instead of hardcoded, so K copies run as K genuinely
// independent Markov chains rather than K identical ones. Build via lattice/src/Makefile.
//
//   generate_gauge --grid Lx.Ly.Lz.Lt --beta 5.6 --seed 1001 \
//                  --StartingType HotStart --Trajectories 200
// Writes NERSC ckpoint_lat.<n> every trajectory into the CURRENT directory -- run each stream
// from its own directory (run/00_generate.sh does this).
#include <Grid/Grid.h>
#include <sstream>
#include <string>

int main(int argc, char **argv) {
  using namespace Grid;
  Grid_init(&argc, &argv);

  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;
  HMCWrapper TheHMC;

  // ---- our additions: beta and seed from the command line ----
  RealD beta = 5.6;
  int seed = 1;
  if (GridCmdOptionExists(argv, argv + argc, "--beta"))
    beta = std::stod(GridCmdOptionPayload(argv, argv + argc, "--beta"));
  if (GridCmdOptionExists(argv, argv + argc, "--seed")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--seed");
    GridCmdOptionInt(s, seed);
  }

  TheHMC.Resources.AddFourDimGrid("gauge");

  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 1;
  CPparams.format = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  // independent RNG seeds per stream, derived from --seed (10 distinct integers)
  std::ostringstream ss, ps;
  for (int i = 0; i < 5; ++i) ss << (seed + i) << (i < 4 ? " " : "");
  for (int i = 0; i < 5; ++i) ps << (seed + 5 + i) << (i < 4 ? " " : "");
  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = ss.str();
  RNGpar.parallel_seeds = ps.str();
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  WilsonGaugeActionR Waction(beta);
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  TheHMC.TheAction.push_back(Level1);

  TheHMC.Parameters.MD.MDsteps = 20;
  TheHMC.Parameters.MD.trajL = 1.0;

  std::cout << GridLogMessage << "generate_gauge: beta=" << beta << " seed=" << seed
            << " serial_seeds='" << RNGpar.serial_seeds
            << "' parallel_seeds='" << RNGpar.parallel_seeds << "'" << std::endl;

  TheHMC.ReadCommandLine(argc, argv);   // --grid, --Trajectories, --StartingType, ...
  TheHMC.Run();

  Grid_finalize();
  return 0;
}
