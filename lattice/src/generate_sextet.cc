// generate_sextet.cc -- dynamical SU(3) N_f=2 sextet (two-index symmetric) Wilson HMC, with CLI
// (beta, mass, seed), for the Eternal Dawn gamma_m gate. This is rung 3: the candidate walking
// ensemble with sextet SEA quarks (the quenched valence scan, run/04, only put sextet quarks in
// the measurement, not the vacuum). Built on Grid's own higher-rep HMC machinery
// (GenericHMCRunnerHirep + WilsonTwoIndexSymmetricFermion + TwoFlavourPseudoFermionAction), exactly
// as tests/hmc/Test_hmc_WilsonTwoIndexSymmetricFermionGauge.cc -- the ONLY changes are that beta,
// mass, the RNG seeds and the integrator cadence are read from the command line instead of being
// hardcoded, so we can tune (beta, mass) and run independent streams. No representation code of ours.
//
//   generate_sextet --grid Lx.Ly.Lz.Lt --beta 5.4 --mass -0.95 --seed 1001 \
//                   --mdsteps 20 --save-interval 5 --StartingType HotStart --Trajectories 2000
// Writes NERSC ckpoint_lat.<n> into the CURRENT directory. Build via lattice/src/Makefile.
//
// (beta, mass) ARE THE PHYSICS KNOBS: tune mass toward the sextet critical point (run/04 brackets
// it) and choose beta in the scaling window; their values decide walks-vs-conformal. trajL=1.
#include <Grid/Grid.h>
#include <sstream>
#include <string>

int main(int argc, char **argv) {
  using namespace Grid;
  Grid_init(&argc, &argv);

  // higher-rep HMC: evolve the fundamental links, with the sextet (two-index symmetric) sea
  typedef Representations<FundamentalRepresentation, TwoIndexSymmetricRepresentation> TheRepresentations;
  typedef GenericHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper;
  typedef WilsonTwoIndexSymmetricImplR FermionImplPolicy;
  typedef WilsonTwoIndexSymmetricFermionD FermionAction;
  typedef typename FermionAction::FermionField FermionField;

  HMCWrapper TheHMC;

  // ---- our additions: beta, mass, seed, integrator steps, checkpoint cadence from the CLI ----
  RealD beta = 5.4;        // bare gauge coupling -- TUNE (scaling window)
  RealD mass = -0.95;      // sextet Wilson mass  -- TUNE toward critical (run/04 brackets m_c)
  int seed = 1;
  int mdsteps = 20;        // leapfrog steps per trajectory; raise ~V^{1/4} / for the fermion force
  int saveInterval = 5;    // write a NERSC config every Nth trajectory
  if (GridCmdOptionExists(argv, argv + argc, "--beta"))
    beta = std::stod(GridCmdOptionPayload(argv, argv + argc, "--beta"));
  if (GridCmdOptionExists(argv, argv + argc, "--mass"))
    mass = std::stod(GridCmdOptionPayload(argv, argv + argc, "--mass"));
  if (GridCmdOptionExists(argv, argv + argc, "--seed")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--seed");
    GridCmdOptionInt(s, seed);
  }
  if (GridCmdOptionExists(argv, argv + argc, "--mdsteps")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--mdsteps");
    GridCmdOptionInt(s, mdsteps);
  }
  if (GridCmdOptionExists(argv, argv + argc, "--save-interval")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--save-interval");
    GridCmdOptionInt(s, saveInterval);
  }

  TheHMC.Resources.AddFourDimGrid("gauge");

  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = saveInterval;
  CPparams.format = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  std::ostringstream ss, ps;
  for (int i = 0; i < 5; ++i) ss << (seed + i) << (i < 4 ? " " : "");
  for (int i = 0; i < 5; ++i) ps << (seed + 5 + i) << (i < 4 ? " " : "");
  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = ss.str();
  RNGpar.parallel_seeds = ps.str();
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // ---- actions: sextet N_f=2 pseudofermion (Level1) + Wilson gauge (Level2, finer) ----
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  TwoIndexSymmetricRepresentation::LatticeField U(GridPtr);     // sextet links (updated by the Hirep runner)

  FermionAction FermOp(U, *GridPtr, *GridRBPtr, mass);
  ConjugateGradient<FermionField> CG(1.0e-8, 2000, false);
  TwoFlavourPseudoFermionAction<FermionImplPolicy> Nf2(FermOp, CG, CG);
  Nf2.is_smeared = false;   // CRITICAL: route the force through the representation chain, not the
                            // (absent) smearing chain. Without this the sea-quark force is zeroed
                            // -- the action is still evaluated (Metropolis is correct) but the MD
                            // sees no fermion force, i.e. the ensemble is quenched-in-disguise.

  WilsonGaugeActionR Waction(beta);

  ActionLevel<HMCWrapper::Field, TheRepresentations> Level1(1);
  Level1.push_back(&Nf2);
  ActionLevel<HMCWrapper::Field, TheRepresentations> Level2(4);  // gauge updated 4x per fermion step
  Level2.push_back(&Waction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  TheHMC.Parameters.MD.MDsteps = mdsteps;
  TheHMC.Parameters.MD.trajL = 1.0;

  std::cout << GridLogMessage << "generate_sextet: SU(3) Nf=2 sextet  beta=" << beta
            << " mass=" << mass << " seed=" << seed << " mdsteps=" << mdsteps
            << " saveInterval=" << saveInterval << std::endl;

  TheHMC.ReadCommandLine(argc, argv);   // --grid, --Trajectories, --Thermalizations, --StartingType, ...
  TheHMC.Run();

  Grid_finalize();
  return 0;
}
