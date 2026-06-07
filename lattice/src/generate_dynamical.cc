// generate_dynamical.cc -- dynamical SU(3) N_f=2 FUNDAMENTAL Wilson HMC, with CLI (beta, mass, seed),
// for the Eternal Dawn torsiton PRODUCTION run. This is the sea the quenched pilot (run/06 on
// pure-gauge configs) lacked: here the fermion determinant is included, so the vacuum is the real,
// fluctuating, chiral-symmetry-breaking condensate -- the "weirder than mean field" vacuum where the
// physical torsiton mass, and any excited rungs (candidate further generations), actually live.
//
// FAITHFUL THEORY (not the dead sextet): the torsiton is the SU(3)-FUNDAMENTAL baryon (colour from
// the Pauli label), so the sea quarks are fundamental -- the standard, simplest dynamical SU(3)
// setup (no higher-rep machinery). N_f=2 (two-flavour pseudofermion) is the cheap QCD-like choice
// that confines and breaks chiral symmetry; N_f is a physics knob. Built on Grid's GenericHMCRunner
// exactly as tests/hmc/Test_hmc_WilsonFermionGauge.cc -- the only changes are CLI (beta, mass, seed,
// integrator cadence) so we can tune and run independent streams.
//
//   generate_dynamical --grid Lx.Ly.Lz.Lt --beta 5.6 --mass -0.5 --seed 1001 \
//                      --mdsteps 20 --save-interval 5 --StartingType HotStart --Trajectories 2000
// Writes NERSC ckpoint_lat.<n> into the CURRENT directory (run each stream in its own dir).
//
// (beta, mass) ARE THE KNOBS: lower mass toward the DYNAMICAL critical point for a light pion (the
// sea suppresses the exceptional configs that walled the quenched scan, so you can go lighter); pick
// beta in the scaling window. Build via lattice/src/Makefile. trajL=1.
#include <Grid/Grid.h>
#include <sstream>
#include <string>

int main(int argc, char **argv) {
  using namespace Grid;
  Grid_init(&argc, &argv);

  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;     // fundamental: no Hirep machinery
  typedef WilsonImplR FermionImplPolicy;
  typedef WilsonFermion<FermionImplPolicy> FermionAction;  // this Grid lacks the WilsonFermionR
                                                           // typedef; instantiate the template
  typedef typename FermionAction::FermionField FermionField;

  HMCWrapper TheHMC;

  // ---- CLI: beta, mass, seed, integrator steps, checkpoint cadence ----
  RealD beta = 5.6;        // bare gauge coupling -- TUNE (scaling window)
  RealD mass = -0.5;       // fundamental Wilson sea-quark mass -- TUNE toward the dynamical critical
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

  std::ostringstream ss, ps;                              // independent seeds per stream
  for (int i = 0; i < 5; ++i) ss << (seed + i) << (i < 4 ? " " : "");
  for (int i = 0; i < 5; ++i) ps << (seed + 5 + i) << (i < 4 ? " " : "");
  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = ss.str();
  RNGpar.parallel_seeds = ps.str();
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // ---- actions: N_f=2 fundamental pseudofermion (Level1) + Wilson gauge (Level2, finer) ----
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  LatticeGaugeField U(GridPtr);                          // fundamental links (updated by the runner)

  FermionAction FermOp(U, *GridPtr, *GridRBPtr, mass);
  ConjugateGradient<FermionField> CG(1.0e-8, 2000, false);
  TwoFlavourPseudoFermionAction<FermionImplPolicy> Nf2(FermOp, CG, CG);
  Nf2.is_smeared = false;                                // no smearing chain (we evolve the bare links)

  WilsonGaugeActionR Waction(beta);

  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Nf2);
  ActionLevel<HMCWrapper::Field> Level2(4);              // gauge updated 4x per fermion step
  Level2.push_back(&Waction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  TheHMC.Parameters.MD.MDsteps = mdsteps;
  TheHMC.Parameters.MD.trajL = 1.0;

  std::cout << GridLogMessage << "generate_dynamical: SU(3) Nf=2 fundamental  beta=" << beta
            << " mass=" << mass << " seed=" << seed << " mdsteps=" << mdsteps
            << " saveInterval=" << saveInterval << std::endl;

  TheHMC.ReadCommandLine(argc, argv);   // --grid, --Trajectories, --Thermalizations, --StartingType, ...
  TheHMC.Run();

  Grid_finalize();
  return 0;
}
