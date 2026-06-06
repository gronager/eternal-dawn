// measure_modenumber.cc -- the low Dirac spectrum of a gauge config (Eternal Dawn, gamma_m gate).
// RUNG 1 of the gamma_m ladder: measure the eigenvalues of the Hermitian Dirac operator on a real
// gauge configuration, using Grid's FUNDAMENTAL Wilson fermion (its own tested operator -- we add
// no representation code here). The program prints the lowest |lambda| (one "EIG <value>" line
// each); the validated Python side (cartasis_sims.lattice: mode_number_from_eigenvalues ->
// anomalous_dimension_from_mode_number) counts nu(M) and fits gamma_m. Keeping the C++ to just the
// spectrum shrinks the debug surface.
//
// VALIDATION FIRST (stand on the analytic answer of rung 0):
//   measure_modenumber --grid 32.32.32.32 --free --mass <m_crit> --Nstop 80
//     -> identity links = free field; the |lambda| must reproduce free_wilson_mode_number(),
//        i.e. nu(M)~M^4, gamma_m~0. Only once that matches do we trust an interacting config:
//   measure_modenumber --grid 32.32.32.64 --config <NERSC cfg> --mass <m>
//
// The eigenvalues come from Grid's ImplicitlyRestartedLanczos on M^dag M (eigenvalues |lambda|^2),
// Chebyshev-accelerated at the low end. The IRL/Chebyshev knobs are CLI flags -- tune convergence
// without recompiling (watch the "Nconv" line; raise --Nm / --cheb-ord if it under-converges).
//
// Build via lattice/src/Makefile (adds this to PROGS).
#include <Grid/Grid.h>
#include <cmath>

using namespace Grid;

static int cli_int(char **a, char **e, const std::string &flag, int def) {
  if (!GridCmdOptionExists(a, e, flag)) return def;
  std::string s = GridCmdOptionPayload(a, e, flag); int v = def; GridCmdOptionInt(s, v); return v;
}
static RealD cli_real(char **a, char **e, const std::string &flag, RealD def) {
  if (!GridCmdOptionExists(a, e, flag)) return def;
  return std::stod(GridCmdOptionPayload(a, e, flag));
}

// Fundamental Wilson fermion, double precision. Bound to the primary template + the explicit
// double impl (WilsonImplD) -- this Grid build dropped the WilsonFermionR/D convenience typedefs,
// and vComplexD matches the SIMD layout the grids below are built with.
typedef WilsonImplD                  FermImpl;
typedef WilsonFermion<FermImpl>      FermionAction;
typedef FermionAction::FermionField  FermionField;

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  // ---- gauge field: identity (--free, the analytic validation) or a NERSC config ----
  LatticeGaugeField Umu(UGrid);
  if (GridCmdOptionExists(argv, argv + argc, "--free")) {
    SU<Nc>::ColdConfiguration(Umu);                 // U = 1 everywhere -> free Dirac operator
    std::cout << GridLogMessage << "measure_modenumber: FREE field (identity links)" << std::endl;
  } else {
    std::string cfg = GridCmdOptionPayload(argv, argv + argc, "--config");
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, cfg);
    std::cout << GridLogMessage << "measure_modenumber: config " << cfg << std::endl;
  }

  // ---- fundamental Wilson Dirac operator (Grid's own, tested) ----
  RealD mass = cli_real(argv, argv + argc, "--mass", -0.5);   // Wilson mass; near-critical for light modes
  FermionAction Dw(Umu, *UGrid, *UrbGrid, mass);

  // Hermitian positive operator M^dag M; its eigenvalues are |lambda|^2
  MdagMLinearOperator<FermionAction, FermionField> HermOp(Dw);

  // ---- IRL knobs (CLI-tunable so convergence is fixed without a recompile) ----
  int Nstop   = cli_int(argv, argv + argc, "--Nstop", 60);    // eigenvalues we want converged
  int Nk      = cli_int(argv, argv + argc, "--Nk", 80);       // working set
  int Nm      = cli_int(argv, argv + argc, "--Nm", 160);      // Krylov dimension (Nm > Nk > Nstop)
  RealD resid = cli_real(argv, argv + argc, "--resid", 1e-8);
  int MaxIt   = cli_int(argv, argv + argc, "--MaxIt", 200);
  // Chebyshev low-mode filter for IRL. Convention (cf. Grid tests/lanczos/Test_dwf_lanczos.cc,
  // which uses Cheby(6e-7, 5.5, 4001)): cheb-lo ~ 0 (bottom of spectrum), cheb-hi = the spectral
  // MAX of M^dag M (so the whole bulk is in-band and suppressed), and cheb-ord must be LARGE
  // (thousands) -- that order is what actually amplifies the low modes. A small order (we had 60)
  // gives no separation and IRL returns junk biased to the high end. For free Wilson the M^dag M
  // max is 8^2 = 64, so cheb-hi defaults just above that.
  RealD chLo  = cli_real(argv, argv + argc, "--cheb-lo", 0.001);
  RealD chHi  = cli_real(argv, argv + argc, "--cheb-hi", 70.0);
  int chOrd   = cli_int(argv, argv + argc, "--cheb-ord", 2000);

  // IRL requires Nm > Nk > Nstop. Enforce it here with a clear message rather than letting Grid
  // trip a cryptic 'k < Nm' assert deep in ImplicitlyRestartedLanczos (e.g. when --Nm is lowered
  // below the default --Nk). eval/evec are sized from Nm just below, so fixing it now is safe.
  if (Nk <= Nstop) { Nk = Nstop + std::max(1, Nstop / 2); std::cout << GridLogMessage
      << "measure_modenumber: Nk <= Nstop; bumping Nk to " << Nk << std::endl; }
  if (Nm <= Nk)    { Nm = Nk + Nstop;                     std::cout << GridLogMessage
      << "measure_modenumber: Nm <= Nk; bumping Nm to " << Nm << std::endl; }

  Chebyshev<FermionField> Cheby(chLo, chHi, chOrd);            // amplify the low end for IRL
  FunctionHermOp<FermionField> AccelOp(Cheby, HermOp);
  PlainHermOp<FermionField> PlainOp(HermOp);
  ImplicitlyRestartedLanczos<FermionField> IRL(AccelOp, PlainOp, Nstop, Nk, Nm, resid, MaxIt);

  std::vector<RealD> eval(Nm);
  std::vector<FermionField> evec(Nm, UGrid);
  FermionField src(UGrid);
  GridParallelRNG RNG(UGrid);
  RNG.SeedFixedIntegers(std::vector<int>({1, 2, 3, 4}));
  gaussian(RNG, src);

  int Nconv = 0;
  IRL.calc(eval, evec, src, Nconv);

  std::cout << GridLogMessage << "measure_modenumber: Nconv=" << Nconv
            << " (eigenvalues of M^dag M; |lambda| = sqrt)" << std::endl;
  for (int i = 0; i < Nconv; ++i)
    std::cout << "EIG " << std::setprecision(12) << std::sqrt(std::abs(eval[i])) << std::endl;

  Grid_finalize();
  return 0;
}
