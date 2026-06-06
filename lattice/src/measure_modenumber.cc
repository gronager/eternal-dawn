// measure_modenumber.cc -- Chebyshev moments of the Dirac mode number (Eternal Dawn, gamma_m gate).
// The gamma_m gate via the Kernel Polynomial Method (KPM), the scalable, solver-free route used in
// the mode-number literature. Instead of computing eigenvalues (the Lanczos route, which needs a
// high-order filter and convergence babysitting) we estimate the Chebyshev moments
//     mu_k = Tr T_k(Xtilde),   Xtilde = (D^dag D)/(xmax/2) - 1   in [-1,1]
// by a stochastic trace -- random sources and D^dag D matvecs ONLY, no eigenvalues, no CG. The
// validated Python side (cartasis_sims.lattice.mode_number_from_chebyshev_moments) combines the
// moments with the Jackson-damped step to give nu(M) for the whole M grid, then gamma_m. One run
// of moments covers every threshold M.
//
// Representation: --rep fundamental (default) or --rep sextet (2-index symmetric, the candidate
// walking sector). Both use Grid's own tested Wilson operator -- fundamental WilsonFermion, or
// WilsonTwoIndexSymmetricFermionD on the sextet links built from the fundamental config by Grid's
// TwoIndexRep::update_representation. No representation code of our own.
//
// VALIDATION FIRST (stand on the analytic answer): run --free and check the moments reproduce the
// free Wilson mode number that cartasis_sims already validates analytically:
//     measure_modenumber --grid 16.16.16.16 --free --mass 0.0 --Nmom 600 --Nnoise 12 > moments.txt
//   --rep sextet on --free must give exactly TWICE the fundamental mode number (6 vs 3 colours,
//   same |lambda|): mu_0 = 24*sites and the same gamma_m. Only then an interacting config:
//     measure_modenumber --grid 32.32.32.64 --config <NERSC cfg> --rep sextet --mass <m> --Nmom 800
//
// The only matvec is Grid's Wilson D^dag D (MdagMLinearOperator), already confirmed correct (its
// spectral ceiling matched the free maximum). The recurrence stays bounded because Xtilde maps the
// spectrum into [-1,1] (|T_k| <= 1), so high Nmom is stable.
//
// Output (stdout): "XMAX <value>" then "MOMENT <k> <mu_k>" for k=0..Nmom. Build via the Makefile.
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
static std::string cli_str(char **a, char **e, const std::string &flag, const std::string &def) {
  if (!GridCmdOptionExists(a, e, flag)) return def;
  return GridCmdOptionPayload(a, e, flag);
}

// KPM Chebyshev moments mu_k = Tr T_k(Xtilde), Xtilde = (D^dag D)/(xmax/2) - 1, by a stochastic
// trace. Templated on the Wilson action so the SAME code serves the fundamental and the sextet
// (two-index symmetric) representations -- only the gauge field type and colour count differ.
// `ncol` is the representation's colour dimension (3 fundamental, 6 sextet); mu_0 = Tr I =
// 4(spin) * ncol * sites is known exactly and used to fix the noise normalisation.
template <class Action>
std::vector<RealD> kpm_moments(typename Action::GaugeField &U,
                               GridCartesian *UGrid, GridRedBlackCartesian *UrbGrid,
                               GridParallelRNG &RNG, RealD mass, int Nmom, int Nnoise,
                               int powIts, RealD &xmax, RealD ncol) {
  typedef typename Action::FermionField Field;
  Action Dw(U, *UGrid, *UrbGrid, mass);
  MdagMLinearOperator<Action, Field> HermOp(Dw);

  // spectral bound: largest eigenvalue of D^dag D via power iteration (must bound the spectrum so
  // Xtilde stays in [-1,1]). Override with --xmax on a known case.
  if (xmax <= 0.0) {
    Field v(UGrid), w(UGrid);
    gaussian(RNG, v);
    v = v * (1.0 / std::sqrt(norm2(v)));
    RealD lam = 0.0;
    for (int i = 0; i < powIts; ++i) {
      HermOp.HermOp(v, w);
      lam = std::sqrt(norm2(w));                     // -> largest eigenvalue as v aligns
      v = w * (1.0 / lam);
    }
    xmax = 1.1 * lam;                                // 10% headroom
    std::cout << GridLogMessage << "measure_modenumber: power-method lambda_max(D^dag D) ~ "
              << lam << std::endl;
  }
  std::cout << "XMAX " << std::setprecision(12) << xmax << std::endl;
  const RealD a = xmax / 2.0;                         // Xtilde = HermOp/a - 1
  const RealD dim = 4.0 * ncol * (RealD)UGrid->gSites();   // Tr I = Ns(4) * ncol * sites

  std::vector<RealD> mu(Nmom + 1, 0.0);
  Field eta(UGrid), T0(UGrid), T1(UGrid), T2(UGrid), tmp(UGrid);
  for (int s = 0; s < Nnoise; ++s) {
    gaussian(RNG, eta);
    T0 = eta;                                         // T_0(Xtilde) eta = eta
    HermOp.HermOp(eta, tmp);                          // T_1 = Xtilde eta = HermOp(eta)/a - eta
    T1 = tmp * (1.0 / a) - eta;
    mu[0] += real(innerProduct(eta, T0));
    mu[1] += real(innerProduct(eta, T1));
    for (int k = 2; k <= Nmom; ++k) {
      HermOp.HermOp(T1, tmp);                         // T_k = 2 Xtilde T_{k-1} - T_{k-2}
      T2 = tmp * (2.0 / a) - T1 * 2.0 - T0;
      mu[k] += real(innerProduct(eta, T2));
      T0 = T1;
      T1 = T2;
    }
    std::cout << GridLogMessage << "measure_modenumber: noise source " << (s + 1) << "/" << Nnoise
              << " done" << std::endl;
  }
  for (auto &m : mu) m /= (RealD)Nnoise;
  const RealD scale = dim / mu[0];                    // force mu_0 == 4*ncol*sites
  for (auto &m : mu) m *= scale;
  return mu;
}

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  // ---- fundamental gauge field: identity (--free, the analytic validation) or a NERSC config ----
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

  RealD mass  = cli_real(argv, argv + argc, "--mass", -0.5);    // Wilson mass; near-critical for light modes
  int Nmom    = cli_int(argv, argv + argc, "--Nmom", 600);      // Chebyshev moments (k = 0..Nmom)
  int Nnoise  = cli_int(argv, argv + argc, "--Nnoise", 12);     // stochastic sources
  RealD xmax  = cli_real(argv, argv + argc, "--xmax", 0.0);     // spectral bound; <=0 => estimate
  int powIts  = cli_int(argv, argv + argc, "--power-iters", 30);
  int seed    = cli_int(argv, argv + argc, "--seed", 1234);
  std::string rep = cli_str(argv, argv + argc, "--rep", "fundamental");

  GridParallelRNG RNG(UGrid);
  RNG.SeedFixedIntegers(std::vector<int>({seed, seed + 1, seed + 2, seed + 3}));

  // ---- dispatch on representation. Fundamental: act on Umu directly. Sextet (2-index symmetric):
  //      build the 6-dim representation links from the fundamental via Grid's TwoIndexRep, then run
  //      the identical KPM on the sextet Wilson D^dag D. On --free the rep links are identity, so
  //      the sextet mode number is exactly 2x the fundamental (6 vs 3 colours) at the same gamma_m.
  std::vector<RealD> mu;
  if (rep == "sextet" || rep == "twoindexsym") {
    std::cout << GridLogMessage << "measure_modenumber: representation = sextet (two-index symmetric, dim "
              << TwoIndexSymmetricRepresentation::Dimension << ")" << std::endl;
    TwoIndexSymmetricRepresentation Rep(UGrid);
    Rep.update_representation(Umu);
    mu = kpm_moments<WilsonTwoIndexSymmetricFermionD>(
        Rep.U, UGrid, UrbGrid, RNG, mass, Nmom, Nnoise, powIts, xmax,
        (RealD)TwoIndexSymmetricRepresentation::Dimension);
  } else {
    std::cout << GridLogMessage << "measure_modenumber: representation = fundamental (dim " << Nc << ")"
              << std::endl;
    mu = kpm_moments<WilsonFermion<WilsonImplD>>(
        Umu, UGrid, UrbGrid, RNG, mass, Nmom, Nnoise, powIts, xmax, (RealD)Nc);
  }

  for (int k = 0; k <= Nmom; ++k)
    std::cout << "MOMENT " << k << " " << std::setprecision(12) << mu[k] << std::endl;

  Grid_finalize();
  return 0;
}
