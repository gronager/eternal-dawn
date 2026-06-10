// measure_decay_constant.cc -- the substrate scale on the lattice: the pion decay constant f_pi and
// the chiral condensate scale v, the missing input for the framework's sharp prediction
//
//     sigma = 2 * pi * v^2        (target L2: ONE substrate sets BOTH the mass scale and the
//                                  confinement string tension).
//
// We already have sqrt(sigma)*a = 0.385 (static potential, run/01) and m_pi*a = 1.269 (spectrum,
// run/06). What is owed is v (equivalently f_pi). This program supplies it on the EXISTING dynamical
// N_f=2 SU(3) FUNDAMENTAL Wilson ensemble.
//
// WHAT & HOW. Reusing measure_baryon's (validated) solver and pion contraction verbatim: read one
// configuration, solve a point-source Wilson-quark propagator S(x;0), and build the zero-momentum
// pseudoscalar-pseudoscalar correlator
//
//     C_PP(t) = sum_x Tr[ gamma5 S(x,t;0) gamma5 S(x,t;0)^dag ] = sum_x Tr[ S(x)^dag S(x) ]
//
// (the SECOND equality is gamma5-hermiticity, S(x;0) = gamma5 S(0;x)^dag gamma5 -- the SAME object
// measure_baryon already prints as C_pi). At zero momentum this is the standard PP correlator whose
// time profile gives m_pi and whose amplitude gives the pseudoscalar overlap.
//
// EXTRACTION (PCAC / pseudoscalar-density route -- robust, needs only S and the bare quark mass). On
// the lattice the analysis (cartasis_sims.lattice.decay_constant) fits the two-sided cosh
//
//     C_PP(t) = (G_PP / (2 m_pi)) * (e^{-m_pi t} + e^{-m_pi (T-t)}),   G_PP = |<0|Pbar|pi>|^2,
//
// for m_pi*a and G_PP, then takes
//
//     f_pi * a = (2 * m_q / m_pi^2) * sqrt(G_PP)      (PCAC, bare/unrenormalised),
//
// with m_q the BARE Wilson quark mass = the input MASS shifted by the critical mass m_crit. This is
// the UNRENORMALISED number: the renormalised f_pi owes Z_A, Z_P and m_crit (each an O(1) factor),
// which are NOT supplied here and are flagged in the output and the analysis.
//
// NORMALISATION CONVENTION. We use the f_pi ~ 92 MeV (chiral F_pi) normalisation; multiply by
// sqrt(2) for the f_pi = sqrt(2) F_pi ~ 130 MeV convention. Stated here and in the analysis so the
// number is never ambiguous.
//
// THE SUBSTRATE SCALE v. The framework identifies v = f_pi (the NJL / chiral-soliton identification --
// an ASSUMPTION of the framework, stated not hidden). So v*a = f_pi*a here. The alternative
// v^2 ~ |<qbar q>| (the chiral condensate) is NOT computed by this program (it would need a
// disconnected/all-to-all trace); it is noted as the cheap-if-available cross-check.
//
//   measure_decay_constant --grid Lx.Ly.Lz.Lt --config <NERSC cfg> --mass <m> [--cg-tol 1e-8] \
//                          [--cg-maxit 30000]
//
// Build with lattice/src/Makefile (uses grid-config). Self-contained: uses the SAME core Grid
// primitives as measure_baryon (WilsonFermionD, ConjugateGradient, Gamma, sliceSum), so it is robust
// across Grid versions. The emitted "PP t C_PP(t)" rows feed run/11_decay_constant.sh -> decay_raw.dat.
#include <Grid/Grid.h>
#include <iomanip>

using namespace Grid;

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  typedef WilsonImplD FImpl;

  Coordinate latt = GridDefaultLatt();
  const int Tdir = Nd - 1;
  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(
      latt, GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  if (!GridCmdOptionExists(argv, argv + argc, "--config") ||
      !GridCmdOptionExists(argv, argv + argc, "--mass")) {
    std::cout << "usage: measure_decay_constant --grid Lx.Ly.Lz.Lt --config <NERSC cfg> --mass <m> "
                 "[--cg-tol 1e-8] [--cg-maxit 30000]" << std::endl;
    Grid_finalize();
    return 1;
  }
  std::string cfg = GridCmdOptionPayload(argv, argv + argc, "--config");
  RealD mass = std::stod(GridCmdOptionPayload(argv, argv + argc, "--mass"));
  RealD cg_tol = 1.0e-8;
  int cg_maxit = 30000;
  if (GridCmdOptionExists(argv, argv + argc, "--cg-tol"))
    cg_tol = std::stod(GridCmdOptionPayload(argv, argv + argc, "--cg-tol"));
  if (GridCmdOptionExists(argv, argv + argc, "--cg-maxit")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--cg-maxit");
    GridCmdOptionInt(s, cg_maxit);
  }

  // ---- load the (dynamical N_f=2 fundamental) configuration ----
  LatticeGaugeFieldD Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, cfg);
  RealD plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << "# config " << cfg << "  plaquette " << std::setprecision(10) << plaq
            << "  mass " << mass << std::endl;
  // CONVENTION reminder in the raw stream: f_pi ~ 92 MeV (F_pi); analysis owes Z_A, Z_P, m_crit.
  std::cout << "# f_pi via PCAC: f_pi*a = (2 m_q / m_pi^2) sqrt(G_PP), m_q = mass - m_crit (m_crit OWED)"
            << std::endl;
  std::cout << "# UNRENORMALISED (Z_A, Z_P owed); v = f_pi is the framework NJL identification; "
               "f_pi~92 MeV (F_pi) convention" << std::endl;

  // ---- Wilson quark, point source at the origin (identity in spin x colour) -- as measure_baryon --
  WilsonFermionD Dw(Umu, *UGrid, *UrbGrid, mass);
  LatticePropagatorD src(UGrid);
  src = Zero();
  {
    typename LatticePropagatorD::vector_object::scalar_object Sid = Zero();
    for (int s = 0; s < Ns; ++s)
      for (int c = 0; c < Nc; ++c) Sid()(s, s)(c, c) = ComplexD(1.0, 0.0);
    Coordinate origin({0, 0, 0, 0});
    pokeSite(Sid, src, origin);
  }

  // ---- solve the 12 propagator columns: (Mdag M) x = Mdag b, then M x = b (as measure_baryon) ----
  LatticePropagatorD prop(UGrid);
  prop = Zero();
  MdagMLinearOperator<WilsonFermionD, LatticeFermionD> HermOp(Dw);
  ConjugateGradient<LatticeFermionD> CG(cg_tol, cg_maxit);
  for (int s = 0; s < Ns; ++s) {
    for (int c = 0; c < Nc; ++c) {
      LatticeFermionD fsrc(UGrid), fmd(UGrid), fsol(UGrid);
      PropToFerm<FImpl>(fsrc, src, s, c);
      Dw.Mdag(fsrc, fmd);
      fsol = Zero();
      CG(HermOp, fmd, fsol);
      FermToProp<FImpl>(prop, fsol, s, c);
    }
  }

  // ---- the zero-momentum PSEUDOSCALAR correlator (exactly measure_baryon's pion) ----
  // C_PP(t) = sum_x Tr[ gamma5 S(x,t;0) gamma5 S(x,t;0)^dag ] = sum_x Tr[ S(x)^dag S(x) ] by
  // gamma5-hermiticity (gamma5 S gamma5 = S(0;x)^dag), the robust one-trace pion anchor.
  LatticeComplexD pion(UGrid);
  pion = trace(adj(prop) * prop);
  std::vector<TComplexD> pion_t;
  sliceSum(pion, pion_t, Tdir);

  // ---- emit the PP correlator (real part): PP  t  C_PP(t)  (feeds decay_raw.dat / decay_constant) ----
  std::cout << "# PP  t  C_PP(t)" << std::endl;
  for (size_t t = 0; t < pion_t.size(); ++t) {
    ComplexD cp = TensorRemove(pion_t[t]);
    std::cout << "PP " << t << " " << std::setprecision(12) << real(cp) << std::endl;
  }

  Grid_finalize();
  return 0;
}
