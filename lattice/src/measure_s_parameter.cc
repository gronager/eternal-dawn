// measure_s_parameter.cc -- a FIRST lattice estimate of the electroweak S parameter (Peskin-Takeuchi)
// for the dynamical N_f=2 SU(3) FUNDAMENTAL Wilson ensemble, via the isovector VECTOR minus
// AXIAL-VECTOR current correlators (Pi_V - Pi_A). This tests the framework's load-bearing assumption
// (Ch. "The Forces from the Field"): does the electroweak-breaking sector clear S < 0.1, or sit in the
// QCD-like graveyard at S ~ 0.25?
//
// HONEST EXPECTATION (baked in here and in the analysis). The FUNDAMENTAL N_f=2 rep is QCD-like (NOT
// near-conformal): this proxy is EXPECTED to land near the QCD value S ~ 0.25. That is INFORMATIVE --
// it characterises the sector and TESTS the assumption -- but it is NOT a validation of S < 0.1, which
// needs the near-conformal / propagating-torsion regime (target L5). Do NOT oversell: this is the
// QCD-like control point, not the walking result.
//
// WHAT & HOW. Reusing measure_decay_constant's (validated) solver and zero-momentum meson machinery
// VERBATIM: read one configuration, solve a point-source Wilson-quark propagator S(x;0), and build the
// CONNECTED isovector correlators. For an isovector interpolator O_Gamma = qbar Gamma q, the connected
// zero-momentum correlator is
//
//   C_Gamma(t) = sum_x Re Tr[ Gamma S(x,t;0) Gammabar gamma5 S(x,t;0)^dag gamma5 ],
//                Gammabar = gamma5 Gamma^dag gamma5.
//
// (The DISCONNECTED pieces cancel for the ISOVECTOR current -- the two flavours enter with opposite
// charge, so the quark-loop trace cancels in the difference -- so only the connected contraction above
// is needed. Stated, not hidden.) We implement two channels, SUMMED over the spatial polarisations
// i = 1,2,3 (x,y,z):
//   * VECTOR        Gamma = gamma_i           -> C_V(t)   (the rho channel)
//   * AXIAL-VECTOR  Gamma = gamma_i gamma5     -> C_A(t)   (the a1 channel)
// The pion channel Gamma = gamma5 (= trace(adj(prop)*prop)) is already in measure_decay_constant; here
// the ONLY change is the gamma insertion, using the SAME Gamma::Algebra access the existing code uses.
//
// For these gammas Gammabar = gamma5 Gamma^dag gamma5 evaluates to: gamma_i for the vector
// (gamma5 gamma_i^dag gamma5 = -(-gamma_i) = ... ) and gamma_i gamma5 for the axial -- but rather than
// hard-code the algebra we apply the definition literally with the Gamma objects, so the contraction is
// correct by construction for any insertion.
//
// EXTRACTION (analysis: cartasis_sims.lattice.s_parameter_proxy). C_V and C_A are cosh-fit (the SAME
// _cosh_fit_pp / correlator_mass helpers as the decay constant) for the rho mass M_V, the a1 mass M_A,
// and the corresponding amplitudes / decay constants F_V, F_A. The lowest-resonance (vector-meson-
// dominance / Weinberg-sum-rule-saturated) S-proxy is
//
//   S = 4 pi ( F_V^2 / M_V^2  -  F_A^2 / M_A^2 ).
//
// Also reported: the direct walking diagnostic M_A / M_V (-> 1 if conformal, ~1.6 QCD-like) and the
// integrated chiral order parameter sum_t (C_V(t) - C_A(t)). HONEST: lowest-resonance saturation;
// currents UNRENORMALISED (Z_V, Z_A owed); single heavy sea mass; fundamental rep expected QCD-like.
//
//   measure_s_parameter --grid Lx.Ly.Lz.Lt --config <NERSC cfg> --mass <m> [--cg-tol 1e-8] \
//                       [--cg-maxit 30000]
//
// Build with lattice/src/Makefile (uses grid-config). Self-contained: uses the SAME core Grid
// primitives as measure_decay_constant (WilsonFermionD, ConjugateGradient, Gamma, sliceSum), so it is
// robust across Grid versions. The emitted "V t C_V(t)" / "A t C_A(t)" rows feed run/12_s_parameter.sh
// -> sparam_raw.dat (mirroring decay_raw.dat).
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
    std::cout << "usage: measure_s_parameter --grid Lx.Ly.Lz.Lt --config <NERSC cfg> --mass <m> "
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
  // CONVENTION/HONESTY reminder in the raw stream: connected isovector V & A; currents UNRENORMALISED.
  std::cout << "# isovector CONNECTED V (gamma_i, rho) and A (gamma_i*gamma5, a1), summed over i=1,2,3"
            << std::endl;
  std::cout << "# S-proxy S = 4 pi (F_V^2/M_V^2 - F_A^2/M_A^2); UNRENORMALISED (Z_V, Z_A owed); "
               "fundamental rep EXPECTED QCD-like (~0.25), NOT walking -- does NOT validate S<0.1 (L5)"
            << std::endl;

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

  // ---- the connected isovector V and A correlators (the only new physics vs measure_decay_constant)
  // For O_Gamma = qbar Gamma q the connected zero-momentum correlator is
  //   C_Gamma(t) = sum_x Re Tr[ Gamma S(x) Gammabar g5 S(x)^dag g5 ],  Gammabar = g5 Gamma^dag g5.
  // We build "propbar" = g5 S^dag g5 once (gamma5-hermiticity image of the backward leg), then for each
  // spatial polarisation i sum the trace with Gamma = gamma_i (vector) and Gamma = gamma_i*gamma5 (axial).
  Gamma g5(Gamma::Algebra::Gamma5);
  Gamma gx(Gamma::Algebra::GammaX), gy(Gamma::Algebra::GammaY), gz(Gamma::Algebra::GammaZ);
  Gamma gi[3] = {gx, gy, gz};

  LatticePropagatorD propbar(UGrid);
  propbar = g5 * adj(prop) * g5;                          // = S(0;x), the gamma5-hermitian backward leg

  LatticeComplexD cv(UGrid), ca(UGrid);
  cv = Zero();
  ca = Zero();
  for (int i = 0; i < 3; ++i) {
    Gamma G = gi[i];                                       // vector insertion gamma_i
    Gamma Gbar = g5 * G * g5;                              // Gammabar = g5 Gamma^dag g5 (g5 g5 = 1 here)
    cv = cv + trace(G * prop * Gbar * propbar);
    // axial insertion gamma_i*gamma5 (a1 channel); Gammabar applied the same way
    Gamma GA = G * g5;
    Gamma GAbar = g5 * GA * g5;
    ca = ca + trace(GA * prop * GAbar * propbar);
  }
  std::vector<TComplexD> cv_t, ca_t;
  sliceSum(cv, cv_t, Tdir);
  sliceSum(ca, ca_t, Tdir);

  // ---- emit the V and A correlators (real parts): feed sparam_raw.dat / s_parameter_proxy ----
  std::cout << "# V  t  C_V(t)   (vector, rho channel)" << std::endl;
  for (size_t t = 0; t < cv_t.size(); ++t) {
    ComplexD c = TensorRemove(cv_t[t]);
    std::cout << "V " << t << " " << std::setprecision(12) << real(c) << std::endl;
  }
  std::cout << "# A  t  C_A(t)   (axial-vector, a1 channel)" << std::endl;
  for (size_t t = 0; t < ca_t.size(); ++t) {
    ComplexD c = TensorRemove(ca_t[t]);
    std::cout << "A " << t << " " << std::setprecision(12) << real(c) << std::endl;
  }

  Grid_finalize();
  return 0;
}
