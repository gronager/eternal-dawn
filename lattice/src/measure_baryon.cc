// measure_baryon.cc -- the torsiton on the lattice: quenched valence spectroscopy (Eternal Dawn L4).
//
// The torsiton is the SU(3)-fundamental BARYON (colour from the Pauli label, Ch. "Forces"): three
// quarks, one of each colour, the same object QCD calls the nucleon. This program "rediscovers the
// bound ground state on the lattice" -- it reads ONE pure-gauge NERSC configuration, solves a
// point-source Wilson-quark propagator, and contracts two correlators:
//
//   * the PION  C_pi(t) = sum_x Tr[ S(x)^dag S(x) ]  -- the robust anchor. It always works (one
//     trace, gamma5-hermiticity), it CALIBRATES the quark mass (m_pi tells us how chiral we are),
//     and via GMOR (m_pi^2 ~ m_q <qbar q>) it is the lattice's window on the chiral condensate --
//     the "far weirder than mean field" vacuum that, if anything, carries the mass hierarchy.
//   * the NUCLEON C_N(t) -- the TORSITON. The two-Wick-contraction baryon two-point function with
//     the diquark C*gamma5 and the positive-parity projector. Its plateau in the effective mass
//     m_eff(t) = ln[C(t)/C(t+1)] is the torsiton rest mass (in lattice units). A clean plateau at
//     m_N > 0 IS the bound ground state, found non-perturbatively.
//
// QUENCHED / VALENCE: the propagator is solved on a PURE-GAUGE config (no dynamical-fermion
// determinant), so this is the cheap pilot -- it reuses the L1 ensemble and needs no new HMC. It
// gets the ground-state baryon mass and stands up the pipeline; dynamical fermions near the chiral
// limit (the real spectrum, and any excited rungs = candidate further generations) come after.
//
//   measure_baryon --grid Lx.Ly.Lz.Lt --config <NERSC cfg> --mass <m> [--cg-tol 1e-8] \
//                  [--cg-maxit 30000]
//
// Build with lattice/src/Makefile (uses grid-config). Self-contained: uses core Grid primitives
// (WilsonFermionD, ConjugateGradient, peekColour, Gamma, sliceSum) so it is robust across Grid
// versions. NOTE on the nucleon contraction: the TWO-term structure, the C*gamma5 diquark, the
// colour epsilon, and the parity projector are standard; the relative sign and the spin-transpose
// placement should be validated once against a known nucleon plateau (a quenched m_N/sqrt(sigma) ~
// 1.8-2.0 at these masses) -- a sign slip shows up as a noisy or wrong-sloped correlator, not a crash.
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
    std::cout << "usage: measure_baryon --grid Lx.Ly.Lz.Lt --config <NERSC cfg> --mass <m> "
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

  // ---- load the pure-gauge configuration ----
  LatticeGaugeFieldD Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, cfg);
  RealD plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << "# config " << cfg << "  plaquette " << std::setprecision(10) << plaq
            << "  mass " << mass << std::endl;

  // ---- Wilson quark, point source at the origin (identity in spin x colour) ----
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

  // ---- solve the 12 propagator columns: (Mdag M) x = Mdag b, then M x = b ----
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

  // ---- (1) the PION correlator: C_pi(t) = sum_x Tr[S(x)^dag S(x)] (gamma5-hermitian PP) ----
  LatticeComplexD pion(UGrid);
  pion = trace(adj(prop) * prop);
  std::vector<TComplexD> pion_t;
  sliceSum(pion, pion_t, Tdir);

  // ---- (2) the NUCLEON (torsiton) correlator ----
  // chi = eps^{abc} (q^a C g5 q^b) q^c, degenerate q = S. Two Wick terms; positive-parity projector
  // P_+ = (1 + gamma_t)/2. Colour epsilon by the 6 signed permutations of (0,1,2).
  Gamma Cg5 = Gamma(Gamma::Algebra::GammaY) * Gamma(Gamma::Algebra::GammaT) *
              Gamma(Gamma::Algebra::Gamma5);
  Gamma gT(Gamma::Algebra::GammaT);
  const int eps[6][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}};
  const double sgn[6] = {+1, +1, +1, -1, -1, -1};

  LatticeComplexD nucl(UGrid);
  nucl = Zero();
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      int a = eps[i][0], b = eps[i][1], c = eps[i][2];
      int ap = eps[j][0], bp = eps[j][1], cp = eps[j][2];
      double s_eps = sgn[i] * sgn[j];
      LatticeSpinMatrixD Saa = peekColour(prop, a, ap);
      LatticeSpinMatrixD Sbb = peekColour(prop, b, bp);
      LatticeSpinMatrixD Scc = peekColour(prop, c, cp);
      // P_+ applied to the spectator block:  Pscc = (Scc + gamma_t Scc)/2
      LatticeSpinMatrixD Pscc = 0.5 * (Scc + gT * Scc);
      // the C*g5 diquark of the (a,b) pair, with a spin transpose on the b-leg. Sbb is a bare
      // LatticeSpinMatrix (colour peeled off), so its only matrix index is spin: plain transpose()
      // is the spin transpose (this Grid's transposeSpin is the index-templated variant).
      LatticeSpinMatrixD Dab = Cg5 * Saa * Cg5 * transpose(Sbb);
      // term 1 (direct):   Tr[P_+ S^cc'] * Tr[ Dab ]
      // term 2 (exchange): Tr[ P_+ S^cc' * (Cg5 (S^bb')^T Cg5 S^aa') ]
      LatticeSpinMatrixD Dex = Cg5 * transpose(Sbb) * Cg5 * Saa;
      nucl = nucl + s_eps * (trace(Pscc) * trace(Dab) + trace(Pscc * Dex));
    }
  }
  std::vector<TComplexD> nucl_t;
  sliceSum(nucl, nucl_t, Tdir);

  // ---- emit both correlators (real parts): t  C_pi(t)  C_N(t) ----
  std::cout << "# t  C_pi(t)  C_N(t)" << std::endl;
  for (size_t t = 0; t < pion_t.size(); ++t) {
    ComplexD cp = TensorRemove(pion_t[t]);
    ComplexD cn = TensorRemove(nucl_t[t]);
    std::cout << t << " " << std::setprecision(12) << real(cp) << " " << real(cn) << std::endl;
  }

  Grid_finalize();
  return 0;
}
