// measure_condensate_3pt.cc -- the GENUINE in-medium condensate profile inside the torsiton:
// the connected nucleon scalar 3-point <N| qbar q(r) | N>, the definitive measurement behind the
// s_T lever (Eternal Dawn, Ch. 11 Generations; the refinement of measure_bag_profile's one-body bag).
//
// WHAT & WHY. measure_bag_profile gives the ONE-body dressed-quark bag rho(r)=Tr[S^dag S]; it rises
// into the lever's window toward chiral but is a proxy. The genuine object is the three-body, in-
// medium chiral condensate seen BY the torsiton: how much the scalar density qbar q is depleted
// (chiral symmetry restored) inside the bag. That is the connected 3-point
//     G3(r,tau) = <N(t_snk) [qbar q](x,tau) Nbar(0)>_conn ,   r=|x|, 0 < tau < t_snk,
// whose plateau in tau, divided by the 2-point C_N(t_snk), is the spatial condensate profile. Its
// width is s_T directly (no one-body dictionary), and its shape settles electron-vs-tau (which term
// of the configurational mass dominates).
//
// METHOD (fixed-sink sequential source). The sequential source is the DERIVATIVE of measure_baryon's
// (validated) nucleon contraction w.r.t. the spectator-quark propagator S^{cc'} -- it opens that line
// for the insertion while the diquark (a,b) lines stay forward. From
//   C_N = sum eps eps' s [ tr(P+ S^cc') tr(Dab) + tr(P+ S^cc' Dex) ],  Dab=Cg5 S^aa' Cg5 (S^bb')^T,
// the derivative d C_N / d S^{cc'}_{rho sigma} is
//   sigma_seq^{cc'}_{rho sigma} = (P+)_{sigma rho} tr(Dab) + (Dex P+)_{sigma rho}
//                               = transpose( P+ * tr(Dab) + Dex * P+ )_{rho sigma}.
// SELF-CHECK (the correctness gate): C_N is LINEAR in S^cc', so sum_{x_snk} sigma_seq . S(x_snk;0)
// MUST reproduce C_N(t_snk) (Euler). We emit that reconstruction (CHK) next to the standard C_N --
// if they match, the sequential SOURCE is right. The sequential SOLVE then uses gamma5-hermiticity:
//   eta = g5 sigma_seq^dag g5 at t_snk,  M Sigma = eta,  G3(x) = Tr[ g5 Sigma(x)^dag g5 * S(x;0) ].
// (Derivation: G3 = sum_xsnk Tr[ sigma_seq(xsnk) g5 S(x;xsnk)^dag g5 S(x;0) ] = Tr[g5 Sigma^dag g5 S].)
//
//   measure_condensate_3pt --grid Lx.Ly.Lz.Lt --config <cfg> --mass <m> [--sink-time t] [--cg-tol 1e-8]
//
// This is the SPECTATOR-line connected scalar density (the diquark-line insertions share the same
// spatial SHAPE by degeneracy; the flavour-summed normalisation is an O(1) factor). Self-contained
// core-Grid primitives, mirroring measure_baryon (contraction) and measure_bag_profile (radius bins).
// VALIDATE FIRST: confirm CHK == C_N(t_snk) before trusting the profile; a mismatch flags a one-line
// sign/transpose fix in sigma_seq.
#include <Grid/Grid.h>
#include <cmath>
#include <iomanip>

using namespace Grid;

// gauge-covariant Gaussian (Wuppertal) smearing on the SPATIAL directions, NORMALISED form (a convex
// average, stable for any alpha) -- ported from measure_baryon_gevp. Smearing the SOURCE widens the
// ground-state overlap so the 3-pt plateaus at smaller t_snk (the point source needs ~8 slices to
// reach the ground state; smearing shortens that). It leaves the self-check intact: sigma_seq and the
// forward prop both use the SAME smeared-source propagator, so sigma_seq . S = C_N still holds.
static void gauss_smear(LatticePropagatorD &prop, const LatticeGaugeFieldD &U, RealD alpha, int Niter) {
  if (Niter <= 0) return;
  RealD norm = 1.0 / (1.0 + 2.0 * (Nd - 1) * alpha);
  for (int n = 0; n < Niter; ++n) {
    LatticePropagatorD nbr(prop.Grid());
    nbr = Zero();
    for (int mu = 0; mu < Nd - 1; ++mu) {                  // spatial only
      LatticeColourMatrixD Umu = PeekIndex<LorentzIndex>(U, mu);
      nbr += Umu * Cshift(prop, mu, +1) + adj(Cshift(Umu, mu, -1)) * Cshift(prop, mu, -1);
    }
    prop = norm * (prop + alpha * nbr);
  }
}

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  typedef WilsonImplD FImpl;
  Coordinate latt = GridDefaultLatt();
  const int Tdir = Nd - 1;
  const int Lt = latt[Tdir];
  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(
      latt, GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  if (!GridCmdOptionExists(argv, argv + argc, "--config") ||
      !GridCmdOptionExists(argv, argv + argc, "--mass")) {
    std::cout << "usage: measure_condensate_3pt --grid L.L.L.T --config <cfg> --mass <m> "
                 "[--sink-time t] [--cg-tol 1e-8] [--smear-iters 0] [--smear-alpha 0.25]" << std::endl;
    Grid_finalize();
    return 1;
  }
  std::string cfg = GridCmdOptionPayload(argv, argv + argc, "--config");
  RealD mass = std::stod(GridCmdOptionPayload(argv, argv + argc, "--mass"));
  int t_snk = Lt / 2;
  if (GridCmdOptionExists(argv, argv + argc, "--sink-time")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--sink-time");
    GridCmdOptionInt(s, t_snk);
  }
  RealD cg_tol = 1.0e-8;
  if (GridCmdOptionExists(argv, argv + argc, "--cg-tol"))
    cg_tol = std::stod(GridCmdOptionPayload(argv, argv + argc, "--cg-tol"));
  int smear_iters = 0;                                  // Wuppertal source-smearing steps (0 = point)
  RealD smear_alpha = 0.25;                             // smearing weight; radius ~ sqrt(2*iters*alpha)
  if (GridCmdOptionExists(argv, argv + argc, "--smear-iters")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--smear-iters");
    GridCmdOptionInt(s, smear_iters);
  }
  if (GridCmdOptionExists(argv, argv + argc, "--smear-alpha"))
    smear_alpha = std::stod(GridCmdOptionPayload(argv, argv + argc, "--smear-alpha"));

  LatticeGaugeFieldD Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, cfg);
  RealD plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << "# config " << cfg << "  plaquette " << std::setprecision(10) << plaq
            << "  mass " << mass << "  sink-time " << t_snk << "  smear-iters " << smear_iters
            << "  smear-alpha " << smear_alpha << std::endl;

  WilsonFermionD Dw(Umu, *UGrid, *UrbGrid, mass);
  MdagMLinearOperator<WilsonFermionD, LatticeFermionD> HermOp(Dw);
  ConjugateGradient<LatticeFermionD> CG(cg_tol, 30000);

  // ---- forward propagator S(x;0) from a (smeared) source at the origin (as measure_baryon) ----
  LatticePropagatorD src(UGrid);
  src = Zero();
  {
    typename LatticePropagatorD::vector_object::scalar_object Sid = Zero();
    for (int s = 0; s < Ns; ++s)
      for (int c = 0; c < Nc; ++c) Sid()(s, s)(c, c) = ComplexD(1.0, 0.0);
    Coordinate origin({0, 0, 0, 0});
    pokeSite(Sid, src, origin);
  }
  gauss_smear(src, Umu, smear_alpha, smear_iters);       // widen the source -> ground state at small t
  LatticePropagatorD prop(UGrid);
  prop = Zero();
  for (int s = 0; s < Ns; ++s)
    for (int c = 0; c < Nc; ++c) {
      LatticeFermionD fsrc(UGrid), fmd(UGrid), fsol(UGrid);
      PropToFerm<FImpl>(fsrc, src, s, c);
      Dw.Mdag(fsrc, fmd);
      fsol = Zero();
      CG(HermOp, fmd, fsol);
      FermToProp<FImpl>(prop, fsol, s, c);
    }

  Gamma Cg5 = Gamma(Gamma::Algebra::GammaY) * Gamma(Gamma::Algebra::GammaT) *
              Gamma(Gamma::Algebra::Gamma5);
  Gamma gT(Gamma::Algebra::GammaT);
  Gamma g5(Gamma::Algebra::Gamma5);
  const int eps[6][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}};
  const double sgn[6] = {+1, +1, +1, -1, -1, -1};

  // identity spin matrix (for P+ = (1 + gamma_t)/2 as a SpinMatrix), built like the point source
  LatticeSpinMatrixD Ihat(UGrid);
  {
    typename LatticeSpinMatrixD::vector_object::scalar_object id = Zero();
    for (int s = 0; s < Ns; ++s) id()(s, s)() = ComplexD(1.0, 0.0);
    Ihat = id;
  }

  // ---- (A) the standard nucleon 2-point C_N(t) (validated contraction, for the self-check) ----
  LatticeComplexD nucl(UGrid);
  nucl = Zero();
  // ---- (B) the sequential source sigma_seq^{cc'} = d C_N / d S^{cc'} (open the spectator c-line) ----
  // grid-initialised accumulators (Grid Lattice has no default ctor, so seed the vector from UGrid)
  std::vector<LatticeSpinMatrixD> acc(Nc * Nc, LatticeSpinMatrixD(UGrid));
  for (auto &A : acc) A = Zero();

  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) {
      int a = eps[i][0], b = eps[i][1], c = eps[i][2];
      int ap = eps[j][0], bp = eps[j][1], cp = eps[j][2];
      double s_eps = sgn[i] * sgn[j];
      LatticeSpinMatrixD Saa = peekColour(prop, a, ap);
      LatticeSpinMatrixD Sbb = peekColour(prop, b, bp);
      LatticeSpinMatrixD Scc = peekColour(prop, c, cp);
      LatticeSpinMatrixD Pscc = 0.5 * (Scc + gT * Scc);
      LatticeSpinMatrixD Dab = Cg5 * Saa * Cg5 * transpose(Sbb);
      LatticeSpinMatrixD Dex = Cg5 * transpose(Sbb) * Cg5 * Saa;
      nucl = nucl + s_eps * (trace(Pscc) * trace(Dab) + trace(Pscc * Dex));
      // sigma_seq^{cc'}_{rho sigma} = transpose( P+ * tr(Dab) + Dex * P+ ),  P+ = 0.5(I + gamma_t)
      LatticeComplexD trDab = trace(Dab);
      LatticeSpinMatrixD Pp_trDab = 0.5 * (Ihat * trDab + (gT * Ihat) * trDab);
      LatticeSpinMatrixD DexPp = 0.5 * (Dex + Dex * gT);
      acc[c * Nc + cp] = acc[c * Nc + cp] + s_eps * transpose(Pp_trDab + DexPp);
    }
  std::vector<TComplexD> nucl_t;
  sliceSum(nucl, nucl_t, Tdir);

  // assemble sigma_seq as a propagator and restrict it to the sink time slice
  LatticePropagatorD seqsrc(UGrid);
  seqsrc = Zero();
  for (int c = 0; c < Nc; ++c)
    for (int cp = 0; cp < Nc; ++cp) pokeColour(seqsrc, acc[c * Nc + cp], c, cp);
  LatticeInteger tcoor(UGrid);
  LatticeCoordinate(tcoor, Tdir);
  LatticePropagatorD zeroP(UGrid);
  zeroP = Zero();
  seqsrc = where(tcoor == Integer(t_snk), seqsrc, zeroP);

  // ---- self-check: sum_x Tr[ sigma_seq(x) . transpose(S(x;0)) ] at t_snk must equal C_N(t_snk) ----
  LatticeComplexD chk(UGrid);
  chk = trace(seqsrc * transpose(prop));
  std::vector<TComplexD> chk_t;
  sliceSum(chk, chk_t, Tdir);

  // ---- sequential solve: eta = g5 sigma_seq^dag g5 ; Sigma = M^{-1} eta ----
  LatticePropagatorD eta(UGrid);
  eta = g5 * adj(seqsrc) * g5;
  LatticePropagatorD Sig(UGrid);
  Sig = Zero();
  for (int s = 0; s < Ns; ++s)
    for (int c = 0; c < Nc; ++c) {
      LatticeFermionD fsrc(UGrid), fmd(UGrid), fsol(UGrid);
      PropToFerm<FImpl>(fsrc, eta, s, c);
      Dw.Mdag(fsrc, fmd);
      fsol = Zero();
      CG(HermOp, fmd, fsol);
      FermToProp<FImpl>(Sig, fsol, s, c);
    }

  // ---- the connected scalar 3-point density: G3(x) = Tr[ g5 Sigma(x)^dag g5 * S(x;0) ] ----
  LatticePropagatorD Sigbar(UGrid);
  Sigbar = g5 * adj(Sig) * g5;
  LatticeComplexD g3(UGrid);
  g3 = trace(Sigbar * prop);

  // ---- spatial radius (nearest periodic image), as measure_bag_profile ----
  LatticeInteger rsq(UGrid);
  rsq = Zero();
  int rsq_max = 0;
  for (int mu = 0; mu < Nd - 1; ++mu) {
    Integer L = latt[mu];
    LatticeInteger coor(UGrid);
    LatticeCoordinate(coor, mu);
    LatticeInteger Lfull(UGrid);
    Lfull = L;
    LatticeInteger two_x(UGrid);
    two_x = coor + coor;
    LatticeInteger dx = where(two_x > Lfull, Lfull - coor, coor);
    rsq = rsq + dx * dx;
    rsq_max += (L / 2) * (L / 2);
  }
  LatticeComplexD ones(UGrid);
  ones = ComplexD(1.0, 0.0);
  LatticeComplexD zeroC(UGrid);
  zeroC = Zero();

  // ---- emit: the self-check vs 2-point, the sum-rule SR(tau)=sum_x G3, and the profile G3(r,tau) ----
  std::cout << "# CHK  t_snk  Re[recon_2pt]  Re[C_N(t_snk)]   (these two MUST agree)" << std::endl;
  std::cout << "CHK " << t_snk << " " << std::setprecision(12)
            << real(TensorRemove(chk_t[t_snk])) << " " << real(TensorRemove(nucl_t[t_snk]))
            << std::endl;

  std::cout << "# C2  t  C_N(t)" << std::endl;
  for (int t = 0; t < Lt; ++t)
    std::cout << "C2 " << t << " " << std::setprecision(12) << real(TensorRemove(nucl_t[t]))
              << std::endl;

  std::vector<TComplexD> sr_t;
  sliceSum(g3, sr_t, Tdir);
  std::cout << "# SR  tau  sum_x G3(tau)   (/ C_N(t_snk) -> scalar charge; flat in the plateau)"
            << std::endl;
  for (int t = 0; t < Lt; ++t)
    std::cout << "SR " << t << " " << std::setprecision(12) << real(TensorRemove(sr_t[t]))
              << std::endl;

  std::cout << "# P3  r2  tau  G3_sum  count" << std::endl;
  for (int bb = 0; bb <= rsq_max; ++bb) {
    LatticeComplexD shell = where(rsq == Integer(bb), g3, zeroC);
    LatticeComplexD cnt = where(rsq == Integer(bb), ones, zeroC);
    std::vector<TComplexD> sh, sc;
    sliceSum(shell, sh, Tdir);
    sliceSum(cnt, sc, Tdir);
    double tot = 0.0;
    for (size_t t = 0; t < sc.size(); ++t) tot += real(TensorRemove(sc[t]));
    if (tot < 0.5) continue;
    for (size_t t = 0; t < sh.size(); ++t)
      std::cout << "P3 " << bb << " " << t << " " << std::setprecision(12)
                << real(TensorRemove(sh[t])) << " "
                << (long)std::llround(real(TensorRemove(sc[t]))) << std::endl;
  }

  Grid_finalize();
  return 0;
}
