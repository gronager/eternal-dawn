// measure_baryon_gevp.cc -- the EXCITED torsiton spectrum: a variational (multi-operator) nucleon
// correlator matrix for the GEVP (Eternal Dawn L4, excited rungs = candidate generations).
//
// THE TEST. The framework's thesis: the torsiton tower has exactly THREE radial rungs (the three
// generations) -- the ground state, and two radial excitations whose densities cross zero once and
// twice -- and the condensate supports the third but NOT a fourth. A single operator (measure_baryon)
// only sees the ground state cleanly. To resolve the EXCITED states you need several operators with
// different spatial overlaps and the generalised eigenvalue problem (GEVP): the eigenvalues of
// C(t) v = lambda(t) C(t0) v decay as e^{-E_n t}, so an N-operator basis resolves N states. If a
// clean third nucleon state appears and a fourth does NOT, that supports the three-generation thesis;
// if the tower keeps going, it does not. Either way the lattice can SAY something.
//
// The operators here are the nucleon interpolator (same Cg5 diquark + parity projector as
// measure_baryon) built from quark fields GAUSSIAN-SMEARED at different radii -- a cheap, standard
// variational basis (wider smearing overlaps more excited radial structure). Output: the symmetric
// N x N matrix C_ij(t) = <N_i(t) Nbar_j(0)>, source-smeared i, sink-smeared j.
//
//   measure_baryon_gevp --grid Lx.Ly.Lz.Lt --config <cfg> --mass <m> \
//                       --smear-iters 0,20,50 [--smear-alpha 0.25] [--cg-tol 1e-8]
// (smear-iters = a list of Wuppertal step counts, one operator each; 0 = point. radius ~
//  sqrt(2*iters*alpha). NORMALISED smearing -- stable for any alpha, unlike the un-normalised form.)
// Self-contained core-Grid primitives (covariant Laplacian via Cshift + links, as the APE smearing
// in measure_potential). The nucleon contraction is measure_baryon's, validated there.
//
// -------------------------------------------------------------------------------------------------
// NODAL BASIS (the radial / "generation tower"):  --basis nodal  (or --nodal)
// -------------------------------------------------------------------------------------------------
// The Gaussian basis above is built from Wuppertal smearings at a few iteration counts. Those
// sources are all NODELESS (each is ~ exp(alpha*Delta) eta, a positive bump that just gets wider);
// they are nearly collinear and overlap the radial EXCITATIONS very weakly, so the GEVP cannot
// resolve the excited rungs. The radial excitations of a baryon (analogues of the hydrogenic
// 1s/2s/3s wavefunctions) carry 0, 1, 2 radial NODES -- their densities cross zero. A basis that
// is to see them must itself contain node-bearing operators.
//
// The lightweight, standard way to inject nodes WITHOUT full distillation / Laplacian-Heaviside
// (Peardon et al.) is covariant-Laplacian POLYNOMIAL ("nodal" / Laguerre-filtered) smearing: take
// the gauge-covariant 3D Laplacian
//     Delta phi(x) = sum_{i spatial} [ U_i(x) phi(x+i) + U_i^dag(x-i) phi(x-i) - 2 phi(x) ]
// (the SAME covariant kernel the Wuppertal step uses -- gauss_smear is one normalised Jacobi step
// of (1 + alpha*Delta_norm)), and build the basis as low-order polynomials of Delta acting on a
// common Gaussian-smeared point source g = Gaussian(eta):
//     phi_0 = g                                 (0 nodes  -- nodeless ground-state-like, "1s")
//     phi_1 = (1 - c1*Delta) g                  (1 node   -- "2s"-like)
//     phi_2 = (1 - c2*Delta + c2p*Delta^2) g    (2 nodes  -- "3s"-like)
// In the continuum a Gaussian g ~ exp(-r^2/2w^2) has Delta g ~ (r^2/w^4 - 3/w^2) g, so each factor
// (1 - c*Delta) flips the sign of the source at a radius ~ sqrt(1/c)*w: a polynomial of degree k in
// Delta acting on a Gaussian is a degree-2k polynomial in r times the Gaussian, i.e. a Laguerre-
// like radial profile L_k(r^2/w^2) exp(-r^2/2w^2) that crosses zero k times -- k radial nodes.
// The c_k below are fixed (from <Delta>, <Delta^2> on the Gaussian, see set_nodal_coeffs) so the
// nodes sit roughly at the Gaussian's first and second standard radii; this is approximate, not a
// diagonalised radial basis.
//
// HONEST CAVEAT (keep the project's honesty bar): this is POLYNOMIAL-FILTERED nodal smearing -- a
// cheap stand-in for the rigorous radial basis (full distillation / LapH, which builds the operators
// from the low Laplacian eigenvectors and gets the node POSITIONS variationally). Here the node
// positions are APPROXIMATE (set by the polynomial coefficients, not optimised), so the excited
// masses carry a basis-systematic on top of the statistical error. It is enough to give the GEVP
// genuine overlap with the 1-node and 2-node sectors -- to SEE whether a third (and not a fourth)
// rung resolves -- but resolving 3 vs 4 generations cleanly remains statistics- and basis-limited.
//
//   measure_baryon_gevp --grid Lx.Ly.Lz.Lt --config <cfg> --mass <m> --basis nodal \
//                       [--nodal-nops 3] [--nodal-iters 30] [--smear-alpha 0.25] [--cg-tol 1e-8]
// (nodal-nops = number of operators = max nodes + 1, default 3 -> {0,1,2} nodes; nodal-iters = the
//  Gaussian width shared by all nodal operators, as a Wuppertal iteration count.)
// Output format is IDENTICAL to the Gaussian path ("i j t Re(C_ij)"), so gevp_spectrum reads it
// unchanged.
#include <Grid/Grid.h>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace Grid;

// gauge-covariant Gaussian (Wuppertal) smearing on the SPATIAL directions, NORMALISED form -- a
// convex average, stable for any alpha>0 (the un-normalised (1+a Lap) blows up once a>1/6):
//   prop <- (prop + alpha * sum_i[U_i(x) P(x+i) + U_i^dag(x-i) P(x-i)]) / (1 + 2*Nspatial*alpha)
// The smearing radius grows as ~sqrt(2*Niter*alpha); Niter is the operator knob (0 = point).
static void gauss_smear(LatticePropagatorD &prop, const LatticeGaugeFieldD &U, RealD alpha, int Niter) {
  if (Niter <= 0) return;
  RealD norm = 1.0 / (1.0 + 2.0 * (Nd - 1) * alpha);
  for (int n = 0; n < Niter; ++n) {
    LatticePropagatorD nbr(prop.Grid());
    nbr = Zero();
    for (int mu = 0; mu < Nd - 1; ++mu) {                      // spatial only
      LatticeColourMatrixD Umu = PeekIndex<LorentzIndex>(U, mu);
      nbr += Umu * Cshift(prop, mu, +1) + adj(Cshift(Umu, mu, -1)) * Cshift(prop, mu, -1);
    }
    prop = norm * (prop + alpha * nbr);
  }
}

// the gauge-covariant 3D (spatial) Laplacian, applied to a propagator:
//   Delta P(x) = sum_{i spatial} [ U_i(x) P(x+i) + U_i^dag(x-i) P(x-i) - 2 P(x) ]
// the SAME covariant kernel the Wuppertal step uses (gauss_smear's neighbour sum is the off-diagonal
// part); here we keep the full Laplacian (with the -2*Nspatial diagonal) so polynomials in Delta have
// the standard node structure. Used to build the nodal (Laguerre-filtered) basis below.
static LatticePropagatorD covariant_laplacian(const LatticePropagatorD &prop,
                                              const LatticeGaugeFieldD &U) {
  LatticePropagatorD out(prop.Grid());
  out = Zero();
  for (int mu = 0; mu < Nd - 1; ++mu) {                       // spatial only
    LatticeColourMatrixD Umu = PeekIndex<LorentzIndex>(U, mu);
    out += Umu * Cshift(prop, mu, +1) + adj(Cshift(Umu, mu, -1)) * Cshift(prop, mu, -1)
           - 2.0 * prop;
  }
  return out;
}

// Coefficients c_k for the degree-k nodal polynomial P_k(Delta) so that P_k(Delta) g carries k
// radial nodes, with g a Gaussian-smeared source. On a Gaussian g ~ exp(-r^2/2w^2) the covariant
// Laplacian acts (in the free/continuum limit) as Delta g ~ (r^2/w^4 - Nspatial/w^2) g, i.e. Delta
// has a characteristic eigenscale lam ~ 1/w^2 on g. Choosing c1 ~ w^2 makes (1 - c1*Delta) g change
// sign once (one node near r ~ w); the degree-2 polynomial (1 - c2 Delta + c2p Delta^2) g is tuned
// to cross zero twice. We estimate the scale w^2 from the Wuppertal iteration count and alpha via the
// smearing-radius relation r^2 ~ 2*Niter*alpha, so w^2 ~ Niter*alpha (the Gaussian variance). The
// resulting node positions are APPROXIMATE (this is filtered, not diagonalised, smearing -- see the
// header caveat); they only need to give the GEVP non-trivial overlap with the 1- and 2-node sectors.
struct NodalCoeffs { RealD c1, c2, c2p; };
static NodalCoeffs set_nodal_coeffs(int gauss_iters, RealD alpha) {
  RealD w2 = std::max(1.0, (RealD)gauss_iters * alpha);       // Gaussian variance ~ Niter*alpha
  NodalCoeffs k;
  // 1 node: zero of (1 - c1*lam) at lam ~ 1/w2  ->  c1 ~ w2 (node near r ~ w).
  k.c1 = w2;
  // 2 nodes: a quadratic in Delta with two positive roots straddling 1/w2, placed at lam1 ~ 0.6/w2
  // and lam2 ~ 1.8/w2; (1 - c2 lam + c2p lam^2) = c2p (lam - lam1)(lam - lam2), normalised to 1 at
  // lam=0:  c2 = (lam1+lam2)/(lam1 lam2),  c2p = 1/(lam1 lam2).
  RealD lam1 = 0.6 / w2, lam2 = 1.8 / w2;
  k.c2 = (lam1 + lam2) / (lam1 * lam2);
  k.c2p = 1.0 / (lam1 * lam2);
  return k;
}

// build the k-node nodal source phi_k = P_k(Delta) g from a Gaussian-smeared source g:
//   k=0: g                              k=1: (1 - c1 Delta) g
//   k=2: (1 - c2 Delta + c2p Delta^2) g
// Delta is the covariant 3D Laplacian above. Reuses the codebase's covariant kernel; only the
// polynomial assembly is new. Normalised to unit Frobenius norm so the GEVP matrix is well-scaled.
static LatticePropagatorD nodal_source(const LatticePropagatorD &g, const LatticeGaugeFieldD &U,
                                       int knode, const NodalCoeffs &k) {
  LatticePropagatorD phi = g;
  if (knode >= 1) {
    LatticePropagatorD Lg = covariant_laplacian(g, U);
    phi = g - k.c1 * Lg;
    if (knode >= 2) {
      LatticePropagatorD L2g = covariant_laplacian(Lg, U);     // Delta^2 g
      phi = g - k.c2 * Lg + k.c2p * L2g;
    }
  }
  RealD nrm = std::sqrt(norm2(phi));
  if (nrm > 0.0) phi = phi * (1.0 / nrm);
  return phi;
}

// the nucleon (torsiton) zero-momentum correlator C(t) for a given propagator (measure_baryon's
// two-Wick contraction, Cg5 diquark, positive-parity projector).
static std::vector<ComplexD> nucleon_corr(const LatticePropagatorD &prop, int Tdir) {
  Gamma Cg5 = Gamma(Gamma::Algebra::GammaY) * Gamma(Gamma::Algebra::GammaT) *
              Gamma(Gamma::Algebra::Gamma5);
  Gamma gT(Gamma::Algebra::GammaT);
  const int eps[6][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}};
  const double sgn[6] = {+1, +1, +1, -1, -1, -1};
  LatticeComplexD nucl(prop.Grid());
  nucl = Zero();
  for (int i = 0; i < 6; ++i) {
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
    }
  }
  std::vector<TComplexD> slice;
  sliceSum(nucl, slice, Tdir);
  std::vector<ComplexD> out(slice.size());
  for (size_t t = 0; t < slice.size(); ++t) out[t] = TensorRemove(slice[t]);
  return out;
}

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
    std::cout << "usage: measure_baryon_gevp --grid L.L.L.T --config <cfg> --mass <m> "
                 "[--basis gauss|nodal] --smear-iters 0,20,50 [--smear-alpha 0.25] [--cg-tol 1e-8]\n"
                 "       nodal basis: [--basis nodal | --nodal] [--nodal-nops 3] [--nodal-iters 30]"
              << std::endl;
    Grid_finalize();
    return 1;
  }
  std::string cfg = GridCmdOptionPayload(argv, argv + argc, "--config");
  RealD mass = std::stod(GridCmdOptionPayload(argv, argv + argc, "--mass"));
  RealD cg_tol = 1.0e-8;
  RealD smear_alpha = 0.25;
  if (GridCmdOptionExists(argv, argv + argc, "--cg-tol"))
    cg_tol = std::stod(GridCmdOptionPayload(argv, argv + argc, "--cg-tol"));
  if (GridCmdOptionExists(argv, argv + argc, "--smear-alpha"))
    smear_alpha = std::stod(GridCmdOptionPayload(argv, argv + argc, "--smear-alpha"));

  // basis selector: "gauss" (default, Wuppertal-iteration operators) or "nodal" (covariant-Laplacian
  // polynomial-filtered operators carrying 0,1,2,... radial nodes -- the radial / generation tower).
  bool nodal = false;
  if (GridCmdOptionExists(argv, argv + argc, "--nodal")) nodal = true;
  if (GridCmdOptionExists(argv, argv + argc, "--basis"))
    nodal = (GridCmdOptionPayload(argv, argv + argc, "--basis") == "nodal");

  std::vector<int> iters{0, 20, 50};                          // smearing-step counts = the operators
  if (GridCmdOptionExists(argv, argv + argc, "--smear-iters")) {
    iters.clear();
    std::stringstream ss(GridCmdOptionPayload(argv, argv + argc, "--smear-iters"));
    std::string tok;
    while (std::getline(ss, tok, ',')) iters.push_back(std::stoi(tok));
  }

  // nodal-basis knobs: number of operators (= max nodes + 1) and the shared Gaussian width (as a
  // Wuppertal iteration count). All nodal operators share one Gaussian g; they differ only in the
  // node-injecting polynomial P_k(Delta). Default 3 operators -> {0,1,2} nodes.
  int nodal_nops = 3;
  int nodal_iters = 30;
  if (GridCmdOptionExists(argv, argv + argc, "--nodal-nops"))
    nodal_nops = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--nodal-nops"));
  if (GridCmdOptionExists(argv, argv + argc, "--nodal-iters"))
    nodal_iters = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--nodal-iters"));

  const int N = nodal ? nodal_nops : (int)iters.size();
  NodalCoeffs ncoeff = set_nodal_coeffs(nodal_iters, smear_alpha);

  LatticeGaugeFieldD Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, cfg);
  RealD plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << "# config " << cfg << "  plaquette " << std::setprecision(10) << plaq
            << "  mass " << mass << "  smear_alpha " << smear_alpha << std::endl;
  if (nodal)
    std::cout << "# basis nodal (covariant-Laplacian polynomial; " << N << " operators = nodes 0.."
              << (N - 1) << ", gauss-iters " << nodal_iters << ", coeffs c1=" << ncoeff.c1
              << " c2=" << ncoeff.c2 << " c2p=" << ncoeff.c2p << ")" << std::endl;
  else
    std::cout << "# basis gauss (Wuppertal-iteration operators)" << std::endl;

  WilsonFermionD Dw(Umu, *UGrid, *UrbGrid, mass);
  MdagMLinearOperator<WilsonFermionD, LatticeFermionD> HermOp(Dw);
  ConjugateGradient<LatticeFermionD> CG(cg_tol, 30000);

  // base point source at the origin
  LatticePropagatorD src0(UGrid);
  src0 = Zero();
  {
    typename LatticePropagatorD::vector_object::scalar_object Sid = Zero();
    for (int s = 0; s < Ns; ++s)
      for (int c = 0; c < Nc; ++c) Sid()(s, s)(c, c) = ComplexD(1.0, 0.0);
    Coordinate origin({0, 0, 0, 0});
    pokeSite(Sid, src0, origin);
  }

  // The shared Gaussian source for the nodal basis (smeared once, then node-filtered per operator).
  LatticePropagatorD gauss_base = src0;
  if (nodal) gauss_smear(gauss_base, Umu, smear_alpha, nodal_iters);

  // build the N smeared SOURCES, then solve a propagator from each.
  //   gauss basis:  src_i = Gaussian_{iters[i]}(eta)
  //   nodal basis:  src_i = P_i(Delta) Gaussian_{nodal_iters}(eta)  (i radial nodes)
  std::vector<LatticePropagatorD> prop_src;
  for (int i = 0; i < N; ++i) {
    LatticePropagatorD src(UGrid);
    if (nodal) {
      src = nodal_source(gauss_base, Umu, i, ncoeff);          // i nodes
    } else {
      src = src0;
      gauss_smear(src, Umu, smear_alpha, iters[i]);
    }
    LatticePropagatorD prop(UGrid);
    prop = Zero();
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
    prop_src.push_back(prop);
    if (nodal)
      std::cout << "# solved source nodal-operator " << i << " (" << i << " nodes)" << std::endl;
    else
      std::cout << "# solved source smearing " << i << " iters " << iters[i] << std::endl;
  }

  // C_ij(t): source i, SINK operator j (apply operator j to the propagator's sink, then contract).
  // The sink operator mirrors the source: a Wuppertal smear (gauss basis) or the same node-filter
  // P_j(Delta) on the Gaussian-smeared propagator (nodal basis).
  std::cout << "# i  j  t  Re(C_ij)" << std::endl;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      LatticePropagatorD prop = prop_src[i];
      if (nodal) {
        gauss_smear(prop, Umu, smear_alpha, nodal_iters);      // shared Gaussian sink width
        prop = nodal_source(prop, Umu, j, ncoeff);             // then j-node filter at the sink
      } else {
        gauss_smear(prop, Umu, smear_alpha, iters[j]);
      }
      std::vector<ComplexD> C = nucleon_corr(prop, Tdir);
      for (size_t t = 0; t < C.size(); ++t)
        std::cout << i << " " << j << " " << t << " " << std::setprecision(12) << real(C[t])
                  << std::endl;
    }
  }

  Grid_finalize();
  return 0;
}
