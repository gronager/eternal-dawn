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
                 "--smear-iters 0,20,50 [--smear-alpha 0.25] [--cg-tol 1e-8]" << std::endl;
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
  std::vector<int> iters{0, 20, 50};                          // smearing-step counts = the operators
  if (GridCmdOptionExists(argv, argv + argc, "--smear-iters")) {
    iters.clear();
    std::stringstream ss(GridCmdOptionPayload(argv, argv + argc, "--smear-iters"));
    std::string tok;
    while (std::getline(ss, tok, ',')) iters.push_back(std::stoi(tok));
  }
  const int N = iters.size();

  LatticeGaugeFieldD Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, cfg);
  RealD plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << "# config " << cfg << "  plaquette " << std::setprecision(10) << plaq
            << "  mass " << mass << "  smear_alpha " << smear_alpha << std::endl;

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

  // for each SOURCE smearing i: smear the source, solve the propagator
  std::vector<LatticePropagatorD> prop_src;
  for (int i = 0; i < N; ++i) {
    LatticePropagatorD src = src0;
    gauss_smear(src, Umu, smear_alpha, iters[i]);
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
    std::cout << "# solved source smearing " << i << " iters " << iters[i] << std::endl;
  }

  // C_ij(t): source i, SINK smearing j (smear the propagator's sink, then contract)
  std::cout << "# i  j  t  Re(C_ij)" << std::endl;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      LatticePropagatorD prop = prop_src[i];
      gauss_smear(prop, Umu, smear_alpha, iters[j]);
      std::vector<ComplexD> C = nucleon_corr(prop, Tdir);
      for (size_t t = 0; t < C.size(); ++t)
        std::cout << i << " " << j << " " << t << " " << std::setprecision(12) << real(C[t])
                  << std::endl;
    }
  }

  Grid_finalize();
  return 0;
}
