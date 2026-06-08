// measure_bag_profile.cc -- the SHAPE of the torsiton's mass-giving core (Eternal Dawn, the s_T lever).
//
// THE MEASUREMENT the generations mechanism owes (Ch. 11, Sec. Generations). The fermion mass
// hierarchy is set by ONE number: the sharpness of the dynamically-massive "bag" the radial rungs
// couple to (the configurational-mass core size s_T). A sharp core spreads the consecutive rungs
// n=1,2,3 across the full charged-lepton span (x3477); a broad core flattens them (x~50). So the
// whole hierarchy reduces to a quantity the lattice can MEASURE: how sharp is the bag.
//
// We measure the gauge-invariant scalar density profile of a dressed quark emitted from a point
// source,
//     rho(x) = Tr_{spin,colour}[ S(x;0)^dag S(x;0) ]            (the pion point-to-point integrand),
// resolved by the spatial distance r=|x| from the source on each Euclidean time slice. rho(r) is the
// spatial spread of the constituent (dynamically massive) quark -- the BAG. Its half-density radius
// R0 and its surface (wall) thickness `a` are the bag's size and sharpness; the sharpness RATIO
// a/R0, in the plateau, is the lattice's read-out of s_T (analysed in cartasis_sims.lattice). The
// spatially-summed profile is exactly the pion correlator C_pi(t)=sum_x rho(x,t), emitted too, so the
// same run calibrates the scale (m_pi a) it is measured in.
//
// rho = Tr[S^dag S] is manifestly gauge invariant (S -> g(x) S g(0)^dag leaves Tr[S^dag S] fixed),
// so NO gauge fixing is needed -- this is the robust, single-solve probe of the bag. It is the ONE-
// body (constituent) profile; the genuine three-body in-medium condensate <N|qbar q(r)|N> is the
// sequential-source 3-point refinement (a labelled next step, not this program).
//
//   measure_bag_profile --grid Lx.Ly.Lz.Lt --config <NERSC cfg> --mass <m> [--cg-tol 1e-8]
//
// Self-contained core-Grid primitives (LatticeCoordinate, where, sliceSum), as measure_baryon.
#include <Grid/Grid.h>
#include <cmath>
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
    std::cout << "usage: measure_bag_profile --grid Lx.Ly.Lz.Lt --config <NERSC cfg> --mass <m> "
                 "[--cg-tol 1e-8]" << std::endl;
    Grid_finalize();
    return 1;
  }
  std::string cfg = GridCmdOptionPayload(argv, argv + argc, "--config");
  RealD mass = std::stod(GridCmdOptionPayload(argv, argv + argc, "--mass"));
  RealD cg_tol = 1.0e-8;
  if (GridCmdOptionExists(argv, argv + argc, "--cg-tol"))
    cg_tol = std::stod(GridCmdOptionPayload(argv, argv + argc, "--cg-tol"));

  // ---- load the configuration, build the Wilson quark ----
  LatticeGaugeFieldD Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, cfg);
  RealD plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu);
  std::cout << "# config " << cfg << "  plaquette " << std::setprecision(10) << plaq
            << "  mass " << mass << "  grid " << latt[0] << "." << latt[1] << "." << latt[2]
            << "." << latt[3] << std::endl;

  // ---- point source at the origin; solve the 12-column propagator (as measure_baryon) ----
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
  LatticePropagatorD prop(UGrid);
  prop = Zero();
  MdagMLinearOperator<WilsonFermionD, LatticeFermionD> HermOp(Dw);
  ConjugateGradient<LatticeFermionD> CG(cg_tol, 30000);
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

  // ---- the gauge-invariant scalar density rho(x) = Tr[S^dag S] (the pion integrand) ----
  LatticeComplexD rho(UGrid);
  rho = trace(adj(prop) * prop);

  // ---- spatial radius squared from the origin, with the nearest periodic image per axis ----
  // nearest image: dx = min(x, L-x), written as 2x>L ? L-x : x to avoid left-scalar / parity issues,
  // and using Lattice-Lattice comparisons (a constant L field) which are the robust Grid idiom.
  LatticeInteger rsq(UGrid);
  rsq = Zero();
  int rsq_max = 0;
  for (int mu = 0; mu < Nd - 1; ++mu) {                  // spatial axes only
    Integer L = latt[mu];
    LatticeInteger coor(UGrid);
    LatticeCoordinate(coor, mu);
    LatticeInteger Lfull(UGrid);
    Lfull = L;                                           // constant field = lattice extent
    LatticeInteger two_x(UGrid);
    two_x = coor + coor;                                 // 2x (avoid scalar*lattice ambiguity)
    LatticeInteger dx = where(two_x > Lfull, Lfull - coor, coor);   // min(x, L-x)
    rsq = rsq + dx * dx;
    rsq_max += (L / 2) * (L / 2);
  }

  // ---- for each r^2 shell: sum rho and count sites on every time slice (one sliceSum per shell) ----
  LatticeComplexD ones(UGrid);
  ones = ComplexD(1.0, 0.0);
  LatticeComplexD zero(UGrid);
  zero = Zero();

  std::cout << "# PROF  r2  t  rho_sum  count   (shell-summed scalar density and site count)"
            << std::endl;
  for (int b = 0; b <= rsq_max; ++b) {
    LatticeComplexD shell_rho = where(rsq == Integer(b), rho, zero);
    LatticeComplexD shell_cnt = where(rsq == Integer(b), ones, zero);
    std::vector<TComplexD> sr, sc;
    sliceSum(shell_rho, sr, Tdir);
    sliceSum(shell_cnt, sc, Tdir);
    // skip empty shells (no lattice site has this r^2): count == 0 on every slice
    double tot = 0.0;
    for (size_t t = 0; t < sc.size(); ++t) tot += real(TensorRemove(sc[t]));
    if (tot < 0.5) continue;
    for (size_t t = 0; t < sr.size(); ++t)
      std::cout << "PROF " << b << " " << t << " " << std::setprecision(12)
                << real(TensorRemove(sr[t])) << " "
                << (long)std::llround(real(TensorRemove(sc[t]))) << std::endl;
  }

  // ---- the spatially-summed profile IS the pion correlator: scale calibration in the same run ----
  std::vector<TComplexD> pion_t;
  sliceSum(rho, pion_t, Tdir);
  std::cout << "# PION  t  C_pi(t)" << std::endl;
  for (size_t t = 0; t < pion_t.size(); ++t)
    std::cout << "PION " << t << " " << std::setprecision(12) << real(TensorRemove(pion_t[t]))
              << std::endl;

  Grid_finalize();
  return 0;
}
