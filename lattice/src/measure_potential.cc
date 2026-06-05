// measure_potential.cc -- static quark potential from timelike Wilson loops (Eternal Dawn L1-L2).
// Reads one NERSC gauge configuration and prints the average timelike Wilson loop W(R,T) for a
// ladder of spatial R and temporal T. The run script loops over configs and averages; the
// Python analyzer (cartasis_sims.lattice) extracts V(R) = ln[W(R,T)/W(R,T+1)] and fits the
// Cornell form V(R)=c-alpha/R+sigma*r for the string tension sigma.
//
// SPATIAL APE SMEARING (--smear N [--smear-alpha A]):
//   Unsmeared Wilson loops have poor overlap with the gluon-string ground state, so the
//   effective potential is contaminated by excited states at the small T where the signal/noise
//   is still usable -- which biases V(R), and the string tension, HIGH. APE-smearing the SPATIAL
//   links (the temporal links are left untouched, so the transfer-matrix/potential interpretation
//   is preserved) boosts the ground-state overlap, so V_eff(R,T) plateaus already at small T.
//   This is the single biggest signal/noise win for the static potential. Default N=0 reproduces
//   the original unsmeared behaviour exactly (the validated path), so this never regresses.
//
//   measure_potential --grid Lx.Ly.Lz.Lt --config <NERSC cfg> [--rmax R] [--tmax T] \
//                     [--smear 20] [--smear-alpha 0.5]
//
// Build with lattice/src/Makefile (uses grid-config). Written against Grid's WilsonLoops API
// (Grid/qcd/utils/WilsonLoops.h): avgPlaquette, avgTimelikeWilsonLoop are static members.
#include <Grid/Grid.h>
#include <iomanip>

using namespace Grid;

// One step of APE smearing on the SPATIAL links only (directions 0..Nd-2; the temporal
// direction Nd-1 is left exactly as-is). For each spatial link U_mu(x) we replace it by the
// projection back onto SU(3) of  (1-alpha) U_mu(x) + (alpha/6) * (sum of the 4 spatial staples).
// Uses only core Grid primitives (PeekIndex/PokeIndex, Cshift, adj, ProjectOnGroup) so it is
// robust across Grid versions.
static void ape_smear_spatial(LatticeGaugeFieldD &U, RealD alpha) {
  GridBase *grid = U.Grid();
  LatticeGaugeFieldD Uin = U;                 // read from the old field, write into U
  const int Tdir = Nd - 1;                    // temporal direction: never smeared
  for (int mu = 0; mu < Nd - 1; ++mu) {       // spatial links only
    LatticeColourMatrixD Umu = PeekIndex<LorentzIndex>(Uin, mu);
    LatticeColourMatrixD staple(grid);
    staple = Zero();
    for (int nu = 0; nu < Nd - 1; ++nu) {     // spatial neighbours only (skip mu and time)
      if (nu == mu) continue;
      LatticeColourMatrixD Unu = PeekIndex<LorentzIndex>(Uin, nu);
      // upper staple: U_nu(x) U_mu(x+nu) U_nu(x+mu)^dag
      staple += Unu * Cshift(Umu, nu, +1) * adj(Cshift(Unu, mu, +1));
      // lower staple: U_nu(x-nu)^dag U_mu(x-nu) U_nu(x-nu+mu)
      LatticeColourMatrixD Unu_dn = Cshift(Unu, nu, -1);
      staple += adj(Unu_dn) * Cshift(Umu, nu, -1) * Cshift(Cshift(Unu, nu, -1), mu, +1);
    }
    LatticeColourMatrixD Vmu = (1.0 - alpha) * Umu + (alpha / 6.0) * staple;
    ProjectOnGroup(Vmu);                      // reunitarise back into SU(3)
    PokeIndex<LorentzIndex>(U, Vmu, mu);
  }
  (void)Tdir;
}

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  typedef PeriodicGimplR Gimpl;                 // standard periodic SU(3) gauge implementation

  Coordinate latt = GridDefaultLatt();
  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(
      latt, GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  Gimpl::GaugeField Umu(UGrid);

  if (!GridCmdOptionExists(argv, argv + argc, "--config")) {
    std::cout << "usage: measure_potential --grid Lx.Ly.Lz.Lt --config <NERSC cfg> "
                 "[--rmax R] [--tmax T] [--smear N] [--smear-alpha A]" << std::endl;
    Grid_finalize();
    return 1;
  }
  std::string cfg = GridCmdOptionPayload(argv, argv + argc, "--config");
  int Rmax = latt[0] / 2, Tmax = latt[Nd - 1] / 2;
  int nsmear = 0;
  RealD smear_alpha = 0.5;
  std::string s;
  if (GridCmdOptionExists(argv, argv + argc, "--rmax")) {
    s = GridCmdOptionPayload(argv, argv + argc, "--rmax"); GridCmdOptionInt(s, Rmax);
  }
  if (GridCmdOptionExists(argv, argv + argc, "--tmax")) {
    s = GridCmdOptionPayload(argv, argv + argc, "--tmax"); GridCmdOptionInt(s, Tmax);
  }
  if (GridCmdOptionExists(argv, argv + argc, "--smear")) {
    s = GridCmdOptionPayload(argv, argv + argc, "--smear"); GridCmdOptionInt(s, nsmear);
  }
  if (GridCmdOptionExists(argv, argv + argc, "--smear-alpha")) {
    smear_alpha = std::stod(GridCmdOptionPayload(argv, argv + argc, "--smear-alpha"));
  }

  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, cfg);

  RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(Umu);
  std::cout << "# config " << cfg << "  plaquette " << std::setprecision(10) << plaq << std::endl;
  if (nsmear > 0) {
    for (int n = 0; n < nsmear; ++n) ape_smear_spatial(Umu, smear_alpha);
    std::cout << "# spatial APE smearing: " << nsmear << " steps, alpha=" << smear_alpha
              << std::endl;
  }
  std::cout << "# R  T  W(R,T)" << std::endl;
  for (int R = 1; R <= Rmax; ++R) {
    for (int T = 1; T <= Tmax; ++T) {
      RealD W = WilsonLoops<Gimpl>::avgTimelikeWilsonLoop(Umu, R, T);
      std::cout << R << " " << T << " " << std::setprecision(12) << W << std::endl;
    }
  }

  Grid_finalize();
  return 0;
}
