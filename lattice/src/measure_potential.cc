// measure_potential.cc -- static quark potential from timelike Wilson loops (Eternal Dawn L1-L2).
// Reads one NERSC gauge configuration and prints the average timelike Wilson loop W(R,T) for a
// ladder of spatial R and temporal T. The run script loops over configs and averages; the
// Python analyzer (cartasis_sims.lattice) extracts V(R) = ln[W(R,T)/W(R,T+1)] and fits the
// Cornell form V(R)=c-alpha/R+sigma*R for the string tension sigma.
//
// Build with lattice/src/Makefile (uses grid-config). Written against Grid's WilsonLoops API
// (Grid/qcd/utils/WilsonLoops.h): avgPlaquette, avgTimelikeWilsonLoop are static members.
#include <Grid/Grid.h>
#include <iomanip>

using namespace Grid;

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  typedef PeriodicGimplR Gimpl;                 // standard periodic SU(3) gauge implementation

  Coordinate latt = GridDefaultLatt();
  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(
      latt, GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  Gimpl::GaugeField Umu(UGrid);

  if (!GridCmdOptionExists(argv, argv + argc, "--config")) {
    std::cout << "usage: measure_potential --grid Lx.Ly.Lz.Lt --config <NERSC cfg> "
                 "[--rmax R] [--tmax T]" << std::endl;
    Grid_finalize();
    return 1;
  }
  std::string cfg = GridCmdOptionPayload(argv, argv + argc, "--config");
  int Rmax = latt[0] / 2, Tmax = latt[Nd - 1] / 2;
  std::string s;
  if (GridCmdOptionExists(argv, argv + argc, "--rmax")) {
    s = GridCmdOptionPayload(argv, argv + argc, "--rmax"); GridCmdOptionInt(s, Rmax);
  }
  if (GridCmdOptionExists(argv, argv + argc, "--tmax")) {
    s = GridCmdOptionPayload(argv, argv + argc, "--tmax"); GridCmdOptionInt(s, Tmax);
  }

  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, cfg);

  RealD plaq = WilsonLoops<Gimpl>::avgPlaquette(Umu);
  std::cout << "# config " << cfg << "  plaquette " << std::setprecision(10) << plaq << std::endl;
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
