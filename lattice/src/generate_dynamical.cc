// generate_dynamical.cc -- dynamical SU(3) N_f=2 FUNDAMENTAL Wilson HMC, with CLI (beta, mass, seed),
// for the Eternal Dawn torsiton PRODUCTION run. This is the sea the quenched pilot (run/06 on
// pure-gauge configs) lacked: here the fermion determinant is included, so the vacuum is the real,
// fluctuating, chiral-symmetry-breaking condensate -- the "weirder than mean field" vacuum where the
// physical torsiton mass, and any excited rungs (candidate further generations), actually live.
//
// FAITHFUL THEORY (not the dead sextet): the torsiton is the SU(3)-FUNDAMENTAL baryon (colour from
// the Pauli label), so the sea quarks are fundamental -- the standard, simplest dynamical SU(3)
// setup. N_f=2 (two-flavour pseudofermion) is the cheap QCD-like choice that confines and breaks
// chiral symmetry. Built on Grid's GenericHMCRunner.
//
//   generate_dynamical --grid Lx.Ly.Lz.Lt --beta 5.6 --mass -0.5 --seed 1001 \
//                      --mdsteps 20 --save-interval 5 --StartingType HotStart --Trajectories 2000
//   [--hasenbusch "-0.2,0.2"]   # OPT-IN accelerator -- see SAFETY below
// Writes NERSC ckpoint_lat.<n> into the CURRENT directory (run each stream in its own dir).
//
// ============================================================================================
// SPEEDUP (for the light-sea Stage-2 run -- the 96 GB GH200 and the 3-day wall clock)
// --------------------------------------------------------------------------------------------
//  (1) EVEN-ODD (red-black) PRECONDITIONING -- DEFAULT, ALWAYS ON.
//      TwoFlavourEvenOddPseudoFermionAction instead of the plain TwoFlavourPseudoFermionAction.
//      This is an EXACT determinant identity: det(M) = det(M_ee) det(M_oo - M_eo M_ee^{-1} M_oe),
//      so it CANNOT change the ensemble -- it only halves the solver's condition number (~2x).
//      The worst it can do if mistuned is lower the acceptance; it cannot bias the physics.
//
//  (2) HASENBUSCH MASS PRECONDITIONING -- OPT-IN via --hasenbusch, OFF by default.
//      Splits the N_f=2 determinant into a telescoping product of mass-ratio factors:
//        det(M(m_sea)^d M(m_sea))
//          = det( M(m_sea)^d M(m_sea) / M(m_1)^d M(m_1) )      <- expensive, but SMOOTH force
//          * det( M(m_1)^d M(m_1)    / M(m_2)^d M(m_2) ) * ...
//          * det( M(m_n)^d M(m_n) )                             <- cheap, heavy, stiff force
//      with the Hasenbusch masses m_1 < m_2 < ... HEAVIER than the sea mass m_sea. Each factor has
//      a smaller MD force, so you take bigger steps / fewer CG iterations -> another ~2-3x. The
//      numerator carries the LIGHTER mass, the denominator the HEAVIER (Grid convention -- VERIFY,
//      see below). The product telescopes back to the exact sea determinant IF the convention is
//      right; if it is BACKWARDS it silently samples the WRONG action and acceptance still looks
//      fine. THEREFORE Hasenbusch is GATED and MUST be cross-validated (next).
//
//  ---- MANDATORY VALIDATION before committing 3 days to this binary ----
//   a) It COMPILES against your Grid (I cannot compile Grid here -- see "VERIFY AT COMPILE").
//   b) SHORT run (~200-400 trajectories), check acceptance ~0.7-0.9 and dH ~ O(1).
//   c) If using --hasenbusch: generate ~20-40 decorrelated configs and CHECK THE PLAQUETTE (and,
//      cheaply, m_pi via run/06) AGREE with the existing plain m=-0.5 ensemble within errors.
//      Even-odd alone (no --hasenbusch) is an exact identity and needs only (a),(b).
//
//  (3) LOW-MODE DEFLATION -- OPT-IN via --deflate N, and ONLY in a binary built -DUSE_DEFLATION.
//      The genuine 96 GB exploit: hold N low eigenvectors of the sea operator resident and use them
//      as a deflated initial guess, cutting CG iterations near the chiral point. EXPERIMENTAL
//      SCAFFOLD -- the Grid eigensolver API is version-specific and the in-run subspace refresh is an
//      owed hook; see the USE_DEFLATION block. Exact (action solve refines to tolerance) and
//      reversible (subspace frozen for the run, recompute on restart) by construction, but the
//      acceleration decays as the gauge field decorrelates until the refresh hook is added. Default
//      OFF and behind the compile flag so it can NEVER break the validated even-odd+Hasenbusch build.
//
//  Still OFF on purpose: mixed-precision CG (single-precision inner solve, ~1.5-2x; safe but the
//  SP-grid wiring is fiddly) -- a clean follow-up once the above are validated.
// ============================================================================================
//
// (beta, mass) ARE THE KNOBS: lower mass toward the DYNAMICAL critical point for a light pion (the
// sea suppresses the exceptional configs that walled the quenched scan, so you can go lighter); pick
// beta in the scaling window. Build via lattice/src/Makefile. trajL=1.
#include <Grid/Grid.h>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#ifdef USE_DEFLATION
// ============================================================================================
// LOW-MODE DEFLATION accelerator (the 96 GB VRAM filler) -- EXPERIMENTAL SCAFFOLD, build with
// -DUSE_DEFLATION only. This is NOT validated blind: the Grid Lanczos / eigensolver API is
// version-specific, and the production-grade win needs the subspace REFRESHED at trajectory
// boundaries (see the staleness note below), which must be hooked into your Grid's HMC runner.
//
// SAFETY (must hold for the ensemble to stay correct):
//   * EXACTNESS: the action (accept/reject) solve must still converge to the tight CG tolerance --
//     deflation only changes the INITIAL GUESS, so a converged solve gives the identical action.
//   * REVERSIBILITY: the deflation subspace must be FROZEN within a trajectory (and across any
//     reversibility check). Refresh ONLY between trajectories. The simple scaffold below freezes
//     the subspace for the WHOLE run (built once on the starting config) -> trivially reversible
//     and exact, but the acceleration DECAYS as the gauge field decorrelates; recompute it on each
//     checkpoint/restart (cheap, automatic) to keep it fresh. Sustained in-run refresh = owed hook.
//
// A self-contained deflated solver: form the low-mode initial guess  psi0 = sum_i (v_i^dag b)/lam_i v_i
// in the stored eigenbasis, then refine with an ordinary CG. Cannot bias physics (CG refines to tol).
// NB: defined inside namespace Grid so OperatorFunction/RealD/ComplexD/Zero/innerProduct resolve here
// (the `using namespace Grid` is local to main()).
namespace Grid {
template <class Field>
class DeflatedCG : public OperatorFunction<Field> {
 public:
  std::vector<Field> &evec;          // eigenvectors of the SAME operator being solved (VRAM resident)
  std::vector<RealD> &eval;          // their eigenvalues
  ConjugateGradient<Field> &cg;      // the refiner (tight tolerance -> exact action)
  DeflatedCG(std::vector<Field> &v, std::vector<RealD> &e, ConjugateGradient<Field> &c)
      : evec(v), eval(e), cg(c) {}
  void operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) {
    psi = Zero();
    for (size_t i = 0; i < evec.size(); ++i) {           // deflated initial guess in the low modes
      ComplexD a = innerProduct(evec[i], src);
      psi = psi + (a / eval[i]) * evec[i];
    }
    cg(Linop, src, psi);                                 // refine to tolerance: result is exact-to-tol
  }
};
}  // namespace Grid
#endif

int main(int argc, char **argv) {
  using namespace Grid;
  Grid_init(&argc, &argv);

  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;     // fundamental: no Hirep machinery
  typedef WilsonImplR FermionImplPolicy;
  typedef WilsonFermion<FermionImplPolicy> FermionAction;  // this Grid lacks the WilsonFermionR
                                                           // typedef; instantiate the template
  typedef typename FermionAction::FermionField FermionField;

  HMCWrapper TheHMC;

  // ---- CLI: beta, mass, seed, integrator steps, checkpoint cadence, hasenbusch ----
  RealD beta = 5.6;        // bare gauge coupling -- TUNE (scaling window)
  RealD mass = -0.5;       // fundamental Wilson sea-quark mass -- TUNE toward the dynamical critical
  int seed = 1;
  int mdsteps = 20;        // leapfrog steps per trajectory; raise ~V^{1/4} / for the fermion force
  int saveInterval = 5;    // write a NERSC config every Nth trajectory
  std::string hasenStr;    // OPT-IN: comma-separated Hasenbusch masses, HEAVIER than `mass`
  if (GridCmdOptionExists(argv, argv + argc, "--beta"))
    beta = std::stod(GridCmdOptionPayload(argv, argv + argc, "--beta"));
  if (GridCmdOptionExists(argv, argv + argc, "--mass"))
    mass = std::stod(GridCmdOptionPayload(argv, argv + argc, "--mass"));
  if (GridCmdOptionExists(argv, argv + argc, "--seed")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--seed");
    GridCmdOptionInt(s, seed);
  }
  if (GridCmdOptionExists(argv, argv + argc, "--mdsteps")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--mdsteps");
    GridCmdOptionInt(s, mdsteps);
  }
  if (GridCmdOptionExists(argv, argv + argc, "--save-interval")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--save-interval");
    GridCmdOptionInt(s, saveInterval);
  }
  if (GridCmdOptionExists(argv, argv + argc, "--hasenbusch"))
    hasenStr = GridCmdOptionPayload(argv, argv + argc, "--hasenbusch");
  int deflateN = 0;        // OPT-IN low-mode deflation (the VRAM filler); needs -DUSE_DEFLATION build
  if (GridCmdOptionExists(argv, argv + argc, "--deflate")) {
    std::string s = GridCmdOptionPayload(argv, argv + argc, "--deflate");
    GridCmdOptionInt(s, deflateN);
  }

  // parse the (optional) Hasenbusch masses: must be HEAVIER (larger) than the sea `mass`, ascending
  std::vector<RealD> hmasses;
  if (!hasenStr.empty()) {
    std::stringstream hs(hasenStr);
    std::string tok;
    while (std::getline(hs, tok, ',')) {
      if (!tok.empty()) hmasses.push_back(std::stod(tok));
    }
    std::sort(hmasses.begin(), hmasses.end());           // ascending: m_1 < m_2 < ... (all > m_sea)
  }

  TheHMC.Resources.AddFourDimGrid("gauge");

  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = saveInterval;
  CPparams.format = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  std::ostringstream ss, ps;                              // independent seeds per stream
  for (int i = 0; i < 5; ++i) ss << (seed + i) << (i < 4 ? " " : "");
  for (int i = 0; i < 5; ++i) ps << (seed + 5 + i) << (i < 4 ? " " : "");
  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = ss.str();
  RNGpar.parallel_seeds = ps.str();
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // ---- actions: N_f=2 fundamental pseudofermion(s) (Level1) + Wilson gauge (Level2, finer) ----
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  LatticeGaugeField U(GridPtr);                          // fundamental links (updated by the runner)

  // One CG per inversion (double precision). Tight stop for the accept/reject Hamiltonian.
  // VERIFY AT COMPILE: ConjugateGradient<FermionField>(tol, maxiter, err_on_no_conv).
  ConjugateGradient<FermionField> CG(1.0e-8, 4000, false);

  // The mass ladder: the SEA mass first, then the (heavier) Hasenbusch masses. The fermion operators
  // own no state we mutate, but they reference U, which the runner evolves; build them on the heap so
  // they outlive this scope while the actions hold pointers to them.
  std::vector<RealD> ladder;
  ladder.push_back(mass);
  for (auto mh : hmasses) ladder.push_back(mh);          // ascending, all heavier than the sea

  std::vector<FermionAction *> fops;
  for (auto m : ladder)
    fops.push_back(new FermionAction(U, *GridPtr, *GridRBPtr, m));

  // ---- solver selection: plain CG, or (opt-in, -DUSE_DEFLATION) low-mode deflated CG ----
  // The action constructors take an OperatorFunction&; default to CG. With deflation we build a
  // low-mode subspace of the SEA operator (fops[0], the stiffest/most expensive solves) and refine
  // through it. (The Hasenbusch HEAVY solves see the sea-mass subspace as a harmless, suboptimal
  // guess -- CG still refines to tolerance, so it stays exact; for full benefit each mass wants its
  // own subspace, an owed refinement.)
  OperatorFunction<FermionField> *slvp = &CG;            // pointer, so deflation can REPOINT it
#ifdef USE_DEFLATION
  std::vector<FermionField> defl_evec;
  std::vector<RealD> defl_eval;
  DeflatedCG<FermionField> *DCG = nullptr;
  if (deflateN > 0) {
    std::cout << GridLogMessage << "deflation: building " << deflateN
              << " low modes of the SEA even-odd M^dag M (VRAM-resident) ..." << std::endl;
    // ---- Implicitly-Restarted Lanczos on the even-odd Schur operator M_ee^dag M_ee (Hermitian PD).
    // VERIFY AT COMPILE -- Grid's IRL ctor signature and the Cheby/Op wrappers vary by version; the
    // TUNE constants (spectral window lo/hi, Chebyshev order, Krylov sizes Nk/Nm) need a few tries.
    // If anything here won't build, this whole block is behind -DUSE_DEFLATION, so the safe binary is
    // unaffected; and if calc() returns 0 modes, slvp stays &CG (graceful fallback to plain CG).
    typedef SchurDiagMooeeOperator<FermionAction, FermionField> SchurOp;
    SchurOp HermOp(*fops[0]);                              // M_ee^dag M_ee on the red-black checkerboard
    const RealD cheb_lo = 1.0e-3, cheb_hi = 64.0;          // TUNE: bracket the M^dag M spectrum
    const int cheb_ord = 60;                               // TUNE: Chebyshev order (higher = sharper)
    Chebyshev<FermionField> cheby(cheb_lo, cheb_hi, cheb_ord);
    FunctionHermOp<FermionField> ChebyOp(cheby, HermOp);   // filtered op (amplifies the low end)
    PlainHermOp<FermionField> PlainOp(HermOp);             // bare op (for the Ritz/eigenvalue pass)
    const int Nstop = deflateN;                            // wanted converged modes
    const int Nk = deflateN + 20;                          // TUNE: working set
    const int Nm = deflateN + 40;                          // TUNE: max Krylov dim (Nm > Nk > Nstop)
    ImplicitlyRestartedLanczos<FermionField> IRL(ChebyOp, PlainOp, Nstop, Nk, Nm, 1.0e-8, 200);
    GridBase *frb = fops[0]->FermionRedBlackGrid();        // the exact checkerboard grid the Schur op solves on
    defl_eval.resize(Nm);
    defl_evec.reserve(Nm);                                 // Grid Lattice is not vector-copy-constructible:
    for (int i = 0; i < Nm; ++i)                           // build each RB-grid field in place (no copy)
      defl_evec.emplace_back(frb);                         // eigenvectors live on the RB checkerboard grid
    FermionField src0(frb);
    GridParallelRNG pRNG(frb);                             // local RNG for the Lanczos start vector
    pRNG.SeedFixedIntegers(std::vector<int>{seed, seed + 1, seed + 2, seed + 3});
    gaussian(pRNG, src0);
    int Nconv = 0;
    IRL.calc(defl_eval, defl_evec, src0, Nconv);
    defl_evec.resize(Nconv);
    defl_eval.resize(Nconv);
    DCG = new DeflatedCG<FermionField>(defl_evec, defl_eval, CG);
    if (!defl_evec.empty()) slvp = DCG;                  // repoint the solver at the deflated one
    std::cout << GridLogMessage << "deflation: " << defl_evec.size()
              << " modes held; subspace FROZEN for the run (recompute on restart to refresh)."
              << std::endl;
  }
#else
  if (deflateN > 0)
    std::cout << GridLogMessage << "WARNING: --deflate " << deflateN
              << " ignored -- rebuild with -DUSE_DEFLATION to enable the deflation accelerator."
              << std::endl;
#endif
  OperatorFunction<FermionField> &slv = *slvp;           // CG by default, the deflated solver if built

  // The pseudofermion actions. EVEN-ODD preconditioned throughout (exact; ~2x).
  //  - no Hasenbusch (ladder size 1): a single EvenOdd two-flavour det(M_sea^d M_sea).
  //  - with Hasenbusch: a telescoping chain of EvenOdd RATIO actions (light numerator / heavy
  //    denominator) capped by a single EvenOdd two-flavour action on the heaviest mass.
  // VERIFY AT COMPILE: class names + ctor signatures in YOUR Grid --
  //   TwoFlavourEvenOddPseudoFermionAction<Impl>(FermOp&, CG&, CG&);
  //   TwoFlavourEvenOddRatioPseudoFermionAction<Impl>(NumOp&, DenOp&, CG&, CG&);   // num=LIGHT
  // (Some Grid versions name the ratio class TwoFlavourEvenOddRatioPseudoFermionAction; if absent,
  //  the non-EO TwoFlavourRatioPseudoFermionAction exists -- but prefer the EO form.)
  std::vector<Action<HMCWrapper::Field> *> pfActions;
  if (ladder.size() == 1) {
    auto *Nf2 = new TwoFlavourEvenOddPseudoFermionAction<FermionImplPolicy>(*fops[0], slv, slv);
    Nf2->is_smeared = false;
    pfActions.push_back(Nf2);
  } else {
    // ratios: det(M_i^d M_i / M_{i+1}^d M_{i+1}), i = 0 .. n-1  (numerator LIGHTER = fops[i])
    for (size_t i = 0; i + 1 < ladder.size(); ++i) {
      auto *ratio = new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(
          *fops[i], *fops[i + 1], slv, slv);            // num = fops[i] (lighter), den = fops[i+1]
      ratio->is_smeared = false;
      pfActions.push_back(ratio);
    }
    // cap: the heaviest mass as a plain EvenOdd two-flavour determinant
    auto *cap = new TwoFlavourEvenOddPseudoFermionAction<FermionImplPolicy>(*fops.back(), slv, slv);
    cap->is_smeared = false;
    pfActions.push_back(cap);
  }

  WilsonGaugeActionR Waction(beta);

  ActionLevel<HMCWrapper::Field> Level1(1);              // fermion level (coarse)
  for (auto *a : pfActions) Level1.push_back(a);
  ActionLevel<HMCWrapper::Field> Level2(4);              // gauge updated 4x per fermion step
  Level2.push_back(&Waction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  TheHMC.Parameters.MD.MDsteps = mdsteps;
  TheHMC.Parameters.MD.trajL = 1.0;

  std::cout << GridLogMessage << "generate_dynamical: SU(3) Nf=2 fundamental  beta=" << beta
            << " mass=" << mass << " seed=" << seed << " mdsteps=" << mdsteps
            << " saveInterval=" << saveInterval << "  even-odd=ON";
  if (hmasses.empty()) {
    std::cout << "  hasenbusch=OFF (exact det, no cross-check needed)";
  } else {
    std::cout << "  hasenbusch=[";
    for (size_t i = 0; i < hmasses.size(); ++i)
      std::cout << hmasses[i] << (i + 1 < hmasses.size() ? "," : "");
    std::cout << "]  <-- CROSS-CHECK plaquette/m_pi vs the plain ensemble before production!";
  }
  std::cout << std::endl;

  TheHMC.ReadCommandLine(argc, argv);   // --grid, --Trajectories, --Thermalizations, --StartingType, ...
  TheHMC.Run();

  for (auto *a : pfActions) delete a;
  for (auto *f : fops) delete f;
#ifdef USE_DEFLATION
  if (DCG) delete DCG;
#endif
  Grid_finalize();
  return 0;
}
