# DESI DR2 vs the parent-accretion dark energy prediction

## The prediction (Chapter 6)

Dark energy is the baby universe's response to its parent black hole's growth, so
the effective equation of state is NOT w=-1 but tracks the parent's accretion
history. A parent that is still accreting but with a saturating fractional growth
rate gives, in CPL terms, **w0 > -1 and wa < 0** -- dynamical dark energy that
was more negative in the past and rises toward today.

## The data (DESI DR2, March 2025; arXiv:2503.14738, 2503.14743)

VERIFY central values against the published tables before quoting.

| dataset | w0 | wa | sigma vs LCDM |
|---|---|---|---|
| DESI+CMB | -0.42 +/- 0.21 | -1.75 +/- 0.58 | 3.1 |
| DESI+CMB+Pantheon+ | -0.838 +/- 0.055 | -0.62 +/- 0.22 | 2.8 |
| DESI+CMB+Union3 | -0.667 +/- 0.088 | -1.09 +/- 0.31 | 3.8 |
| DESI+CMB+DESY5 | -0.752 +/- 0.057 | -0.86 +/- 0.22 | 4.2 |

**Every combination has w0 > -1 and wa < 0** -- all four sit in the framework's
predicted wedge, and LCDM (-1, 0) sits at its corner, disfavoured at 2.8-4.2 sigma.
The *sign pattern* the framework predicted is exactly what DESI sees.

## The honest caveat: phantom crossing

DESI's CPL fit implies w(z) < -1 (phantom) in the past, crossing to w > -1 today.
The minimal toy model here, w(a) = -1 + eps*p0*a^q with eps,p0,q > 0, **cannot go
phantom** -- it stays >= -1 at all z. So:

- it matches the *present-day* values (w0, wa) by construction (eps, q tuned), and
  reproduces the *sign* and the rise of w toward today;
- it does NOT reproduce the phantom past (w < -1 for z > ~0.5).

A phantom phase needs more than a monotonically growing parent -- e.g. a parent
whose accretion *accelerated* over the relevant epoch, or membrane back-reaction
that effectively gives w < -1 (super-negative pressure as the membrane area grows
faster than volume). Whether the persistent Cartasis interface naturally does this
is open and is the right next calculation. Caveat stated, not hidden.

## Verdict

A genuine, pre-registered *sign* prediction (w0 > -1, wa < 0) that matches the
headline DESI DR2 result and disfavours LCDM -- without a cosmological constant,
as a consequence of the parent black hole growing. The amplitude and the phantom
crossing are not yet derived (eps is a free membrane coupling; phantom needs
accelerating accretion or membrane back-reaction). Encouraging and falsifiable:
if future data drive (w0, wa) back to (-1, 0), the parent-accretion picture fails.

## To compute next

1. Derive the membrane back-reaction equation of state from the junction
   conditions (does growing membrane area give effective w<-1?).
2. Map a physical parent accretion history Mdot(t) (e.g. Eddington-limited then
   saturating) to w(z) and refit -- replace the ansatz with dynamics.
3. Fold dark/baryon (f~1/6) and BH-growth (M_crit) cross-checks into the same
   data chapter.
