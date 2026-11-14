# Complete global-bound example

This example is part of release 1.2.0. It uses new normalized rational data and
independent coefficients, with two reservoirs, one run-of-river plant, and 8/16
period transport delays. The 24-period hydraulic anchor is explicitly feasible;
it is not a newly optimized MILP. All controls and volumes have finite local
boxes at exactly r=m=1/50. Costs use synthetic monetary units (SMU).

Verify every saved mathematical certificate, without NumPy, SciPy or a solver:

```bash
python3 -S verify_global_bounds_demo.py
```

Rebuild the three outer LPs and one certified inner LP, then independently
verify their new evidence (requires the pinned NumPy/SciPy dependencies):

```bash
python3 reproduce_global_bounds.py
```

`python3 reproduce_all.py` also replays the retained examples, for 232 LP solves
in total. No real-bank source, licensed solver or network access is required.

## Supplied evidence

- `input.json`: all rational coefficients, profiles, delays, technical bounds,
  local radii, initial storage, water values, failure price and hydraulic anchor.
- `certificate.json`: every retained upper plane and its expansion point; all
  three complete rational LP matrices, finite bounds and objective constants;
  inequality/equality multipliers, nonzero stationarity residuals, exact box
  corrections, and the two complete reference-feasible upper-point vectors.
- `../global_bounds_demo.py`: signed-Hessian upper envelopes, three accumulating
  outer rounds, a full-Taylor certified inner LP and exact hydraulic realization.
- `../verify_global_bounds_demo.py`: separately written, standard-library-only
  verifier, importing neither the generator nor its Hessian construction.

## Why the bounds hold

For each actual upper plane, the verifier expands `upper - production` exactly.
It maps the entire physical box to [0,1]^3 and computes tensor-product Bernstein
coefficients. Bernstein basis functions are nonnegative and sum to one, so a
nonnegative lower coefficient bound proves nonnegativity everywhere on the box.
Rational bisection is available to cover a box by certified subboxes. A fixed
maximum depth rejects an inconclusive proof; sampled points never authorize
acceptance. All 216 planes in this example pass directly on their full boxes.
The adversarial tests also cover a positive polynomial needing subdivision and
an interior violation missed by corner tests.

The verifier independently reconstructs the physical hydraulics, all finite
boxes, production rows, linear objective and restored constant from the input.
For an outer LP `min C+c'z`, `Az<=b`, `Ez=h`, `l<=z<=u`, any nonnegative lambda
and unrestricted nu give a lower bound

    L = C - lambda'b - nu'h + sum_j min(r_j*l_j, r_j*u_j),
    r = c + A'lambda + E'nu.

No stationarity or numerical dual feasibility is assumed. Each quantity is a
Fraction, and every retained residual and bound is checked exactly. All three
stored residual vectors are nonzero; ignoring their contribution is forbidden.
The best of the three valid lower bounds remains valid.

Both upper points are checked for exact hydraulic equalities, technical/local
bounds, nonnegative demand service, and support by the expanded nonlinear law.
They have identical hydraulic controls: scheduled service retains a safe
commitment, while final service uses all available generation up to demand. The
same objective, including its constant, is evaluated for both.

The resulting stored interval has lower endpoint 1564.7725536159999 SMU and
final upper endpoint 1565.7827293243256 SMU (these decimal values are illustrative;
the JSON contains exact rationals). The final suboptimality cap is at most
1.010176 SMU after upward rounding. The scheduled cap is at most 1.863030 SMU.
These concern the continuous reference optimum inside the declared boxes.
Normalization uses the all-demand failure penalty, 4860 SMU: it is neither an
optimal-operating-cost percentage nor a fraction of maximum attainable saving.

The six tests reject a corrupt upper plane, dual lower value, restored constant,
nonfinite box, inexact hydraulic point and unsupported service. This is a fixed
verification example developed before publication, with no statistical sample
or global-MILP optimality claim. It illustrates the procedure without recreating
the restricted twelve-case real-bank cost intervals.
