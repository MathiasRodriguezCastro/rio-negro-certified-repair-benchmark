# Synthetic certified hydropower repair benchmark v1.2.0

This MIT-licensed standalone package accompanies *Certified Affine Production
Bounds for Local Hydropower Schedule Repair in the Río Negro Cascade*. Code,
inputs and initial MILP witnesses are independently written or synthetically
generated. No Río Negro coefficients, real-bank scenarios, inherited native
classes, network weights or proprietary solvers are distributed here.
Third-party dependencies retain their own licenses.

## New in v1.2.0: complete global-bound certificate

The additional [24-period example](global_bounds_example/README.md) provides all
input data, upper production planes, three rational outer LPs, dual residual
corrections and both exact feasible reference points. A separately written
verifier proves each upper plane on its entire box using rational Bernstein
coefficients, reconstructs the model, and checks the dual and hydraulic witness.
It runs with the Python standard library alone:

```bash
python3 -S verify_global_bounds_demo.py
python3 reproduce_global_bounds.py  # repeats 4 open-source LPs
python3 reproduce_all.py            # all 232 LPs, including the retained 228
```

All v1.1.0 scientific files and earlier release archives are retained unchanged.
This synthetic example demonstrates the global-bounding procedure; it does not
redistribute or independently recreate the restricted real-bank inputs. Its
optimum is the continuous reference optimum within fixed local boxes. Normalized
gap caps use the all-demand failure penalty, not optimal operating cost.

## Retained v1.1.0 replay

```bash
python3 -m pip install -r requirements.txt
python3 reproduce_full_state.py
```

Python 3.10 or later, NumPy and SciPy/HiGHS suffice. Replay needs no network or
external project directory. The new entry point retains the original 192-LP
replay and adds 36 full-state LPs: 228 solves, 5,184 independently checked rational
whole-box production rows, and all 12 complete initial MILP witnesses. To check
saved outputs without rerunning optimization:

```bash
python3 reproduce_full_state.py --audit-only
```

The original `benchmark.py`, `reproduce.py`, test data and reference results from
v1.0.1 remain unchanged. `python3 reproduce.py` still performs the original
fixed-storage replay. The v1.0.1 public release and archive remain immutable.
`FULL_STATE_REPRODUCTION_CHECKS.json` documents the complete new isolated replay;
`REPRODUCTION_CHECKS.json` documents the original sub-replay. Costs and exact
certificates are compared; elapsed times are observations that vary on replay.

## New in v1.1.0: Taylor including storage

`full_state.py` forms `px*x + py*y + pv*v + k`, including the ordinary storage
Taylor derivative. It certifies this prescribed form over exactly the same
boxes and initial MILPs as v1.0.1. This derivative itself is standard calculus;
the point is retaining uniform lower-bound validity for the more informative
form. The error remains convex in turbine flow at fixed other coordinates.
At each flow endpoint it separates into quadratics in spill and storage.
Endpoint and clipped stationary-point evaluation therefore uses at most 8, 12,
and 12 candidates for the three synthetic plant curvature patterns. A generic
family with both concave spill and concave storage sections could require 18;
we do not claim twelve for all coefficients.

The original production monotonicity does not imply error monotonicity once
`pv` is included. `interior_volume_witness.json` gives an independent synthetic
Baygorria-type curvature example in which both volume endpoints have negative
error but the interior maximum is positive. Omitting that state candidate would
produce an invalid certificate. It contains synthetic coefficients, not an
additional real Río Negro trajectory or coefficient set.

Three complete full-state certificates use identical forms, anchors and boxes:
componentwise Hessian remainder, centered cubic subdivision (4 by 4 by 2), and
exact rational candidates. Shifts round upward and the final exported intercept
rounds downward; an independent completed-square calculation checks the actual
rounded rows, including `pv`. `results/full_state/metrics.csv` also reports paired
scheduled and final economic differences from the original fixed-storage
candidate LP. No configuration is selected from these comparisons. The cases
were already inspected; these are retrospective verification examples, not new
independent statistical validation or an out-of-distribution test.

## Model and data

All flows, storage, generation, time and cost are normalized synthetic quantities,
not a calibration of a real power system. Three reservoirs obey

`v[i,t] = v[i,t-1] + 0.2*(a[i,t] + upstream[i,t] - x[i,t] - y[i,t])`.

Upstream flow has a one-period delay to reservoir 2 and a two-period delay to
reservoir 3; releases before the horizon are zero. The 24-period limits are
`0 <= x <= 1`, `0 <= y <= 0.8`, and `2 <= v <= 6`. For each plant,

`P(x,y,v) = a*x + b*x*v + c*x*v^2 + d*x*(x+y) + e*x*(x+y)^2`.

The three coefficient tuples are `(1.25,0.12,-0.006,-0.08,0)`,
`(1.9,-0.08,0.003,-0.07,0)`, and `(1.45,0.10,-0.004,-0.08,0.006)`.
These are deliberately different from the real benchmark. The relevant state
derivative has a fixed sign and `F_xx >= 0`, `F_yy <= 0` over the technical box.
Exact comparisons treat every stored binary64 coefficient as its rational value.

The input generator uses seed 20261113, with disjoint offsets for three
development and twelve test cases. No observed hydrology or demand is used.
The initial MILP uses two head bands `[2,4]`, `[4,6]`, their midpoints 3 and 5,
binary band selectors, split turbine flows, and a coarse head-only production
model. Its approximation can over-credit supply. Complete selectors and split
flows are included so the original MILP witness can be checked independently.
Its saved convergence/gap is an observation, not an exact integer optimality
proof. To regenerate initial MILPs in a separate copy, remove that copy's
`data/test` directory and run `python3 benchmark.py generate --split test`.
The published reference repair uses the supplied frozen MILP witnesses. A new
MILP run can choose a different tied incumbent on another solver/platform;
regenerating such starts creates a new benchmark realization rather than a
requirement for replaying the supplied reference.

## Certificates, decisions and coordination

All four LP variants use the same initial flows, reference volumes and explicit
boxes (`r=0.08`, `m=0.04`). They differ only in the production shift: Taylor
(zero shift), derivative-based bound, 4-by-4 centered interval bound, or exact
finite candidates. Each mathematical certificate protects a complete box; an
LP status is additionally checked with a normalized point tolerance of `1e-7`.
Signed credits are accounting variables. Demand service is bounded by aggregate
credit. Available physical production can be freely curtailed in this toy model.

`scheduled_cost` values returned service. `final_cost` evaluates maximum reference
availability at the same hydraulic point. Their difference is explicit; final
credit is not used to hide an unsupported raw promise. Taylor outputs are audited
and retained even when their scheduled service is unsupported.

The six coordination problems change service prices, reverse those prices,
increase terminal-water values, add a terminal-storage floor, add eight-period
service floors, or add aggregate-flow ramps and selected-hour service caps.
The original production rows and boxes remain unchanged. The development-only
implementation probes are described in the paper's supplementary material.
Settings are frozen before the twelve test cases are evaluated.

During the implementation audit after the first small replay, the mixed-partial
expression was corrected to `P_xy = d + 4*e*x + 2*e*y`, with an independent exact
finite-difference Hessian test. The previous exported bounds had also passed
full-box rational verification. All small comparisons are nevertheless regenerated
with the corrected expression, using identical data, initial MILP starts and
parameters. The protocol records both code versions. These twelve cases are
verification examples, not a claim of an untouched prospective test after that
correction. No historical real-bank campaign is modified.

For each valid candidate baseline, cached and rebuilt certificates are compared
on the same six problems, alternating their order by case and variant. Both
paths use fresh LP calls, so the comparison isolates certificate reconstruction,
not a warm-start advantage. Rebuilt rows must have exactly the same hashes.
The 72 cached coordination solves and 72 counterfactual rebuilt solves are
reported separately. Different coordination objectives are not claims of better
original economic cost. These exact-rational Python construction times are not
performance estimates for the native Río Negro implementation.

`test_certificates.py` additionally checks nonquadratic/nonpolynomial concave
slices, a flat interval of stationary maxima, affine/degenerate cases, multiple
monotone states, and counterexamples when the derivative hypotheses are omitted.
The runnable hydropower kernel handles the explicitly stated polynomial; the
general theorem does not promise closed-form roots for every smooth function.

## Files and limits

- `data/`: generated scenarios and full initial MILP witnesses.
- `results/test/`: every returned point, affine row, solver status and timer.
- `protocol.json`: frozen code/configuration hashes and experiment design.
- `benchmark.py`: generation, exact bounds, LP/MILP models and coordination runs.
- `reproduce.py`: independent full-box and point checks plus optional LP replay.
- `MANIFEST.sha256`: package inventory; no external files are required.

This is an independently distributable verification example, not a replacement
for the 10,000-case real-bank evaluation or a claim of operational realism.
Public releases are hosted at https://github.com/MathiasRodriguezCastro/rio-negro-certified-repair-benchmark/releases . A persistent DOI remains pending.
The official dependency interfaces are documented at
[SciPy linprog](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html)
and [SciPy milp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.milp.html).
