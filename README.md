# Synthetic certified hydropower repair benchmark

This self-contained example accompanies *Certified Affine Production Bounds for
Local Hydropower Schedule Repair in the Río Negro Cascade*. Everything in this
directory is newly written or synthetically generated. The code, documentation
and generated data in this directory are distributed under the included MIT
license. This license does **not** apply to any historical Río Negro material.
Third-party dependencies retain their own licenses and are installed separately.

There are no measured profiles, Río Negro coefficients, inherited native sources,
network weights, licensed solvers or external project imports. Python 3.10.12,
NumPy 2.2.6 and SciPy 1.15.3 were used. Python 3.10 or later is required. HiGHS
comes through SciPy; no separate commercial installation is needed.

```bash
python3 -m pip install -r requirements.txt
python3 reproduce.py
```

The default replay checks every saved input and MILP witness, independently
verifies every certified affine row on its full box with exact rational
arithmetic, audits all returned points, and reruns the 192 small LPs from the
supplied initial dispatches. It compares costs and mathematical certificates;
elapsed times are observations and are not expected to reproduce exactly.
For a saved-output audit without new optimization, use:

```bash
python3 reproduce.py --audit-only
```

`REPRODUCTION_CHECKS.json` states what was actually reconstructed. Replay uses
only this directory, Python, NumPy, SciPy and the standard library. There is no
network access during replay. Installing dependencies is a separate step.

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
No public hosting location or persistent DOI is asserted by creating this package.
The official dependency interfaces are documented at
[SciPy linprog](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html)
and [SciPy milp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.milp.html).
