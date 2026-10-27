# Plan: "How many replicates do we need?" screen simulator

## Goal
Given a proposed screen (compounds, doses, targets, plate format) and assay statistics, tell the user the replicates per condition needed to reach a target power, and show the resulting plate layout. One Flask app, backend and frontend in a single file (`app.py`).

## Inputs
| Group | Field | Notes |
|---|---|---|
| Design | plate format | 96 (8x12), 144 (9x16, assumed), 384 (16x24) |
| | compounds, doses/compound, targets | full factorial: compound x dose x target |
| | vehicle control wells/plate | on every plate |
| | skip edge wells | mitigates edge effect |
| Stats | minimum effect (delta), well SD | signal units |
| | alpha, target power | Bonferroni optional (over compound x dose comparisons) |
| Control | replicate override, seed | blank override = auto; seed fixes randomisation |

## Method
1. **Conditions** = compounds x doses x targets.
2. **Usable wells** per plate = total - edge wells (if skipped); treatment wells per plate = usable - control wells.
3. For n = 2..100 replicates: plates = ceil(conditions x n / treatment wells per plate); control n = plates x controls per plate.
4. **Power** of a two-sided t-test, treatment (n) vs pooled vehicle (n_ctrl), exact noncentral t:
   `ncp = delta / (sd * sqrt(1/n + 1/n_ctrl))`, df = n + n_ctrl - 2.
5. Required n = smallest n with power >= target (alpha divided by compounds x doses if Bonferroni).
6. **Layout**: shuffle all treatment wells (seeded), split across plates, add controls to every plate, shuffle within plate.

## Outputs (dashboard)
- KPIs: required replicates, replicates used, achieved power, conditions, treatment/control wells, plates, adjusted alpha.
- Power-vs-replicates curve with target line and the n in use.
- Per-plate tabs: well grid coloured by compound, darker = higher dose, grey = vehicle, dashed = skipped; hover shows compound/dose/target.

## API
`POST /api/plan` (JSON of the inputs) -> plan + curve + layout; 400 with `{error}` on bad input. `GET /` serves the page, which recomputes on every input change (200 ms debounce).

## Phase 2 design from data (simulate -> estimate -> redesign)
Closes the loop: the power calc above needs delta and SD, which in practice come from a previous screen. Simulated data stands in for that screen and uses the same pipe as real data.

1. **Simulate** a screen on the current layout. Per well:
   `y = baseline + plate_effect + edge_effect + top * dose^h / (EC50^h + dose^h) + noise`
   - truth table per compound x target: active or not, EC50, Hill slope h, max effect (inactive = 0); hit rate is an input.
   - nuisance: plate offset ~ N(0, sd_plate), optional edge gradient, noise ~ N(0, sd_well) (optionally multiplicative/CV).
   - seeded, so a run reproduces. Returns a long table: plate, well, compound, dose, target, value.
2. **Ingest**: `POST /api/simulate` returns that table; `POST /api/upload` accepts the same columns as CSV, so real data replaces simulated data with no other change.
3. **Estimate** (numpy/scipy only): well SD from vehicle wells (pooled and per plate, giving a plate variance component), per-compound x target hit call (vehicle-corrected, multiplicity-adjusted), 4-parameter Hill fit for EC50 and max effect with CI.
4. **Redesign phase 2** from the estimates:
   - carry forward only hits (plus controls); drop inactive compounds so plates are spent on what matters;
   - re-centre the dose series around each hit's EC50 (e.g. half-log steps spanning 0.1x-10x);
   - set delta to the smallest effect worth resolving (default: the weakest carried hit's estimated effect) and SD to the estimated SD, then re-run the power calc for the replicate count;
   - include the plate variance component by raising the control count or blocking by plate if it is large;
   - emit the phase 2 layout through the existing plate view.
5. **Validate** the planner: repeat simulate -> estimate -> redesign N times (Monte Carlo) and report hit recovery, false-positive rate and the spread of fitted EC50 against truth. Doubles as a check of the analytic power.

New UI: a "Simulate screen" button, then panels for estimates (SD, plate variance, hit table, fitted curves) and a "Design phase 2" button that fills the form with the derived inputs and shows the new layout beside phase 1.

Build order: (a) simulate + long table, (b) estimators, (c) phase 2 derivation into the existing form, (d) CSV upload, (e) Monte Carlo validation.

## Known limits / next steps
- Pairwise t-test vs control, not Dunnett; wells assumed independent with equal SD.
- Randomisation is unblocked: no plate/row/column effects, no true DOE designs (fractional factorial, response surface).
- Layout is random, not optimised (e.g. no balancing of compounds across plates, no edge-effect modelling).
- Next, if needed: Dunnett power, balanced block design, CSV export of the layout (power now also checked by Monte Carlo in the phase 2 section).

## Run
```
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python app.py   # http://localhost:5000
```
