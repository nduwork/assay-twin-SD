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

## Known limits / next steps
- Pairwise t-test vs control, not Dunnett; wells assumed independent with equal SD.
- Randomisation is unblocked: no plate/row/column effects, no true DOE designs (fractional factorial, response surface).
- Layout is random, not optimised (e.g. no balancing of compounds across plates, no edge-effect modelling).
- Next, if needed: Dunnett or simulation-based power, plate-effect variance component, balanced block design, CSV export of the layout.

## Run
```
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python app.py   # http://localhost:5000
```
