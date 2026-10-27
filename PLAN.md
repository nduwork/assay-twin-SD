# Plan: "How many replicates do we need?" screen simulator

## Goal
Given a proposed screen (compounds, doses, targets, plate format) and assay statistics, tell the user the replicates per condition needed to reach a target power, and show the resulting plate layout. One Flask app, backend and frontend in a single file (`app.py`).

## Inputs
| Group | Field | Notes |
|---|---|---|
| Design | plate format | 96 (8x12), 144 (9x16, assumed), 384 (16x24) |
| | compounds, doses/compound, targets | full factorial: compound x dose x target |
| | top concentration, dilution factor | dose series in µM, default 10 µM top, 1:3 down |
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

## Follow-up dose-response design (drug-discovery flow)
Supersedes the earlier "Phase 2" sketch. The service takes a primary screen and designs the 96-well dose-response follow-up for the top leads.

### Built so far (generic calculator, `app.py`)
A dose-response simulator on the generic calculator's factorial layout; reusable for (a) and (g) of the build order below, but not yet the single-dose primary screen (% activity, Z', QC) this flow describes.
- `y = baseline + plate_offset + edge + top * x^h / (EC50^h + x^h) + noise`, x from top concentration and dilution factor.
  Truth per compound x target: active with p = hit rate (0.25); EC50 log-uniform half a log beyond each end of the tested range; hill U(0.7, 2); top U(0.5, 2) x delta, increasing only.
  Nuisance: plate offset SD (default 0.5 x well SD), optional edge gradient (halves per ring inward), additive or multiplicative (CV = SD / baseline) noise. Own seed.
- `POST /api/simulate`: long table (plate, well, compound, target, dose, conc_uM, value) + truth table + plate offsets; page shows a value heatmap on the plate grid and downloads CSV.
- `POST /api/mc`: Monte Carlo of the pooled t-test vs vehicle at 12 replicate counts: power at exactly delta (matches analytic within MC error when plate SD = 0, drops from 0.83 to 0.72 at the defaults with plate offsets), hit recovery over true actives, false-positive rate.

### Artifact dashboard (`dashboard.html`)
Hosting is a concern, so the follow-up design is also a static artifact page: everything runs in the browser (4PL fit, Fisher-information power, blocked layout, downloads), no server.
- **Input:** the simulator's long table (`sample_screen.csv`): `plate, well, compound, target, dose, conc_uM, value`. Vehicle wells are `compound = vehicle` (conc 0). A single-dose screen (`compound`, `value`) also works.
- **Ingest:** each plate's vehicle mean is subtracted; well SD comes from within-plate vehicle spread, plate-offset SD from the spread of plate vehicle means. Top concentration, dose count and dilution factor are filled from the data. With several targets the user picks one.
- **Leads:** ranked by effect at the top dose. Where a compound has 4 or more doses, top/EC50/Hill come from a Levenberg-Marquardt fit (bottom fixed at 0 after vehicle subtraction); otherwise from the one-dose assumption. On the example file the fitted top and EC50 of the five leads are within about 30% of the simulator's truth.
- **Design and outputs:** as in the sections above, with a smallest efficacy difference to resolve (readout units) as the success criterion. Visible inputs are limited to top leads, plates, dose points, top concentration, that difference and target power; the rest sit under Advanced settings.
- **Not yet in the page:** Monte Carlo confirmation of the power numbers (the Flask `/api/mc` does this for the generic calculator).

### Data we have and collect
- **Primary screen (one dose, often n=1-2):** compound id, concentration, raw signal, % activity normalised to plate controls, plate, well, per-plate QC (Z', control CV, signal window). Hit calls come from % activity, z-score or robust MAD.
- **Follow-up (to be collected):** 8-12 point dose-response per lead, replicates across plates, fitted with a 4-parameter logistic (4PL): `y = bottom + (top - bottom) / (1 + (EC50/x)^hill)`. Efficacy is `top` (Emax), potency is EC50.
- **Confounders to design against:** plate/batch offset, edge effect, row/column drift, time-on-deck and dispense-order effects, tip carryover, DMSO tolerance, assay-interference compounds (fluorescent, aggregators), cell passage/day. Each is something the layout or the power calc must account for, not just ignore.

### User flow
1. **Simulate a primary screen** (hit rate, activity distribution, plate noise are inputs), download as CSV. Uploading a real CSV with the same columns replaces it.
2. **Set constraints:** number of top leads, plates available (96-well), controls per plate, dose points, top concentration and dilution factor, liquid handler limits (guess: max dispenses per run, max plates per run, minutes per plate), alpha, target power.
3. **Service proposes the design** and the dashboard shows it (below).

### Success criterion
The design is acceptable when there is enough power (default 80%, alpha 0.05 after Bonferroni) to distinguish the efficacy of the leads from each other. The binding comparison is the closest pair of leads in the ranking (the smallest expected |d|), since if that pair is resolved the wider gaps are too. The budget solver picks replicates to meet this; if the plates or handler limits cannot, the dashboard reports the achievable power and the smallest difference that is resolvable instead.

### Design logic
1. **Lead selection:** rank by % activity after QC (drop plates with Z' < 0.5 and flagged interferers); take the top N.
2. **Priors for the 4PL:** the single dose gives only top-ish information, so assume hill = 1, bottom = control level, top from the screen activity, and EC50 uncertain over about 3 logs. The dose series is a half-log ladder spanning that range so the design is robust to a wrong EC50 guess.
3. **Noise:** SD from vehicle/positive wells and replicate spread; plate variance as its own component.
4. **Budget:** wells per lead = doses x replicates; plates needed = ceil(total wells / usable wells per plate). The replicate count is the largest that fits the plates and handler limits, or the smallest that reaches target power. If the two disagree the dashboard flags the design infeasible and shows what to relax.
5. **Efficacy power (primary goal = tell compounds apart on efficacy):** for each pair of leads, the effect size is `d = (top_i - top_j) / sigma`. Variance of each fitted top comes from the 4PL Jacobian, `Var(theta) = sigma^2 (J'J)^-1`, evaluated at the planned doses and replicates. Test is a Wald z with Bonferroni over the pairs. Reported: power for each pair, the minimum detectable d at target power, and the expected d for the simulated leads. EC50 differences are reported as a secondary result (on the log scale).
6. **Check by simulation:** repeat simulate -> fit 4PL -> test many times (Monte Carlo), which confirms the analytic power, catches fits that fail to converge, and gives hit recovery and EC50 error against the simulated truth.
7. **Layout (96 wells):** plate is the block. Each lead's full dose series appears on every plate it is on, with replicates spread across plates so plate effects cancel in the comparison. Vehicle (max) and positive (min) controls are on every plate, split between left and right columns to expose drift. Edge wells hold buffer or vehicle only. Positions are randomised within plate (seeded); the dispense order runs low to high dose within a compound to limit carryover.
8. **Confounder panel (interactive):** toggles for edge effect, plate offset SD, drift and carryover. Each re-runs the power calc with inflated variance or the simulation with that effect on, and shows how much power it costs and which layout choice protects against it.

### Dashboard
- Inputs panel and the constraint check (wells, plates, handler time against limits).
- Lead table: simulated screen activity, planned doses, replicates, expected d, power.
- Power: pairwise power matrix, power vs replicates, minimum detectable d.
- Fitted 4PL curves with confidence bands from the simulation.
- Plate view: per-plate tabs of the 96-well layout (compound, dose, control, buffer), coloured as before, plus a distribution summary (wells per lead and per dose per plate).
- Downloads: layout CSV (plate, well, compound, concentration, role, dispense order), a liquid handler worklist CSV (source, destination, volume; format to be matched to the instrument), and a design summary (JSON, and the screen CSV).

### Endpoints
`POST /api/screen` (simulate primary screen, returns CSV rows); `POST /api/design` (screen rows + constraints -> leads, doses, replicates, power, layout); `GET /api/download/<layout|worklist|summary>`. The existing `POST /api/plan` stays for the generic replicate calculator.

### Build order
(a) screen simulator + CSV, (b) 4PL fit (scipy `curve_fit`) and Jacobian power, (c) lead selection and budget solver, (d) blocked layout + downloads, (e) dashboard wiring, (f) confounder toggles, (g) Monte Carlo check.

### Calls I made (to confirm or overturn)
- "Enough power to distinguish the efficacy of compounds" (confirmed) is taken as Emax differences between leads, checked on the closest-ranked pair, EC50 secondary.
- Wald test with Bonferroni, backed by Monte Carlo, rather than a full mixed model.
- Liquid handler limits are a plain count/time budget, not a specific instrument's worklist format.
- 96-well only for this flow (the generic calculator still supports 144 and 384).

## Known limits / next steps
- Pairwise t-test vs control, not Dunnett; wells assumed independent with equal SD.
- Randomisation is unblocked: no plate/row/column effects, no true DOE designs (fractional factorial, response surface).
- Layout is random, not optimised (e.g. no balancing of compounds across plates, no edge-effect modelling).
- Next, if needed: Dunnett power, balanced block design, CSV export of the layout (power now also checked by Monte Carlo in the phase 2 section).

## Run
```
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python app.py   # runs the __main__ checks, then http://localhost:5000
# Windows: .venv\Scripts\python app.py
```
