"""Replicate-number planner: DOE screen proposal + plate layout + simulated screen. Run: python app.py"""
import math
import random

import numpy as np
from flask import Flask, jsonify, render_template_string, request
from scipy.stats import nct, t

app = Flask(__name__)

# rows x cols
PLATES = {"96": (8, 12), "144": (9, 16), "384": (16, 24)}


def power(n, n_ctrl, delta, sd, alpha):
    """Two-sided t-test, treatment (n) vs pooled vehicle control (n_ctrl)."""
    df = n + n_ctrl - 2
    ncp = delta / (sd * math.sqrt(1 / n + 1 / n_ctrl))
    tc = t.ppf(1 - alpha / 2, df)
    # lower tail via symmetry, P(T < -tc | ncp) = P(T > tc | -ncp): nct.cdf(-tc) returns NaN for some df/ncp
    return float(nct.sf(tc, df, ncp) + nct.sf(tc, df, -ncp))


def layout(conds, n, plates, rows, cols, edge, ctrl_pp, seed):
    """Randomised block fill: shuffle treatment wells, spread over plates, controls on every plate."""
    wells = [(r, c) for r in range(rows) for c in range(cols)]
    if edge:
        wells = [(r, c) for r, c in wells if 0 < r < rows - 1 and 0 < c < cols - 1]
    rng = random.Random(seed)
    items = [i for i in range(len(conds)) for _ in range(n)]
    rng.shuffle(items)
    per = len(wells) - ctrl_pp
    out = []
    for p in range(plates):
        chunk = items[p * per:(p + 1) * per] + [-1] * ctrl_pp  # -1 = vehicle control
        rng.shuffle(chunk)
        out.append([{"r": r, "c": c, "k": k} for (r, c), k in zip(wells, chunk)])
    return out


def concs(a):
    """Dilution series in µM, lowest first: top / dil^(doses-1) ... top."""
    top, dil, d = float(a["top"]), float(a["dil"]), int(a["doses"])
    if top <= 0 or dil <= 1:
        raise ValueError("top concentration must be > 0 and dilution factor > 1")
    return [top / dil ** (d - 1 - i) for i in range(d)]


def plan(a):
    rows, cols = PLATES[a["plate"]]
    comps, doses, targets = int(a["compounds"]), int(a["doses"]), int(a["targets"])
    cs = concs(a)
    # order target > compound > dose, so condition k belongs to compound x target pair k // doses
    conds = [{"compound": f"C{i + 1}", "dose": d + 1, "conc": cs[d], "target": f"T{k + 1}"}
             for k in range(targets) for i in range(comps) for d in range(doses)]
    ctrl_pp = int(a["ctrl_pp"])
    usable = rows * cols - (2 * (rows + cols) - 4 if a["edge"] else 0)
    if usable <= ctrl_pp:
        raise ValueError("no usable wells left for treatments")
    # Bonferroni across all treatment-vs-control comparisons (one family per target)
    k = comps * doses
    alpha = float(a["alpha"]) / (k if a["bonferroni"] else 1)
    curve, n_req = [], None
    for n in range(2, 101):
        plates = math.ceil(len(conds) * n / (usable - ctrl_pp))
        pw = power(n, plates * ctrl_pp, float(a["delta"]), float(a["sd"]), alpha)
        curve.append({"n": n, "power": pw})
        if n_req is None and pw >= float(a["power"]):
            n_req = n
    n_use = int(a["n"]) if a.get("n") else (n_req or 100)
    plates = math.ceil(len(conds) * n_use / (usable - ctrl_pp))
    return {
        "conds": conds, "n_required": n_req, "n_used": n_use, "curve": curve,
        "alpha_adj": alpha, "plates": plates, "usable": usable,
        "treat_wells": len(conds) * n_use, "ctrl_wells": plates * ctrl_pp,
        "rows": rows, "cols": cols,
        "achieved_power": curve[min(n_use, 100) - 2]["power"],
        "layout": layout(conds, n_use, plates, rows, cols, a["edge"], ctrl_pp, int(a["seed"])),
    }


def sim_params(a):
    sd = float(a["sd"])
    sd_plate = a.get("sd_plate")
    return {"base": float(a.get("baseline") or 100), "sd": sd, "delta": float(a["delta"]),
            "sd_plate": 0.5 * sd if sd_plate in (None, "") else float(sd_plate),
            "cv": bool(a.get("cv")), "edge_grad": float(a.get("edge_grad") or 0),
            "hit_rate": float(a.get("hit_rate", 0.25))}


def draw_truth(rng, size, lo, hi, delta, hit_rate):
    """Per compound x target: active flag, EC50 log-uniform from half a log below to half a log above
    the tested range (so some hits are partial or invisible), Hill slope, top (0 if inactive)."""
    active = rng.random(size) < hit_rate
    ec50 = 10 ** rng.uniform(math.log10(lo) - 0.5, math.log10(hi) + 0.5, size)
    hill = rng.uniform(0.7, 2, size)
    top = np.where(active, rng.uniform(0.5, 2, size) * delta, 0.0)
    return active, ec50, hill, top


def hill_resp(conc, ec50, h, top):
    return top * conc ** h / (ec50 ** h + conc ** h)


def edge_effect(r, c, rows, cols, grad):
    """grad at the outer ring, halving with each ring inward."""
    return grad * 0.5 ** np.minimum(np.minimum(r, c), np.minimum(rows - 1 - r, cols - 1 - c))


def observe(mean, rng, p):
    """Well noise on the mean signal: additive N(0, sd), or multiplicative with CV = sd / baseline."""
    e = rng.standard_normal(np.shape(mean))
    return mean * (1 + e * p["sd"] / p["base"]) if p["cv"] else mean + e * p["sd"]


def simulate(a, pl):
    """One seeded screen on the plan's layout: long table, truth table, per-plate values for the heatmap."""
    p, rng = sim_params(a), np.random.default_rng(int(a["sim_seed"]))
    doses, conds = int(a["doses"]), pl["conds"]
    cs = np.array([c["conc"] for c in conds])
    active, ec50, hill, top = draw_truth(rng, len(conds) // doses, cs.min(), cs.max(), p["delta"], p["hit_rate"])
    plate_off = rng.normal(0, p["sd_plate"], pl["plates"])
    pi, r, c, k = (np.array(v) for v in zip(*[(i, w["r"], w["c"], w["k"])
                                             for i, plate in enumerate(pl["layout"]) for w in plate]))
    j = np.maximum(k, 0) // doses  # compound x target pair (masked out for vehicle)
    resp = np.where(k >= 0, hill_resp(np.where(k >= 0, cs[k], 0.0), ec50[j], hill[j], top[j]), 0.0)
    y = observe(p["base"] + plate_off[pi] + edge_effect(r, c, pl["rows"], pl["cols"], p["edge_grad"]) + resp, rng, p)
    vehicle = {"compound": "vehicle", "target": "", "dose": 0, "conc_uM": 0.0}
    rows = [{"plate": int(i) + 1, "well": f"{chr(65 + rr)}{cc + 1}",
             **({"compound": conds[kk]["compound"], "target": conds[kk]["target"], "dose": conds[kk]["dose"],
                 "conc_uM": conds[kk]["conc"]} if kk >= 0 else vehicle),
             "value": round(float(v), 3)} for i, rr, cc, kk, v in zip(pi, r, c, k, y)]
    truth = [{"compound": conds[q * doses]["compound"], "target": conds[q * doses]["target"],
              "active": bool(active[q]), "ec50": float(ec50[q]), "hill": float(hill[q]), "top": float(top[q])}
             for q in range(len(active))]
    split = np.cumsum([len(plate) for plate in pl["layout"]])[:-1]
    return {"rows": rows, "truth": truth, "plate_offsets": plate_off.round(3).tolist(),
            "values": [v.round(3).tolist() for v in np.split(y, split)]}


def mc(a, pl, ns, reps):
    """Monte Carlo of the analysis the power calc assumes: each condition vs pooled vehicle, pooled t-test
    at the adjusted alpha, with the simulator's nuisance (plate offsets, edge, noise model). Per n returns
    power for a probe with effect exactly delta (checks the analytic curve), and with a fresh truth table
    per rep, hit recovery (true active pairs with >= 1 significant dose) and false-positive rate.
    Wells land on random plates/positions rather than the exact layout."""
    p, rng = sim_params(a), np.random.default_rng(int(a["sim_seed"]))
    rows, cols, ctrl_pp, doses = pl["rows"], pl["cols"], int(a["ctrl_pp"]), int(a["doses"])
    rr, cc = np.array([(r, c) for r in range(rows) for c in range(cols)
                       if not a["edge"] or (0 < r < rows - 1 and 0 < c < cols - 1)]).T
    ev = edge_effect(rr, cc, rows, cols, p["edge_grad"])
    cs = np.array([c["conc"] for c in pl["conds"]])
    K, j, out = len(cs), np.arange(len(cs)) // doses, []
    for n in ns:
        P = math.ceil(K * n / (pl["usable"] - ctrl_pp))
        m = P * ctrl_pp
        off = rng.normal(0, p["sd_plate"], (reps, P))

        def wells(pidx):  # baseline + plate offset + edge at a random usable position
            return (p["base"] + np.take_along_axis(off, pidx.reshape(reps, -1), 1).reshape(pidx.shape)
                    + ev[rng.integers(0, len(ev), pidx.shape)])
        ctrl = observe(wells(np.tile(np.repeat(np.arange(P), ctrl_pp), (reps, 1))), rng, p)
        active, ec50, hill, top = draw_truth(rng, (reps, K // doses), cs.min(), cs.max(), p["delta"], p["hit_rate"])
        resp = np.concatenate([hill_resp(cs, ec50[:, j], hill[:, j], top[:, j]), np.full((reps, 1), p["delta"])], 1)
        y = observe(wells(rng.integers(0, P, (reps, K + 1, n))) + resp[..., None], rng, p)
        sp2 = ((n - 1) * y.var(2, ddof=1) + (m - 1) * ctrl.var(1, ddof=1, keepdims=True)) / (n + m - 2)
        tstat = (y.mean(2) - ctrl.mean(1, keepdims=True)) / np.sqrt(sp2 * (1 / n + 1 / m))
        sig = np.abs(tstat) > t.ppf(1 - pl["alpha_adj"] / 2, n + m - 2)
        hit = sig[:, :-1].reshape(reps, -1, doses).any(2)
        out.append({"n": int(n), "power": float(sig[:, -1].mean()),
                    "recovery": float(hit[active].mean()) if active.any() else None,
                    "fpr": float(hit[~active].mean()) if (~active).any() else None})
    return out


def api(fn):
    try:
        a = request.get_json()
        return jsonify(fn(a, plan(a)))
    except (ValueError, KeyError, ZeroDivisionError) as e:
        return jsonify(error=str(e)), 400


@app.post("/api/plan")
def api_plan():
    return api(lambda a, pl: pl)


@app.post("/api/simulate")
def api_simulate():
    return api(simulate)


@app.post("/api/mc")
def api_mc():
    def run(a, pl):
        xm = min(100, max(10, 2 * (pl["n_required"] or 50)))  # same x range as the page's power chart
        reps = min(5000, int(a.get("reps") or 300))
        return {"reps": reps, "points": mc(a, pl, sorted({int(round(x)) for x in np.linspace(2, xm, 12)}), reps)}
    return api(run)


@app.get("/")
def index():
    return render_template_string(PAGE)


PAGE = """<!doctype html><meta charset=utf-8><meta name=viewport content="width=device-width,initial-scale=1">
<title>Replicate Planner</title>
<style>
:root{--bg:#fff;--fg:#1c1e21;--mut:#6b7280;--card:#f4f5f7;--acc:#2563eb}
@media(prefers-color-scheme:dark){:root{--bg:#15171a;--fg:#e6e8eb;--mut:#9aa1ab;--card:#1f2226;--acc:#60a5fa}}
body{margin:0;font:14px system-ui;background:var(--bg);color:var(--fg);display:grid;grid-template-columns:260px 1fr;min-height:100vh}
@media(max-width:800px){body{grid-template-columns:1fr}}
aside{background:var(--card);padding:16px}main{padding:16px;min-width:0}
label{display:block;margin:8px 0 2px;color:var(--mut)}input,select{width:100%;box-sizing:border-box;padding:4px;font:inherit}
.kpi{display:flex;gap:12px;flex-wrap:wrap}.kpi div{background:var(--card);padding:8px 14px;border-radius:6px}
.kpi b{display:block;font-size:22px}
#plate{display:grid;gap:2px;margin:8px 0;max-width:900px}
#plate div{aspect-ratio:1;border-radius:50%;font-size:9px;display:flex;align-items:center;justify-content:center;color:#000}
button.tab{margin:2px;padding:3px 8px}button.on{background:var(--acc);color:#fff}
svg{width:100%;max-width:500px;background:var(--card);border-radius:6px}
#err{color:#dc2626}.mut{color:var(--mut)}table{border-collapse:collapse;font-size:12px}td,th{padding:2px 8px;text-align:left}
details{margin:4px 0;max-height:300px;overflow:auto}.lg span{display:inline-block;margin:2px 8px 2px 0}.lg i{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:3px}
</style>
<aside><h3>Experiment</h3><form id=f oninput=go()>
<label>Plate<select name=plate><option>96<option>144<option>384</select></label>
<label>Compounds<input name=compounds type=number value=8 min=1></label>
<label>Doses / compound<input name=doses type=number value=4 min=1></label>
<label>Top concentration (µM)<input name=top type=number value=10 step=any></label>
<label>Dilution factor (1:x)<input name=dil type=number value=3 step=any></label>
<label>Targets (assays/cell lines)<input name=targets type=number value=2 min=1></label>
<label>Vehicle ctrl wells / plate<input name=ctrl_pp type=number value=8 min=1></label>
<label><input type=checkbox name=edge checked style="width:auto"> Skip edge wells</label>
<label><input type=checkbox name=bonferroni checked style="width:auto"> Bonferroni (per target)</label>
<h3>Statistics</h3>
<label>Min. effect (Δ, signal units)<input name=delta type=number value=20 step=any></label>
<label>Well-to-well SD<input name=sd type=number value=15 step=any></label>
<label>α<input name=alpha type=number value=0.05 step=any></label>
<label>Target power<input name=power type=number value=0.8 step=any></label>
<label>Override replicates (blank = auto)<input name=n type=number min=2></label>
<label>Randomisation seed<input name=seed type=number value=1></label>
<h3>Simulation</h3>
<label>Hit rate (active compound × target)<input name=hit_rate type=number value=0.25 step=any min=0 max=1></label>
<label>Baseline signal<input name=baseline type=number value=100 step=any></label>
<label>Plate offset SD (blank = 0.5 × well SD)<input name=sd_plate type=number step=any min=0></label>
<label><input type=checkbox name=cv style="width:auto"> Multiplicative noise (CV = SD / baseline)</label>
<label>Edge gradient, outer ring (0 = off)<input name=edge_grad type=number value=0 step=any></label>
<label>Simulation seed<input name=sim_seed type=number value=1></label>
<label>Monte Carlo reps<input name=reps type=number value=300 min=50 max=5000></label>
</form></aside>
<main><h2>How many replicates do we need?</h2><div id=err></div>
<div class=kpi id=kpi></div>
<h3>Power vs replicates <button type=button onclick=runMC()>Run Monte Carlo</button> <span id=mcs></span></h3>
<svg id=pc viewBox="0 0 500 220"></svg><div class=lg id=pcl></div>
<h3>Simulated screen <button type=button onclick="view='values';sim()">Simulate screen</button>
 <button type=button id=dl onclick=csv() disabled>Download CSV</button></h3><div id=simsum></div>
<h3>Plate <span id=views></span> <span id=tabs></span></h3><div class=lg id=lg></div><div id=plate></div></main>
<script>
const f=document.getElementById('f'),$=id=>document.getElementById(id);let timer,cur=0,D,S,M,view='layout';
const fmt=x=>+(+x).toPrecision(2);
function args(){const o={};for(const e of f.elements){if(!e.name)continue;
 o[e.name]=e.type=='checkbox'?e.checked:e.value}return o}
async function post(u){const r=await fetch(u,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(args())});
 const j=await r.json();$('err').textContent=j.error||'';return j.error?null:j}
function go(){clearTimeout(timer);timer=setTimeout(async()=>{
 const j=await post('/api/plan');if(!j)return;D=j;M=null;cur=Math.min(cur,D.layout.length-1);
 if(S)sim();else draw()},200)}
async function sim(){const j=await post('/api/simulate');if(!j)return;S=j;
 const v=S.values.flat();S.lo=Math.min(...v);S.hi=Math.max(...v);draw()}
async function runMC(){$('mcs').textContent='running…';const j=await post('/api/mc');$('mcs').textContent='';if(j){M=j;draw()}}
function csv(){const h=['plate','well','compound','target','dose','conc_uM','value'];
 const s=[h.join(','),...S.rows.map(r=>h.map(k=>r[k]).join(','))].join('\\n');
 const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([s],{type:'text/csv'}));
 a.download=`screen_seed${args().sim_seed}.csv`;a.click()}
function col(k){if(k<0)return'#9ca3af';const c=D.conds[k],ci=+c.compound.slice(1),nd=+args().doses;
 return `hsl(${(ci*47)%360},65%,${80-55*(c.dose-1)/Math.max(1,nd-1)}%)`}
function heat(v){return `hsl(25,85%,${92-62*(v-S.lo)/((S.hi-S.lo)||1)}%)`}
function draw(){
 $('kpi').innerHTML=[['Replicates / condition',D.n_required??'>100'],['Using',D.n_used],
  ['Achieved power',(D.achieved_power*100).toFixed(0)+'%'],['Conditions',D.conds.length],
  ['Treatment wells',D.treat_wells],['Control wells',D.ctrl_wells],['Plates',D.plates],['Adj. α',D.alpha_adj.toExponential(1)]]
  .map(([a,b])=>`<div>${a}<b>${b}</b></div>`).join('');
 // power curve
 const pts=D.curve.filter(p=>p.n<=Math.min(100,Math.max(10,(D.n_required||50)*2)));const xm=pts.at(-1).n;
 const X=n=>40+(n-2)/(xm-2)*440,Y=p=>200-p*180;const tg=+args().power;
 $('pc').innerHTML=`<polyline fill=none stroke="var(--acc)" stroke-width=2 points="${pts.map(p=>X(p.n)+','+Y(p.power)).join(' ')}"/>
 <line x1=40 x2=480 y1=${Y(tg)} y2=${Y(tg)} stroke=#dc2626 stroke-dasharray=4 />
 <line x1=${X(D.n_used)} x2=${X(D.n_used)} y1=20 y2=200 stroke="var(--mut)" stroke-dasharray=2 />
 <text x=40 y=214 font-size=10 fill=currentColor>replicates 2…${xm}</text><text x=4 y=24 font-size=10 fill=currentColor>1.0</text><text x=4 y=200 font-size=10 fill=currentColor>0</text>
 ${M?M.points.filter(p=>p.n<=xm).map(p=>`<circle cx=${X(p.n)} cy=${Y(p.power)} r=3.5 fill="var(--acc)"><title>n=${p.n}: MC power at Δ ${(p.power*100).toFixed(1)}%</title></circle>`
  +(p.recovery==null?'':`<circle cx=${X(p.n)} cy=${Y(p.recovery)} r=3.5 fill=#ea580c><title>n=${p.n}: hit recovery ${(p.recovery*100).toFixed(1)}%</title></circle>`)
  +(p.fpr==null?'':`<rect x=${X(p.n)-3} y=${Y(p.fpr)-3} width=6 height=6 fill=#16a34a><title>n=${p.n}: false-positive rate ${(p.fpr*100).toFixed(1)}%</title></rect>`)).join(''):''}`;
 $('pcl').innerHTML='<span><i style="background:var(--acc);border-radius:0;height:2px"></i>analytic power at Δ</span><span><i style="background:#dc2626;border-radius:0;height:2px"></i>target</span>'
  +(M?`<span><i style="background:var(--acc)"></i>MC power at Δ</span><span><i style="background:#ea580c"></i>MC hit recovery (true actives)</span><span><i style="background:#16a34a;border-radius:0"></i>MC false-positive rate</span><span>${M.reps} reps / point, includes plate offset + edge + noise model</span>`:'');
 // simulated screen summary + truth table
 $('dl').disabled=!S;
 if(S){const act=S.truth.filter(q=>q.active),cs=D.conds.map(c=>c.conc),lo=Math.min(...cs),hi=Math.max(...cs);
  const out=act.filter(q=>q.ec50<lo||q.ec50>hi).length;
  $('simsum').innerHTML=`<p>${S.rows.length} wells · ${act.length} of ${S.truth.length} compound × target pairs active · ${out} with EC50 outside the tested ${fmt(lo)}–${fmt(hi)} µM range · plate offsets ${S.plate_offsets.map(x=>x.toFixed(1)).join(', ')}</p>
  <details><summary>Truth table</summary><table><tr><th>compound<th>target<th>active<th>EC50 µM<th>Hill<th>top</tr>${S.truth.map(q=>`<tr><td>${q.compound}<td>${q.target}<td>${q.active?'yes':''}<td>${q.active?fmt(q.ec50):''}<td>${q.active?q.hill.toFixed(2):''}<td>${q.top.toFixed(1)}</tr>`).join('')}</table></details>`}
 else $('simsum').innerHTML='<p class=mut>Not simulated yet.</p>';
 const hv=S&&view=='values';
 $('views').innerHTML=S?['layout','values'].map(v=>`<button type=button class="tab ${v==view?'on':''}" onclick="view='${v}';draw()">${v=='values'?'simulated values':'layout'}</button>`).join(''):'';
 $('tabs').innerHTML=D.layout.map((_,i)=>`<button type=button class="tab ${i==cur?'on':''}" onclick="cur=${i};draw()">${i+1}</button>`).join('');
 const seen=new Set(D.conds.map(c=>c.compound));
 $('lg').innerHTML=hv?`<span>${S.lo.toFixed(0)}</span><i style="width:160px;height:10px;border-radius:0;background:linear-gradient(90deg,${heat(S.lo)},${heat(S.hi)})"></i><span>${S.hi.toFixed(0)} signal</span>`
  :[...seen].map(c=>`<span><i style="background:hsl(${(+c.slice(1)*47)%360},65%,55%)"></i>${c}</span>`).join('')+'<span><i style="background:#9ca3af"></i>vehicle</span><span>darker = higher dose</span>';
 const g=$('plate');g.style.gridTemplateColumns=`repeat(${D.cols+1},1fr)`;
 const m={};D.layout[cur].forEach((w,j)=>m[w.r+','+w.c]=[w.k,j]);
 let h='<div></div>'+[...Array(D.cols).keys()].map(c=>`<div>${c+1}</div>`).join('');
 for(let r=0;r<D.rows;r++){h+=`<div>${String.fromCharCode(65+r)}</div>`;
  for(let c=0;c<D.cols;c++){if(!m[r+','+c]){h+='<div style="border:1px dashed var(--mut)"></div>';continue}
   const [k,j]=m[r+','+c],q=k<0?null:D.conds[k],v=hv?S.values[cur]?.[j]:undefined;
   const t=q?`${q.compound} ${fmt(q.conc)} µM / ${q.target}`:'vehicle control';
   h+=v===undefined?`<div style="background:${col(k)}" title="${t}">${q?q.target:'V'}</div>`
    :`<div style="background:${heat(v)};color:${(v-S.lo)/((S.hi-S.lo)||1)>.6?'#fff':'#000'}" title="${t}: ${v.toFixed(1)}">${q?'':'V'}</div>`}}
 g.innerHTML=h}
go()
</script>"""

if __name__ == "__main__":
    a = dict(plate="96", compounds=8, doses=4, targets=2, ctrl_pp=8, edge=True, bonferroni=True,
             delta=20, sd=15, alpha=0.05, power=0.8, n="", seed=1, top=10, dil=3,
             hit_rate=0.25, baseline=100, sd_plate="", cv=False, edge_grad=0, sim_seed=1)
    p = plan(a)
    assert p["n_required"] and p["achieved_power"] >= 0.8 and p["plates"] >= 1
    assert all(sum(w["k"] >= 0 for w in pl) <= p["usable"] for pl in p["layout"])
    assert sum(w["k"] >= 0 for pl in p["layout"] for w in pl) == p["treat_wells"]
    assert abs(p["conds"][3]["conc"] - 10) < 1e-9 and abs(p["conds"][0]["conc"] - 10 / 27) < 1e-9

    # simulation: seeded, one row per filled well, EC50 within half a log of the tested range
    s = simulate(a, p)
    assert s == simulate(a, p) and s != simulate(dict(a, sim_seed=2), p)
    assert len(s["rows"]) == sum(map(len, p["layout"])) == sum(map(len, s["values"]))
    assert all(10 / 27 / 10 ** 0.5 <= q["ec50"] <= 10 * 10 ** 0.5 for q in s["truth"])
    assert all((q["top"] > 0) == q["active"] for q in s["truth"])
    # no plate effect: vehicle SD ~ well SD, inactive compounds flat at baseline
    b = dict(a, plate="384", compounds=40, ctrl_pp=32, n=10, sd_plate=0)
    pb = plan(b)
    sb = simulate(b, pb)
    veh = np.array([r["value"] for r in sb["rows"] if r["compound"] == "vehicle"])
    assert abs(veh.std(ddof=1) / 15 - 1) < 0.1 and abs(veh.mean() - 100) < 2
    off = {(q["compound"], q["target"]) for q in sb["truth"] if not q["active"]}
    flat = np.array([r["value"] for r in sb["rows"] if (r["compound"], r["target"]) in off])
    assert abs(flat.mean() - 100) < 1
    assert abs(np.mean([q["active"] for q in sb["truth"]]) - 0.25) < 0.1
    # multiplicative noise with CV = sd / baseline gives the same vehicle SD
    vc = np.array([r["value"] for r in simulate(dict(b, cv=True), pb)["rows"] if r["compound"] == "vehicle"])
    assert abs(vc.std(ddof=1) / 15 - 1) < 0.1
    # edge gradient: outer ring sits ~grad above the interior
    e = dict(b, edge=False, edge_grad=40, hit_rate=0)
    ring = [(min(ord(r["well"][0]) - 65, int(r["well"][1:]) - 1, 15 - (ord(r["well"][0]) - 65),
                 23 - (int(r["well"][1:]) - 1)), r["value"]) for r in simulate(e, plan(e))["rows"]]
    outer = np.mean([v for d, v in ring if d == 0]) - np.mean([v for d, v in ring if d >= 3])
    assert 30 < outer < 45, outer

    # Monte Carlo under the analytic assumptions reproduces the analytic power at n_required
    m = mc(dict(a, sd_plate=0), p, [p["n_required"]], 4000)[0]
    assert abs(m["power"] - p["achieved_power"]) < 0.03, (m, p["achieved_power"])
    print("ok", p["n_required"], p["plates"], "analytic", round(p["achieved_power"], 3), "mc", m,
          "with plate SD", mc(a, p, [p["n_required"]], 4000)[0])
    app.run(debug=True, port=5000)
