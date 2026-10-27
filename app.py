"""Replicate-number planner: DOE screen proposal + plate layout. Run: python app.py"""
import math
import random

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
    return float(1 - nct.cdf(tc, df, ncp) + nct.cdf(-tc, df, ncp))


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


def plan(a):
    rows, cols = PLATES[a["plate"]]
    comps, doses, targets = int(a["compounds"]), int(a["doses"]), int(a["targets"])
    conds = [{"compound": f"C{i + 1}", "dose": d + 1, "target": f"T{k + 1}"}
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


@app.post("/api/plan")
def api_plan():
    try:
        return jsonify(plan(request.get_json()))
    except (ValueError, KeyError, ZeroDivisionError) as e:
        return jsonify(error=str(e)), 400


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
#err{color:#dc2626}.lg span{display:inline-block;margin:2px 8px 2px 0}.lg i{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:3px}
</style>
<aside><h3>Experiment</h3><form id=f oninput=go()>
<label>Plate<select name=plate><option>96<option>144<option>384</select></label>
<label>Compounds<input name=compounds type=number value=8 min=1></label>
<label>Doses / compound<input name=doses type=number value=4 min=1></label>
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
</form></aside>
<main><h2>How many replicates do we need?</h2><div id=err></div>
<div class=kpi id=kpi></div>
<h3>Power vs replicates</h3><svg id=pc viewBox="0 0 500 220"></svg>
<h3>Plate layout <span id=tabs></span></h3><div class=lg id=lg></div><div id=plate></div></main>
<script>
const f=document.getElementById('f'),$=id=>document.getElementById(id);let timer,cur=0,D;
function args(){const o={};for(const e of f.elements){if(!e.name)continue;
 o[e.name]=e.type=='checkbox'?e.checked:e.value}return o}
function go(){clearTimeout(timer);timer=setTimeout(async()=>{
 const r=await fetch('/api/plan',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(args())});
 const j=await r.json();$('err').textContent=j.error||'';if(j.error)return;D=j;cur=0;draw()},200)}
function col(k){if(k<0)return'#9ca3af';const c=D.conds[k],ci=+c.compound.slice(1),nd=+args().doses;
 return `hsl(${(ci*47)%360},65%,${80-55*(c.dose-1)/Math.max(1,nd-1)}%)`}
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
 <text x=40 y=214 font-size=10 fill=currentColor>replicates 2…${xm}</text><text x=4 y=24 font-size=10 fill=currentColor>1.0</text><text x=4 y=200 font-size=10 fill=currentColor>0</text>`;
 $('tabs').innerHTML=D.layout.map((_,i)=>`<button type=button class="tab ${i==cur?'on':''}" onclick="cur=${i};draw()">${i+1}</button>`).join('');
 const nc=D.conds.length,seen=new Set(D.conds.map(c=>c.compound));
 $('lg').innerHTML=[...seen].map(c=>`<span><i style="background:hsl(${(+c.slice(1)*47)%360},65%,55%)"></i>${c}</span>`).join('')+'<span><i style="background:#9ca3af"></i>vehicle</span><span>darker = higher dose</span>';
 const g=$('plate');g.style.gridTemplateColumns=`repeat(${D.cols+1},1fr)`;
 const m={};D.layout[cur].forEach(w=>m[w.r+','+w.c]=w.k);
 let h='<div></div>'+[...Array(D.cols).keys()].map(c=>`<div>${c+1}</div>`).join('');
 for(let r=0;r<D.rows;r++){h+=`<div>${String.fromCharCode(65+r)}</div>`;
  for(let c=0;c<D.cols;c++){const k=m[r+','+c];
   if(k===undefined){h+='<div style="border:1px dashed var(--mut)"></div>';continue}
   const q=k<0?null:D.conds[k];
   h+=`<div style="background:${col(k)}" title="${q?`${q.compound} dose ${q.dose} / ${q.target}`:'vehicle control'}">${q?q.target:'V'}</div>`}}
 g.innerHTML=h}
go()
</script>"""

if __name__ == "__main__":
    a = dict(plate="96", compounds=8, doses=4, targets=2, ctrl_pp=8, edge=True, bonferroni=True,
             delta=20, sd=15, alpha=0.05, power=0.8, n="", seed=1)
    p = plan(a)
    assert p["n_required"] and p["achieved_power"] >= 0.8 and p["plates"] >= 1
    assert all(sum(w["k"] >= 0 for w in pl) <= p["usable"] for pl in p["layout"])
    assert sum(w["k"] >= 0 for pl in p["layout"] for w in pl) == p["treat_wells"]
    print("ok", p["n_required"], p["plates"])
    app.run(debug=True, port=5000)
