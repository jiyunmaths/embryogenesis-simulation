/* Live views of actual solver snapshots; rendering never advances the model. */
"use strict";
const $ = id => document.getElementById(id);
const groups = [
  ["Run & resolution", ["seed", "steps", "max_cells", "save_every", "grid", "extent", "dt", "interface_width"]],
  ["Activator–inhibitor signaling", ["signaling", "signal_beta", "signal_da", "signal_dh", "signal_partition_noise", "signal_fate_gain", "graph_contact_cutoff"]],
  ["Division & cytokinesis", ["division_interval", "cycle_jitter", "division_orientation", "axis_degeneracy", "cytokinesis_duration", "ring_strength", "neck_threshold", "division_overlap_tolerance"]],
  ["Cell identity", ["differentiation", "competence_cells", "fate_rate", "fate_threshold", "neighbor_inhibition", "exposure_bias", "fate_noise", "partition_noise"]],
  ["Mechanical feedback", ["feedback", "surface_tension", "volume_stiffness", "repulsion", "adhesion", "fate_adhesion", "fate_tension"]],
  ["Apical–basal polarity", ["polarity_enabled", "polarity_rate", "polarity_alignment", "polarity_decay", "polarity_tension"]]
];
const labels = {seed:"Random seed", steps:"Total steps", max_cells:"Cell limit", save_every:"Record every (steps)", grid:"Grid per axis", extent:"Domain half-width", dt:"Time step", interface_width:"Interface width", signaling:"Enable signaling", signal_beta:"Inhibitor reaction rate", signal_da:"Activator coupling", signal_dh:"Inhibitor coupling", signal_partition_noise:"Signal partition noise", signal_fate_gain:"Signal → identity gain", graph_contact_cutoff:"Contact cutoff", division_interval:"Division interval", cycle_jitter:"Cycle jitter", division_orientation:"Spindle orientation", axis_degeneracy:"Axis degeneracy", cytokinesis_duration:"Cytokinesis duration", ring_strength:"Ring strength", neck_threshold:"Neck threshold", division_overlap_tolerance:"Overlap tolerance", differentiation:"Enable differentiation", competence_cells:"Competence cell count", fate_rate:"Identity response rate", fate_threshold:"Identity threshold", neighbor_inhibition:"Neighbor inhibition", exposure_bias:"Exposure bias", fate_noise:"Identity noise", partition_noise:"Identity partition noise", feedback:"Enable mechanical feedback", surface_tension:"Surface tension", volume_stiffness:"Volume stiffness", repulsion:"Cell repulsion", adhesion:"Cell adhesion", fate_adhesion:"Identity adhesion contrast", fate_tension:"Identity tension contrast", polarity_enabled:"Enable polarity", polarity_rate:"Exposure response", polarity_alignment:"Neighbor alignment", polarity_decay:"Polarity decay", polarity_tension:"Directional tension"};
const help = {
  seed:"Reset with the same settings and seed reproduces the same simulation.",
  save_every:"Snapshot spacing, not solver speed. Pause and completion also record a frame.",
  grid:"Dense grid along each of the three axes. Higher resolution costs more memory and time.",
  dt:"Changing dt also changes simulated duration (steps × dt). Numerical stability constraints are checked by the server.",
  interface_width:"Must be at least 0.6 grid spacings. Adjust together with grid and extent.",
  signal_da:"Dimensionless exchange rate on the normalized contact graph, not a physical diffusion coefficient.",
  signal_dh:"Dimensionless exchange rate on the normalized contact graph, not a physical diffusion coefficient.",
  division_orientation:"Shape aligns the spindle with the longest cell axis; isotropic is a comparison control.",
  fate_threshold:"Labels for the continuous regulatory state. These do not establish irreversible commitment.",
  feedback:"Enable identity-dependent adhesion/tension and polarity-dependent tension. Disabling this retains polarity dynamics but removes its mechanical effect.",
  competence_cells:"Minimum population at which the identity switch responds."
};
const palette = {teal:"#209d91", gold:"#d89f38", blue:"#70a9df", coral:"#ee957f", gray:"#9caab0"};
let schema, state, frames = [], selected = 0, followLive = true, busy = false, connected = false;
let requestEpoch = 0, camera = {yaw:.65, pitch:.4, zoom:1}, paintPending = false;
let fields = new Map();
const number = (value, digits=2) => Number.isFinite(value) ? value.toFixed(digits) : "—";
const small = value => Number.isFinite(value) ? (value !== 0 && Math.abs(value) < .01 ? value.toExponential(1) : number(value, 3)) : "—";

function buildParameters() {
  $("parameters").replaceChildren();
  groups.forEach(([title, names], index) => {
    const group = document.createElement("details"); group.className = "parameter-group"; group.open = index === 0;
    const summary = document.createElement("summary"); summary.textContent = title; group.append(summary);
    const grid = document.createElement("div"); grid.className = "parameter-grid"; group.append(grid);
    names.forEach(name => {
      const type = schema.types[name], wrapper = document.createElement("div");
      wrapper.className = "parameter-field" + (type === "boolean" ? " boolean wide" : type === "string" ? " wide" : "");
      const label = document.createElement("label"); label.htmlFor = `param-${name}`; label.textContent = labels[name] || name;
      const input = document.createElement(type === "string" ? "select" : "input"); input.id = `param-${name}`; input.name = name;
      input.title = help[name] || name.replaceAll("_", " ");
      if (type === "string") (schema.choices[name] || []).forEach(value => {const option = document.createElement("option"); option.value=value; option.textContent=value === "shape" ? "Longest cell axis" : "Isotropic (control)"; input.append(option);});
      else if (type === "boolean") input.type="checkbox";
      else {input.type="number"; input.step=type === "integer" ? "1" : "any"; input.required=true; input.min="0"; if (Array.isArray(schema.limits[name])) {input.min=schema.limits[name][0]; input.max=schema.limits[name][1];}}
      input.addEventListener("input", updateControls); wrapper.append(label,input); grid.append(wrapper); fields.set(name,input);
    });
    $("parameters").append(group);
  });
  $("parameters").addEventListener("submit", event => event.preventDefault());
}
function fillParameters(config) {
  for (const [name, input] of fields) {if(schema.types[name] === "boolean") input.checked=config[name]; else input.value=config[name];}
  updateControls();
}
function readParameters() {
  const config = {};
  for (const [name, input] of fields) {
    const type = schema.types[name]; let value = type === "boolean" ? input.checked : type === "string" ? input.value : Number(input.value);
    const valid = type === "boolean" || type === "string" || (input.value !== "" && Number.isFinite(value) && input.checkValidity() && (type !== "integer" || Number.isInteger(value)));
    input.setAttribute("aria-invalid", String(!valid)); if (!valid) return null;
    config[name] = value;
  }
  return config;
}
function sameConfig(a,b) {return !!a && !!b && Object.keys(a).every(key => a[key] === b[key]);}
function updateControls() {
  const config = schema ? readParameters() : null, dirty = state && config && !sameConfig(config,state.config);
  const unavailable = !connected || !state || busy, running = state?.state === "running";
  for (const input of fields.values()) input.disabled=unavailable || running;
  $("run").disabled = unavailable || !config || running || ["completed","error"].includes(state?.state) || (state?.state === "paused" && dirty);
  $("run-label").textContent = state?.state === "paused" ? "Resume" : "Run";
  $("pause").disabled = unavailable || !running;
  $("reset").disabled = unavailable || !config;
  $("defaults").disabled = unavailable || running;
  $("download").disabled = !frames.length;
  $("parameter-note").classList.toggle("dirty", !!dirty);
  $("parameter-note").textContent = !config && schema ? "Enter valid numbers before starting or resetting." : dirty ? (state.state === "ready" ? "Settings changed. Run or Reset will apply them." : "Settings changed. Reset applies them and returns to a zygote.") : running ? "The solver is running. Pause to edit parameters or inspect its current state." : "Parameters apply when you start a run. Reset returns to a single cell using these settings.";
}
async function api(path, body) {
  const response = await fetch(path, body === undefined ? {cache:"no-store"} : {method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(body)});
  const data = await response.json(); if (!response.ok) throw new Error(data.error || `Request failed (${response.status})`); return data;
}
function showError(message) {$("error").hidden=!message; $("error").textContent=message || "";}
function acceptState(next) {
  if (state && next.generation < state.generation) return;
  const fresh = !state || next.generation !== state.generation;
  const oldSelected = frames[selected]?.revision;
  if (fresh) {frames=[]; selected=0; followLive=true; showError(null);}
  const known = new Set(frames.map(frame=>frame.revision));
  for (const frame of next.frames) if(!known.has(frame.revision)) {frames.push(frame); known.add(frame.revision);}
  frames = frames.filter(frame => frame.revision >= next.history_start_revision).sort((a,b)=>a.revision-b.revision);
  state=next; connected=true;
  if (fresh) {fillParameters(next.config); const extent=next.config.extent; $("cut").min=-extent; $("cut").max=extent; $("cut").step=extent/100; $("cut").value=extent;}
  if (followLive) selected=Math.max(0,frames.length-1);
  else selected=Math.max(0,frames.findIndex(frame=>frame.revision===oldSelected));
  if (next.error) showError(next.error);
  renderStatus(); updateControls(); schedulePaint();
}
async function command(action) {
  const config = readParameters(); if (action !== "pause" && !config) return;
  busy=true; ++requestEpoch; updateControls(); showError(null);
  $("connection").textContent=action === "pause" ? "Finishing current step…" : action === "reset" ? "Preparing zygote…" : "Starting solver…";
  try {const next=await api(`/api/${action}`, action === "pause" ? {} : {config}); acceptState(next);}
  catch(error) {showError(error.message);}
  finally {busy=false; ++requestEpoch; updateControls(); if(state) renderStatus();}
}
async function poll() {
  const epoch=requestEpoch;
  try {
    if (!busy) {
      let next=await api(`/api/state?after=${state?.revision ?? -1}`);
      if (epoch !== requestEpoch) return;
      // A revision can be reused after Reset; fetch complete new history then.
      if(state && next.generation !== state.generation && !next.frames.some(frame=>frame.revision===next.history_start_revision)) next=await api("/api/state?after=-1");
      if (epoch === requestEpoch) acceptState(next);
    }
  } catch(error) {
    if (epoch !== requestEpoch) return;
    connected=false; $("connection").textContent="Disconnected · retrying"; $("state").className="state-pill error"; $("state").replaceChildren(document.createTextNode("Disconnected")); updateControls();
  } finally {setTimeout(poll, state?.state === "running" ? 400 : 900);}
}
function renderStatus() {
  const names={ready:"Ready", running:"Running", paused:"Paused", completed:"Completed", error:"Solver error"};
  $("state").className=`state-pill ${state.state}`; $("state").replaceChildren(document.createElement("i"),document.createTextNode(names[state.state] || state.state));
  $("run-detail").textContent=`Step ${state.step.toLocaleString()} / ${state.total_steps.toLocaleString()} · ${number(state.step*state.config.dt)} simulated time`;
  $("connection").textContent=`Local session · ${number(state.elapsed_seconds,1)} s compute`;
  $("progress").style.width=`${Math.min(100,100*state.step/state.total_steps)}%`;
  $("history-note").textContent=`${frames.length} retained frames · recorded every ${state.config.save_every} steps${state.history_start_revision > 0 ? " · earlier frames pruned" : ""}. Timeline inspection does not rewind the solver.`;
}
function schedulePaint() {if (!paintPending) {paintPending=true; requestAnimationFrame(()=>{paintPending=false; renderFrame();});}}
function renderFrame() {
  const frame=frames[selected]; if(!frame || !state) return;
  const m=frame.metrics; $("scene-empty").hidden=true;
  $("metric-cells").textContent=m.cells; $("metric-cap").textContent=`/ ${state.config.max_cells} cells`;
  $("metric-division").textContent=m.dividing_cells ? `${m.dividing_cells} cell${m.dividing_cells === 1 ? "" : "s"} in cytokinesis` : m.cells === 1 ? "One zygote. An open possibility." : "Progressive cleavage · deformable cells";
  $("metric-time").textContent=number(m.time); $("metric-step").textContent=`Displayed frame: step ${frame.step}`;
  $("metric-signal").textContent=small(m.activator_std); $("metric-modes").textContent=state.config.signaling ? `${m.unstable_graph_modes} growing modes on this contact graph` : "Signaling disabled";
  $("metric-shape").textContent=number(m.axis_ratio);
  $("scene-caption").textContent=`${m.cells === 1 ? "A single zygote" : `${m.cells} cells, a shared geometry`} · t = ${number(m.time)}`;
  $("view-mode").textContent=followLive ? "LATEST FRAME" : "HISTORICAL FRAME";
  $("live").classList.toggle("active",followLive); $("live").setAttribute("aria-pressed",String(followLive));
  $("timeline").max=Math.max(0,frames.length-1); $("timeline").value=selected;
  $("frame-info").textContent=`Frame ${selected+1} / ${frames.length}`;
  $("volume-error").textContent=`${number(100*m.relative_volume_error)}%`; $("polarity-value").textContent=number(m.mean_polarity); $("uncommitted").textContent=m.uncommitted;
  const issues=[];
  if(m.max_cell_volume_error > .05) issues.push(`Largest cell volume error: ${number(100*m.max_cell_volume_error,1)}%.`);
  if(m.boundary_occupancy > .01) issues.push("Cells are approaching the computational boundary.");
  if(m.min_radius_grid_cells < 4) issues.push("Smallest cells are poorly resolved on this grid.");
  if(m.overdue_divisions > 0) issues.push(`${m.overdue_divisions} division(s) delayed by the mechanical criteria.`);
  if(m.clipped_fraction > .01) issues.push("Substantial phase-field clipping; inspect time step and resolution.");
  $("diagnostic-title").textContent=issues.length ? "Inspect numerical limitations" : "Inspect before interpreting";
  $("diagnostic-title").closest("section").classList.toggle("warning",!!issues.length);
  $("diagnostic-note").textContent=issues.length ? issues.join(" ") : "Cell identity labels indicate regulatory tendencies, not proven commitment. Graph growth rates describe a frozen, linearized system.";
  drawScene(frame); drawCharts();
}
function canvasContext(canvas) {
  const rect=canvas.getBoundingClientRect(), ratio=Math.min(window.devicePixelRatio || 1,2);
  const width=Math.max(1,Math.round(rect.width*ratio)), height=Math.max(1,Math.round(rect.height*ratio));
  if(canvas.width !== width || canvas.height !== height) {canvas.width=width; canvas.height=height;}
  const ctx=canvas.getContext("2d"); ctx.setTransform(ratio,0,0,ratio,0,0); ctx.clearRect(0,0,rect.width,rect.height);
  return {ctx,width:rect.width,height:rect.height};
}
function rotate(point) {
  const [x,y,z]=point, cy=Math.cos(camera.yaw), sy=Math.sin(camera.yaw), cp=Math.cos(camera.pitch), sp=Math.sin(camera.pitch);
  const xx=cy*x+sy*z, zz=-sy*x+cy*z; return [xx,cp*y-sp*zz,sp*y+cp*zz];
}
function mix(a,b,t) {return a.map((v,i)=>Math.round(v*(1-t)+b[i]*t));}
function cellColor(cell) {
  const mode=$("color").value;
  if(mode === "lineage") {const h=(cell.id*137.508)%360; return `hsl(${h},57%,69%)`;}
  const value=mode === "fate" ? cell.fate : Math.tanh(Math.log(Math.max(cell[mode],1e-12)));
  return `rgb(${mix([204,219,204],value < 0 ? [95,158,215] : [238,143,117],Math.min(1,Math.abs(value))).join(",")})`;
}
function drawScene(frame) {
  const {ctx,width,height}=canvasContext($("view")), extent=state.config.extent;
  const scale=Math.min(width,height)*.36*camera.zoom*1.6/extent, cut=Number($("cut").value);
  const project = point => {const p=rotate(point); return [width/2+p[0]*scale,height/2-p[1]*scale,p[2]];};
  // A subtle coordinate grid provides depth without prescribing embryo shape.
  ctx.strokeStyle="rgba(135,177,170,.09)"; ctx.lineWidth=1;
  for(let k=-4;k<=4;k++) {
    for(const pair of [[[k*extent/4,-extent*.7,-extent],[k*extent/4,-extent*.7,extent]],[[-extent,-extent*.7,k*extent/4],[extent,-extent*.7,k*extent/4]]]) {const a=project(pair[0]),b=project(pair[1]); ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();}
  }
  if($("contacts").checked && frame.graph?.weights) {
    const cellsById=new Map(frame.cells.map(cell=>[cell.id,cell]));
    const ordered=(frame.graph.ids || frame.cells.map(cell=>cell.id)).map(id=>cellsById.get(id));
    ctx.strokeStyle="rgba(174,210,200,.38)";ctx.lineWidth=1;
    ordered.forEach((cell,i)=>{if(!cell || cell.center[2]>cut) return; for(let j=i+1;j<ordered.length;j++){const other=ordered[j];if(!other || other.center[2]>cut || !(frame.graph.weights[i]?.[j]>0))continue;const a=project(cell.center),b=project(other.center);ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();}});
  }
  const points=[];
  for(const cell of frame.cells) {const color=cellColor(cell); for(const point of cell.points) if(point[2]<=cut) {const p=project(point);points.push({x:p[0],y:p[1],z:p[2],color});}}
  points.sort((a,b)=>a.z-b.z);
  const radius=Math.max(1.2,Math.min(8,scale*(2*extent/(state.config.grid-1))*.58));
  for(const p of points) {ctx.globalAlpha=.78+.2*Math.max(0,Math.min(1,(p.z/extent+1)/2));ctx.fillStyle=p.color;ctx.beginPath();ctx.arc(p.x,p.y,radius,0,Math.PI*2);ctx.fill();}
  ctx.globalAlpha=1;
  if($("arrows").checked) for(const cell of frame.cells) {
    if(cell.center[2]>cut || !cell.polarity || Math.hypot(...cell.polarity)<.02)continue;
    const a=project(cell.center),b=project(cell.center.map((v,i)=>v+.3*cell.polarity[i])); const angle=Math.atan2(b[1]-a[1],b[0]-a[0]);
    ctx.strokeStyle="#f2f7d6";ctx.fillStyle="#f2f7d6";ctx.lineWidth=1.5;ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();ctx.beginPath();ctx.moveTo(b[0],b[1]);ctx.lineTo(b[0]-7*Math.cos(angle-.45),b[1]-7*Math.sin(angle-.45));ctx.lineTo(b[0]-7*Math.cos(angle+.45),b[1]-7*Math.sin(angle+.45));ctx.closePath();ctx.fill();
  }
  ctx.font="10px system-ui";
  [["x",[1,0,0],"#df9a88"],["y",[0,1,0],"#8fbca6"],["z",[0,0,1],"#8bacd7"]].forEach(([label,axis,color])=>{const p=rotate(axis);ctx.strokeStyle=color;ctx.fillStyle=color;ctx.beginPath();ctx.moveTo(width-50,height-48);ctx.lineTo(width-50+p[0]*22,height-48-p[1]*22);ctx.stroke();ctx.fillText(label,width-53+p[0]*31,height-45-p[1]*31);});
  $("cut-value").textContent=cut >= extent-.001 ? "Full embryo" : `z ≤ ${number(cut)}`;
  const mode=$("color").value;
  const keys=mode === "fate" ? [[palette.blue,"Identity B (−)"],["#ccdbcc","Undecided"],[palette.coral,"Identity A (+)"]] : mode === "lineage" ? [["#add5aa","Distinct color per cell ID"]] : [[palette.blue,"Below 1"],["#ccdbcc","Equilibrium: 1"],[palette.coral,"Above 1"]];
  setLegend($("legend"),keys);
}
function setLegend(element,keys) {element.replaceChildren(...keys.map(([color,label])=>{const span=document.createElement("span"),swatch=document.createElement("i");swatch.className="key";swatch.style.background=color;span.append(swatch,document.createTextNode(label));return span;}));}
function drawChart(canvas,series,options={}) {
  const {ctx,width,height}=canvasContext(canvas), left=43,right=12,top=12,bottom=27;
  const w=width-left-right,h=height-top-bottom, times=frames.map(f=>f.metrics.time), start=times[0] || 0, end=Math.max(start+state.config.dt,times.at(-1) || 0);
  const values=series.flatMap(s=>frames.map(f=>s.value(f.metrics))).filter(Number.isFinite);
  let ymin=options.min ?? Math.min(0,...values), ymax=Math.max(options.max ?? 0,...values);
  if (ymax-ymin < 1e-12) ymax=ymin+(options.tiny ? .001 : 1);
  else ymax+=(ymax-ymin)*.12;
  const x=t=>left+(t-start)/(end-start)*w,y=v=>top+h-(v-ymin)/(ymax-ymin)*h;
  ctx.font="10px system-ui";ctx.lineWidth=1;
  for(let i=0;i<4;i++){const value=ymin+(ymax-ymin)*i/3,yy=y(value);ctx.strokeStyle="#e8edeb";ctx.beginPath();ctx.moveTo(left,yy);ctx.lineTo(width-right,yy);ctx.stroke();ctx.fillStyle="#83918c";ctx.textAlign="right";ctx.fillText(Math.abs(value)>0 && Math.abs(value)<.01 ? value.toExponential(0) : number(value, ymax-ymin>5 ? 0 : 2),left-7,yy+3);}
  ctx.textAlign="left";ctx.fillText(number(start,1),left,height-8);ctx.textAlign="right";ctx.fillText(`${number(end,1)}  t`,width-right,height-8);
  for(const s of series) {ctx.strokeStyle=s.color;ctx.lineWidth=2;ctx.beginPath();let active=false;frames.forEach(f=>{const value=s.value(f.metrics);if(!Number.isFinite(value)){active=false;return;} if(!active)ctx.moveTo(x(f.metrics.time),y(value));else ctx.lineTo(x(f.metrics.time),y(value));active=true;});ctx.stroke();const f=frames[selected],value=s.value(f.metrics);if(Number.isFinite(value)){ctx.fillStyle=s.color;ctx.beginPath();ctx.arc(x(f.metrics.time),y(value),3,0,Math.PI*2);ctx.fill();}}
  const marker=x(frames[selected].metrics.time);ctx.strokeStyle="#83918c66";ctx.setLineDash([3,4]);ctx.beginPath();ctx.moveTo(marker,top);ctx.lineTo(marker,top+h);ctx.stroke();ctx.setLineDash([]);
}
function drawCharts() {
  drawChart($("signals-chart"),[{color:palette.teal,value:m=>m.activator_std},{color:palette.gold,value:m=>m.inhibitor_std}],{min:0,tiny:true});
  const mode=$("organization").value; let series, title, options={min:0};
  if(mode === "shape") {title="Shape through time";series=[{label:"Axis ratio",color:palette.teal,value:m=>m.axis_ratio}];options={min:1};}
  else if(mode === "volume") {title="Volume conservation";series=[{label:"Total error (%)",color:palette.teal,value:m=>100*m.relative_volume_error},{label:"Largest cell error (%)",color:palette.gold,value:m=>100*m.max_cell_volume_error}];options={};}
  else {title="Identity through time";series=[{label:"Identity A",color:palette.coral,value:m=>m.fate_a},{label:"Identity B",color:palette.blue,value:m=>m.fate_b},{label:"Undecided",color:palette.gray,value:m=>m.uncommitted}];}
  $("organization-title").textContent=title; setLegend($("organization-key"),series.map(s=>[s.color,s.label]));drawChart($("organization-chart"),series,options);
}

$("run").addEventListener("click",()=>command("run"));$("pause").addEventListener("click",()=>command("pause"));$("reset").addEventListener("click",()=>command("reset"));
$("defaults").addEventListener("click",()=>fillParameters(schema.defaults));
$("about").addEventListener("click",()=>$("about-dialog").showModal());
for(const id of ["close-about","about-done"]) $(id).addEventListener("click",()=>$("about-dialog").close());
$("timeline").addEventListener("input",()=>{followLive=false;selected=Number($("timeline").value);schedulePaint();});
$("live").addEventListener("click",()=>{followLive=true;selected=Math.max(0,frames.length-1);schedulePaint();});
for(const id of ["color","arrows","contacts","cut","organization"]) $(id).addEventListener("input",schedulePaint);
$("camera").addEventListener("click",()=>{camera={yaw:.65,pitch:.4,zoom:1};schedulePaint();});
let drag=null;
$("view").addEventListener("pointerdown",event=>{drag={id:event.pointerId,x:event.clientX,y:event.clientY};$("view").setPointerCapture(event.pointerId);});
$("view").addEventListener("pointermove",event=>{if(!drag || drag.id !== event.pointerId)return;camera.yaw+=(event.clientX-drag.x)*.008;camera.pitch=Math.max(-1.5,Math.min(1.5,camera.pitch+(event.clientY-drag.y)*.008));drag.x=event.clientX;drag.y=event.clientY;schedulePaint();});
for(const name of ["pointerup","pointercancel","lostpointercapture"]) $("view").addEventListener(name,()=>{drag=null;});
$("view").addEventListener("wheel",event=>{event.preventDefault();camera.zoom=Math.max(.4,Math.min(3,camera.zoom*Math.exp(-event.deltaY*.001)));schedulePaint();},{passive:false});
new ResizeObserver(schedulePaint).observe($("view"));new ResizeObserver(schedulePaint).observe($("signals-chart"));
$("download").addEventListener("click",()=>{
  const payload={format:"embryo-dashboard-frames-v1",config:state.config,generation:state.generation,history_start_revision:state.history_start_revision,notes:"Retained visualization frames only; not a restart checkpoint. Earlier frames may have been pruned.",frames};
  const url=URL.createObjectURL(new Blob([JSON.stringify(payload)],{type:"application/json"}));const link=document.createElement("a");link.href=url;link.download=`embryo-seed-${state.config.seed}-step-${frames.at(-1).step}.json`;link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
});
async function initialize() {
  try {schema=await api("/api/schema");buildParameters();fillParameters(schema.defaults);await poll();}
  catch(error) {showError(`Cannot initialize dashboard: ${error.message}. Retrying…`);setTimeout(initialize,2000);}
}
initialize();
