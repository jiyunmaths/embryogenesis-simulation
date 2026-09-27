/* Recorded continuum fields. Playback and rendering never advance a solver. */
"use strict";
(() => {
  const el = id => document.getElementById(id);
  const fmt = (v, digits=3) => Number.isFinite(v) ? (Math.abs(v)>0 && Math.abs(v)<.001 ? v.toExponential(2) : v.toFixed(digits)) : "—";
  const percent = v => Number.isFinite(v) ? `${fmt(100*v)}%` : "Undefined (uniform field)";
  const colors = ["#209d91", "#d89f38", "#839ca6", "#b77dba"];
  const ramp = [[53,42,135],[36,108,170],[27,158,145],[140,190,85],[245,223,77]];
  let metadata=null, data=null, row=null, frame=0, epoch=0, loading=false;
  let playing=false, timer=null, visible=false, paintPending=false;
  let camera={yaw:.65,pitch:.45,zoom:1}, drag=null, lookup=null;
  const checks = {
    initial_fields_agree_under_restriction:"Same physical initial field across meshes",
    positive_concentrations_all_steps:"Positive concentrations at every step",
    discrete_reaction_amount_balance:"Reaction-adjusted amount balance",
    persistent_fine_mesh_contrast:"Persistent finest-mesh contrast",
    fine_mesh_late_pattern_stationary:"Late field change ≤ 1%",
    spatial_final_difference_decreases:"Spatial difference decreases on refinement",
    fine_spatial_final_error_below_5_percent:"Finest spatial difference ≤ 5%",
    fine_spatial_final_correlation_above_0_99:"Finest activator correlation ≥ 0.99",
    temporal_error_decreases:"Temporal difference decreases on refinement",
    fine_temporal_error_below_0_5_percent:"Finer temporal difference ≤ 0.5%",
    equal_diffusivity_control_loses_contrast:"Equal-diffusivity control loses contrast"
  };
  async function get(url) {
    const response=await fetch(url,{cache:"no-store"});
    const value=await response.json();
    if(!response.ok) throw new Error(value.error || `Request failed (${response.status})`);
    return value;
  }
  function message(text, error=false) {
    el("benchmark-message").hidden=!text;
    el("benchmark-message").textContent=text || "";
    el("benchmark-message").classList.toggle("error",error);
  }
  function stop() {
    playing=false;
    if(timer!==null) clearTimeout(timer);
    timer=null;
    el("benchmark-play").textContent="▶ Play";
  }
  function tick() {
    if(!playing || !visible || !data) return;
    if(frame>=data.times.length-1) {stop();return;}
    frame++; paint();
    if(frame===data.times.length-1) stop();
    else timer=setTimeout(tick,900);
  }
  function controls() {
    for(const id of ["benchmark-play","benchmark-reset","benchmark-timeline","benchmark-slice","benchmark-species","benchmark-view","benchmark-axis","benchmark-camera"]) el(id).disabled=loading || !data;
    el("benchmark-loading").hidden=!loading;
  }
  async function loadRun(id) {
    stop(); const token=++epoch;
    const previousTime=data?.times[frame] ?? 0;
    loading=true; controls(); message(null);
    try {
      const next=await get(`/api/benchmark/run?id=${encodeURIComponent(id)}`);
      if(token!==epoch) return;
      const nextRow=metadata.runs.find(r=>r.id===id);
      if(!nextRow || next.id!==id) throw new Error("Unexpected benchmark run response");
      data=next; row=nextRow;
      frame=data.times.reduce((best,t,i)=>Math.abs(t-previousTime)<Math.abs(data.times[best]-previousTime) ? i : best,0);
      lookup=new Int32Array(data.n**3); lookup.fill(-1);
      data.indices.forEach((ijk,index)=>{lookup[(ijk[0]*data.n+ijk[1])*data.n+ijk[2]]=index;});
      const oldMax=Number(el("benchmark-slice").max), fraction=oldMax>0 ? Number(el("benchmark-slice").value)/oldMax : 1;
      el("benchmark-slice").max=data.n-1;
      el("benchmark-slice").value=Math.round(fraction*(data.n-1));
      el("benchmark-content").hidden=false;
    } catch(error) {
      if(token!==epoch) return;
      // Keep a failed selection from being mistaken for the previous fields.
      data=null; row=null;
      el("benchmark-content").hidden=true;
      message(`Cannot load recorded fields: ${error.message}. Switch to Live embryo and back to retry.`,true);
      metadata=null;
    } finally {
      if(token===epoch) {loading=false; controls(); paint();}
    }
  }
  async function initialize() {
    if(metadata || loading) return;
    loading=true; controls(); message("Loading the nonlinear experiment…");
    try {
      const result=await get("/api/benchmark");
      if(!result.available) {
        message(`${result.message}\n${result.command}\n${result.hint}`,true);
        return;
      }
      metadata=result;
      el("benchmark-run").replaceChildren();
      for(const run of result.runs) {
        const option=document.createElement("option"); option.value=run.id; option.textContent=run.label;
        el("benchmark-run").append(option);
      }
      el("benchmark-run").value=result.default_run;
      el("benchmark-source").textContent=`Source: ${result.source}`;
      el("benchmark-checks").replaceChildren();
      const passed=Object.values(result.checks).filter(v=>v===true).length;
      el("benchmark-check-summary").textContent=`${passed} / ${Object.keys(result.checks).length} acceptance checks passed`;
      for(const [key,value] of Object.entries(result.checks)) {
        const item=document.createElement("li"); item.textContent=`${value ? "✓ Pass" : "✕ Fail"} · ${checks[key] || key.replaceAll("_"," ")}`;
        item.className=value ? "passed" : "failed"; el("benchmark-checks").append(item);
      }
      const spatial=result.spatial_comparisons.at(-1), temporal=result.temporal_comparisons.at(-1);
      el("benchmark-summary").textContent=`Final finest-pair spatial difference: ${percent(spatial.final.relative_joint_rms)}. Finer time-pair maximum sampled difference: ${percent(temporal.maximum_relative_joint_rms)}. Evidence applies to this perturbation and this fixed domain.`;
      await loadRun(result.default_run);
    } catch(error) {metadata=null;message(`Cannot read benchmark: ${error.message}. Switch workspaces to retry.`,true);}
    finally {loading=false;controls();}
  }
  function workspace(which) {
    visible=which==="benchmark";
    el("embryo-workspace").hidden=visible; el("benchmark-workspace").hidden=!visible;
    for(const name of ["embryo","benchmark"]) {
      const active=name===which; el(`tab-${name}`).classList.toggle("active",active);
      el(`tab-${name}`).setAttribute("aria-pressed",String(active));
    }
    if(visible) {initialize();paint();}
    else {stop();if(typeof schedulePaint==="function") schedulePaint();}
  }
  function paint() {
    if(paintPending) return;
    paintPending=true;
    requestAnimationFrame(()=>{paintPending=false;if(visible && data && !loading) render();});
  }
  function context(canvas) {
    const rect=canvas.getBoundingClientRect(), ratio=Math.min(window.devicePixelRatio || 1,2);
    const w=Math.max(1,Math.round(rect.width*ratio)),h=Math.max(1,Math.round(rect.height*ratio));
    if(canvas.width!==w || canvas.height!==h){canvas.width=w;canvas.height=h;}
    const ctx=canvas.getContext("2d");ctx.setTransform(ratio,0,0,ratio,0,0);ctx.clearRect(0,0,rect.width,rect.height);
    return {ctx,width:rect.width,height:rect.height};
  }
  function color(value, range) {
    const v=Math.max(0,Math.min(1,(value-range[0])/Math.max(range[1]-range[0],1e-12)))*(ramp.length-1);
    const i=Math.min(ramp.length-2,Math.floor(v)),t=v-i;
    return `rgb(${ramp[i].map((x,j)=>Math.round(x*(1-t)+ramp[i+1][j]*t)).join(",")})`;
  }
  function rotated(point) {
    const center=data.edges.map(e=>(e[0]+e.at(-1))/2),length=data.edges[0].at(-1)-data.edges[0][0];
    const [x,y,z]=point.map((v,i)=>(v-center[i])/length);
    const cy=Math.cos(camera.yaw),sy=Math.sin(camera.yaw),cp=Math.cos(camera.pitch),sp=Math.sin(camera.pitch);
    const xx=cy*x+sy*z,zz=-sy*x+cy*z;
    return [xx,cp*y-sp*zz,sp*y+cp*zz];
  }
  function drawField(values, range, axis, cut) {
    const {ctx,width,height}=context(el("benchmark-canvas"));
    const axes=[0,1,2].filter(i=>i!==axis), names=["X","Y","Z"];
    if(el("benchmark-view").value==="slice") {
      const [u,v]=axes,eu=data.edges[u],ev=data.edges[v];
      const scale=Math.min((width-100)/(eu.at(-1)-eu[0]),(height-140)/(ev.at(-1)-ev[0]));
      const left=(width-(eu.at(-1)-eu[0])*scale)/2, bottom=height/2+(ev.at(-1)-ev[0])*scale/2+15;
      data.indices.forEach((ijk,i)=>{if(ijk[axis]!==cut)return;
        const x=left+(eu[ijk[u]]-eu[0])*scale,y=bottom-(ev[ijk[v]+1]-ev[0])*scale;
        ctx.fillStyle=color(values[i],range);
        ctx.fillRect(x,y,(eu[ijk[u]+1]-eu[ijk[u]])*scale+.2,(ev[ijk[v]+1]-ev[ijk[v]])*scale+.2);
      });
      ctx.fillStyle="#abc4bc";ctx.font="11px system-ui";ctx.textAlign="center";
      ctx.fillText(`${names[u]} [${fmt(eu[0],1)}, ${fmt(eu.at(-1),1)}]`,width/2,bottom+20);
      ctx.textAlign="left";ctx.fillText(names[v],left-20,bottom-(ev.at(-1)-ev[0])*scale/2);
      return;
    }
    const scale=Math.min(width,height)*.64*camera.zoom;
    const project=p=>{const q=rotated(p);return [width/2+q[0]*scale,height/2-q[1]*scale+18,q[2]];};
    const faces=[];
    function kept(ijk) {
      return ijk[axis]<=cut && ijk.every(v=>v>=0 && v<data.n) && lookup[(ijk[0]*data.n+ijk[1])*data.n+ijk[2]]>=0;
    }
    data.indices.forEach((ijk,index)=>{
      if(ijk[axis]>cut)return;
      for(let d=0;d<3;d++) for(const side of [-1,1]) {
        const neighbor=ijk.slice();neighbor[d]+=side;if(kept(neighbor))continue;
        const other=[0,1,2].filter(a=>a!==d);
        const points=[[0,0],[1,0],[1,1],[0,1]].map(bits=>{
          const p=ijk.map((v,a)=>data.edges[a][v]);
          p[d]=data.edges[d][ijk[d]+(side===1 ? 1 : 0)];
          other.forEach((a,j)=>{p[a]=data.edges[a][ijk[a]+bits[j]];});return project(p);
        });
        faces.push({points,z:points.reduce((sum,p)=>sum+p[2],0)/4,color:color(values[index],range)});
      }
    });
    faces.sort((a,b)=>a.z-b.z);
    for(const face of faces) {
      ctx.fillStyle=face.color;ctx.strokeStyle=face.color;ctx.lineWidth=.45;
      ctx.beginPath();face.points.forEach((p,i)=>i===0 ? ctx.moveTo(p[0],p[1]) : ctx.lineTo(p[0],p[1]));
      ctx.closePath();ctx.fill();ctx.stroke();
    }
    // Axes share the field's camera, making rotation and cut direction legible.
    const origin=data.edges.map(e=>(e[0]+e.at(-1))/2),r0=rotated(origin);
    ctx.font="10px system-ui";ctx.textAlign="center";
    for(let axisIndex=0;axisIndex<3;axisIndex++) {
      const p=origin.slice();p[axisIndex]+=.18*(data.edges[0].at(-1)-data.edges[0][0]);const r=rotated(p);
      const x=55+(r[0]-r0[0])*180,y=height-70-(r[1]-r0[1])*180;
      ctx.strokeStyle="#b2cbc2";ctx.fillStyle="#b2cbc2";ctx.beginPath();ctx.moveTo(55,height-70);ctx.lineTo(x,y);ctx.stroke();ctx.fillText(names[axisIndex],x,y-6);
    }
  }
  function legend(id,series) {
    const container=el(id);container.replaceChildren();
    for(const s of series){const span=document.createElement("span"),dot=document.createElement("i");dot.style.background=s.color;span.append(dot,document.createTextNode(s.name));container.append(span);}
  }
  function chart(id,series,time,late=false) {
    const {ctx,width,height}=context(el(id)),left=51,right=13,top=13,bottom=27,w=width-left-right,h=height-top-bottom;
    const end=metadata.parameters.duration;
    let max=0;for(const s of series)for(const p of s.points)if(Number.isFinite(p[1]))max=Math.max(max,p[1]);
    max=max>1e-12 ? max*1.12 : .001;
    const x=t=>left+t/end*w,y=v=>top+h-v/max*h;
    if(late){ctx.fillStyle="#eff2ef";ctx.fillRect(x(row.late_window_start),top,x(end)-x(row.late_window_start),h);}
    ctx.font="10px system-ui";
    for(let i=0;i<4;i++) {const v=max*i/3,yy=y(v);ctx.strokeStyle="#e8edeb";ctx.beginPath();ctx.moveTo(left,yy);ctx.lineTo(width-right,yy);ctx.stroke();ctx.fillStyle="#738779";ctx.textAlign="right";ctx.fillText(fmt(v,max>2 ? 1 : 3),left-7,yy+3);}
    ctx.textAlign="left";ctx.fillText("0",left,height-7);ctx.textAlign="right";ctx.fillText(`${fmt(end,1)} t`,width-right,height-7);
    for(const s of series){ctx.strokeStyle=s.color;ctx.fillStyle=s.color;ctx.lineWidth=2;ctx.beginPath();let active=false;
      for(const [t,v] of s.points){if(!Number.isFinite(v)){active=false;continue;}if(active)ctx.lineTo(x(t),y(v));else ctx.moveTo(x(t),y(v));active=true;}ctx.stroke();
      if(s.dots)for(const [t,v] of s.points){if(!Number.isFinite(v))continue;ctx.beginPath();ctx.arc(x(t),y(v),2.5,0,Math.PI*2);ctx.fill();}
    }
    ctx.strokeStyle="#748e8390";ctx.setLineDash([3,4]);ctx.beginPath();ctx.moveTo(x(time),top);ctx.lineTo(x(time),top+h);ctx.stroke();ctx.setLineDash([]);
  }
  function render() {
    const time=data.times[frame],species=el("benchmark-species").value,axis=Number(el("benchmark-axis").value),cut=Number(el("benchmark-slice").value);
    const range=metadata.ranges[species];
    el("benchmark-count").textContent=row.compartments.toLocaleString();
    el("benchmark-resolution").textContent=`n=${row.n} · control volumes, not biological cells`;
    el("benchmark-time").textContent=fmt(time,2);el("benchmark-dt").textContent=`Solver dt=${row.dt} · β=${row.beta} · Da=${row.da}, Dh=${row.dh}`;
    // History may not include every field time for non-default sample intervals.
    const history=row.history.find(r=>Math.abs(r.time-time)<1e-9);
    el("benchmark-contrast").textContent=history ? fmt(history.std_a) : "Not sampled";
    el("benchmark-late").textContent=Number.isFinite(row.maximum_late_relative_field_change) ? percent(row.maximum_late_relative_field_change) : "—";
    el("benchmark-window").textContent=`Full-run measurement · t=${row.late_window_start}–${row.duration}`;
    el("benchmark-caption").textContent=`${row.kind==="control" ? "Equal-diffusivity control" : "Nonlinear signaling"} · t=${fmt(time,2)} · n=${row.n}`;
    el("benchmark-legend-title").textContent=species==="a" ? "Activator concentration" : "Inhibitor concentration";
    el("benchmark-color-min").textContent=fmt(range[0]);el("benchmark-color-max").textContent=fmt(range[1]);
    el("benchmark-timeline").max=data.times.length-1;el("benchmark-timeline").value=frame;
    el("benchmark-frame").textContent=`Frame ${frame+1} / ${data.times.length}`;
    const edge=data.edges[axis];
    el("benchmark-slice-label").textContent=el("benchmark-view").value==="slice" ? `${["X","Y","Z"][axis]} ∈ [${fmt(edge[cut])}, ${fmt(edge[cut+1])}]` : `${["X","Y","Z"][axis]} ≤ ${fmt(edge[cut+1])}`;
    el("benchmark-scene-help").textContent=el("benchmark-view").value==="slice" ? "Cell-average concentrations · unequal compartment widths preserved" : "Drag to orbit · scroll to zoom · slice to inspect the interior";
    drawField(data[species][frame],range,axis,cut);
    const control=metadata.runs.find(r=>r.kind==="control");
    const series=[{name:"Activator σ",color:colors[0],points:row.history.map(r=>[r.time,r.std_a])},
      {name:"Inhibitor σ",color:colors[1],points:row.history.map(r=>[r.time,r.std_h])},
      {name:`Equal diffusion · n=${control.n} · activator σ`,color:colors[2],points:control.history.map(r=>[r.time,r.std_a])}];
    legend("benchmark-signal-key",series);chart("benchmark-signals",series,time,true);
    const mode=el("benchmark-comparison").value;
    const comparisons=metadata[`${mode}_comparisons`].map((r,i)=>({name:mode==="spatial" ? `n=${r.coarse_n} vs ${r.fine_n}` : `dt=${metadata.parameters.temporal_steps[i]} vs ${metadata.parameters.temporal_steps[i+1]}`,
      color:colors[i%colors.length],dots:true,points:r.snapshots.map(p=>[p.time,Number.isFinite(p.relative_joint_rms) ? 100*p.relative_joint_rms : null])}));
    el("benchmark-comparison-title").textContent=mode==="spatial" ? "Spatial agreement (%)" : "Time-step agreement (%)";
    legend("benchmark-comparison-key",comparisons);chart("benchmark-errors",comparisons,time);
  }
  el("tab-embryo").addEventListener("click",()=>workspace("embryo"));
  el("tab-benchmark").addEventListener("click",()=>workspace("benchmark"));
  el("benchmark-run").addEventListener("change",()=>loadRun(el("benchmark-run").value));
  el("benchmark-play").addEventListener("click",()=>{
    if(!data || loading)return;
    if(playing){stop();return;}
    if(frame===data.times.length-1)frame=0;
    playing=true;el("benchmark-play").textContent="Ⅱ Pause";paint();timer=setTimeout(tick,900);
  });
  el("benchmark-reset").addEventListener("click",()=>{stop();frame=0;paint();});
  el("benchmark-timeline").addEventListener("input",()=>{stop();frame=Number(el("benchmark-timeline").value);paint();});
  for(const id of ["benchmark-species","benchmark-view","benchmark-axis","benchmark-slice","benchmark-comparison"])el(id).addEventListener("input",paint);
  el("benchmark-camera").addEventListener("click",()=>{camera={yaw:.65,pitch:.45,zoom:1};paint();});
  const canvas=el("benchmark-canvas");
  canvas.addEventListener("pointerdown",event=>{if(el("benchmark-view").value!=="3d")return;drag={id:event.pointerId,x:event.clientX,y:event.clientY};canvas.setPointerCapture(event.pointerId);});
  canvas.addEventListener("pointermove",event=>{if(!drag || event.pointerId!==drag.id)return;camera.yaw+=(event.clientX-drag.x)*.008;camera.pitch=Math.max(-1.5,Math.min(1.5,camera.pitch+(event.clientY-drag.y)*.008));drag.x=event.clientX;drag.y=event.clientY;paint();});
  for(const event of ["pointerup","pointercancel","lostpointercapture"])canvas.addEventListener(event,()=>{drag=null;});
  canvas.addEventListener("wheel",event=>{if(el("benchmark-view").value!=="3d")return;event.preventDefault();camera.zoom=Math.max(.4,Math.min(2.5,camera.zoom*Math.exp(-event.deltaY*.001)));paint();},{passive:false});
  new ResizeObserver(paint).observe(canvas);new ResizeObserver(paint).observe(el("benchmark-signals"));
  controls();
})();
