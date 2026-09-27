/* Real playback handlers against recorded-field fixtures, including fetch races. */
const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const path=require('node:path');
const root=path.resolve(__dirname,'..');
const source=fs.readFileSync(path.join(root,'embryo/benchmark.js'),'utf8');
const html=fs.readFileSync(path.join(root,'embryo/dashboard.html'),'utf8');
class Element {
  constructor(){this.value='';this.max='31';this.hidden=false;this.disabled=false;this.textContent='';this.handlers={};this.children=[];this.attributes={};this.style={};this.classList={toggle(){}};}
  append(...children){this.children.push(...children);}
  replaceChildren(...children){this.children=children;}
  setAttribute(name,value){this.attributes[name]=value;}
  addEventListener(event,fn){(this.handlers[event] ||= []).push(fn);}
  dispatch(event){return Promise.all((this.handlers[event] || []).map(fn=>fn({preventDefault(){}})));}
  getBoundingClientRect(){return {width:700,height:300};}
  getContext(){return new Proxy({}, {get:(target,key)=>target[key] || (()=>{})});}
}
const times=[0,20,40,60,75,80,90,100];
function run(id,n=4,kind='spatial'){
  return {id,n,kind,label:id,dt:.05,beta:2,da:.02,dh:kind==='control' ? .02 : .4,compartments:1,
    late_window_start:75,duration:100,maximum_late_relative_field_change:kind==='control' ? null : .004,
    history:times.map(time=>({time,std_a:kind==='control' ? 0 : time*.006,std_h:time*.002}))};
}
function fields(id){return {id,n:4,times,edges:[[0,.25,.5,.75,1],[0,.25,.5,.75,1],[0,.25,.5,.75,1]],indices:[[0,0,0]],a:times.map(t=>[1+t*.01]),h:times.map(t=>[1+t*.005])};}
function metadata(){return {available:true,source:'fixture',default_run:'n4',runs:[run('n4'),run('n8',8),run('stable_control',4,'control')],ranges:{a:[.03,3.1],h:[.5,2]},checks:{positive_concentrations_all_steps:true,persistent_fine_mesh_contrast:false},parameters:{duration:100,temporal_steps:[.05,.025,.0125]},spatial_comparisons:[{coarse_n:4,fine_n:8,final:{relative_joint_rms:.01},snapshots:times.map(time=>({time,relative_joint_rms:time*.0001}))}],temporal_comparisons:[{maximum_relative_joint_rms:.001,snapshots:times.map(time=>({time,relative_joint_rms:time*.00001}))}]};}
const settle=()=>new Promise(resolve=>setImmediate(resolve));
async function client(meta=metadata()){
  const elements=new Map([...html.matchAll(/id="([^"]+)"/g)].map(m=>[m[1],new Element()]));
  for(const [id,value] of Object.entries({'benchmark-species':'a','benchmark-view':'3d','benchmark-axis':'2','benchmark-slice':'31','benchmark-comparison':'spatial'}))elements.get(id).value=value;
  const calls=[],paints=[],timers=new Map();let serial=0,responder=null;
  const context=vm.createContext({document:{getElementById:id=>elements.get(id),createElement:()=>new Element(),createTextNode:text=>({textContent:text})},window:{devicePixelRatio:1},console,
    requestAnimationFrame:fn=>paints.push(fn),ResizeObserver:class{observe(){}},
    setTimeout:fn=>{timers.set(++serial,fn);return serial;},clearTimeout:id=>timers.delete(id),
    fetch:async url=>{calls.push(url);if(responder){const result=await responder(url);if(result)return result;}return {ok:true,json:async()=>url==='/api/benchmark' ? meta : fields(decodeURIComponent(url.split('=')[1]))};}
  });
  vm.runInContext(source,context);
  const paint=()=>{while(paints.length)paints.shift()();};
  return {elements,calls,timers,paint,respond:fn=>{responder=fn;},
    async click(id){await elements.get(id).dispatch('click');await settle();paint();},
    async edit(id,value,event='input'){elements.get(id).value=String(value);await elements.get(id).dispatch(event);await settle();paint();},
    tick(){const item=timers.entries().next().value;if(item){timers.delete(item[0]);item[1]();paint();}}
  };
}

test('benchmark loads lazily and renders actual saved times, checks, and fixed scales',async()=>{
  const c=await client();assert.equal(c.calls.length,0);
  await c.click('tab-benchmark');
  assert.equal(c.elements.get('embryo-workspace').hidden,true);
  assert.equal(c.elements.get('benchmark-content').hidden,false);
  assert.equal(c.elements.get('benchmark-time').textContent,'0.00');
  assert.equal(c.elements.get('benchmark-check-summary').textContent,'1 / 2 acceptance checks passed');
  assert.equal(c.elements.get('benchmark-color-max').textContent,'3.100');
  await c.edit('benchmark-timeline',4);
  assert.equal(c.elements.get('benchmark-time').textContent,'75.00');
  assert.equal(c.elements.get('benchmark-color-max').textContent,'3.100');
  assert.ok(c.calls.every(url=>url.startsWith('/api/benchmark')));
});

test('play, pause, reset and workspace switching never call solver mutation routes',async()=>{
  const c=await client();await c.click('tab-benchmark');await c.click('benchmark-play');c.tick();
  assert.equal(c.elements.get('benchmark-time').textContent,'20.00');
  await c.click('benchmark-play');assert.equal(c.timers.size,0);
  await c.click('benchmark-play');await c.click('tab-embryo');assert.equal(c.timers.size,0);
  await c.click('tab-benchmark');await c.click('benchmark-reset');
  assert.equal(c.elements.get('benchmark-time').textContent,'0.00');
  assert.equal(c.calls.length,2);
});

test('playback stops at final saved frame and starts over when played again',async()=>{
  const c=await client();await c.click('tab-benchmark');await c.edit('benchmark-timeline',6);
  await c.click('benchmark-play');c.tick();
  assert.equal(c.elements.get('benchmark-time').textContent,'100.00');
  assert.equal(c.elements.get('benchmark-play').textContent,'▶ Play');
  assert.equal(c.timers.size,0);
  await c.click('benchmark-play');assert.equal(c.elements.get('benchmark-time').textContent,'0.00');
});

test('switching runs preserves physical snapshot time and shared chemical color scale',async()=>{
  const c=await client();await c.click('tab-benchmark');await c.edit('benchmark-timeline',5);
  await c.edit('benchmark-run','stable_control','change');
  assert.equal(c.elements.get('benchmark-time').textContent,'80.00');
  assert.equal(c.elements.get('benchmark-contrast').textContent,'0.000');
  assert.equal(c.elements.get('benchmark-late').textContent,'—');
  assert.equal(c.elements.get('benchmark-color-max').textContent,'3.100');
  await c.edit('benchmark-species','h');
  assert.equal(c.elements.get('benchmark-color-max').textContent,'2.000');
  assert.match(c.elements.get('benchmark-legend-title').textContent,/Inhibitor/);
});

test('cross-sections and temporal comparisons use their selected controls',async()=>{
  const c=await client();await c.click('tab-benchmark');await c.edit('benchmark-view','slice');
  await c.edit('benchmark-axis','0');await c.edit('benchmark-slice',1);
  assert.equal(c.elements.get('benchmark-slice-label').textContent,'X ∈ [0.250, 0.500]');
  await c.edit('benchmark-comparison','temporal');
  assert.equal(c.elements.get('benchmark-comparison-title').textContent,'Time-step agreement (%)');
});

test('a delayed run response cannot replace a newer selection',async()=>{
  const c=await client();await c.click('tab-benchmark');let release;
  c.respond(url=>url.endsWith('id=n8') ? new Promise(resolve=>{release=()=>resolve({ok:true,json:async()=>fields('n8')});}) : null);
  const pending=c.edit('benchmark-run','n8','change');await settle();
  await c.edit('benchmark-run','stable_control','change');release();await pending;
  assert.match(c.elements.get('benchmark-caption').textContent,/Equal-diffusivity control/);
  assert.equal(c.elements.get('benchmark-run').value,'stable_control');
});

test('missing benchmark shows actionable instructions with playback disabled',async()=>{
  const c=await client({available:false,message:'No completed benchmark.',command:'python -m embryo.nonlinear',hint:'Use --benchmark PATH'});
  await c.click('tab-benchmark');
  assert.match(c.elements.get('benchmark-message').textContent,/python -m embryo.nonlinear/);
  assert.equal(c.elements.get('benchmark-play').disabled,true);
});

test('failed field loads hide stale data and allow a fresh load after switching workspaces',async()=>{
  const c=await client();await c.click('tab-benchmark');
  c.respond(url=>url.endsWith('id=n8') ? {ok:false,json:async()=>({error:'broken fields'})} : null);
  await c.edit('benchmark-run','n8','change');
  assert.equal(c.elements.get('benchmark-content').hidden,true);
  assert.equal(c.elements.get('benchmark-play').disabled,true);
  assert.match(c.elements.get('benchmark-message').textContent,/broken fields/);
  c.respond(()=>null);await c.click('tab-embryo');await c.click('tab-benchmark');
  assert.equal(c.elements.get('benchmark-content').hidden,false);
});
