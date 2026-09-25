/* Client state-machine tests without a browser or third-party DOM dependency.
 * These exercise real event handlers against API fixtures; layout is not tested.
 * Run with: node --test tests/test_dashboard_client.cjs
 */
const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const {execFileSync} = require('node:child_process');
const root = path.resolve(__dirname, '..');
const schema = JSON.parse(execFileSync(process.env.PYTHON || 'python', ['-c',
  'import json; from embryo.dashboard import schema; print(json.dumps(schema()))'],
  {cwd:root, env:{...process.env, OPENBLAS_NUM_THREADS:'1'}}));
const source = fs.readFileSync(path.join(root,'embryo/dashboard.js'),'utf8');
const html = fs.readFileSync(path.join(root,'embryo/dashboard.html'),'utf8');

class Element {
  constructor(tag='div') {
    this.tagName=tag; this.value=''; this.checked=false; this.disabled=false;
    this.children=[]; this.handlers={}; this.style={}; this.attributes={};
    this.classList={toggle(){}}; this.textContent='';
  }
  append(...children) {this.children.push(...children);}
  replaceChildren(...children) {this.children=children;}
  setAttribute(name,value) {this.attributes[name]=value;}
  addEventListener(name,handler) {(this.handlers[name] ||= []).push(handler);}
  dispatch(name) {return Promise.all((this.handlers[name] || []).map(handler=>handler({preventDefault(){}})));}
  checkValidity() {return !(this.min !== undefined && Number(this.value)<Number(this.min)) && !(this.max !== undefined && Number(this.value)>Number(this.max));}
  closest() {return this;}
  getBoundingClientRect() {return {width:600,height:this.tagName==='canvas' ? 200 : 30};}
  getContext() {return new Proxy({}, {get:(target,key)=>target[key] || (()=>{})});}
  showModal() {this.open=true;}
  close() {this.open=false;}
}
function frame(revision=0,step=0) {
  return {revision,step,metrics:{time:step*.015,cells:1,dividing_cells:0,activator_std:0,inhibitor_std:0,unstable_graph_modes:0,axis_ratio:1,relative_volume_error:0,max_cell_volume_error:0,mean_polarity:0,uncommitted:1,fate_a:0,fate_b:0,min_radius_grid_cells:8},cells:[{id:0,fate:0,activator:1,inhibitor:1,center:[0,0,0],polarity:[0,0,0],points:[[0,0,0]]}],graph:{weights:[[0]],ids:[0]}};
}
function snapshot(changes={}) {
  return {state:'ready',generation:0,revision:0,step:0,total_steps:1000,config:{...schema.defaults},error:null,elapsed_seconds:0,history_start_revision:0,frames:[frame()],...changes};
}
const settle = () => new Promise(resolve=>setImmediate(resolve));
async function client() {
  const elements=new Map([...html.matchAll(/id="([^"]+)"/g)].map(match=>[match[1],new Element(match[1].includes('chart') || match[1]==='view' ? 'canvas' : 'div')]));
  elements.get('color').value='fate';elements.get('organization').value='fates';elements.get('arrows').checked=true;
  const calls=[], timers=[], paints=[];
  let server=snapshot(), responder;
  const context=vm.createContext({
    document:{getElementById:id=>elements.get(id),createElement:tag=>new Element(tag),createTextNode:text=>({textContent:text})},
    window:{devicePixelRatio:1},console, setTimeout:fn=>timers.push(fn),
    requestAnimationFrame:fn=>paints.push(fn),ResizeObserver:class{observe(){}},
    fetch:async(url,options)=>{
      calls.push({url,options});
      if(responder){const result=await responder(url,options);if(result)return result;}
      return {ok:true,json:async()=>url==='/api/schema' ? schema : server};
    }
  });
  vm.runInContext(source,context); await settle();
  const run=expression=>vm.runInContext(expression,context);
  const flushPaint=()=>{while(paints.length)paints.shift()();};
  flushPaint();
  return {elements,calls,timers,run,flushPaint,
    setServer:value=>{server=value;}, respond:fn=>{responder=fn;},
    async apply(value){server=value;context.fixture=value;run('acceptState(fixture)');flushPaint();},
    async click(id){await elements.get(id).dispatch('click');await settle();flushPaint();},
    edit(name,value){const input=run(`fields.get(${JSON.stringify(name)})`);if(typeof value==='boolean')input.checked=value;else input.value=String(value);input.dispatch('input');},
  };
}

test('initial zygote exposes every Config parameter and enables Run',async()=>{
  const c=await client();
  assert.deepEqual(JSON.parse(c.run('JSON.stringify(readParameters())')),schema.defaults);
  assert.equal(c.run('fields.size'),Object.keys(schema.defaults).length);
  assert.equal(c.elements.get('run').disabled,false);
  assert.equal(c.elements.get('pause').disabled,true);
  assert.equal(c.elements.get('metric-cells').textContent,1);
  assert.equal(c.elements.get('scene-empty').hidden,true);
});

test('Run sends edits and running locks the form; Pause resumes the same config',async()=>{
  const c=await client();c.edit('seed',31);
  c.respond((url,options)=>{
    if(url==='/api/run') return {ok:true,json:async()=>snapshot({state:'running',generation:1,config:{...schema.defaults,seed:JSON.parse(options.body).config.seed}})};
    if(url==='/api/pause')return {ok:true,json:async()=>snapshot({state:'paused',generation:1,step:13,revision:1,config:{...schema.defaults,seed:31},frames:[frame(),frame(1,13)]})};
  });
  await c.click('run');
  assert.equal(JSON.parse(c.calls.at(-1).options.body).config.seed,31);
  assert.equal(c.elements.get('pause').disabled,false);
  assert.equal(c.run('fields.get("seed").disabled'),true);
  await c.click('pause');
  assert.equal(c.elements.get('run-label').textContent,'Resume');
  assert.equal(c.elements.get('metric-step').textContent,'Displayed frame: step 13');
  assert.equal(c.run('fields.get("seed").disabled'),false);
  assert.equal(c.run('state.config.seed'),31);
});

test('edits after pause require Reset and reset replaces history with a zygote',async()=>{
  const c=await client();await c.apply(snapshot({state:'paused',step:20,revision:1,frames:[frame(),frame(1,20)]}));
  c.edit('seed',44);assert.equal(c.elements.get('run').disabled,true);
  assert.match(c.elements.get('parameter-note').textContent,/Reset/);
  c.respond((url,options)=>url==='/api/reset' ? {ok:true,json:async()=>snapshot({generation:1,config:JSON.parse(options.body).config})} : null);
  await c.click('reset');
  assert.equal(c.run('frames.length'),1);assert.equal(c.run('state.step'),0);
  assert.equal(c.run('state.config.seed'),44);assert.equal(c.elements.get('run').disabled,false);
});

test('invalid numeric edits disable execution; backend validation errors preserve current run',async()=>{
  const c=await client();c.edit('steps','');assert.equal(c.elements.get('run').disabled,true);assert.equal(c.elements.get('reset').disabled,true);
  c.edit('steps',1000);c.edit('dt',1);
  c.respond(url=>url==='/api/run' ? {ok:false,json:async()=>({error:'unstable time step'})} : null);
  await c.click('run');
  assert.equal(c.elements.get('error').textContent,'unstable time step');
  assert.equal(c.run('state.step'),0);assert.equal(c.run('state.config.dt'),.015);
});

test('timeline inspection stays on the chosen frame as live frames arrive',async()=>{
  const c=await client();await c.apply(snapshot({state:'running',step:20,revision:1,frames:[frame(),frame(1,20)]}));
  c.elements.get('timeline').value=0;await c.elements.get('timeline').dispatch('input');c.flushPaint();
  await c.apply(snapshot({state:'running',step:40,revision:2,frames:[frame(2,40)]}));
  assert.equal(c.elements.get('metric-step').textContent,'Displayed frame: step 0');
  assert.equal(c.elements.get('view-mode').textContent,'HISTORICAL FRAME');
  assert.match(c.elements.get('run-detail').textContent,/Step 40/);
  await c.click('live');assert.equal(c.elements.get('metric-step').textContent,'Displayed frame: step 40');
});

test('pruning drops old client frames and marks incomplete retained history',async()=>{
  const c=await client();await c.apply(snapshot({revision:2,step:40,frames:[frame(1,20),frame(2,40)],history_start_revision:1}));
  assert.equal(c.run('frames.length'),2);assert.equal(c.run('frames[0].revision'),1);
  assert.match(c.elements.get('history-note').textContent,/pruned/);
});

test('polling a reset with a reused revision fetches the new generation history',async()=>{
  const c=await client();let count=0;
  c.respond(url=>{
    if(url.startsWith('/api/state')){count++;return {ok:true,json:async()=>snapshot({generation:1,config:{...schema.defaults,seed:42},frames:count===1 ? [] : [frame()]})};}
  });
  await c.run('poll()');c.flushPaint();
  assert.equal(count,2);assert.equal(c.calls.at(-1).url,'/api/state?after=-1');
  assert.equal(c.run('frames.length'),1);assert.equal(c.run('state.config.seed'),42);
});

test('a stale poll cannot undo a Reset response',async()=>{
  const c=await client();let release;
  c.respond(url=>{
    if(url.startsWith('/api/state'))return new Promise(resolve=>{release=()=>resolve({ok:true,json:async()=>snapshot({state:'running',step:40,revision:2,frames:[frame(2,40)]})});});
    if(url==='/api/reset')return {ok:true,json:async()=>snapshot({generation:1})};
  });
  const pending=c.run('poll()');await settle();await c.click('reset');release();await pending;c.flushPaint();
  assert.equal(c.run('state.generation'),1);assert.equal(c.run('state.step'),0);assert.equal(c.run('state.state'),'ready');
});

test('network disconnect disables mutations and reconnection restores them',async()=>{
  const c=await client();c.respond(()=>{throw new Error('offline');});await c.run('poll()');
  assert.equal(c.elements.get('run').disabled,true);assert.match(c.elements.get('connection').textContent,/Disconnected/);
  c.respond(()=>null);await c.run('poll()');assert.equal(c.elements.get('run').disabled,false);
});

test('a stale failed poll cannot disconnect a newer successful Reset',async()=>{
  const c=await client();let release;
  c.respond(url=>{
    if(url.startsWith('/api/state'))return new Promise((resolve,reject)=>{release=()=>reject(new Error('old request failed'));});
    if(url==='/api/reset')return {ok:true,json:async()=>snapshot({generation:1})};
  });
  const pending=c.run('poll()');await settle();await c.click('reset');release();await pending;
  assert.equal(c.elements.get('run').disabled,false);assert.equal(c.run('connected'),true);
});

test('a reset in another tab clears the previous solver error',async()=>{
  const c=await client();await c.apply(snapshot({state:'error',error:'numerical failure'}));
  assert.equal(c.elements.get('error').hidden,false);
  await c.apply(snapshot({generation:1}));
  assert.equal(c.elements.get('error').hidden,true);
});
