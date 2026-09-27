const {test}=require('node:test');
const assert=require('node:assert/strict');
const {draw,clipZ,rgb}=require('../embryo/cell_surface.js');

test('cutaway clips crossing triangles at the actual plane',()=>{
  const polygon=clipZ([[0,0,0],[1,0,0],[0,1,1]],.5);
  assert.equal(polygon.length,4);
  assert(polygon.every(p=>p[2]<=.5));
  assert.deepEqual(polygon.filter(p=>p[2]===.5),[[.5,.5,.5],[0,.5,.5]]);
  assert.equal(clipZ([[0,0,1],[1,0,1],[0,1,1]],0).length,0);
});
test('RGB and lineage HSL colors are retained',()=>{
  assert.deepEqual(rgb('rgb(255,0,128)'),[1,0,128/255]);
  assert.deepEqual(rgb('hsl(120 100% 50%)'),[0,1,0]);
});
test('without WebGL, real mesh triangles render instead of dots',()=>{
  global.document={createElement:()=>({getContext:()=>null})};
  let polygons=0,dots=0;
  const ctx={beginPath(){},moveTo(){},lineTo(){},closePath(){polygons++;},fill(){},arc(){dots++;}};
  const cell={points:[[0,0,0]],mesh:{vertices:[[0,0,0],[1,0,0],[0,1,0]],faces:[[0,1,2]],normals:[[0,0,1],[0,0,1],[0,0,1]],closed:false}};
  const options={cells:[cell],mode:'surfaces',color:()=> 'rgb(200,100,50)',cut:2,rotate:p=>p,project:p=>p};
  const result=draw(ctx,options);
  assert.equal(result.backend,'canvas');assert.equal(result.open,1);assert.equal(polygons,1);assert.equal(dots,0);
  draw(ctx,{...options,mode:'points'});assert.equal(dots,1);
  const legacy=draw(ctx,{...options,cells:[{points:[[0,0,0]]}]});
  assert.equal(legacy.legacy,1);assert.equal(dots,2);
});
