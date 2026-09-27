/* Shared offline/live isosurface renderer. Rendering never advances the solver. */
(function(root){
  'use strict';
  const contexts=new WeakMap();
  function rgb(css){
    const v=css.match(/[-+]?\d*\.?\d+/g).map(Number);
    if(css.startsWith('rgb'))return v.slice(0,3).map(x=>x/255);
    const h=((v[0]%360)+360)%360/60,s=v[1]/100,l=v[2]/100;
    const c=(1-Math.abs(2*l-1))*s,x=c*(1-Math.abs(h%2-1)),m=l-c/2;
    const values=[[c,x,0],[x,c,0],[0,c,x],[0,x,c],[x,0,c],[c,0,x]][Math.floor(h)%6];
    return values.map(a=>a+m);
  }
  function clipZ(vertices,cut){
    const result=[];
    for(let i=0;i<vertices.length;i++){
      const a=vertices[i],b=vertices[(i+1)%vertices.length],ia=a[2]<=cut,ib=b[2]<=cut;
      if(ia)result.push(a);
      if(ia!==ib){const t=(cut-a[2])/(b[2]-a[2]);result.push(a.map((v,j)=>v+t*(b[j]-v)));}
    }
    return result;
  }
  function create(ctx){
    const canvas=document.createElement('canvas');
    const gl=canvas.getContext('webgl',{alpha:true,antialias:true,preserveDrawingBuffer:true,premultipliedAlpha:false});
    if(!gl)throw Error('WebGL unavailable');
    const vertex=`attribute vec3 position; attribute vec3 normal;
      uniform vec4 rotation; uniform vec3 projection;
      varying vec3 n; varying float worldZ;
      vec3 turn(vec3 p){float z=-rotation.y*p.x+rotation.x*p.z;
        return vec3(rotation.x*p.x+rotation.y*p.z,rotation.z*p.y-rotation.w*z,rotation.w*p.y+rotation.z*z);}
      void main(){vec3 p=turn(position);gl_Position=vec4(p.x*projection.x,p.y*projection.y,-p.z*projection.z,1.0);
        n=turn(normal);worldZ=position.z;}`;
    const fragment=`precision highp float; uniform vec3 color; uniform float cut;
      varying vec3 n; varying float worldZ;
      void main(){if(worldZ>cut)discard;
        vec3 normal=normalize(n)*(gl_FrontFacing?1.0:-1.0);
        vec3 light=normalize(vec3(-0.4,0.7,1.0));
        float diffuse=max(dot(normal,light),0.0);
        float spec=pow(max(dot(normal,normalize(light+vec3(0.0,0.0,1.0))),0.0),32.0);
        gl_FragColor=vec4(color*(0.35+0.65*diffuse)+vec3(0.13)*spec,1.0);}`;
    const program=gl.createProgram();
    for(const [type,source] of [[gl.VERTEX_SHADER,vertex],[gl.FRAGMENT_SHADER,fragment]]){
      const shader=gl.createShader(type);gl.shaderSource(shader,source);gl.compileShader(shader);
      if(!gl.getShaderParameter(shader,gl.COMPILE_STATUS))throw Error('Surface shader failed');
      gl.attachShader(program,shader);gl.deleteShader(shader);
    }
    gl.linkProgram(program);if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw Error('Surface program failed');
    const locations={};for(const name of ['rotation','projection','color','cut'])locations[name]=gl.getUniformLocation(program,name);
    locations.position=gl.getAttribLocation(program,'position');locations.normal=gl.getAttribLocation(program,'normal');
    const state={canvas,gl,program,locations,cells:null,buffers:[]};
    canvas.addEventListener('webglcontextlost',e=>e.preventDefault());
    canvas.addEventListener('webglcontextrestored',()=>contexts.delete(ctx));
    return state;
  }
  function gpu(ctx,o){
    let state=contexts.get(ctx);
    if(state===false)return false;
    if(!state){try{state=create(ctx);contexts.set(ctx,state);}catch(error){contexts.set(ctx,false);return false;}}
    const {gl,canvas,program,locations:l}=state;
    if(gl.isContextLost())return false;
    const w=Math.max(1,Math.round(o.width*o.dpr)),h=Math.max(1,Math.round(o.height*o.dpr));
    if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;}
    gl.viewport(0,0,w,h);gl.clearColor(0,0,0,0);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST);gl.disable(gl.CULL_FACE);gl.disable(gl.BLEND);gl.useProgram(program);
    gl.uniform4f(l.rotation,Math.cos(o.yaw),Math.sin(o.yaw),Math.cos(o.pitch),Math.sin(o.pitch));
    gl.uniform3f(l.projection,2*o.scale/o.width,2*o.scale/o.height,1/(4*o.extent));gl.uniform1f(l.cut,o.cut);
    if(state.cells!==o.cells){
      for(const record of state.buffers)gl.deleteBuffer(record.buffer);
      state.buffers=[];state.cells=o.cells;
      for(const cell of o.cells){
        const m=cell.mesh;if(!m?.faces?.length)continue;
        const data=new Float32Array(m.faces.length*18);let k=0;
        for(const face of m.faces)for(const index of face){data.set(m.vertices[index],k);data.set(m.normals[index],k+3);k+=6;}
        const buffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,data,gl.STATIC_DRAW);
        state.buffers.push({cell,buffer,count:m.faces.length*3});
      }
    }
    gl.enableVertexAttribArray(l.position);gl.enableVertexAttribArray(l.normal);
    for(const record of state.buffers){
      gl.bindBuffer(gl.ARRAY_BUFFER,record.buffer);gl.vertexAttribPointer(l.position,3,gl.FLOAT,false,24,0);gl.vertexAttribPointer(l.normal,3,gl.FLOAT,false,24,12);
      gl.uniform3fv(l.color,rgb(o.color(record.cell)));gl.drawArrays(gl.TRIANGLES,0,record.count);
    }
    ctx.drawImage(canvas,0,0,o.width,o.height);return true;
  }
  function software(ctx,o){
    const triangles=[];
    for(const cell of o.cells){
      const m=cell.mesh;if(!m?.faces?.length)continue;
      const base=rgb(o.color(cell));
      for(const face of m.faces){
        const polygon=clipZ(face.map(i=>m.vertices[i]),o.cut);if(polygon.length<3)continue;
        const points=polygon.map(o.project),normal=o.rotate([0,1,2].map(j=>face.reduce((sum,i)=>sum+m.normals[i][j],0)/3));
        const norm=Math.hypot(...normal)||1;
        const diffuse=Math.max(0,(-.4*normal[0]+.7*normal[1]+normal[2])/(Math.sqrt(1.65)*norm));
        triangles.push({points,z:points.reduce((sum,p)=>sum+p[2],0)/points.length,
          color:`rgb(${base.map(c=>Math.round(255*c*(.35+.65*diffuse))).join(',')})`});
      }
    }
    triangles.sort((a,b)=>a.z-b.z);
    for(const t of triangles){ctx.fillStyle=t.color;ctx.beginPath();ctx.moveTo(t.points[0][0],t.points[0][1]);for(const p of t.points.slice(1))ctx.lineTo(p[0],p[1]);ctx.closePath();ctx.fill();}
  }
  function draw(ctx,o){
    const meshes=o.cells.filter(c=>c.mesh?.faces?.length),surface=o.mode!=='points';
    let backend='points';
    if(surface&&meshes.length){backend=gpu(ctx,o)?'webgl':'canvas';if(backend==='canvas')software(ctx,o);}
    const points=[];
    for(const cell of o.cells){if(surface&&cell.mesh)continue;
      for(const v of cell.points||[])if(v[2]<=o.cut)points.push({p:o.project(v),color:o.color(cell)});
    }
    points.sort((a,b)=>a.p[2]-b.p[2]);ctx.globalAlpha=1;
    for(const {p,color} of points){ctx.fillStyle=color;ctx.beginPath();ctx.arc(p[0],p[1],o.radius||1.5,0,2*Math.PI);ctx.fill();}
    return {backend,legacy:o.cells.filter(c=>!c.mesh).length,open:meshes.filter(c=>!c.mesh.closed).length};
  }
  root.CellSurface={draw,clipZ,rgb};
  if(typeof module!=='undefined')module.exports=root.CellSurface;
})(globalThis);
