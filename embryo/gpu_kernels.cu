// Stream-aware tensor-pointer entry points; no allocations or host field copies.
#include "cuda_spatial_bench.cu"

__global__ void exclusion(const float* phi, float* out, int cells, int size) {
 int v=blockIdx.x*blockDim.x+threadIdx.x;if(v>=size)return;
 float total=0;for(int c=0;c<cells;c++){float p=phi[c*size+v];total+=p*p;}
 for(int c=0;c<cells;c++){float p=phi[c*size+v];out[c*size+v]=total-p*p;}
}
__global__ void clip_quality(float* phi, int* quality, int count) {
 __shared__ int clips[256],invalid[256];int t=threadIdx.x,i=blockIdx.x*256+t;
 int a=0,b=0;if(i<count){float x=phi[i];a=(x<0||x>1);b=!isfinite(x);phi[i]=fminf(1.f,fmaxf(0.f,x));}
 clips[t]=a;invalid[t]=b;__syncthreads();
 for(int s=128;s;s>>=1){if(t<s){clips[t]+=clips[t+s];invalid[t]+=invalid[t+s];}__syncthreads();}
 if(!t){atomicAdd(quality,clips[0]);atomicAdd(quality+1,invalid[0]);}
}
__global__ void polarity_update(const double* p,const double* cue,const double* a,
 const double* aligned,double* result,int cells,double dt,double rate,double alignment,double decay) {
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=cells)return;
 double norm=0;for(int d=0;d<3;d++)norm+=p[3*i+d]*p[3*i+d];
 double drive=2*a[i]/(1+a[i]),r[3],mag=0;
 for(int d=0;d<3;d++){int j=3*i+d;double drift=rate*drive*cue[j]+alignment*aligned[j];drift-=(decay+norm)*p[j];r[d]=p[j]+dt*drift;mag+=r[d]*r[d];}
 mag=fmax(1.,sqrt(mag));for(int d=0;d<3;d++)result[3*i+d]=r[d]/mag;
}
__global__ void union_partial(const float* occupied,double* partial,int n,int chunks,double dx,double extent) {
 __shared__ double sums[10][256];int t=threadIdx.x,v=blockIdx.x*256+t,size=n*n*n;
 double val[10]={};if(v<size){int i=v/(n*n),j=(v/n)%n,k=v%n;double w=fminf(occupied[v],1.f);
 double x=(i+.5)*dx-extent,y=(j+.5)*dx-extent,z=(k+.5)*dx-extent;
 val[0]=w;val[1]=w*x;val[2]=w*y;val[3]=w*z;
 val[4]=w*x*x;val[5]=w*y*y;val[6]=w*z*z;val[7]=w*x*y;val[8]=w*x*z;val[9]=w*y*z;}
 for(int q=0;q<10;q++)sums[q][t]=val[q];__syncthreads();
 for(int s=128;s;s>>=1){if(t<s)for(int q=0;q<10;q++)sums[q][t]+=sums[q][t+s];__syncthreads();}
 if(!t)for(int q=0;q<10;q++)partial[q*chunks+blockIdx.x]=sums[q][0];
}

extern "C" int gpu_arrays(const float* p,float* h,float* shell,double* s2,float* occupied,
 float* ex,int cells,int n,void* stream,int device) {
 try{check(cudaSetDevice(device));int size=n*n*n,count=cells*size;auto s=(cudaStream_t)stream;
 spatial_arrays<<<(count+255)/256,256,0,s>>>(p,h,shell,s2,cells,size);
 occupied_array<<<(size+255)/256,256,0,s>>>(h,occupied,cells,size);
 exclusion<<<(size+255)/256,256,0,s>>>(p,ex,cells,size);check(cudaGetLastError());return 0;
 }catch(const std::exception& e){error_text=e.what();return 1;}
}
extern "C" int gpu_geometry(const float* p,const float* h,const float* shell,const float* occupied,
 double* part,double* geom,int cells,int n,double dx,double extent,void* stream,int device) {
 try{check(cudaSetDevice(device));int chunks=(n*n*n+255)/256;auto s=(cudaStream_t)stream;
 partial_geometry<<<dim3(chunks,cells),256,0,s>>>(p,h,shell,occupied,part,n,chunks,dx,extent);
 finish_geometry<<<cells,256,0,s>>>(part,geom,chunks,dx);check(cudaGetLastError());return 0;
 }catch(const std::exception& e){error_text=e.what();return 1;}
}
extern "C" int gpu_force(const float* p,const float* ex,const double* att,const double* centers,
 const double* polarity,const double* tensions,const double* vf,double* gamma,float* out,int* quality,
 int cells,int n,double dx,double extent,double width,double contrast,double repulsion,double dt,void* stream,int device) {
 try{check(cudaSetDevice(device));int count=cells*n*n*n,blocks=(count+255)/256;auto s=(cudaStream_t)stream;
 check(cudaMemsetAsync(quality,0,2*sizeof(int),s));
 gamma_kernel<double><<<blocks,256,0,s>>>(gamma,centers,polarity,tensions,cells,n,dx,extent,width,contrast);
 force_kernel<double><<<blocks,256,0,s>>>(p,ex,att,gamma,vf,out,cells,n,dx,width,repulsion,dt);
 clip_quality<<<blocks,256,0,s>>>(out,quality,count);check(cudaGetLastError());return 0;
 }catch(const std::exception& e){error_text=e.what();return 1;}
}
extern "C" int gpu_polarity(const double* p,const double* cue,const double* a,const double* aligned,
 double* out,int cells,double dt,double rate,double alignment,double decay,void* stream,int device) {
 try{check(cudaSetDevice(device));polarity_update<<<(cells+255)/256,256,0,(cudaStream_t)stream>>>(p,cue,a,aligned,out,cells,dt,rate,alignment,decay);check(cudaGetLastError());return 0;
 }catch(const std::exception& e){error_text=e.what();return 1;}
}
extern "C" int gpu_union(const float* occupied,double* partial,int n,double dx,double extent,void* stream,int device) {
 try{check(cudaSetDevice(device));int chunks=(n*n*n+255)/256;union_partial<<<chunks,256,0,(cudaStream_t)stream>>>(occupied,partial,n,chunks,dx,extent);check(cudaGetLastError());return 0;
 }catch(const std::exception& e){error_text=e.what();return 1;}
}
