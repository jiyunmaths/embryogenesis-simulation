// Numerical control only: retain sub-float32 phase increments in a double carry.
// Existing field/geometry/force arithmetic stays unchanged; never replace the
// accepted resident kernel. Zero carry has the original one-step update.
#include "gpu_kernels.cu"

__global__ void precision_force_kernel(const float* phi,const float* excluded,
 const double* attract,const double* gamma,const double* vf,float* out,
 double* carry,int* rounding,int cells,int n,double dx,double width,
 double repulsion,double dt,int use_carry) {
 __shared__ int active[256],stuck[256];int t=threadIdx.x;
 size_t index=(size_t)blockIdx.x*blockDim.x+t,plane=(size_t)n*n,size=plane*n;
 int a=0,b=0;
 if(index<size*cells){
  int cell=index/size;size_t v=index%size;int i=v/plane,j=(v/n)%n,k=v%n;
  double flux=0,dx2=dx*dx;
  if(i<n-1)flux+=.5*(gamma[index]+gamma[index+plane])*double(phi[index+plane]-phi[index])/dx2;
  if(i>0)flux-=.5*(gamma[index-plane]+gamma[index])*double(phi[index]-phi[index-plane])/dx2;
  if(j<n-1)flux+=.5*(gamma[index]+gamma[index+n])*double(phi[index+n]-phi[index])/dx2;
  if(j>0)flux-=.5*(gamma[index-n]+gamma[index])*double(phi[index]-phi[index-n])/dx2;
  if(k<n-1)flux+=.5*(gamma[index]+gamma[index+1])*double(phi[index+1]-phi[index])/dx2;
  if(k>0)flux-=.5*(gamma[index-1]+gamma[index])*double(phi[index]-phi[index-1])/dx2;
  float value=phi[index],one_minus=1.f-value;
  float derivative=((2.f*value)*one_minus)*(1.f-2.f*value);
  double force=width*width*flux-gamma[index]*double(derivative);
  force+=((vf[cell]*6.)*double(value))*double(one_minus);
  float repel=(float(repulsion)*value)*excluded[index];force-=double(repel);
  force+=double(derivative)*attract[index];double increment=dt*force;
  float original=float(double(value)+increment);
  a=(increment!=0.);b=(a&&original==value);
  if(use_carry){
   double corrected=double(value)+(carry[index]+increment);float rounded=float(corrected);
   carry[index]=(rounded<0.f||rounded>1.f)?0.:corrected-double(rounded);
   out[index]=rounded;
  }else{out[index]=original;carry[index]=0.;}
 }
 active[t]=a;stuck[t]=b;__syncthreads();
 for(int s=128;s;s>>=1){if(t<s){active[t]+=active[t+s];stuck[t]+=stuck[t+s];}__syncthreads();}
 if(!t){atomicAdd(rounding,active[0]);atomicAdd(rounding+1,stuck[0]);}
}

extern "C" int precision_force(const float* p,const float* ex,const double* att,
 const double* centers,const double* polarity,const double* tensions,const double* vf,
 double* gamma,float* out,int* quality,double* carry,int* rounding,
 int cells,int n,int use_carry,double dx,double extent,double width,
 double contrast,double repulsion,double dt,void* stream,int device) {
 try{
  check(cudaSetDevice(device));int count=cells*n*n*n,blocks=(count+255)/256;auto s=(cudaStream_t)stream;
  check(cudaMemsetAsync(quality,0,2*sizeof(int),s));check(cudaMemsetAsync(rounding,0,2*sizeof(int),s));
  gamma_kernel<double><<<blocks,256,0,s>>>(gamma,centers,polarity,tensions,cells,n,dx,extent,width,contrast);
  precision_force_kernel<<<blocks,256,0,s>>>(p,ex,att,gamma,vf,out,carry,rounding,cells,n,dx,width,repulsion,dt,use_carry);
  clip_quality<<<blocks,256,0,s>>>(out,quality,count);check(cudaGetLastError());return 0;
 }catch(const std::exception& e){error_text=e.what();return 1;}
}
