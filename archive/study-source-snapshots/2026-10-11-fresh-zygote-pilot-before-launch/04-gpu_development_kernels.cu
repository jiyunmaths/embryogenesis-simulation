// Isolated developmental extension. No accepted mature source is modified.
#include "gpu_precision_kernels.cu"

__global__ void development_force_kernel(const float* phi,const float* excluded,
 const double* attract,const double* gamma,const double* vf,const double* centers,
 const double* events,float* out,double* carry,int* rounding,int cells,int n,
 double dx,double extent,double width,double repulsion,double dt,double strength) {
 __shared__ int active[256],stuck[256];int t=threadIdx.x;
 size_t index=(size_t)blockIdx.x*256+t,plane=(size_t)n*n,size=plane*n;
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
  double f=width*width*flux-gamma[index]*double(derivative);
  f+=((vf[cell]*6.)*double(value))*double(one_minus);
  float repel=(float(repulsion)*value)*excluded[index];f-=double(repel);
  f+=double(derivative)*attract[index];
  const double* ev=events+7*cell;
  if(ev[0]!=0){
   double x=(i+.5)*dx-extent-centers[3*cell];
   double y=(j+.5)*dx-extent-centers[3*cell+1];
   double z=(k+.5)*dx-extent-centers[3*cell+2];
   double axial=ev[1]*x+ev[2]*y+ev[3]*z;
   double radial=sqrt(fmax(x*x+y*y+z*z-axial*axial,0.));
   double band=exp(-.5*pow((axial-ev[4])/width,2.));
   double outside=.5*(1.+tanh((radial-ev[5]*(1.-ev[6]))/width));
   f-=strength*ev[6]*band*outside*6.*double(value)*double(one_minus);
  }
  double inc=dt*f;float original=float(double(value)+inc);
  a=(inc!=0.);b=(a&&original==value);
  double corrected=double(value)+(carry[index]+inc);float rounded=float(corrected);
  carry[index]=(rounded<0.f||rounded>1.f)?0.:corrected-double(rounded);
  out[index]=rounded;
 }
 active[t]=a;stuck[t]=b;__syncthreads();
 for(int s=128;s;s>>=1){if(t<s){active[t]+=active[t+s];stuck[t]+=stuck[t+s];}__syncthreads();}
 if(!t){atomicAdd(rounding,active[0]);atomicAdd(rounding+1,stuck[0]);}
}

extern "C" int development_force(const float* p,const float* ex,const double* att,
 const double* centers,const double* polarity,const double* tensions,const double* vf,
 double* gamma,float* out,int* quality,double* carry,int* rounding,const double* events,
 int cells,int n,double dx,double extent,double width,double contrast,
 double repulsion,double dt,double strength,void* stream,int device) {
 try{
  check(cudaSetDevice(device));int count=cells*n*n*n,blocks=(count+255)/256;auto s=(cudaStream_t)stream;
  check(cudaMemsetAsync(quality,0,2*sizeof(int),s));check(cudaMemsetAsync(rounding,0,2*sizeof(int),s));
  gamma_kernel<double><<<blocks,256,0,s>>>(gamma,centers,polarity,tensions,cells,n,dx,extent,width,contrast);
  development_force_kernel<<<blocks,256,0,s>>>(p,ex,att,gamma,vf,centers,events,out,carry,rounding,cells,n,dx,extent,width,repulsion,dt,strength);
  clip_quality<<<blocks,256,0,s>>>(out,quality,count);check(cudaGetLastError());return 0;
 }catch(const std::exception& e){error_text=e.what();return 1;}
}
