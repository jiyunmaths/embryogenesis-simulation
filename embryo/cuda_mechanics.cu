// Experimental benchmark only. Preserve the CPU operation order in FP64 mode.
#include <cuda_runtime.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
static thread_local std::string error_text;
static void check(cudaError_t code) {
 if(code!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(code));
}
struct Buffers {
 std::vector<void*> allocations;
 template<class T> T* allocate(size_t count) {
  T* result;check(cudaMalloc((void**)&result,count*sizeof(T)));allocations.push_back(result);return result;
 }
 template<class T> T* upload(const T* source,size_t count) {
  T* result=allocate<T>(count);check(cudaMemcpy(result,source,count*sizeof(T),cudaMemcpyHostToDevice));return result;
 }
 ~Buffers(){for(void* p:allocations)cudaFree(p);}
};
struct Events {
 cudaEvent_t start=nullptr,stop=nullptr;
 Events(){check(cudaEventCreate(&start));check(cudaEventCreate(&stop));}
 ~Events(){if(start)cudaEventDestroy(start);if(stop)cudaEventDestroy(stop);}
};
template<class T> __global__ void gamma_kernel(T* gamma,const double* centers,
 const double* polarity,const double* tensions,int cells,int n,T dx,T extent,T width,T contrast) {
 size_t index=(size_t)blockIdx.x*blockDim.x+threadIdx.x,size=(size_t)n*n*n;
 if(index>=size*cells)return;
 int cell=index/size;size_t v=index%size;int i=v/(n*n),j=(v/n)%n,k=v%n;
 T x=(T(i)+T(.5))*dx-extent-T(centers[3*cell]);
 T y=(T(j)+T(.5))*dx-extent-T(centers[3*cell+1]);
 T z=(T(k)+T(.5))*dx-extent-T(centers[3*cell+2]);
 T r=sqrt(x*x+y*y+z*z+width*width);
 gamma[index]=T(tensions[cell])*(T(1)-contrast*(T(polarity[3*cell])*(x/r)+T(polarity[3*cell+1])*(y/r)+T(polarity[3*cell+2])*(z/r)));
}
template<class T> __global__ void force_kernel(const float* phi,const float* excluded,
 const double* attract,const T* gamma,const double* vf,float* out,int cells,int n,
 T dx,T width,T repulsion,T dt) {
 size_t index=(size_t)blockIdx.x*blockDim.x+threadIdx.x,plane=(size_t)n*n,size=plane*n;
 if(index>=size*cells)return;
 int cell=index/size;size_t v=index%size;int i=v/plane,j=(v/n)%n,k=v%n;
 T flux=0,dx2=dx*dx;
 if(i<n-1)flux+=T(.5)*(gamma[index]+gamma[index+plane])*T(phi[index+plane]-phi[index])/dx2;
 if(i>0)flux-=T(.5)*(gamma[index-plane]+gamma[index])*T(phi[index]-phi[index-plane])/dx2;
 if(j<n-1)flux+=T(.5)*(gamma[index]+gamma[index+n])*T(phi[index+n]-phi[index])/dx2;
 if(j>0)flux-=T(.5)*(gamma[index-n]+gamma[index])*T(phi[index]-phi[index-n])/dx2;
 if(k<n-1)flux+=T(.5)*(gamma[index]+gamma[index+1])*T(phi[index+1]-phi[index])/dx2;
 if(k>0)flux-=T(.5)*(gamma[index-1]+gamma[index])*T(phi[index]-phi[index-1])/dx2;
 float value=phi[index],one_minus=1.f-value;
 float derivative=((2.f*value)*one_minus)*(1.f-2.f*value);
 T force=width*width*flux-gamma[index]*T(derivative);
 force+=((T(vf[cell])*T(6))*T(value))*T(one_minus);
 float repel=(float(repulsion)*value)*excluded[index];force-=T(repel);
 force+=T(derivative)*T(attract[index]);out[index]=float(T(value)+dt*force);
}
extern "C" const char* cuda_mechanics_error(){return error_text.c_str();}
extern "C" int cuda_mechanics_device(char* name,int length) {
 try {check(cudaSetDevice(0));cudaDeviceProp p;check(cudaGetDeviceProperties(&p,0));
 snprintf(name,length,"%s (sm_%d%d)",p.name,p.major,p.minor);return 0;
 }catch(const std::exception& e){error_text=e.what();return 1;}
}
extern "C" int cuda_mechanics_update(const float* phi,const float* excluded,
 const double* attract,const double* centers,const double* polarity,const double* tensions,
 const double* vf,float* out,int cells,int n,double dx,double extent,double width,
 double contrast,double repulsion,double dt,int precision,int repeats,float* milliseconds) {
 try {
 check(cudaSetDevice(0));Buffers memory;size_t count=(size_t)cells*n*n*n;
 auto f=memory.upload(phi,count),ex=memory.upload(excluded,count);
 auto a=memory.upload(attract,count),c=memory.upload(centers,(size_t)cells*3);
 auto p=memory.upload(polarity,(size_t)cells*3),t=memory.upload(tensions,cells),v=memory.upload(vf,cells);
 auto result=memory.allocate<float>(count);int blocks=(count+255)/256;
 double* gd=nullptr;float* gf=nullptr;
 if(precision==64)gd=memory.allocate<double>(count);else gf=memory.allocate<float>(count);
 auto launch=[&]() {
  if(precision==64) {
   gamma_kernel<double><<<blocks,256>>>(gd,c,p,t,cells,n,dx,extent,width,contrast);
   force_kernel<double><<<blocks,256>>>(f,ex,a,gd,v,result,cells,n,dx,width,repulsion,dt);
  }else {
   gamma_kernel<float><<<blocks,256>>>(gf,c,p,t,cells,n,dx,extent,width,contrast);
   force_kernel<float><<<blocks,256>>>(f,ex,a,gf,v,result,cells,n,dx,width,repulsion,dt);
  }
  check(cudaGetLastError());
 };
 // Warm first; repetitions use the same input, not an evolving trajectory.
 launch();check(cudaDeviceSynchronize());Events events;
 check(cudaEventRecord(events.start));for(int r=0;r<repeats;r++)launch();
 check(cudaEventRecord(events.stop));check(cudaEventSynchronize(events.stop));
 check(cudaEventElapsedTime(milliseconds,events.start,events.stop));*milliseconds/=repeats;
 check(cudaMemcpy(out,result,count*sizeof(float),cudaMemcpyDeviceToHost));return 0;
 }catch(const std::exception& e){error_text=e.what();return 1;}
}
