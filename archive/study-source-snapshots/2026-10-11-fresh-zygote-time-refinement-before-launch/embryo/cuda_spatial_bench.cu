// Standalone spatial-operation benchmark. No live simulation backend changes.
#include "cuda_mechanics.cu"
#include <cublas_v2.h>
#include <cstring>
static void blascheck(cublasStatus_t s){if(s!=CUBLAS_STATUS_SUCCESS)throw std::runtime_error("cuBLAS failure "+std::to_string(s));}
struct Blas {cublasHandle_t h;Blas(){blascheck(cublasCreate(&h));}~Blas(){cublasDestroy(h);}};
__global__ void spatial_arrays(const float* phi,float* h,float* shell,double* shell2,int cells,int size) {
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=size*cells)return;
 float p=phi[i],s=p*(1.f-p);h[i]=p*p*(3.f-2.f*p);shell[i]=s;shell2[i]=double(s*s);
}
__global__ void occupied_array(const float* h,float* occupied,int cells,int size) {
 int v=blockIdx.x*blockDim.x+threadIdx.x;if(v>=size)return;
 float sum=0;for(int c=0;c<cells;c++)sum+=h[c*size+v];occupied[v]=sum;
}
__global__ void partial_geometry(const float* phi,const float* h,const float* shell,
 const float* occupied,double* partial,int n,int chunks,double dx,double extent) {
 __shared__ double sums[8][256];int t=threadIdx.x,cell=blockIdx.y,chunk=blockIdx.x;
 int size=n*n*n,plane=n*n,v=chunk*256+t,index=cell*size+v;
 double val[8]={};
 if(v<size){
  int i=v/plane,j=(v/n)%n,k=v%n;
  val[0]=h[index];val[1]=val[0]*((i+.5)*dx-extent);
  val[2]=val[0]*((j+.5)*dx-extent);val[3]=val[0]*((k+.5)*dx-extent);val[4]=shell[index];
  float d=float(dx),gx=(phi[index+(i<n-1?plane:0)]-phi[index-(i>0?plane:0)])/((i==0||i==n-1)?d:2.f*d);
  float gy=(phi[index+(j<n-1?n:0)]-phi[index-(j>0?n:0)])/((j==0||j==n-1)?d:2.f*d);
  float gz=(phi[index+(k<n-1?1:0)]-phi[index-(k>0?1:0)])/((k==0||k==n-1)?d:2.f*d);
  float length=fmaxf(sqrtf(gx*gx+gy*gy+gz*gz),1.e-10f);
  float blocked=fminf(1.f,fmaxf(0.f,2.f*(occupied[v]-h[index])));
  float weight=shell[index]*(1.f-blocked);
  val[5]=(-gx/length)*weight;val[6]=(-gy/length)*weight;val[7]=(-gz/length)*weight;
 }
 for(int q=0;q<8;q++)sums[q][t]=val[q];__syncthreads();
 for(int s=128;s;s>>=1){if(t<s)for(int q=0;q<8;q++)sums[q][t]+=sums[q][t+s];__syncthreads();}
 if(!t)for(int q=0;q<8;q++)partial[(cell*8+q)*chunks+chunk]=sums[q][0];
}
__global__ void finish_geometry(const double* partial,double* result,int chunks,double dx) {
 __shared__ double sums[8][256];int cell=blockIdx.x,t=threadIdx.x;
 for(int q=0;q<8;q++){double x=0;for(int j=t;j<chunks;j+=256)x+=partial[(cell*8+q)*chunks+j];sums[q][t]=x;}
 __syncthreads();for(int s=128;s;s>>=1){if(t<s)for(int q=0;q<8;q++)sums[q][t]+=sums[q][t+s];__syncthreads();}
 if(!t){double h=sums[0][0],den=fmax(sums[4][0],1.e-12);
 result[cell*7]=h*dx*dx*dx;for(int q=1;q<4;q++)result[cell*7+q]=sums[q][0]/h;
 double x=sums[5][0]/den,y=sums[6][0]/den,z=sums[7][0]/den;
 if(sqrt(x*x+y*y+z*z)<1.e-6)x=y=z=0;
 result[cell*7+4]=x;result[cell*7+5]=y;result[cell*7+6]=z;}
}
extern "C" int spatial_benchmark(const float* phi,const double* adhesion,int cells,int n,
 double dx,double extent,int repeats,double* geometry,float* contact,double* attraction,float* times) {
 try {check(cudaSetDevice(0));Buffers b;Blas blas;Events e;int size=n*n*n,count=cells*size,chunks=(size+255)/256;
 auto p=b.upload(phi,count);auto a=b.upload(adhesion,cells*cells);
 auto h=b.allocate<float>(count),shell=b.allocate<float>(count),occupied=b.allocate<float>(size);
 auto s2=b.allocate<double>(count),part=b.allocate<double>(cells*8*chunks),out=b.allocate<double>(cells*7);
 auto c=b.allocate<float>(cells*cells);auto att=b.allocate<double>(count);
 float one=1,zero=0;double done=1,dzero=0;
 auto op=[&](int k){
  if(k==0){spatial_arrays<<<(count+255)/256,256>>>(p,h,shell,s2,cells,size);occupied_array<<<chunks,256>>>(h,occupied,cells,size);}
  if(k==1)blascheck(cublasSgemm(blas.h,CUBLAS_OP_T,CUBLAS_OP_N,cells,cells,size,&one,shell,size,shell,size,&zero,c,cells));
  if(k==2)blascheck(cublasDgemm(blas.h,CUBLAS_OP_N,CUBLAS_OP_N,size,cells,cells,&done,s2,size,a,cells,&dzero,att,size));
  if(k==3){partial_geometry<<<dim3(chunks,cells),256>>>(p,h,shell,occupied,part,n,chunks,dx,extent);finish_geometry<<<cells,256>>>(part,out,chunks,dx);}
  check(cudaGetLastError());
 };
 for(int k=0;k<4;k++)op(k);check(cudaDeviceSynchronize());
 for(int k=0;k<4;k++){check(cudaEventRecord(e.start));for(int r=0;r<repeats;r++)op(k);
 check(cudaEventRecord(e.stop));check(cudaEventSynchronize(e.stop));check(cudaEventElapsedTime(times+k,e.start,e.stop));times[k]/=repeats;}
 check(cudaMemcpy(geometry,out,cells*7*sizeof(double),cudaMemcpyDeviceToHost));
 check(cudaMemcpy(contact,c,cells*cells*sizeof(float),cudaMemcpyDeviceToHost));
 check(cudaMemcpy(attraction,att,count*sizeof(double),cudaMemcpyDeviceToHost));return 0;
 }catch(const std::exception& e){error_text=e.what();return 1;}
}
extern "C" int transfer_benchmark(int megabytes,int repeats,int pinned,float* times) {
 void* host=nullptr;void* device=nullptr;
 try {check(cudaSetDevice(0));size_t size=(size_t)megabytes*1024*1024;
 if(pinned)check(cudaMallocHost(&host,size));else {host=malloc(size);if(!host)throw std::runtime_error("host allocation failed");}
 memset(host,37,size);check(cudaMalloc(&device,size));Events e;
 check(cudaMemcpy(device,host,size,cudaMemcpyHostToDevice));
 for(int direction=0;direction<2;direction++){
 check(cudaEventRecord(e.start));for(int i=0;i<repeats;i++){
 if(!direction)check(cudaMemcpy(device,host,size,cudaMemcpyHostToDevice));else check(cudaMemcpy(host,device,size,cudaMemcpyDeviceToHost));}
 check(cudaEventRecord(e.stop));check(cudaEventSynchronize(e.stop));check(cudaEventElapsedTime(times+direction,e.start,e.stop));times[direction]/=repeats;}
 if(((unsigned char*)host)[0]!=37||((unsigned char*)host)[size-1]!=37)throw std::runtime_error("copy mismatch");
 cudaFree(device);if(pinned)cudaFreeHost(host);else free(host);return 0;
 }catch(const std::exception& e){if(device)cudaFree(device);if(host){if(pinned)cudaFreeHost(host);else free(host);}error_text=e.what();return 1;}
}
