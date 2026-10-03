#include <cmath>
#include <vector>
#include <cstddef>
#include <omp.h>
extern "C" void mechanical_update(const float* phi,const float* excluded,
 const double* attract,const double* centers,const double* polarity,
 const double* tensions,const double* vf,float* result,
 int cells,int n,double dx,double extent,double width,double contrast,
 double repulsion,double dt,int threads) {
 const size_t plane=(size_t)n*n,size=plane*n;const double dx2=dx*dx,w2=width*width;
 #pragma omp parallel for schedule(static) num_threads(threads)
 for(int cell=0;cell<cells;cell++) {
  std::vector<double> gamma(size);const float* f=phi+cell*size;
  const double* p=polarity+3*cell;const double* center=centers+3*cell;
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int k=0;k<n;k++) {
   size_t v=(size_t)i*plane+(size_t)j*n+k;
   double x=(i+.5)*dx-extent-center[0],y=(j+.5)*dx-extent-center[1],z=(k+.5)*dx-extent-center[2];
   double r=std::sqrt(x*x+y*y+z*z+w2);
   gamma[v]=tensions[cell]*(1.-contrast*(p[0]*(x/r)+p[1]*(y/r)+p[2]*(z/r)));
  }
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int k=0;k<n;k++) {
   size_t v=(size_t)i*plane+(size_t)j*n+k,index=cell*size+v;
   double flux=0.;
   if(i<n-1)flux+=.5*(gamma[v]+gamma[v+plane])*(double)(f[v+plane]-f[v])/dx2;
   if(i>0)flux-=.5*(gamma[v-plane]+gamma[v])*(double)(f[v]-f[v-plane])/dx2;
   if(j<n-1)flux+=.5*(gamma[v]+gamma[v+n])*(double)(f[v+n]-f[v])/dx2;
   if(j>0)flux-=.5*(gamma[v-n]+gamma[v])*(double)(f[v]-f[v-n])/dx2;
   if(k<n-1)flux+=.5*(gamma[v]+gamma[v+1])*(double)(f[v+1]-f[v])/dx2;
   if(k>0)flux-=.5*(gamma[v-1]+gamma[v])*(double)(f[v]-f[v-1])/dx2;
   float value=f[v],one_minus=1.f-value;
   float derivative=((2.f*value)*one_minus)*(1.f-2.f*value);
   double force=w2*flux-gamma[v]*(double)derivative;
   force+=((vf[cell]*6.)*(double)value)*(double)one_minus;
   float repel=((float)repulsion*value)*excluded[index];
   force-=(double)repel;
   force+=(double)derivative*attract[index];
   result[index]=(float)((double)value+dt*force);
  }
 }
}
extern "C" void surface_weights(const float *phi,const float *h,const float *occupied,
 float *out,int cells,int n,float dx,int threads) {
 const size_t plane=(size_t)n*n,size=plane*n;
 #pragma omp parallel for schedule(static) num_threads(threads)
 for(int cell=0;cell<cells;cell++) {
  const float *f=phi+cell*size,*occ=h+cell*size;
  float *ox=out+cell*3*size,*oy=ox+size,*oz=oy+size;
  for(int i=0;i<n;i++) {
   size_t im=i>0?plane:0,ip=i<n-1?plane:0;float dx0=(i==0||i==n-1)?dx:2.f*dx;
   for(int j=0;j<n;j++) {
    size_t jm=j>0?n:0,jp=j<n-1?n:0;float dx1=(j==0||j==n-1)?dx:2.f*dx;
    for(int k=0;k<n;k++) {
     size_t v=(size_t)i*plane+(size_t)j*n+k;
     float gx=(f[v+ip]-f[v-im])/dx0,gy=(f[v+jp]-f[v-jm])/dx1;
     float gz=k==0?(f[v+1]-f[v])/dx:k==n-1?(f[v]-f[v-1])/dx:(f[v+1]-f[v-1])/(2.f*dx);
     float length=std::sqrt(gx*gx+gy*gy+gz*gz);length=std::fmax(length,1.e-10f);
     float blocked=2.f*(occupied[v]-occ[v]);blocked=std::fmin(1.f,std::fmax(0.f,blocked));
     float weight=(f[v]*(1.f-f[v]))*(1.f-blocked);
     ox[v]=(-gx/length)*weight;oy[v]=(-gy/length)*weight;oz[v]=(-gz/length)*weight;
    }
   }
  }
 }
}
