/* Fused float32 surface-cue arithmetic with float64 accumulation.
   This differs in reduction order from NumPy and requires trajectory validation. */
#include <math.h>
#include <stddef.h>
void surface_cue(const float *phi,const float *h,const float *occupied,
                 const float *shell_sums,double *out,int cells,int n,float dx) {
    size_t plane=(size_t)n*n, size=plane*n;
    for(int cell=0;cell<cells;cell++) {
        const float *f=phi+cell*size,*occ=h+cell*size;
        double sums[3]={0.,0.,0.};
        for(int i=0;i<n;i++) {
            size_t im=i>0?plane:0,ip=i<n-1?plane:0;
            float dx0=(i==0 || i==n-1)?dx:2.f*dx;
            for(int j=0;j<n;j++) {
                size_t jm=j>0?n:0,jp=j<n-1?n:0;
                float dx1=(j==0 || j==n-1)?dx:2.f*dx;
                for(int k=0;k<n;k++) {
                    size_t v=(size_t)i*plane+(size_t)j*n+k;
                    float gx=(f[v+ip]-f[v-im])/dx0;
                    float gy=(f[v+jp]-f[v-jm])/dx1;
                    float gz=k==0?(f[v+1]-f[v])/dx:k==n-1?(f[v]-f[v-1])/dx:(f[v+1]-f[v-1])/(2.f*dx);
                    float length=sqrtf(gx*gx+gy*gy+gz*gz);
                    length=fmaxf(length,1.e-10f);
                    float blocked=2.f*(occupied[v]-occ[v]);blocked=fminf(1.f,fmaxf(0.f,blocked));
                    float shell=f[v]*(1.f-f[v]);float weight=shell*(1.f-blocked);
                    sums[0]+=(double)((-gx/length)*weight);
                    sums[1]+=(double)((-gy/length)*weight);
                    sums[2]+=(double)((-gz/length)*weight);
                }
            }
        }
        double denom=fmax((double)shell_sums[cell],1.e-12);
        double norm=0.;for(int a=0;a<3;a++){out[cell*3+a]=sums[a]/denom;norm+=out[cell*3+a]*out[cell*3+a];}
        if(sqrt(norm)<1.e-6)for(int a=0;a<3;a++)out[cell*3+a]=0.;
    }
}
