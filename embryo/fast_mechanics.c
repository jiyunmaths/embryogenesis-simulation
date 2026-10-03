/* Optional polar mechanics kernel. No fast-math; zero normal boundary flux. */
#include <math.h>
#include <stddef.h>
void polar_flux(const float *f, double *gamma, double *flux, int n,
                double dx, double extent, const double *center,
                const double *p, double base, double contrast, double width) {
    const size_t plane=(size_t)n*n, size=plane*n;
    for (int i=0;i<n;i++) for(int j=0;j<n;j++) for(int k=0;k<n;k++) {
        size_t v=(size_t)i*plane+(size_t)j*n+k;
        double x=((double)i+.5)*dx-extent-center[0];
        double y=((double)j+.5)*dx-extent-center[1];
        double z=((double)k+.5)*dx-extent-center[2];
        double r=sqrt(x*x+y*y+z*z+width*width);
        gamma[v]=base*(1.-contrast*(p[0]*(x/r)+p[1]*(y/r)+p[2]*(z/r)));
    }
    double dx2=dx*dx;
    for (size_t v=0;v<size;v++) {
        int i=v/plane,j=(v/n)%n,k=v%n;
        size_t offsets[3]={plane,(size_t)n,1};int indices[3]={i,j,k};
        double total=0.;
        for(int axis=0;axis<3;axis++) {
            size_t o=offsets[axis];
            if(indices[axis]<n-1) {
                float difference=f[v+o]-f[v];
                total+=.5*(gamma[v]+gamma[v+o])*(double)difference/dx2;
            }
            if(indices[axis]>0) {
                float difference=f[v]-f[v-o];
                total-=.5*(gamma[v-o]+gamma[v])*(double)difference/dx2;
            }
        }
        flux[v]=total;
    }
}
