#include <cmath>



#ifndef DIMV
#define  DIMV               3                               // R^3 for vector arrangment of particles in the matrix
#endif 

#define ONESIXTH    0.166666666666666666666667




void grad(double *r, 
          int **tri,
          double *dx, 
          double *dy,
          double *dz, 
          double *d,
          int **edges,
          int **op_vertex,
          double *eql, 
          double P, 
          double *k_spr_eff, 
          double *lambda_eff,
          int nedges, 
          int N_tri, 
          int N,
          double *gf);

void calcdr(double *r, 
            double *dx,
            double *dy,
            double *dz, 
            double *d,
            int **edges,
            int nedges);





/*****************************************************************************************/
/*                                   graduniv_byu                                        */
/* This function performs the operation [(I- (v_n)*(v_n)^T)^T]*u / |v|, where v_n is the */
/* normalized vector v and u is other vector. This expresio arises after extracting the  */
/* gradient of a unitary vector and is necessary to calculate the gradient of the        */
/* bending energy                                                                        */
/*****************************************************************************************/
inline void graduniv_byu(double *v_uni, double normv, double *u, double *out_v)
{
    double x,y,z;
    double a,b,c;
    x = v_uni[0];
    y = v_uni[1];
    z = v_uni[2];
    
    a = u[0];
    b = u[1];
    c = u[2];

    out_v[0] = (1 - x*x)*a - x*y*b - x*z*c;
    out_v[1] = (1 - y*y)*b - x*y*a - y*z*c;
    out_v[2] = (1 - z*z)*c - x*z*a - y*z*b;
    
    out_v[0] /= normv;
    out_v[1] /= normv;
    out_v[2] /= normv;
}

