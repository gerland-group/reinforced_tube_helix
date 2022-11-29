#ifndef _TENSOR3_
#define _TENSOR3_

#include <math.h>




void dvdiff3(double *v1, double *v2, double *v_out);
void dmidvec3(double *v1, double *v2, double *v_out);
void dvscal3(double *v, double A);
void dassignv3(double *v_in, double *v_out, double s);
void dvec_list3(double **list, int n, double *v);
void dvec_1D_list3(double *list, int n, double *v);



/****************************************************************/
/*                             ddot3                            */
/* Calculates the dot product between vectors v1 and v2         */
/****************************************************************/
inline double ddot3(double *v1, double *v2)
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}



/****************************************************************/
/*                          dcross3                             */
/* Return a vector of the cross product between v and u         */
/****************************************************************/
inline void dcross3(double *v, double *u, double *v_out)
{
    v_out[0] = v[1]*u[2] - v[2]*u[1];
    v_out[1] = v[2]*u[0] - v[0]*u[2];
    v_out[2] = v[0]*u[1] - v[1]*u[0]; 
}

/****************************************************************/
/*                              dnorm3                          */
/* Calculate the norm of the vector v                           */
/****************************************************************/
inline double dnorm3(double *v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}



#endif
