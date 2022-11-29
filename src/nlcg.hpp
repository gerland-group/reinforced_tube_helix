#ifndef _NLCG_
#define _NLCG_

#include <stdio.h>
#include <string.h> // form memcpy

#define     DIMV             3                                      /* Number of dimensions per particles (R^3) */



// Non-Linear Conjugate Gradient (NLCG) parameters
#define     MAX_ITS_NLCG_FACT              10                      /* Number of CG steps is a factor of the number of particles */                  
#define     MAX_ITS_LINE_SEARCH            20                      /* Number of line search iterations*/
#define     TOLERANCE_CRIT_LINE_SEARCH     1.0e-3
#define     TOLERANCE_CRIT_NLCG            1.0e-15                 /* Error tolerance of the force in the CG minimisation */
#define     SIGMA_0                        0.05                    /* Factor of the line search (by the average rest length of the bonds)  */   
#define     MEAN_CHECK      5                                      /* Number of iterations the gradient should not have changed to scape */                              

                              
                              
// Function declarations
void nlcg(double *x, 
          int **tri, 
          int **edges,
          int **op_edges,
          double *eql,
          double P, 
          double *lambda_eff,
          double K_SPR, 
          double *k_spr_eff,
          int N_tri,
          int N); 

double mean_sqforce_vertex(double *r, int N);
void movedir(double *v, double *d, double alpha, int N);
double ddotn(double *v1, double *v2, int N);
void changesign(double *v, int N);

double dmean(double *v, int N);


#endif
