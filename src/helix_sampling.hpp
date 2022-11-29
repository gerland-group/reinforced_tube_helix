#include <math.h>
#include <limits> //Necessary for numeric_limits
 

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef MAX_DOUBLE
#define MAX_DOUBLE std::numeric_limits<double>::max()
#endif

#ifndef DIMV
#define DIMV    3
#endif


int *helix_sampling(double *r_coords, 
                    int **tri, 
                    double R,
                    double L_body,
                    double turns,
                    double l0, 
                    int N_tri,
                    int N_tot, 
                    int *n_helix);

void find_cylinder_params(double *r, 
                          int Np, 
                          int N_cap, 
                          double *R, 
                          double *L_body);
