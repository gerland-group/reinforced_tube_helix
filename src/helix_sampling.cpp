#include "helix_sampling.hpp"
#include "meshes.hpp"
#include "allvec.hpp"



/******************************************************************************/
/*                                  helix_sampling                            */
/* Finds the particles tracing a helix along the cylinder with coordinates    */
/* r_coords.                                                                  */
/******************************************************************************/
int *helix_sampling(double *r_coords, 
                    int **tri, 
                    double R,
                    double L_body,
                    double turns,
                    double l0, 
                    int N_tri,
                    int N_tot, 
                    int *n_helix)
{
    
    int p;                  // Index of the last vertex inciroirated to the helix
    int Q;                  // Coordination number for the current particle in loop
    int vnn;                // First ring neighbourhood of one of the particles
    
    int *ids_helix;         // Indexes of the particles that define the helix
    int *cn;                // Coordination number of every particle
    int **nnlist;           // First ring neighbourhood of every particle
    
    double t;               // "Time" parameter of the helix
    double s;               // Contour length parametrization
    double dl;              // Effective seperation contour length of the helix
    double L;               // Inverse pitch of the helix
    double a,b,c;           // x,y,z distances between the coordinates and the theoretical helix           
    double chim, tchi;      // Distance between the coordinates and the theoretical helix           
    
    
    L = L_body/(2*PI*turns);
    dl = l0*sqrt(3)/2;
    *n_helix = (int)round(2*PI*turns*sqrt(R*R + L*L)/dl)+1;
    ids_helix = ivector((*n_helix));
    
    cn = coordnum(tri, N_tri, N_tot);
    nnlist = gen_nnlist(tri, cn, N_tri, N_tot);

    // Starting particle
    a = R-r_coords[0];
    b = r_coords[1];
    c = r_coords[2];
    chim = a*a + b*b + c*c;
    ids_helix[0] = 0;
    
    for(int i=1; i<N_tot; i++) {
        a = R-r_coords[i*DIMV];
        b = r_coords[i*DIMV + 1];
        c = r_coords[i*DIMV + 2];
        tchi = a*a + b*b + c*c;
        if(tchi < chim) {
            chim = tchi;
            ids_helix[0] = i;
        }
    }
    
    // Samples the particles nearer the theoretical helix
    int z_pre = 0;
    int n=0;
    while(n < *(n_helix)-1)
    {
    
        // Nearest neigthbours of the last introduced vertex
        p = ids_helix[n]; // Pick the first particle
        Q = cn[p];
        vnn = nnlist[p][0];
        
        n++;
        ids_helix[n] = vnn;
        
        s = n*dl;
        t = s/sqrt(R*R+L*L);
        a = R*cos(t) - r_coords[vnn*DIMV];
        b = R*sin(t) - r_coords[vnn*DIMV+1];
        c = L*t - r_coords[vnn*DIMV+2];
        chim = a*a + b*b + c*c;
            
                
        // Avoid selecting the first particle if is at a lower height
        if(r_coords[vnn*DIMV + 2] < z_pre)
            chim = MAX_DOUBLE;
        
        for(int q=1; q<Q; q++) {
            vnn = nnlist[p][q];
            if(r_coords[vnn*DIMV+2] >= z_pre) {
                    
                a = R*cos(t) - r_coords[vnn*DIMV];
                b = R*sin(t) - r_coords[vnn*DIMV+1];
                c = L*t - r_coords[vnn*DIMV+2];
                tchi = a*a + b*b + c*c;
                
                if(tchi<chim) {
                    ids_helix[n] = vnn;
                    chim = tchi;
                }
                
            }
            
        }
        vnn = ids_helix[n];
        z_pre = r_coords[vnn*DIMV+2];
    } 
    
    free_ivector(cn);
    free_imatrix(nnlist, N_tot);
    
    return ids_helix;
}



/******************************************************************************/
/*                          find_cylinder_params                              */
/* Determines the basis properties of the cylinder knowing that the number    */
/* of veretexes of the caps and that tey are located at the end of the        */
/* coordinates vector r.                                                      */
/******************************************************************************/
void find_cylinder_params(double *r, 
                          int Np, 
                          int N_cap, 
                          double *R, 
                          double *L_body)
{
    
    int N_eff;
    double Rx, Ry;
    N_eff = Np - 2*N_cap;
    Rx=0;
    Ry=0;
    *L_body = 0;
    for(int i=0; i<N_eff; i++){
        if(r[i*DIMV] > Rx)
            Rx = r[i*DIMV];
        if(r[i*DIMV+1] > Ry)
            Ry = r[i*DIMV+1];
        if(r[i*DIMV+2] > *L_body )
            *L_body = r[i*DIMV+2];
    }
    
    *R = 0.5*(Rx+Ry);
}

