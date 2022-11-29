/****************************************************************************

    nlcg: Non-linear conjugate gradient implementation considering pseudocodes 
          and formulas from Shewchuk J.R. (1994), Nocedal J. and Wrigth S. (2000)
           applied to triangulated surfaces.
              
    Copyright (C) 2020  Cesar L. Pastrana

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

****************************************************************************/



#include "nlcg.hpp"
#include "grad.hpp"
#include "allvec.hpp"





/*****************************************************************************************/
/*                                       nlcg                                            */
/* Minimise an objective function with a provided gradient following the non-linear      */
/* conjugate gradient of N particles and  initial coordinates provided in a vector form  */
/* x = [x1, y1, z1, x2, y2, z2, ..., xN, yN, zN].                                        */
/* This implementation is based on the description provided in: An Introduction to the   */
/* Non-linear conjugate gradient method  without the agonizing pain. Shewchuk J.R., 1994 */
/*****************************************************************************************/
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
          int N) 
{
    
    int nedges;                         /* Number of edges in the mesh as determined from the triangulation */
    double L0_CHAR;                     /* Characteristic distance between bonds */

    double delta_d;                     /* Square modulus of the search direction before line search. Used for tolerance crit. */
    double delta0;
    double delta_new, delta_o;          /* Square norm of the gradient before and after a conjugate gradient iteration */
    double delta_mid;                   /* Dot producto between the current configuration gradient and the former */
    
    double eta, eta_prv;                /* Relations for the line search */
    double alpha;                       /* Scaling factor for the line search minimisation */
    double beta;                        /* Scaling factor for the nlcg minimisation (Polak-Riviere)*/
    double mfvf;                        /* Mean square force per particle. Use to scale properties */

    double *dcg;                        /* The search direction for the conjugate gradient */
    double *r;                          /* Gradient for the set of particles */ 
    double *rpi;                        /* Gradient for the set of particles of the previus interaction */ 
    double *dx, *dy, *dz, *d;           /* Distance and distance components between particles */
    
    nedges = 3*N_tri/2;
    
    dcg = dvector(N*DIMV); 
    r   = dvector(N*DIMV); 
    rpi = dvector(N*DIMV); 
    
    dx  = dvector(nedges); 
    dy  = dvector(nedges); 
    dz  = dvector(nedges); 
     d  = dvector(nedges); 
     
     
    // 1. INTIAL CALCULATIONS 
    // 1. a) For the potential (calculate eq. distances of springs)
    calcdr(x, dx, dy, dz,  d, edges, nedges);
    L0_CHAR = dmean(eql, nedges);
    mfvf = 1/(K_SPR*K_SPR*L0_CHAR*L0_CHAR*N);
	    
    // 1. b) For the NLCG
    grad(x, tri, dx, dy, dz, d, edges, op_edges, eql, P, k_spr_eff, lambda_eff, nedges, N_tri, N, r);
    
    changesign(r, DIMV*N);
    memcpy(dcg, r, DIMV*N*sizeof(double));
    memcpy(rpi, r, DIMV*N*sizeof(double));
    
    delta_new =  ddotn(r, dcg, DIMV*N);
    delta0 = delta_new;

    //*****************************************************************************//
    /////////////////     2. NON LINEAR CONJUGATE GRADIENT     //////////////////////
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    int k = 0;              // Controls the reset of conjugate directions
    int itcg = 0;
    int its_meanchk=0;                      // Counts the subset of iteration between 0 and MEAN_CHECK-1 (escape condition)
    double sqgrad_vertex[MEAN_CHECK];       // Squared force at a given NLCG iteration (escape condition)
    while(itcg < MAX_ITS_NLCG_FACT*N*DIMV){
        
        // Line search (Secant method)
        int j = 0;
        
        delta_d = ddotn(dcg, dcg, DIMV*N);
        alpha  = -L0_CHAR*SIGMA_0;
        movedir(x, dcg, SIGMA_0, DIMV*N);
        grad(x, tri, dx, dy, dz, d, edges, op_edges, eql, P, k_spr_eff, lambda_eff, nedges, N_tri, N, r);
        eta_prv =  ddotn(r, dcg, DIMV*N);
        movedir(x, dcg, -SIGMA_0, DIMV*N); // Resets to the original position

        do {
            grad(x, tri, dx, dy, dz, d, edges, op_edges, eql, P, k_spr_eff, lambda_eff, nedges, N_tri, N, r);
            eta = ddotn(r, dcg, DIMV*N);

            if(ddotn(r,r,DIMV*N) == 0.0 || (eta_prv - eta) == 0.0) break;
            alpha *= eta/(eta_prv - eta); 
            movedir(x, dcg, alpha, DIMV*N);
            eta_prv = eta;
            j++;
        } while(j < MAX_ITS_LINE_SEARCH && alpha*alpha*delta_d > TOLERANCE_CRIT_LINE_SEARCH);
        
    
        grad(x, tri, dx, dy, dz, d, edges, op_edges, eql, P, k_spr_eff, lambda_eff, nedges, N_tri, N, r);
        changesign(r ,DIMV*N);
        
        // Polak-Riviere search direction 
        delta_o = delta_new;
        delta_new =  ddotn(r, r, DIMV*N);    if( delta_new  == 0.0) break;
        delta_mid = ddotn(r, rpi, DIMV*N);
        memcpy(rpi, r, DIMV*N*sizeof(double));
       
        beta = (delta_new - delta_mid)/delta_o;
        
        k++;
        if(k == N*DIMV || beta <= 0) {
            beta = 0;
            k=0;
        }
       
        //Update search direction
        for(int i=0; i<DIMV*N; i++)
            dcg[i] = r[i] + beta*dcg[i];
        
        /*        
        sqgrad_vertex[its_meanchk] = mean_sqforce_vertex(r, N);
        if(itcg>2*MEAN_CHECK) {
            if( dmean(sqgrad_vertex, MEAN_CHECK) < TOLERANCE_CRIT_NLCG && 
                sqgrad_vertex[its_meanchk] < TOLERANCE_CRIT_NLCG)
                    break;
        }
        its_meanchk++;
        if(its_meanchk==MEAN_CHECK) { its_meanchk=0; }
        */
        itcg++;
        
    }
    
    /*-----------------------------------------------------------------------------*/
    
    
    /****** Free memory ******/
    
    /* Free vectors */   
    free_dvector(dx);
    free_dvector(dy);
    free_dvector(dz);
    free_dvector(d);
    
    free_dvector(r);
    free_dvector(rpi);
    free_dvector(dcg);
    
}



/*****************************************************************************************/
/*                                    mean_sqforce_vertex                                */
/* Calculates an average square forces per vertex                                        */
/*****************************************************************************************/
double mean_sqforce_vertex(double *r, int N)
{
    double mean_fsq=0;    
    #pragma omp parallel for reduction(+:mean_fsq) num_threads(2)
    for(int i=0; i<N; i++){
        double fx,fy,fz;
        fx = r[i*DIMV];
        fy = r[i*DIMV+1];
        fz = r[i*DIMV+2];
        mean_fsq += fx*fx + fy*fy + fz*fz;
    }
    mean_fsq /= N;
    return mean_fsq;
}

/*****************************************************************************************/
/*                                      movedir                                          */
/* Moves the vector v in the direction given by the vector d and scaled by the factor    */
/* alpha                                                                                 */ 
/*****************************************************************************************/
void movedir(double *v, double *d, double alpha, int N)
{
    for(int i=0; i<N; i++) 
        v[i] += alpha*d[i];
}



/*****************************************************************************************/
/*                                      ddotn                                            */
/* Calculates the dot product of a vector v1 with a vector v2 with dimensions N          */
/*****************************************************************************************/
double ddotn(double *v1, double *v2, int N)
{
    double inp = 0;
    for(int i=0; i<N; i++)
        inp += v1[i]*v2[i];
    return inp;
}


/*****************************************************************************************/
/*                                      changesing                                       */
/* Change the sign of the components in the vector v with size N.                        */
/*****************************************************************************************/
void changesign(double *v, int N)
{
    for(int i=0;i<N;i++)
        v[i] = -v[i];
}


/*****************************************************************************************/
/*                                      dmean                                            */
/* Calculate the mean of the input vector v elements from 0 to the element N-1           */
/*****************************************************************************************/
double dmean(double *v, int N)
{
    double m = 0;
    for(int i=0; i<N; i++)
        m += v[i];
    return m/N;
    
}

