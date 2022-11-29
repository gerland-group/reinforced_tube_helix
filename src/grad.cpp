#include "grad.hpp"
#include "allvec.hpp"
#include "tensor3.hpp"


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
          double *gf)
{
    
    int v1, v2, v3, v4;             // Indexes of the vertexes of the the triangle or the two edges and its oposed vertexes
    double x1, x2, x3;              // x component of the position for every vertex in triangle
    double y1, y2, y3;              // y component of the position for every vertex in triangle
    double z1, z2, z3;              // z component of the position for every vertex in triangle
    double f, fx, fy, fz;           // -forces of the springs (total and components)
    
    double r1[3],r2[3],             // Coordinates of the vectors of the edge and oposes vertexes 
           r3[3],r4[3];             
    double r21[3],r31[3],           // Vector connecting all the vertexes defining the pair of triangles
           r24[3],r34[3],r23[3];
           
    double nI[3], nJ[3];            // Normal vectors for a couple of triangles sharing an edge
    double norm_nI, norm_nJ;        // Norms of the normal vectors
    double gnIbnJ[3],gnJbnI[3];     // Result of the operation (I - n*nT)*u/|n|  
    double f_bend[3];               // Total (or partial) -force contribution for each vertex
    
    
    zeros_dvector(gf, N*DIMV);
    
    //*******************************************************************//
    /////////               1. Pressure gradient               ////////////
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    
    for(int t=0; t<N_tri; t++) {
        v1 = tri[t][0];
        v2 = tri[t][1];
        v3 = tri[t][2];
        
        x1 = r[v1*DIMV];
        x2 = r[v2*DIMV];
        x3 = r[v3*DIMV];
        
        y1 = r[v1*DIMV + 1];
        y2 = r[v2*DIMV + 1];
        y3 = r[v3*DIMV + 1];
        
        z1 = r[v1*DIMV + 2];
        z2 = r[v2*DIMV + 2];
        z3 = r[v3*DIMV + 2];
        
        // Vertex 1
        gf[v1*DIMV]     += P*ONESIXTH*(y2*z3 - y3*z2);
        gf[v1*DIMV + 1] += P*ONESIXTH*(x3*z2 - x2*z3);
        gf[v1*DIMV + 2] += P*ONESIXTH*(x2*y3 - y2*x3);
        
        // Vertex 2
        gf[v2*DIMV]     += P*ONESIXTH*(y3*z1 - y1*z3);
        gf[v2*DIMV + 1] += P*ONESIXTH*(z3*x1 - x3*z1);
        gf[v2*DIMV + 2] += P*ONESIXTH*(x3*y1 - y3*x1);
        
        // Vertex 3
        gf[v3*DIMV]     += P*ONESIXTH*(y1*z2 - y2*z1);
        gf[v3*DIMV + 1] += P*ONESIXTH*(x2*z1 - x1*z2);
        gf[v3*DIMV + 2] += P*ONESIXTH*(x1*y2 - x2*y1);
    }
    //--------------------------------------------------------------//
    
    
    //First, calculate distances between edges in the current configuration
    calcdr(r, dx, dy, dz, d, edges, nedges);
    for(int i=0; i<nedges; i++) { 
        
        //*********************************************************************//
        /////////               2. Stretching gradient               ////////////
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
        v1 = edges[i][0];
        v2 = edges[i][1];
        
        f = k_spr_eff[i]*(d[i] - eql[i])/d[i];
        fx = f*dx[i];
        fy = f*dy[i];
        fz = f*dz[i];
        
        gf[v1*DIMV]     += fx;
        gf[v1*DIMV + 1] += fy;
        gf[v1*DIMV + 2] += fz;
        
        gf[v2*DIMV]     -= fx;
        gf[v2*DIMV + 1] -= fy;
        gf[v2*DIMV + 2] -= fz;
        //-------------------------------------------------------------------//
        
        //*******************************************************************//
        /////////               3. Bending gradient               ////////////
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
        
        v1 = op_vertex[i][0];
        v4 = op_vertex[i][1];
        v2 = edges[i][0];
        v3 = edges[i][1];
        
        dvec_1D_list3(r, v1, r1);
        dvec_1D_list3(r, v2, r2);
        dvec_1D_list3(r, v3, r3);
        dvec_1D_list3(r, v4, r4);
        
        dvdiff3(r3, r1, r31);
        dvdiff3(r2, r1, r21);
        dvdiff3(r2, r4, r24);
        dvdiff3(r3, r4, r34);
        dvdiff3(r2, r3, r23);
        
        // 3.2 Calculate normal vector (cross product and normalization to unity)
        dcross3(r31, r21, nI);
        dcross3(r24, r34, nJ);
        
        norm_nI = dnorm3(nI);
        norm_nJ = dnorm3(nJ);
        
        dvscal3(nI, norm_nI);
        dvscal3(nJ, norm_nJ);
        
        graduniv_byu(nI, norm_nI, nJ, gnIbnJ);
        graduniv_byu(nJ, norm_nJ, nI, gnJbnI);
        
        
        // 3.3 Gradient on each vertex
        // Vertex 1
        dcross3(r23, gnIbnJ, f_bend);
        gf[v1*DIMV]     += lambda_eff[i]*f_bend[0];
        gf[v1*DIMV + 1] += lambda_eff[i]*f_bend[1];
        gf[v1*DIMV + 2] += lambda_eff[i]*f_bend[2];
        
        // Vertex 2
        dcross3(r31, gnIbnJ, f_bend); 
        gf[v2*DIMV]     += lambda_eff[i]*f_bend[0];
        gf[v2*DIMV + 1] += lambda_eff[i]*f_bend[1];
        gf[v2*DIMV + 2] += lambda_eff[i]*f_bend[2];
        
        dcross3(r34, gnJbnI, f_bend); 
        gf[v2*DIMV]     -= lambda_eff[i]*f_bend[0]; // (Minus sign because r34 insted of r43)
        gf[v2*DIMV + 1] -= lambda_eff[i]*f_bend[1];
        gf[v2*DIMV + 2] -= lambda_eff[i]*f_bend[2];
        
        
        // Vertex 3
        dcross3(r21, gnIbnJ, f_bend); 
        gf[v3*DIMV]     -= lambda_eff[i]*f_bend[0]; // (Minus sign because r21 insted of r12)
        gf[v3*DIMV + 1] -= lambda_eff[i]*f_bend[1];
        gf[v3*DIMV + 2] -= lambda_eff[i]*f_bend[2];
        
        dcross3(r24, gnJbnI, f_bend); 
        gf[v3*DIMV]     += lambda_eff[i]*f_bend[0];
        gf[v3*DIMV + 1] += lambda_eff[i]*f_bend[1];
        gf[v3*DIMV + 2] += lambda_eff[i]*f_bend[2];
        
        // Vertex 4
        dcross3(r23, gnJbnI, f_bend); 
        gf[v4*DIMV]     -= lambda_eff[i]*f_bend[0]; // (Minus sign because r32 insted of r23)
        gf[v4*DIMV + 1] -= lambda_eff[i]*f_bend[1];
        gf[v4*DIMV + 2] -= lambda_eff[i]*f_bend[2];
    
        //-------------------------------------------------------------------//        
        
    }
}


/*****************************************************************************************/
/*                                      calcdr                                           */
/* Calculates the distance components and the absolute values of the distances between   */
/* every edge of the mesh                                                                */
/*****************************************************************************************/
void calcdr(double *r, 
            double *dx, 
            double *dy, 
            double *dz, 
            double *d,
            int **edges, 
            int nedges)
{
    
    #pragma omp parallel for num_threads(2)
    for(int i=0; i<nedges; i++) {
        
        int v1, v2;
        v1 = edges[i][0];
        v2 = edges[i][1];
        
        dx[i] = r[v1*DIMV] - r[v2*DIMV];
        dy[i] = r[v1*DIMV + 1] - r[v2*DIMV + 1];
        dz[i] = r[v1*DIMV + 2] - r[v2*DIMV + 2];
        d[i]  = sqrt(dx[i]*dx[i]  + dy[i]*dy[i] + dz[i]*dz[i]);
        
    }
    
}
