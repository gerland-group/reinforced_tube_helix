#include "curvatures.hpp"
#include "tensor3.hpp"


/************************************************************************************************/
/*                                         A_mix                                                */
/*                                                                                              */
/* Calculates the Area mix function (Fig4 in Meyer et al), ie. a Voronoi safe for obtuse        */
/* triangles, on a given triangle of a given particle of interest. Takes as input the           */
/* coordinates of the particles r, the triangle list tri, the particle of interest and the      */
/* triangle id.                                                                                 */
/************************************************************************************************/
double A_mix(double *r, int **tri, int p, int trid)
{
    int v1, v2, v3;
        
    double S = 0;
    double tt[3];
    tt[0] = tri[trid][0];
    tt[1] = tri[trid][1];
    tt[2] = tri[trid][2];
    
    
    // Organize the edges to position of the particle p  on the first position
    if(tt[1] == p) {
        tt[0] = tt[0] + tt[1];
        tt[1] = tt[0] - tt[1];
        tt[0] = tt[0] - tt[1];
    }
    else if (tt[2] == p) {
        tt[0] = tt[0] + tt[2];
        tt[2] = tt[0] - tt[2];
        tt[0] = tt[0] - tt[2];
    }
    
    double r1[3],   r2[3],  r3[3];
    double r21[3], r31[3],  r32[3];
    double tv[3];
    
    double norm_r21, norm_r31, norm_r32;
    
    v1 = tt[0];
    v2 = tt[1];
    v3 = tt[2];
    
    dvec_1D_list3(r, v1, r1);
    dvec_1D_list3(r, v2, r2);
    dvec_1D_list3(r, v3, r3);
    
    dvdiff3(r2,r1, r21);
    dvdiff3(r3,r1, r31);
    dvdiff3(r3,r2, r32);
    
    norm_r21 = dnorm3(r21);
    norm_r31 = dnorm3(r31);
    norm_r32 = dnorm3(r32);

    
    // Determine the Voronoi area or the corrected area if there are obtuse angles
    double cos_alpha;
    double cos_beta, beta;
    double cos_gamma, gamma;
    
    cos_alpha = ddot3(r21, r31)/(norm_r21*norm_r31);
    
    r21[0] *= -1; r21[1] *= -1; r21[2] *= -1;
    cos_beta  = ddot3(r21, r32)/(norm_r21*norm_r32);  
    
    r31[0] *= -1; r31[1] *= -1; r31[2] *= -1;
    r32[0] *= -1; r32[1] *= -1; r32[2] *= -1;
    cos_gamma = ddot3(r32, r31)/(norm_r32*norm_r31);
    
    if( cos_alpha < 0 || cos_beta < 0 || cos_gamma < 0 ) {
        r21[0] *= -1; r21[1] *= -1; r21[2] *= -1;
        dcross3(r21, r32, tv);
        S = dnorm3(tv)/2.0;
        
        // Obtuse angles
        if(cos_alpha < 0)
            S /= 2.0;
        else
            S /= 4.0;
        
    } else { // Voronoi safe
        beta = acos(cos_beta);
        gamma = acos(cos_gamma);
        S = (norm_r21*norm_r21/tan(gamma) + norm_r31*norm_r31/tan(beta) )/8.0;
    }

    return S;
}


