# Cesar L. Pastrana, 2020
import numpy as np
import subprocess
import sys


def det_ncaps():
    """
        det_ncaps:
        Explores the log file to determine the number of particles per cap
    """
    n_caps = 0
    with open('./out_mesh/spherocylinder.log') as f:
        for line in f:
            id=line.split('=')
            if id[0].replace(" ", "")== "[N_PER_CAP]":
                n_caps = int(id[1].replace(" ", ""))
                break
         
    return n_caps



if __name__ == "__main__":

    # Checks
    if len(sys.argv) != 4:
        print("Incorrect arguments. Execute as genspherocyl radius length n_particles")
        sys.exit(0)
        
    r = np.float64(sys.argv[1])
    L = np.float64(sys.argv[2])
    N = np.float64(sys.argv[3])
    
    
    # Generate object
    h = np.sqrt(N*L/(np.sqrt(3.0)*r*np.pi))
    l = 2.0*L/(h*np.sqrt(3))
    w = 2.0*np.pi*r/l
    
    h=int(round(h))
    w=int(round(w))

    print("1. Generate particles")
    subprocess.run(["./bin/genspherocyl", 
                    str(w),
                    str(h),
                    str(l), 
                    str(0),
                    str(0)])
    
    
    # Rotate and offset
    r = np.loadtxt('./out_mesh/init_coords.dat', delimiter='\t',)
    print("2. Rotate and offset")
    
    R_rot = [ [1,0,0], 
              [0,0,-1],
              [0,1,0]]
    R_rot = np.array(R_rot)
    r_rot = np.matmul(R_rot,r.T).T
    
    n_caps = det_ncaps()
    m = min(r_rot[0:-2*n_caps,2])
    r_rot[:,2]-=m
    
    
    # Produce mesh
    print("3. Meshing")
    from meshgen.meshgen import gentrimesh, plot_tri
    tri, n_tri = gentrimesh(r_rot)
    plot_tri(r_rot, tri)
    
    
    # Saving
    print("4. Saving")
    np.savetxt('./out_mesh/init_coords.dat', r_rot, fmt='%f', delimiter='\t')
    np.savetxt('./out_mesh/mesh.dat', tri, fmt='%d', delimiter='\t')
    
