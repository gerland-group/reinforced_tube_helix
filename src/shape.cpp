/*******************************************************************************

    shape:   Non-linear conjugate gradient minimisation of the energy function
             of a pressurized triangularized closed surface
             
    Copyright (C) 2020, Cesar L. Pastrana

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
    
    
    Compile with:
    
    g++ shape.cpp nlcg.cpp grad.cpp meshes.cpp allvec.cpp tensor3.cpp 
        helix_sampling.cpp cpputils.cpp 
        -O2 -lm  -std=c++11 -o ../bin/shape

*******************************************************************************/

#include "shape.hpp"
#include "nlcg.hpp"
#include "meshes.hpp"
#include "curvatures.hpp"
#include "helix_sampling.hpp"
#include "allvec.hpp"
#include "cpputils.hpp"



int main()
{
    
    dvector_vl r_mtx;                       // Coordinates xyz of each particle (from file) [nm]        
    ivector_vl tri_vl;                      // Array of triangles defining the surface in vector libraru form
    
    int Np;                                 // Number of vertexex/particles in the mesh
    int N_tri;                              // Number of faces
   
    double P;                               // Pressure parameter (value loaded from file) [atm -> pN nm]
    double DP;                              // Variation of pressure (value loaded from file)  [atm -> pN nm]
    double KMF;                             // Correction factor of the mesh (value loaded from file)
    double KAPPA;                           // Bending stiffness (obtained from E,h and Poisson ratio) [pN nm]
    double LAMBDA;                          // Effective bending stiffness in the triangular mesh [pN nm]
    double E_3D;                            // 3D Young modulus (value loaded from file) [MPa -> pN/nm^2]
    double h;                               // Shell thickness  (value loaded from file) [nm]
    double Y_2D;                            // 2D Young modulus [N/m -> pN/nm]
    double K_SPR;                           // Stretching stiffness [pN/nm]
    double K_SPR_FACTOR;                    // Factor of K_SPR to the regions of interest (value loaded from file)
    double REST_LEN_FACTOR;                 // Factor of the rest length along the helix
    double TURNS;                           // Number of turns for the helix (value loaded from file)
    int N_CAP;                              // Number of particles of each caps (value loaded from file)
    int LAYERS;                             // Concentric rings for the helix (value loaded from file)
    
    
    //************************************************************************************//
    // 1. Load data                                                                       //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

    // 1.1  Load the file with the initial coordiantes of the particles
    ifstream infile_coords(INIT_COORDS_FNAME);
    string line;
    int l=0;                                 // Line control for data loading
    if(infile_coords.is_open()) {
        while (getline(infile_coords, line))
        {
            double val;
            stringstream ss(line);
            r_mtx.push_back(vector<double>(0));

            while (ss >> val)
                r_mtx[l].push_back(val);
            
            ++l;
        }
        infile_coords.close();
    } else {
        cout << "Error while opening the coordinates file " << INIT_COORDS_FNAME << ". Exiting." << endl;
        exit(EXIT_FAILURE);
    }
        
    // 1.2  Load the file with the mesh data
    ifstream infile_tri(IN_MESH_FNAME);
    l = 0;
    if((infile_tri.is_open()))
    {
        while (getline(infile_tri, line))
        {
            int val;
            stringstream ss(line);
            tri_vl.push_back(vector<int>(0));

            while (ss >> val)
                tri_vl[l].push_back(val);
            
            ++l;
        }
        infile_tri.close();
    } else {
        cout << "Error while opening the mesh file " << IN_MESH_FNAME << ". Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
    
    // 1.3 Load the simulation parameters values 
    DP = 0;
    P = 0;
    E_3D = 30;
    h = 2;
    K_SPR = 50;
    KMF = 2/sqrt(3);
    LAMBDA = E_3D*h*h*h*KMF; 
    K_SPR_FACTOR = 1;
    REST_LEN_FACTOR = 1;
    TURNS = 0;
    LAYERS = 0;
    N_CAP = 0;
    
    // std::ifstream is RAII, i.e. no need to call close
    ifstream infile_params(PARAMETERS_FILE);
    if (infile_params.is_open())
    {
        while (getline(infile_params, line)) 
        {
            line.erase(remove_if(line.begin(), line.end(), ::isspace),line.end());
            
            // Commments with the sharp symbol
            if (line[0] == '#' || line.empty()) continue; 

            auto delim_pos = line.find("=");
            auto tag = line.substr(0, delim_pos);
            auto val = line.substr(delim_pos + 1);

            //Custom coding
            if (tag == "PRESSURE") P = -CONVATM*(std::stod(val));       // Negative
            else if (tag == "DP") DP = -CONVATM*(std::stod(val));       // Negative        
            else if (tag == "KAPPA_MESH_FACTOR") KMF = std::stod(val);
            else if (tag == "3D_YOUNG_MODULUS") E_3D = std::stod(val);
            else if (tag == "SHELL_THICKNESS") h = std::stod(val);
            else if (tag == "K_SPR_FACTOR") K_SPR_FACTOR = std::stod(val);
            else if (tag == "REST_LEN_FACTOR") REST_LEN_FACTOR = std::stod(val);
            else if (tag == "HELIX_TURNS") TURNS = std::stod(val);
            else if (tag == "N_CAP") N_CAP = std::stoi(val);
            else if (tag == "LAYERS") LAYERS = std::stoi(val);
        }
        
        // In principle not necessary beucase ifstream self cleans
        infile_params.close();
    }
    else 
    {
        cout << "Error while opening the parameters file " << PARAMETERS_FILE << ". Exiting." << endl;
        exit(EXIT_FAILURE);
    }
    
   
    
    // 1.4  Calculates basic properties and change the structure of the allocated arrays
    Np = r_mtx.size(); 
    N_tri = tri_vl.size();
    Y_2D = E_3D*h;
    K_SPR = Y_2D*YOUNG_MESH_FACTOR;
    KAPPA = Y_2D*h*h/(12.0*(1.0-1.0/9.0));
    LAMBDA = KAPPA*KMF;
    

    // Copy the set of particles to a vector form r = [x1, y1, z1, x2, y2, z2,...,  xN, yN, zN]
    double *r;
    r = dvector(Np*DIMV);
    
    for(int i = 0; i<Np; i++) {
        r[i*DIMV] = r_mtx[i][0];
        r[i*DIMV + 1] = r_mtx[i][1];
        r[i*DIMV + 2] = r_mtx[i][2];
    }
    r_mtx.clear();
    
    
    // The mesh is transformed from vector to a "new" object
    int **tri;
    tri = imatrix(N_tri,3);
    for(int t = 0; t<N_tri; t++) {
        tri[t][0] = tri_vl[t][0];
        tri[t][1] = tri_vl[t][1];
        tri[t][2] = tri_vl[t][2];
    }
    tri_vl.clear();
    
    
    print_gui(Np, N_tri, P, KAPPA, E_3D, h, K_SPR, 
              K_SPR_FACTOR, REST_LEN_FACTOR, LAYERS, TURNS);
   
    //------------------------------------------------------------------------------------//
    
    //************************************************************************************//
    // 3. Extract parameters and topology of the mesh                                     //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//   
    
    cout << endl << " - Analysis mesh properties... ";
    
    int nedges;                             // Number of edges (determined from the triangulation)
    int **edges;                            // Array with the pair of edges
    int **op_edges;                         // Pair of vertexes defining the oposed vertexes to each edge
    int *cn;                                // Coordination number of each vertex
    int **nnlist;                           // Nearest(first-ring) neigthbours list for each vertex

    double *eql;                            // Rest length of the springs
    double l0;                              // Average distance between vertexes
    
    // Edges properties
    nedges = 3*N_tri/2;
    edges = find_edges(tri, N_tri);
    op_edges = find_oposed_vertexes(tri, N_tri, edges, nedges);
    cn = coordnum(tri, N_tri, Np);
    nnlist = gen_nnlist(tri, cn, N_tri, Np);
    
    // Equilibrium length
    eql = dvector(nedges);
    calc_eql(r, edges, nedges, eql);
    l0 = avedge_length(r, edges, nedges);
    
    //------------------------------------------------------------------------------------//
    
    
    
    //************************************************************************************//
    // 4. Helix sampling                                                                  //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    cout << "Done" << endl << " - Helix determination... "; 
    
    // Calculate parameters of the helix and find particles following the helical path
    int n_helix;                    // Particles main helix;
    double R;                       // Radius helix
    double L_body;                  // Length of the main body
    int *ids_helix;
    
    
    find_cylinder_params(r, Np, N_CAP, &R, &L_body);
    ids_helix = helix_sampling(r, tri, R, L_body, TURNS, l0, N_tri, Np, &n_helix);
    
    // Attribute an specific stretching stiffness (considering patches)
    int nvtx_helix;                 // Particles involved on the helix
    int nedges_helix;               // Total number of vertex in the helix
    int patch_size;                 // Temp size of the patch for the current particle in loop
    int *mod_kstr;                  // Vertex selected with a different stretching
    int *patch_vtx;                 // Vertexes in a patch of a given particle
    double *k_spr_eff;              // Stretching stiffness for every edge modified for helix domains
    double *eql_eff;                // Rest length of the helix domains
    double *lambda_eff;             // Effective bending stiffness of the reinforced regions

    // Default assignment of stretching stiffness and rest lengths
    k_spr_eff = dvector(nedges);
    eql_eff = dvector(nedges);
    lambda_eff = dvector(nedges);
    for(int i=0; i<nedges; i++) {
        k_spr_eff[i] = K_SPR;
        lambda_eff[i] = LAMBDA;
        eql_eff[i] = eql[i];
    }
    
    
    // Factorize every edge belonging to the patch around the helix particles   
    mod_kstr = ivector(Np);
    zeros_ivector(mod_kstr, Np);
    patch_vtx = ivector(nedges);
    patch_size = find_patch(0, nnlist, cn, patch_vtx, LAYERS);
    
    int v;
    for(int i=0; i<n_helix; i++) {
        
        v = ids_helix[i];
        patch_size = find_patch(v, nnlist, cn, patch_vtx, LAYERS);
        
        for(int j=0; j<patch_size; j++) {
            v = patch_vtx[j];
            if(v<Np-2*N_CAP) { // Discard the caps
                mod_kstr[v] = 1;
                for(int e=0; e<nedges; e++) {
                    if(edges[e][0] == v || edges[e][1] == v)  {
                        k_spr_eff[e] = K_SPR_FACTOR*K_SPR; 
                        lambda_eff[e] = K_SPR_FACTOR*LAMBDA; 
                        eql_eff[e] = eql[e]*REST_LEN_FACTOR;
                    }
                }
            }
        }
    }
    
    
    // Determine the total number of edges and particles involved on the helix reinforcement
    nedges_helix=0;
    for(int i=0; i<nedges; i++) {
        if(k_spr_eff[i] != K_SPR)
            nedges_helix++;
    }
    
    nvtx_helix = 0;
    for(int i=0;i<Np; i++){
        if(mod_kstr[i] == 1)
            nvtx_helix++;
    }
    
    
    
    // Determine the total area, the reinforced area and the area exluding the caps from Voronoi areas
    double S_helix, S_total, S_total_no_caps;
    double *voronoi_areas;
    voronoi_areas = dvector(Np);
    zeros_dvector(voronoi_areas, Np);
    S_total = 0; 
    S_helix = 0;
    S_total_no_caps = 0;
    
    for(int i=0; i<Np; i++) {
        voronoi_areas[i] = 0;
        for(int t=0; t<N_tri; t++){
            if(tri[t][0] == i || tri[t][1] == i || tri[t][2] == i) 
                voronoi_areas[i] += A_mix(r, tri, i, t);;
        }
        
        if(mod_kstr[i] == 1)
            S_helix += voronoi_areas[i];
        
        if(i<Np-2*N_CAP)
            S_total_no_caps += voronoi_areas[i];

        S_total += voronoi_areas[i];
    }
    
    free_dvector(voronoi_areas);
    free_ivector(cn);
    free_ivector(patch_vtx); 
    free_imatrix(nnlist, Np);

    
    //------------------------------------------------------------------------------------//
    
    
    
    //************************************************************************************//
    // 5. Non-Linear Conjugate Gradient Minimisation                                      //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    cout << "Done" << endl << " - MINIMISATION RUNNING... "; 
    double miniml0;
    double exec_time;
    time_t now, start, end;   
    char *currt;
    
    now = time(0);
    currt = ctime(&now);
    time(&start); 
    
    string full_path_coords;  
    string press_state;
    string full_path_press;   
    ofstream fid_pres_list;
    fid_pres_list.open(PRESS_LIST_FNAME, ofstream::out);
    fid_pres_list << "step" << "\t"  << "pressure(atm)" << endl;
    
    //  Minimisation in ramp 
    int file_count=0;
    double tp=0;
    if(P < 0) {
        while(tp >= P-DP) {
            nlcg(r, tri, edges, op_edges, eql_eff, tp, lambda_eff, K_SPR, k_spr_eff, N_tri, Np);  
            
            press_state = int2str(file_count);
            full_path_press = strcat(PRESS_PATH, press_state);
            full_path_press = strcat(full_path_press, PRESS_FNAME); 
            save_coordinates(r, Np, full_path_press);
            fid_pres_list << file_count << "\t"  << tp/CONVATM << endl;
            
            tp += DP;
            file_count++;
        }
        nlcg(r, tri, edges, op_edges, eql_eff, P, lambda_eff, K_SPR, k_spr_eff, N_tri, Np);    
        
        press_state = int2str(file_count);
        full_path_press = strcat(PRESS_PATH, press_state);
        full_path_press = strcat(full_path_press, PRESS_FNAME); 
        save_coordinates(r, Np, full_path_press);
        fid_pres_list << file_count << "\t"  << P << endl;
        
    } else if (P==0) {
        nlcg(r, tri, edges, op_edges, eql_eff, P, lambda_eff, K_SPR, k_spr_eff, N_tri, Np);  
    } else {
        while(tp <= P-DP) {
            nlcg(r, tri, edges, op_edges, eql_eff, tp, lambda_eff, K_SPR, k_spr_eff, N_tri, Np);  
            
            press_state = int2str(file_count);
            full_path_press = strcat(PRESS_PATH, press_state);
            full_path_press = strcat(full_path_press, PRESS_FNAME); 
            save_coordinates(r, Np, full_path_press);
            fid_pres_list << file_count << "\t"  << tp << endl;
            
            file_count++;
            tp += DP;
        }
        nlcg(r, tri, edges, op_edges, eql_eff, P, lambda_eff, K_SPR, k_spr_eff, N_tri, Np);  
        
        press_state = int2str(file_count);
        full_path_press = strcat(PRESS_PATH, press_state);
        full_path_press = strcat(full_path_press, PRESS_FNAME); 
        save_coordinates(r, Np, full_path_press);
        fid_pres_list << file_count << "\t"  << P << endl;
               
    }
    fid_pres_list.close();
    

    
    
    miniml0 = avedge_length(r, edges, nedges);
    
    time(&end); 
    exec_time = double(end-start);
    
    
    //-------------------------------------------------------------------------------------//
    
    
    //************************************************************************************//
    // 6. Save data                                                                       //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    
    // Save the final minimised coordiantes
    save_coordinates(r, Np, MCMIN_COORDS_FNAME);
    
    
    // Save the vertexes of the main helix
    ofstream fid_main_helix_vertex;
    fid_main_helix_vertex.open (MAIN_HELIX_VTXS_FILES, ofstream::out);
    for(int i=0; i<n_helix; i++)
        fid_main_helix_vertex << ids_helix[i] << endl;
    fid_main_helix_vertex.close();
    
    
    // Save the vertexes of the helix domain
    ofstream fid_helix_vertex;
    fid_helix_vertex.open (HELIX_VTXS_FILES, ofstream::out);
    for(int i=0; i<Np; i++) {
        if(mod_kstr[i] == 1)
            fid_helix_vertex << i << endl;
    }
    fid_helix_vertex.close();
    
        
    // Save LOG file with simulation parameters.
    ofstream fid_logfile;
    fid_logfile.open (LOG_FILE, ofstream::out);
    fid_logfile << "[NLCG_VERSION] = " << 5.1 << endl;
    fid_logfile << "[DATE] = " << currt << endl;
    fid_logfile << "[EXEC_TIME](s) = " << exec_time << endl;
    fid_logfile << "[#PARTICLES] = " << Np << endl;
    fid_logfile << "[#FACES] = " << N_tri << endl;
    fid_logfile << "[#EDGES] = " << nedges << endl;
    fid_logfile << "[CHARACTERISTIC_REST_LENGTH]_0(nm) = " << l0 << endl;
    fid_logfile << "[CHARACTERISTIC_REST_LENGTH]_minim(nm) = " << miniml0 << endl;
    fid_logfile << "[PRESSURE](atm) = " << -P/CONVATM << endl;
    fid_logfile << "[DP](atm) = " << -DP/CONVATM << endl;
    fid_logfile << "[BENDING_STIFNESS](pN nm) = " << KAPPA << endl;
    fid_logfile << "[BENDING_STIFNESS_MESH](pN nm) = " << LAMBDA << endl;
    fid_logfile << "[3D_YOUNG_MODULUS](MPa) = " << E_3D << endl;
    fid_logfile << "[2D_YOUNG_MODULUS](pN/nm) = " << Y_2D << endl;
    fid_logfile << "[SHELL_THICKNESS](pN/nm) = " << h << endl;
    fid_logfile << "[STRETCHING_STIFFNESS](pN/nm) = " << K_SPR << endl;
    fid_logfile << "[STRETCHING_STIFFNESS_HELIX](pN/nm) = " << K_SPR*K_SPR_FACTOR << endl;
    fid_logfile << "[REST_LEN_FACTOR] = " << REST_LEN_FACTOR << endl;
    fid_logfile << "[HELIX_TURNS] = " << TURNS << endl;
    fid_logfile << "[HELIX_RADIUS](nm) = " << R << endl;
    fid_logfile << "[HELIX_PITCH](nm/rad) = " << L_body/(2*PI*TURNS) << endl;
    fid_logfile << "[#HELIX_LAYERS] = " << LAYERS << endl;
    fid_logfile << "[#PARTICLES_HELIX] = " << nvtx_helix << endl;
    fid_logfile << "[#EDGES_HELIX] = " << nedges_helix << endl; 
    fid_logfile << "[#PARTICLES_PER_CAP] = " << N_CAP << endl;
    fid_logfile << "[BODY_LENGTH](nm) = " << L_body << endl; 
    fid_logfile << "[TOTAL_AREA](nm^2) = " << S_total << endl;
    fid_logfile << "[TOTAL_AREA_NO_CAPS](nm^2) = " << S_total_no_caps << endl;
    fid_logfile << "[AREA_HELIX](nm^2) = " << S_helix << endl;
    fid_logfile.close();
    
    
    cout << "Done" <<  endl << "  (execution time: " << exec_time << " seconds)" << endl << endl;
    
    
    //-------------------------------------------------------------------------------------//
    
    
    // Free memmory
    free_ivector(ids_helix);
    free_ivector(mod_kstr);
    free_dvector(r);  
    free_dvector(k_spr_eff);
    free_dvector(eql_eff);
    free_dvector(lambda_eff);
    free_imatrix(tri, N_tri);
    free_imatrix(edges, nedges);

    
    return 0;
}



/*****************************************************************************************/
/*                                      calc_eql                                         */
/* This function calculates the rest length of every spring and sets on eql. Used to     */
/* determine the rest length of the springs                                              */
/*****************************************************************************************/
void calc_eql(double *r, int **edges, int nedges, double *eql)
{
    int v1,v2;
    double x,y,z;
    for(int i=0; i<nedges; i++) {
        v1 = edges[i][0];
        v2 = edges[i][1];
        x = r[v1*DIMV] - r[v2*DIMV];
        y = r[v1*DIMV + 1] - r[v2*DIMV + 1];
        z = r[v1*DIMV + 2] - r[v2*DIMV + 2];
        eql[i] = sqrt(x*x + y*y + z*z);
    }
}


/*****************************************************************************************/
/*                                  avedge_length                                        */
/* Calculates average rest length of the bonds in the configuration r                    */
/*****************************************************************************************/
double avedge_length(double *r, int **edges, int nedges)
{
    int v1, v2;
    double x, y, z;
    double L = 0;
    for(int i=0; i<nedges; i++) {
        v1 = edges[i][0];
        v2 = edges[i][1];
        x = r[v1*DIMV] - r[v2*DIMV];
        y = r[v1*DIMV + 1] - r[v2*DIMV + 1];
        z = r[v1*DIMV + 2] - r[v2*DIMV + 2];
        L += sqrt(x*x + y*y + z*z);
    }
    return L/nedges;
}



/*****************************************************************************************/
/*                                      print_gui                                        */
/* Prints the initial screen of indicating the parameters used                           */
/*****************************************************************************************/
void print_gui(int N, int T, double P_atm, double Kappa, double E3D, 
               double h, double K_SPR, double K_fact, double REST_LEN_FACTOR, 
               int layers, double turns)
{
    
    cout << endl << "P R E S S U R I Z E D   M E S H : N L C G   M I N I M I S A T I O N" << endl; 
    cout << "---------------------------------------------------------------------" << endl << endl;
    cout << " - Particles: " <<  N << endl;
    cout << " - Triangles: " << T << endl;
    cout << " - Pressure (p): " <<  P_atm/CONVATM << " atm"<< endl;
    cout << " - Bending stiffness (Kappa): " <<  Kappa << " pN nm"<< endl;
    cout << " - 3D Young modulus (E): " <<  E3D << " MPa"<< endl;
    cout << " - 2D Young modulus (Y): " <<  E3D*h << " pN/nm"<< endl;
    cout << " - Shell thickness (h): " <<  h << " nm"<< endl;
    cout << " - Stretching stiffness (K_SPR): " <<  K_SPR << " pN/nm" << endl;
    cout << " - Stretching stiffness helix (K_SPR_FACTOR): " <<  K_fact << endl;
    cout << " - Rest length bonds helix (REST_LEN_FACTOR): " <<  REST_LEN_FACTOR << endl;
    cout << " - Helix turns: " <<  turns << endl;
    cout << " - Layers: " <<  layers << endl;
   
}




/*****************************************************************************************/
/*                                  save_coordinates                                     */
/* Save the coordinates in a friendly matrix arrangment                                  */
/*****************************************************************************************/
void save_coordinates(double *r, int Np, string full_path)
{
    ofstream fid_coords;
    fid_coords.open (full_path, ofstream::out);
    for(int i = 0; i<Np; i++) 
        fid_coords << r[i*DIMV] << "\t" << r[i*DIMV + 1] << "\t" << r[i*DIMV + 2] << endl;
    fid_coords.close();
}

