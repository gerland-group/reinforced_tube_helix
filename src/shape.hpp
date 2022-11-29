#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream> 
#include <algorithm>
#include <string>
#include <ctime>
#include <vector>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define  DIMV               3                               // R^3 for vector arrangment of particles in the matrix
        
#define  CONVATM            0.101325                        // Conversion factor between atm (input file) and pN/nm^2
#define  YOUNG_MESH_FACTOR  0.866025403784439               // Conversion factor between 2D Young modulus and spring constant


/***************************************** FILE PROPERTIES  *******************************************/
#define INIT_COORDS_FNAME       "init_coords.dat"                //  Path of the initial file to load coardinates
#define MCMIN_COORDS_FNAME      "minim_coords.dat"               //  Minimised coordiantes after NLCG
#define IN_MESH_FNAME           "mesh.dat"                       //  File name to save trajectories (local directory)
#define PARAMETERS_FILE         "params.conf"                    //  Parameters file
#define PRESS_PATH              "./press/"                         // Folder to save pressurisation snapshots
#define PRESS_FNAME             "_press.dat"                       // File name to save pressurisation snapshots
#define MAIN_HELIX_VTXS_FILES   "./helix/main_helix_vertexes.dat"  //  File with the vertexes involved in the main helix
#define HELIX_VTXS_FILES        "./helix/helix_vertexes.dat"        //  File with the vertexes involved in the helix (patches and main)
#define PRESS_LIST_FNAME        "./summary/pressures.dat"                     //  List of pressures used
#define LOG_FILE                "./summary/output.log"              //  Log file name with simulation properties (local directory)



/*********************************   Namespaces and type definitions   *******************************/

// Standard namespace
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::ofstream; 
using std::ifstream; 
using std::string;
using std::getline;
using std::stringstream;
using std::nothrow; // Do not send exception but a null pointer in new

typedef vector<vector<double>> dvector_vl;
typedef vector<vector<int>> ivector_vl;



/***********************************    FUNCTION DECLARATIONS    ************************************/

double avedge_length(double *r, int **edges, int nedges);
void calc_eql(double *r, int **edges, int nedges, double *eql);

void print_gui(int N, int T, double P_atm, double Kappa, double E3D, 
               double h, double K_SPR, double K_fact, double REST_LEN_FACTOR, 
               int layers, double turns);

void save_coordinates(double *r, int Np, string ull_path);
