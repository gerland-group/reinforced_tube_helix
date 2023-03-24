# Pressurized spherocylinders with helical reinforcement

Cesar L. Pastrana, 2022 (c)

This software was used in the manuscript "Pressure-induced shape-shifting of helical bacteria" by Pastrana C.L, Qiu L., Armon S., Gerland U., and Amir A.
[Soft Matter](https://pubs.rsc.org/en/content/articlelanding/2023/SM/D2SM01044E, 19(12) 2023: 2224-2230).

## Description
The surface is described using a discrete bead and spring model. The total energy of the system at a pressure $p$ is expressed as:

$$ E = \frac{k_s}{2}\sum_{\langle i,j\rangle} (r_{ij} - r_{ij}^0)^2 +  k_b\sum_{\langle\alpha,\beta\rangle}(1 - \hat{n}_{\alpha}\cdot\hat{n}_{\beta}) -pV$$

where the first sum runs over all pairs of connected notes $i,j$ constituting an edge of the mesh and $r_{ij} = |r_i - r_j|$; and the second sum is over pair of triangles $\alpha,\beta$ sharing an edge, $\hat{n}$ are their respective normal vectors and $k_b$ is a bending stiffness. For a detailed description of the relation between the discrete $k_x$ and the continum variables (Young modulus, thickness and Poisson ratio) see Soung H.S. and Nelson D, R. \([Phys. Rev. A 38,1005--1018, 1988](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.38.1005)\).

In this work, the relax configuration used as input is a spherocylinder. A helical domain with the desired properties is reinforced to be nearly undeformable.

The non-linear conjugate gradient algorithm is used to find the minimal energy configuration.


## Compilation
A Makefile is provided. To build, simply execute:

```
make
```
in the main folder. This will generate the compiled code `shape` in the `bin` folder.

## Execution
The simulation is launched by simply executing `./shape` in the `bin` folder.
Pressurisation proceeds in a step-wise fashion from $0$ to a target pressure $p$ atm, steps of $\Delta p$ provided by the user. 


#### Input
The program requires as input the following files:

- `init_coords.dat`: Array of $N\times 3$ with the $N$ vertices coordinates defining the surface in the relax configuration. The units of the input coordiantes coordinates are nm. 
- `mesh.dat`: Array of $T\times 3$ with the $T$ triangles defining the connectivity of the triangulated surface. 
- `params.conf`: This file contains the simulation parameters. The meaning of each parameter is indicated in commented blocks in the file.
  The main parameters are the target pressure $p$ (`PRESSURE`), the step $\Delta p$ (`DP`) and the reinforcement factor $K$ (`K_SPR_FACTOR`).

For the appropiate determination of the helical region to be reinforced, the spherocilinder should have the main body oriented along the $z$ axis. Moreover, the lowest part of the cylinder should be at the $z$ = 0. The coordinates of the hemispherical caps are located at the end of the `init_coords.dat` file and the number of particles per cap needs to be specified in the variable `N_CAP` of the file `params.conf`.  Spherocyilinders with custom length and radius can be be generated using the accesory Python script in `src/genspherocyl/genspherocyl.py`.


#### Output
The following output files are generated after execution of `shape`:

- `./minim_coords.dat`:  Coordinate at the target pressure.
- `./press/X_press.dat`: Coordinates for each step of pressurization.
- `./summary/pressures.dat`: Pressures in each step.
- `./summary/output.dat`: Summaary of mechanical and mesh parameters employed.
- `./helix/main_helix_vertexes.dat`: Indices of the verticies defining the main helix (layer zero)
- `./helix/helix_vertexes.dat`: Every index involved in the reinforced area


