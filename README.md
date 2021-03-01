## The parallelizable Lotka-Volterra Metacommunity assembly Model (pLVMCM)
##### Jacob O'Sullivan
j.osullivan@qmul.ac.uk | j.osullivan@zoho.com

## Compilaton instructions:

### Install dependencies/build container image

List of dependencies (most recent compatible version):
- Armaillo (most recent)
- Boost (most recent)
- Sundials 2.7.0
- CMake 3.5.1
- MPI (chosen distribution)

To run software in a container environment, build singularity container from defintion file `/LVMCM_src/LVMCM/LVMCM_mpi.def` (debootstrap may be required).
For details see Singularity docs (https://singularity.lbl.gov/).

#### Libarary install/build

##### Boost

For most recent release:

```
sudo apt-get install libboost-all-dev
```

It is recommended to install Boost before Armadillo

##### Armadillo

If not present already, install LAPACK, BLAS and cmake
ARPACK is required for sparse matrix support
```
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
sudo apt-get install cmake
sudo apt-get install libarpack++2-dev
```

Download tar.xz file here: http://arma.sourceforge.net/download.html
Navigate to download location and extract, generate make files and build

```
tar -xf armadillo-10.2.0.tar.xz
cmake -DALLOW_FLEXIBLAS_LINUX=ON . # allow linking to BLAS libraries LINUX
```

See output - check that dependencies have been located

```
sudo make install
```

#### SUNDIALS 2.7

Download 2.7 distribution here: https://computing.llnl.gov/projects/sundials/sundials-software.
Navigate to download location:

```
tar -xf sundials-2.7.0.tar.gz
cd sundials-2.7.0
mkdir instdir
mkdir builddir
cd builddir
cmake .. # installdir defaults to /use/local unless CMAKE_INSTALL_PREFIX is set
make
sudo make install
```

### Buid executable - GCC compiler

Navigate to the subdirectory `LVMCM_src/LVMCM` and run the following:

```
rm -rf build; mkdir build; cd build
cmake ..
make VERBOSE=1
```

### Buid executable - Intel compiler

Intel compilers and dependencies are currently installed natively on the Apocrita cluster (btx718 user space only).
If attempting to compile in alternative user space, I recommend using GCC as described above unless all dependencies are available.
Note that using Intel 2020 requires a different CMakeLists.txt file.
This is currently saved in the LVMCM repo under the filename CMakeLists_int.txt.
To compile using Intel, first run the following `mv CMakeLists.txt CMakeLists_gcc.txt; mv CMakeLists)int.txt CMakeLists.txt`.
Intel compilers produce an improvement in run time performance therefore if available these should be used.

From within btx718 userspace, the following commands load the necessary private modules:

```
# load modules including cmake
module load use.own
module load intel intelmpi dependencies cmake
```

If running from an alternative user space, note the location of the Armadillo and Sundials install files and ensure cmake is pointing to the relevant directories. In the case that depenencies are installed in the folder `$HOME/install` (as in btx718 userspace), the following will build the exectuable:

```
rm -rf build\_intel; mkdir build_intel; cd build_intel
cmake .. \
-DARMADILLO\_LIBRARY=$HOME/install/armadillo/lib64/libarmadillo.so \
-DSUNDIALS\_NVECSERIAL_LIB:FILEPATH=$HOME/install/sundials2.7/lib/libsundials_nvecserial.so \
-DSUNDIALS\_KINSOL_LIB:FILEPATH=$HOME/install/sundials2.7/lib/libsundials_kinsol.so \
-DSUNDIALS\_CVODE_LIB:FILEPATH=$HOME/install/sundials2.7/lib/libsundials_cvode.so \
-DCMAKE\_BUILD_TYPE=Release
make VERBOSE=1
```

This will create and exectuable file called `LVMCM` which can be run to initialize or load a metacommunity simulation.

#### Lotka-volterra dynamics

##### Competitive metacommunity

$\frac{dB_{ix}}{dt}=B_{ix}\biggl(R_{ix}-\sum_{j=1}^{S}{C_{ij}B_{jx}}\biggr)-eB_{ix}(t)+\sum_{y\in \mathcal{N}(x)}{d(x,y)B_{iy}}$

with $d(x,y)=-e\ \text{for}\ x=y$ and $\frac{e}{k_y}\text{exp}(-\delta_{xy}\ell^{-1})$ otherwise ($k$ is the degree of the node $y$ and $\delta_{xy}$ the Euclidean distance between $x$ and $y$).


##### Biparite metacommunity

$\frac{dB_{ix}^{p}}{dt}=B_{ix}^{p}\biggl(R_{ix}-\sum_{j=1}^{S^{p}}{C_{ij}B_{jx}^{p}}-\sum_{l=1}^{S^{c}}{A_{il}B_{lx}^{c}}\biggr)-eB_{ix}^{p}(t)+\sum_{y\in \mathcal{N}(x)}{d(x,y)B_{iy}^{p}}$

$\frac{dB_{kx}^{c}}{dt}=\epsilon B_{kx}^{c}\biggl(\sum_{j=1}^{S^{p}}A_{kj}B_{jk}^{c}-\rho\biggr)-e B_{kx}^{c}+\sum_{y\in \mathcal{N}(x)}{d(x,y)B_{ky}^{c}}$

where $\rho$ is a respiration term, $\epsilon$ the assimilation efficiency and $d(x,y)$ defined as above.

#### Sampling of abiotic/biotic values

##### Spatial structure:
- The cartesian coordinates are sampled either from a 2D uniform distribution in the range $[0,\sqrt{N}]$ or, if the argument `-Lattice T` is passed to the executable, from a 2D lattice of area $N$.
In the latter case a square number of nodes is required.
- For random graphs, edges are allocated using the [Gabriel](https://en.wikipedia.org/wiki/Gabriel_graph) algorithm or, if the argument `-Gab F` is passed to the executable, using a complete graph.

##### Species ecological traits:
- Competition coefficients are sampled from a discrete distribution with program arguments `-c_ij X Y` representing the intesity and probability of non-zero competitive interaction coefficients, or, if the program argument `-Discr F` is passed to the exectuable, from a [beta distribution](https://en.wikipedia.org/wiki/Beta_distribution) (characterized by two shape parameters).
- Trohpic links are sampled from a log normal distribution characterized by two parameters, the standard deviation of the _normal_ distribution set using program argument `-F` and a scaling parameter set using `-alpha`.
- The emigration rate and dispersal length, fixed for all species, are set using the argument `-dispL X Y`. If the emigration rate is set to -1.0, each species will be allocated a unique emigration rate in the range $[0,1]$.

##### Environmental modelling
- The environment is either modelled implicitly 'through the eyes of the species' or explicitly. In the latter case, selected by passing the program argument `-envVar X`. For X>0, X explicit environmental distributions are generated and species are allocated environmental tolerance coefficents which define the impact of a given environmental variable on their growth rate.

#### List of program arguments

In the parameters below defaults are listed in square brackets

Input parameters:
- `-a`: trophic link scaling parameter [0.0]
- `-b`: the full bath of the biomass matrix for a model to be imported [{} (empty string)]
- `-c`: 2x competition parameter [0.5, 0.5]
- `-d`: emigration rate and dispersal length parameters [0.1, 0.1]
- `-dd`: path to storage for iterative domain decomposed solution (testing purposes) [{} (empty string)]
- `-dn`: select dispersal matrix normalisation - 0: effort normalised; 1: degree normalised; 2: passive dispersal [0]
- `-e`: the number of explicitly modelled environmental variables [0]
- `-f`: output folder [{} (empty string)]
- `-g`: cap on regional diversity (use with caution, if greater that regional limits, assembly with not complete) [0]
- `-i`: maximum number of invasions [0]
- `-id`: job ID used for automatic checkpointing ["NA"]
- `-is`: invasion size (is*S+1 invaders generated in each iteration of the assembly algorithm) [0.05]
- `-n`: number of nodes [4]
- `-o`: 3x output variables used in generating file names - a key parameter, experiment name, replicate number [0, DEFAULT, 0]
- `-p`: spatial autocorrelation of the environment [1]
- `-pp`: probability of invading a producer species, set to <1.0 for bipartite model [0.0]
- `-r`: consumer respiration rate parameter [0.0]
- `-s`: trophic link distribution parameter [0.0]
- `-si`: Schwartz iteration number and window size [2, 100]
- `-st`: standard deviation of environmental fluctuations
- `-t`: number of timesteps simulated between each invasion [500]
- `-v`: the standard deviation of the environment/growth rate distribution [0.01]
- `-x`: path to spatial newtork if not random generated (should be Armadillo matrix raw_ascii) [{} (empty string)]

Switches
- `-C`: set to `F` to switch off interspecific interactions between producer species
- `-D`: set to `F` to select continuous (beta) distribution in competitive overlap coefficients
- `-G`: if set to `F`, complete graph selected
- `-O`: output switch, if set to `F` will not write to file
- `-R`: if set to `F`, a spatial network modelled using a 2D lattice
- `-S`: if set to `T`, model with regularly write to file
- `-SC`: if set to `T`, select symmetric competition model
- `-Z`: if set to 2: all random seeds are set to 1 for reproducibility; if set to 1: only random seeds used in generating enviroment set to 1

Perturbation experiments/analysis objects
- `-F`: fragmentation/conservation area experiment
- `-H`: harvesting experiment
- `-K`: interative node/edge removal experiment
- `-L`: long distance distpersal experiment
- `-T`: generate long timeseries trajectory
- `-W`: warming experiment

#### Test implementation

After compilation, navigate to the build folder and the following command:

```
./LVMCM -o 1 testComp 1 -n 16 -p 1 -i 100 -d 0.2 1.0 -c 0.3 0.3 -v 0.1 -t 1000 -Z 2 -O F"
```

This will assemble a competitive metacommunity of 16 nodes but will not write to file.

#### Simulation output

Once `invMax` is reached, and in the case `-O F` is _not_ passed to the executable, a folder called SimulationData will be generated as a subdirectory in `/<path>/<to>/<output_directory>/`.
The various model matrices will be stored in this folder with a file path which records the number of nodes, the experiment name, the date, the requested invasions, a key parameter and the replicate number.
For example, the assembly above will output the following matrices in the folder `/<path>/<to>/<output_directory>/SimulationData/N=32/testAssembly/<date>/`

- `<date>_testAssembly(1000)1bMat0.mat`: matrix of biomasses representing final state of the metacommunity
- `<date>_testAssembly(1000)1bMat_src0.mat`: disrete matrix with 1 corresponding to source populations, -1 to sink populations and 0 to biomasses below the detection threshold of 10^-4 biomass units.
- `<date>_testAssembly(1000)1bMat_c0.mat`: matrix of biomasses consumer species
- `<date>_testAssembly(1000)1bMat_c_src0.mat`: source-sink allocation consumer species
- `<date>_testAssembly(1000)1rMat0.mat`: matrix of local intrinsic growth rates
- `<date>_testAssembly(1000)1S0.mat`: column vector recording the regional species richness as a function of each invasion event
- `<date>_testAssembly(1000)1cMat0.mat`: matrix of competitive overlap coefficients
- `<date>_testAssembly(1000)1network0.mat`: cartesian coordinate of the spatial network
- `<date>_testAssembly(1000)1dMat_n0.mat`: dispersal matrix generated as in EcoLetts paper
- `<date>_testAssembly(1000)1environ0.mat`: explicitly modelled environmental distributions
- `<date>_testAssembly(1000)1tMat0.mat`: species environmental tolerances
- `<date>_testAssembly(1000)1params0.mat`: a list of model parameters

#### Example assemblies

With the repo LVMCM_src cloned into the home directory, enter an R session and run through the script LVMCM_src/LVMCM/RCode/LVMCM_example.R. Note the package dependencies in the scripts LVMCM_src/LVMCM/RCode/plotting_functions.R and LVMCM_src/LVMCM/RCode/read_results.R will need to be installed.
