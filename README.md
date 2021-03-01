## The parallelizable Lotka-Volterra Metacommunity assembly Model (pLVMCM)
##### Jacob O'Sullivan
j.osullivan@qmul.ac.uk | j.osullivan@zoho.com

# Contents

- [Download](x)
- [Compilation](x)
  - [Dependencies](x)
    - [Boost](x)
    - [Armadillo](x)
    - [Sundials](x)
  - [Build executable](x)
- [Model elements](x)
  - [Lotka-Volterra dynamics](x)
  - [Sampling of abiotic/biotic values](x)
    - [Spatial structure](x)
    - [Species ecological traits](x)
    - [Environmental modelling](x)
- [List of program arguments](x)
- [Simulation output](x)
- [Test implementation](x)
- [The emergence of autonomous turnover](x)


# Download

Clone this repository into your home directory using the command:

```
git clone https://github.com/jacobosullivan/LVMCM_src.git
```

The software does not need to be in the home directory to run, however data are saved to the home directory by default and the accompanying R scripts assume this to be the case.

# Compilation

## Dependencies
MacOS or Window users... install Linux

List of dependencies (most recent compatible version):
- Armadillo (most recent)
- Boost (most recent)
- Sundials 2.7.0
- CMake (most recent)

### Boost

For most recent release:

```
sudo apt-get install libboost-all-dev
```

It is recommended to install Boost before Armadillo

### Armadillo

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

### Sundials

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

## Build executable

Navigate to the subdirectory `LVMCM_src/LVMCM` and run the following:

```
rm -rf build; mkdir build; cd build
cmake ..
make VERBOSE=1
```

This will create and executable file called `LVMCM` which can be run to initialize or load a metacommunity simulation.

# Model elements

## Lotka-volterra dynamics

For detail on the dynamical system used to model species biomasses temporal evolution, see published material

## Sampling of abiotic/biotic values

### Spatial structure
- The cartesian coordinates are sampled either from a 2D uniform distribution or, if the argument `-R F` is passed to the executable, from a 2D lattice.
In the latter case a square number of nodes is required.
- For random graphs, edges are allocated using the [Gabriel](https://en.wikipedia.org/wiki/Gabriel_graph) algorithm or, if the argument `-G F` is passed to the executable, using a complete graph.

### Species ecological traits
- The distribution from which competitive coefficients is sample is controlled using the program argument `-D` which can take the following integer values:
  - 0: Discrete distribution
  - 1: Beta distribution
  - 2: Discretized beta distribution
  - 3: Normal distribution
  - 4: Discretized normal distribution
- For 0, 1, and 3, two shape parameters are passed using the argument `-c X Y`.
X is the probability of non-zero interactions (0), the first shape parameter (1), or the mean of the distribution (3).
Y is the value of non-zero interactions (0), the second shape parameter (1), or the variance of the distribution (3).
- For 2 and 4 and addition parameter is passed after the argument `-D`, e.g. `-D 4 0.3` which defines the probability of non-zero interactions for the discretized continuous distribution.
- Trophic links are sampled from a log normal distribution characterized by two parameters, the standard deviation of the normal distribution set using program argument `-F` and a scaling parameter set using `-a`.
- The emigration rate and dispersal length, fixed for all species, are set using the argument `-d X Y`. If the emigration rate is set to -1.0, each species will be allocated a unique emigration rate in the range 0 to 1.

### Environmental modelling
- The environment is either modelled implicitly 'through the eyes of the species' or explicitly. In the latter case, selected by passing the program argument `-e X`. For X>0, X explicit environmental distributions are generated and species are allocated environmental tolerance coefficients which define the impact of a given environmental variable on their growth rate.

# List of program arguments
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
- `-id`: job ID used for automatic check pointing ["NA"]
- `-is`: invasion size (is*S+1 invaders generated in each iteration of the assembly algorithm) [0.05]
- `-n`: number of nodes [4]
- `-o`: 3x output variables used in generating file names - a key parameter, experiment name, replicate number [0, DEFAULT, 0]
- `-p`: spatial autocorrelation of the environment [1]
- `-pp`: probability of invading a producer species, set to <1.0 for bipartite model [0.0]
- `-r`: consumer respiration rate parameter [0.0]
- `-s`: trophic link distribution parameter [0.0]
- `-sc`: path to data for scaling of local interaction matrix [{} (empty string)]
- `-si`: Schwartz iteration number and window size [2, 100]
- `-sk`: environmental sensitivity shape parameter controlling the fundamental niche width [vec {0.0}]
- `-st`: standard deviation of environmental fluctuations
- `-t`: number of time steps simulated between each invasion [500]
- `-v`: the standard deviation of the environment/growth rate distribution [0.01]
- `-x`: path to spatial network if not random generated (should be Armadillo matrix raw_ascii) [{} (empty string)]

Switches
- `-C`: set to `F` to switch off interspecific interactions between producer species
- `-D`: set to `F` to select continuous (beta) distribution in competitive overlap coefficients
- `-G`: if set to `F`, complete graph selected
- `-O`: output switch, if set to `F` will not write to file
- `-R`: if set to `F`, a spatial network modelled using a 2D lattice
- `-S`: if set to `T`, model with regularly write to file
- `-SC`: if set to `T`, select symmetric competition model
- `-SS`: if set to `T`, will generate source-sink matrix without assembling
- `-Z`: if set to 2: all random seeds are set to 1 for reproducibility; if set to 1: only random seeds used in generating enviroment set to 1

Perturbation experiments/analysis objects
- `-B`: compute the spatio-temporal beta-diversity **after every invasion**
- `-F`: fragmentation/conservation area experiment
- `-H`: harvesting experiment
- `-K`: interative node/edge removal experiment
- `-L`: long distance dispersal experiment
- `-T`: generate long time series trajectory
- `-W`: warming experiment

# Simulation output

Once `invMax` is reached, and in the case `-O F` is _not_ passed to the executable, a folder called SimulationData will be generated as a subdirectory by default in the folder `~/LVMCM_src/SimulationData`, however if the argument `-f /<path>/<to>/<output>` is passed, alternative storage locations can be requested.
The various model matrices will be stored in this folder with a file path which records the number of nodes, the experiment name, the date, the requested invasions, a key parameter and the replicate number.
For example, the assembly above will output the following matrices in the folder `/<path>/<to>/<output_directory>/SimulationData/N=32/testAssembly/<date>/`

- `<date>_testAssembly(1000)1bMat0.mat`: matrix of biomasses representing final state of the metacommunity
- `<date>_testAssembly(1000)1bMat_src0.mat`: discrete matrix with 1 corresponding to source populations, -1 to sink populations and 0 to biomasses below the detection threshold of 10^-4 biomass units.
- `<date>_testAssembly(1000)1bMat_c0.mat`: matrix of biomasses consumer species
- `<date>_testAssembly(1000)1bMat_c_src0.mat`: source-sink allocation consumer species
- `<date>_testAssembly(1000)1rMat0.mat`: matrix of local intrinsic growth rates
- `<date>_testAssembly(1000)1S0.mat`: column vector recording the regional species richness as a function of each invasion event
- `<date>_testAssembly(1000)1cMat0.mat`: matrix of competitive overlap coefficients
- `<date>_testAssembly(1000)1network0.mat`: cartesian coordinate of the spatial network
- `<date>_testAssembly(1000)1dMat_n0.mat`: dispersal matrix
- `<date>_testAssembly(1000)1envMat0.mat`: explicitly modelled environmental distributions
- `<date>_testAssembly(1000)1tMat0.mat`: species environmental tolerances
- `<date>_testAssembly(1000)1params0.mat`: a list of model parameters


# Test implementation

After compilation, navigate to the build folder and the following command:

```
./LVMCM -o 1 testComp 1 -n 16 -p 1 -i 100 -d 0.2 1.0 -c 0.3 0.3 -v 0.1 -t 1000 -Z 2 -O F"
```

This will assemble a competitive metacommunity of 16 nodes but will not write to file. In case of error messages, check all dependencies have been installed.

# Example assemblies: The emergence of autonomous turnover

In figure S4 of O'Sullivan et al. (2021) Intrinsic ecological dynamics drive biodiversity turnover in model metacommunities we show how autonomous turnover occurs at the onset of local structural instability for various random matrices.
Here we demonstrate how the LVMCM can be used to show this emergent phenomenon for three test cases, with A_ij sampled from
1. A discrete distribution with c1=c2=0.3
2. A discretized beta distribution with connectance 0.4 and mean/variance equal to 1. above
3. A discretized normal distribution with connectance 0.4 and mean/variance equal to 1. above

All spatial parameters are as in the paper.

Run the following commands from the directory `~/LVMCM_src/LVMCM/build`

```
PARFILE=$HOME/LVMCM_src/LVMCM/parFiles/autonomous_turnover_example/autonomous_turnover_example_pars.txt
./LVMCM $(sed -n "1p" $PARFILE) # assembly 1
./LVMCM $(sed -n "2p" $PARFILE) # assembly 2
./LVMCM $(sed -n "3p" $PARFILE) # assembly 3
```

These simulations will several hours to complete depending on the system. For the example models in the parameter file included with this repository random seeds have been fixed and corresponding output has been included in the folder `~/LVMCM_src/SimulationData`
To reproduce the result shown in Figure S4, run the assemblies or explore the example data using the R script `~/LVMCM_src/LVMCM/RCode/autonomous_turnover_example.R`
