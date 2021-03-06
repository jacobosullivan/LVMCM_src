## The Lotka-Volterra Metacommunity assembly Model (LVMCM)
##### Jacob O'Sullivan
j.osullivan@qmul.ac.uk | j.osullivan@zoho.com

A highly flexible framework for simulating spatially structured ecological communities in heterogeneous landscapes written in C++

<p align="center">
<img width="450" height="310" src="https://github.com/jacobosullivan/LVMCM_src/blob/master/icon.png?raw=true">
</p>

# Contents

- [Download](https://github.com/jacobosullivan/LVMCM_src#download)
- [Compilation](https://github.com/jacobosullivan/LVMCM_src#compilation)
  - [Dependencies](https://github.com/jacobosullivan/LVMCM_src#dependencies)
    - [Boost](https://github.com/jacobosullivan/LVMCM_src#boost)
    - [Armadillo](https://github.com/jacobosullivan/LVMCM_src#armadillo)
    - [Sundials](https://github.com/jacobosullivan/LVMCM_src#sundials)
  - [Build executable](https://github.com/jacobosullivan/LVMCM_src#build-executable)
- [Model elements](https://github.com/jacobosullivan/LVMCM_src#model-elements)
  - [Lotka-Volterra dynamics](https://github.com/jacobosullivan/LVMCM_src#lotka-volterra-dynamics)
  - [Sampling of abiotic/biotic values](https://github.com/jacobosullivan/LVMCM_src#sampling-of-abioticbiotic-values)
    - [Spatial structure](https://github.com/jacobosullivan/LVMCM_src#spatial-structure)
    - [Species ecological traits](https://github.com/jacobosullivan/LVMCM_src#species-ecological-traits)
    - [Environmental modelling](https://github.com/jacobosullivan/LVMCM_src#environmental-modelling)
- [List of program arguments](https://github.com/jacobosullivan/LVMCM_src#list-of-program-arguments)
- [Simulation output](https://github.com/jacobosullivan/LVMCM_src#simulation-output)
- [Test implementation](https://github.com/jacobosullivan/LVMCM_src#test-implementation)
- [Multi-threading](https://github.com/jacobosullivan/LVMCM_src#multi-threading)
- [The emergence of autonomous turnover](https://github.com/jacobosullivan/LVMCM_src#example-assemblies-the-emergence-of-autonomous-turnover)


# Download

Clone this repository into your home directory using the command:

```
git clone https://github.com/jacobosullivan/LVMCM_src.git
```

The software does not need to be in the home directory to run. The output directory can be selected at run-time, but the default location is `~/LVMCM_src/SimulationData` and the accompanying R script, which processes data from the example assemblies, assumes this to be the case.

# Compilation

## Dependencies
Note, this software was developed on Ubuntu 18.04-20.04. While compatibility with MacOS is likely, it is not guaranteed. At present we do not have a Windows distribution.

List of dependencies (most recent compatible version):
- [Boost](https://www.boost.org/) (most recent 1.75)
- [Armadillo](http://arma.sourceforge.net/) (most recent, 10.2)
- [CMake](https://cmake.org/) (most recent 3.2)
- [Sundials](https://computing.llnl.gov/projects/sundials/) **2.7.0** (work to update ODE implementation is on-going)

### Boost

For most recent release:

```
sudo apt-get install libboost-all-dev
```

It is recommended to install Boost before Armadillo

### Armadillo

If not present already, install LAPACK, BLAS and cmake.
ARPACK is also required for sparse matrix support
```
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
sudo apt-get install cmake
sudo apt-get install libarpack++2-dev
```

Download tar.xz file [here](http://arma.sourceforge.net/download.html)
Navigate to download location and extract, generate make files and build:

```
tar -xf armadillo-10.2.0.tar.xz
cd armadillo-10.2.0
cmake -DALLOW_FLEXIBLAS_LINUX=ON . # allow linking to BLAS libraries LINUX
```

See output - check that dependencies have been located

```
sudo make install
```

### Sundials

Download 2.7 distribution [here](https://computing.llnl.gov/projects/sundials/sundials-software).
Navigate to download location:

```
tar -xf sundials-2.7.0.tar.gz
cd sundials-2.7.0
mkdir builddir
cd builddir
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/ ..
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

For details on the dynamical system used to model species biomass dynamics, see published material:
- O'Sullivan et al. (2019) [Metacommunity‐scale biodiversity regulation and the self‐organised emergence of macroecological patterns](https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.13294)
- O'Sullivan et al. (2021) [Intrinsic ecological dynamics drive biodiversity turnover in model metacommunities](https://www.biorxiv.org/content/10.1101/2020.05.22.110262v1)

## Sampling of abiotic/biotic values

### Spatial structure
- The Cartesian coordinates are sampled either from a 2D uniform distribution or, if the argument `-R F` is passed to the executable, from a 2D lattice.
In the latter case a square number of nodes is required.
- For random graphs, edges are allocated using the [Gabriel](https://en.wikipedia.org/wiki/Gabriel_graph) algorithm or, if the argument `-G F` is passed to the executable, using a complete graph.

### Species ecological traits
- The distribution from which competitive coefficients are sampled is controlled using the program argument `-D` which can take the following integer values:
  - 0: Discrete distribution
  - 1: Beta distribution
  - 2: Discretized beta distribution
  - 3: Normal distribution
  - 4: Discretized normal distribution
- For 0, 1, and 3, two shape parameters are passed using the argument `-c X Y`.
X is the probability of non-zero interactions (0), the first shape parameter (1), or the mean of the distribution (3).
Y is the value of non-zero interactions (0), the second shape parameter (1), or the variance of the distribution (3).
- For 2 and 4 and additional parameter is passed after the argument `-D`, e.g. `-D 4 0.3` which defines the probability of non-zero interactions for the discretized continuous distribution.
- Trophic links are sampled from a log normal distribution characterized by two parameters, the standard deviation of the normal distribution used to generate log-normal trophic interaction coefficients is set using program argument `-s` and a scaling parameter set using `-a`.
- The emigration rate (X) and dispersal length (Y), fixed for all species, are set using the argument `-d X Y`. If the emigration rate is set to -1.0, each species will be allocated a unique emigration rate in the range 0 to 1.

### Environmental modelling
- The environment is either modelled implicitly 'through the eyes of the species' or explicitly. In the latter case, selected by passing the program argument `-e X`. For X>0, X explicit environmental distributions are generated and species are allocated environmental tolerance coefficients which define the impact of a given environmental variable on their growth rate.

# List of program arguments
In the parameters below defaults are listed in square brackets

Input parameters:
- `-a`: trophic link scaling parameter [0.0]
- `-b`: the full path of the biomass matrix for a model to be imported [{} (empty string)]
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

Once `invMax` is reached, and in the case `-O F` is not passed to the executable, a folder called SimulationData will be generated as a subdirectory by default in the folder `~/LVMCM_src/`, however if the argument `-f /<path>/<to>/<output>` is passed, alternative storage locations can be requested.
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

Note, not all of the above are generated in each assembly. For example, competitive metacommunities do not include consumer biomass matrices. Only matrices with non-zero dimensions during simulation are output.

# Test implementation

After compilation, navigate to the build folder and the following command:

```
./LVMCM -o 1 testComp 1 -n 16 -p 1 -i 100 -d 0.2 1.0 -c 0.3 0.3 -v 0.1 -t 1000 -Z 2 -O F
```

This will assemble a competitive metacommunity of 16 nodes but will not write to file. In case of error messages, check all dependencies have been properly installed.

# Multi-threading

The Armadillo linear algebra library allows multi-threaded, shared memory parallelisation of matrix operations via [LAPACK, BLAS and high speed replacements (e.g. OpenBLAS, MKL)](http://arma.sourceforge.net/faq.html#dependencies). These should automatically be located by the CMake installer when building source code, however if stored in a non-standard location this may need to be explicitly defined when compiling. Note that the all creation threads and allocation of work is done automatically by the library and you may not find multi-threading is activated for small system (matrix) sizes.

To set the number of OpenBLAS threads run the following command prior to execution of the model sofware:

```
export OPENBLAS_NUM_THREADS=X
```

with `X` the desired (maximum) number of threads.

# Example assemblies: The emergence of autonomous turnover

In figure S4 of O'Sullivan et al. (2021) "Intrinsic ecological dynamics drive biodiversity turnover in model metacommunities" we show how autonomous turnover occurs at the onset of local structural instability for various random matrices.
Here we demonstrate how the LVMCM can be used to show this emergent phenomenon for three test cases, with A_ij sampled from
1. A discrete distribution with c1=c2=0.3
2. A discretized beta distribution with connectance 0.4 and mean/variance equal to 1. above
3. A discretized normal distribution with connectance 0.4 and mean/variance equal to 1. above

All spatial parameters are as in the paper.

The program arguments required to assemble these three metacommunity models are included in the file `/LVMCM_src/LVMCM/parFiles/autonomous_turnover_example/autonomous_turnover_example_pars.txt`.

Run the following commands from the directory `~/LVMCM_src/LVMCM/build`

```
PARFILE=$HOME/LVMCM_src/LVMCM/parFiles/autonomous_turnover_example/autonomous_turnover_example_pars.txt
./LVMCM $(sed -n "1p" $PARFILE) # assembly 1
./LVMCM $(sed -n "2p" $PARFILE) # assembly 2
./LVMCM $(sed -n "3p" $PARFILE) # assembly 3
```

Each simulation will take several hours to complete depending on the system and the three models will generate around 2.6GB of simulation data.

Note the example assemblies explore the rate of [autonomous turnover](https://www.biorxiv.org/content/10.1101/2020.05.22.110262v1) as a function of local diversity saturation. In the example assemblies, we show how the emergence of autonomous turnover coincides with the onset of local diversity limits via ecological structural instability. This requires generating and storing time series (at low temporal resolution for efficiency) after each successive iteration of the invasion algorithm (experiment selected via the argument `-B` included in the parameter file). Such experiments require large amounts of storage. Single assemblies recorded at the end state only and without time series typically require ~10MB of uncompressed storage, though this is of course sensitive to system size.

To run these simulations in parallel on a HPC system compatible with the Oracle (Sun) Grid Engine system submit the file `~/LVMCM_src/LVMCM/shFiles/autonomous_turnover_example.sh` to the scheduler with the `qsub` command or equivalent.

The data generated by these simulations has been stored as a ~1GB tar archive on [figshare.com](https://figshare.com/articles/dataset/Intrinsic_ecological_dynamics_drive_biodiversity_turnover_in_model_metacommunities_Supporting_data/14139644). To download run the following command:

```
wget -O ~/LVMCM_src/autonomous_turnover_example.tar.gz --no-check-certificate https://ndownloader.figshare.com/files/26663276
cd ~/LVMCM_src/
tar -xf autonomous_turnover_example.tar.gz
```

Once the assemblies are complete *or* the example data has downloaded and extracted, explore the outcome using the R script `~/LVMCM_src/LVMCM/RCode/autonomous_turnover_example.R`
