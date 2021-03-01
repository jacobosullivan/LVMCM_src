# Compile and implementation instructions Intel 2020

Compilation - select (rename) intel compatible CMakeLists.txt

```
# load modules including cmake
module load use.own
module load intel intelmpi dependencies cmake

rm -rf build_intel; mkdir build_intel; cd build_intel
cmake .. \
-DARMADILLO_LIBRARY=$HOME/install/armadillo/lib64/libarmadillo.so \
-DSUNDIALS_NVECSERIAL_LIB:FILEPATH=$HOME/install/sundials2.7/lib/libsundials_nvecserial.so \
-DSUNDIALS_KINSOL_LIB:FILEPATH=$HOME/install/sundials2.7/lib/libsundials_kinsol.so \
-DSUNDIALS_CVODE_LIB:FILEPATH=$HOME/install/sundials2.7/lib/libsundials_cvode.so \
-DCMAKE_BUILD_TYPE=Release

make VERBOSE=1
```

implementation using matrix multiplication performance test main.cpp

```
# load modules
module load use.own
module load intel intelmpi dependencies

## non parallel version
./LVMCM -Test 0 10 -tMax 100 -O F

## parallel verson
export I_MPI_HYDRA_TOPOLIB=ipl # This doesn't always seem necessary!
mpirun -n 2 ./LVMCM -Test 0 10 -tMax 100 -MPI T -Bisec 1 -Schwartz 2 100 -O F
```
