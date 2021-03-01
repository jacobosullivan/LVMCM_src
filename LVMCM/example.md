Example implementations of the LVMCM

Assembly:
```
mpirun -np 2 /home/jack/gitClones/LVMCM_src/LVMCM/buildLocal/LVMCM -o 2 Temp 1 -N 32 -Phi 1 -invMax 1000 -disp 0.01 0.5 -c_ij 0.5 0.5 -var_e 0.01 -g 1.0 -z F -O F
```

Warming:

```
mpirun -np 1 ./LVMCM -new F -bMat "/home/jack/Dropbox/SimulationData/N=128/tempGrad_experiment/2020-2-18/2020-2-18(10000)0bMat0.mat"
```

