
# Compilation notes

 - Chombo must be compiled and $CHOMBO_HOME must be set to (chombo_root)/lib
 - For compilation for the Baden-Wurttemburg cluster, set ```cxxoptflags = -march=native```
 - For optimizations, add also ```-g3``` to ```cxxoptflags```

## Optional Flags
To Use AH Finder, the ```XTRACPPFLAGS``` must include

```
                    -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -DUSE_AHFINDER
``` 
and ```syslibflags``` must also include 
```
                    -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc
```



# Parameter Files 

 - There are a number of parameter files in *params/* folder that use different flags
 - The one we've currently been using as our 'production' file is *params_Dynamic_HighRes_AH.txt*
        - This uses 8 refinement levels with a computational box size of 128 and activates the AH locator routine




# Simulation Optimization

- For a computational box size of L=128, we need 8 refinement levels and a global grid scaling of *grid_scaling=2.0*
- A simulation size of L=224 N=128, we need 8 refinement levels and a global grid scaling of *grid_scaling=2.0* 
        
    - This ensures the Horizon is sufficiently resolved, in accordance with Katy's recommendations for sufficient resolution

 
 - We can determine the number of cells in the simulation by doing a quick test run using the exact parameter file you want. 
            Then use VISIT to execute a query on the NumZones. This returns the number of cells. Divide this number by the min box size
            to determine the number of boxes. Then determine the job parameters, maintaining the relation $$(\text{num boxes}) = (\text{num nodes}) * (\text{num cpus per node}) / (\text{num cpus per process})$$

    - Example: With the results above, we have 1280 boxes. With 16 tasks per node and 8 cpus per task, we need 40 nodes to have approximately
                                1 box per task






# General Notes on Matter Field and Initial Conditions

index of spatial Proca vector is DOWN
index of electric field vector is UP
shift vector is index UP

advection terms are the gradient projected along the shift: $\beta^j \partial_j x_i$


Current initial condition of Proca field:
- $X_x = A*e^{-r/r0} * \frac{1}{\sqrt{\gamma}} = A*\chi^{3/2}*e^{-r/r0}$
- $X_y = 0$
- $X_z = 0$
- $\phi = 0$
- $E^i = 0$

Kerr black metric is in form of quasi-isotropic coordinates

The first few seconds of the simulation consistute the *relaxation time* for the conformal factor. We use the ```ChiRelaxation``` method to minimize violations of the Hamiltonian constraint. This is accomplished by the modification of the ```evaluateRHS``` method of the level class. Before a simulation time of $t_{\chi}$, we evolve the formula $$\partial_t \chi = \chi*H*s$$, where $s$ is the relaxation speed. Note: This differs from the implementation of ```ChiRelaxation``` is base GRChombo. In the original source code, they divided by $\chi$, instead of multiplying it. This fails for Kerr in Quasi-Isotropic coordinates since the conformal factor goes to zero near the singularity, hence the flow equation becomes unusable. 

