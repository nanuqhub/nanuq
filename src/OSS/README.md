# OSS: Ocean Surface State

Stuff that deals with the liquid ocean surface state.

The ocean surface state is needed for two reasons:
 * A/ as a bottom "boundary condition" for the sea-ice model
 * B/ to compute air-sea fluxes over open (liquid) ocean
 

### Sources

 * `ossprs.F90` (sucessor of `NEMO/sbcssm.F90`) obtains a prescribed surface state of the ocean read into netCDF files
 
 * `osscpl.F90` (a greatly simplified version of `NEMO/sbccpl.F90`) is the OASIS interface that obtains the surface state of the ocean from the ocean model that NANUQ is coupled to, and, in return, sends the surface fluxes required by the 3D ocean model as surface boundary conditions (regardless of the presence of sea-ice or not above)
 
 
