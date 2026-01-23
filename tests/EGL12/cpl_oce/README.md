## Playing with the EGL12 (East-Greenland) ocean/sea-ice coupled test-case 

As opposed to HUDSON12, this configuration uses prescribed lateral boundary conditions (BDYs) for both NANUQ (sea-ice) and OCE (NEMO 3D ocean component) !

### The `nemo.exe`to use

You have compiled the `nemo.exe` executable of NEMO 4.2.2 against OASIS MCP for the OCE (3D liquid ocean) component only (no SI3)!
To do so, the content of the `cfgs/EGL12/cpp_EGL12.cpp` file is restricted to the following pre-compilation keys:

    bld::tool::fppkeys   key_xios key_qco key_oasis3

Similarly, the line you have to add in the `cfgs/ref_cfgs.txt`, relevant to the EGL12 config, must only read: 

    EGL12 OCE

Prior to compilation, copy the corrected version of `sbccpl.F90` found under `tests/fix_nemo/4.2.2/` in NANUQ to the `cfgs/EGL12/MY_SRC` of your NEMO distribution (it will override the current `src/OCE/SBC/sbccpl.F90` of NEMO at compilation time).

Use the same compiler and libraries (netCDF, OASIS and XIOS) as those used for compiling `nanuq.exe` for the coupled setup, technically it is recommended to use the same `arch_<HOST>_OA3.fcm` file as that used for NANUQ (see the main README).


### The `nanuq.exe` to use

The NANUQ executable has been compiled in generic coupled mode against OASIS MCP (see main README file), it should therefore be there:

```../../../cfgs/generic_cpl_oce/BLD/bin/nanuq.exe```


### Preparing the production directory

You have downloaded and extracted the `INPUT_NANUQ_DISTRIB.tgz` archive `<SOMEWHERE>`.

You have created a directory `ERA5_Arctic` under `<SOMEWHERE>` in which you have downloaded all the netCDF files, for the years of interest, found [here](https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/SASIP/model-forcings/atmo_forcing/ERA5_Arctic/catalog.html): 
It is an Arctic extraction of the 1-hourly ERA5 surface atmospheric state, 1 netCDF file per atmospheric variable.

Make a copy or create symbolic links, in here, of all the netCDF files found under the 3 following directories:
* `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/EGL12/`
* `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/EGL12/cpl_OCE/`
* `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/EGL12/BDY/`
* `<SOMEWHERE>/ERA5_Arctic`


Copy or create a symbolic link of the `nemo.exe` executable (see above).

Copy or create a symbolic link of the `nanuq.exe` executable (see above).


The `prepare_prod_dir.sh` can do all this for you, all you have to do is set the full path to the `nanuq.exe` and `nemo.exe`, set the `DIR_NC_IN` variable to the full path to `INPUT_NANUQ_DISTRIB` directory, and the `FATM_DIR` variable to the full path to `ERA5_Arctic`.

Now you can launch the simulation on `N` cores.


### Launching coupled ocean--sea-ice simulation with XIOS in server mode

Here on `N=32` cores:

```mpirun -n 2 ./xios_server.exe : -n 20 ./nemo.exe : -n 10 ./nanuq.exe``


Follow the progression of the simulation:

```tail -f nanuq.output``` for NANUQ
```tail -f ocean.output``` for NEMO

or 

```tail -f nanuq.stat``` for NANUQ
```tail -f run.stat```   for NEMO


### Cleaning the current directory

Use the `clean.sh` with care!
