## Playing with the HUDSON12 ocean/sea-ice coupled test-case


### The `nemo.exe`to use

You have compiled the `nemo.exe` executable of NEMO 4.2.2 against OASIS MCP for the OCE (3D liquid ocean) component only (no SI3)!
To do so, the content of the `cfgs/HUDSON12/cpp_HUDSON12.cpp` file is restricted to the following pre-compilation keys:

    bld::tool::fppkeys   key_xios key_qco key_oasis3

Similarly, the line you have to add in the `cfgs/ref_cfgs.txt`, relevant to the HUDSON12 config, must only read: 

    HUDSON12 OCE

Prior to compilation, copy the corrected version of `sbccpl.F90` found under `tests/fix_nemo/4.2.2/` in NANUQ to the `cfgs/HUDSON12/MY_SRC` of your NEMO distribution (it will override the current `src/OCE/SBC/sbccpl.F90` of NEMO at compilation time).

Use the same compiler and libraries (netCDF, OASIS and XIOS) as those used for compiling `nanuq.exe` for the coupled setup, technically it is recommended to use the same `arch_<HOST>_OA3.fcm` file as that used for NANUQ (see the main README).


### The `nanuq.exe` to use

The NANUQ executable has been compiled in generic coupled mode against OASIS MCP (see main README file), it should therefore be there:

```../../../cfgs/generic_cpl_oce/BLD/bin/nanuq.exe```


### Preparing the production directory

You have downloaded and extracted the `INPUT_NANUQ_DISTRIB.tgz` archive `<SOMEWHERE>`.

You have created a directory `ERA5-HUDSON12` under `<SOMEWHERE>` in which you have downloaded all the netCDF files found [here](https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/SASIP/model-forcings/nanuq/ERA5-HUDSON12/catalog.html) (16 GiB!): 
It is the 1-hourly ERA5 surface atmospheric forcing for 1997 interpolated on the HUDSON12 domain, 1 netCDF file per atmospheric variable.

Make a copy or create symbolic links, in here, of all the netCDF files found under the 3 following directories:
* `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/HUDSON12/`
* `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/HUDSON12/cpl_OCE/`
* `<SOMEWHERE>/ERA5-HUDSON12`

Adjust the host-specific value of paths such as `DIR_NC_IN` & `HUDSON12_FATM_DIR` into the file located under `<NANUQ_REPO>/tests/paths_nanuq_data.bash`; `INPUT_NANUQ_DISTRIB` should provide the full path to the `INPUT_NANUQ_DISTRIB` you have downloaded.

Copy or create a symbolic link of the `nemo.exe` executable (see above).

Copy or create a symbolic link of the `nanuq.exe` executable (see above).

The `prepare_prod_dir.sh` can do all this for you, example:

`./prepare_prod_dir.sh` prepares the current directory for the run...

Now you can launch the simulation on `N` cores.

### Using XIOS in library mode

Here on `N=30` cores:

```mpirun -n 16 ./nemo.exe : -n 14 ./nanuq.exe```


Follow the progression of the simulation:

```tail -f nanuq.output``` for NANUQ
```tail -f ocean.output``` for NEMO

or 

```tail -f nanuq.stat``` for NANUQ
```tail -f run.stat```   for NEMO


### Same, but using XIOS in "detached" mode (server)

Copy or create a symbolic link of the `xios_server.exe` in the current directory.

Edit the `iodef.xml` to have `true` instead of `false` for the `"using_server"` variable:

```<variable id="using_server" type="bool">true</variable>```

You can launch the simulation on `N` & `M` cores for NEMO & NANUQ and `P` cores for the XIOS server:

```mpirun -n <N> ./nemo.exe: -n <M> ./nanuq.exe: -n <P> ./xios_server.exe```


### Cleaning the current directory

Use the `clean.sh` with care!
