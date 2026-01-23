## Playing with the EGL12 (East-Greenland) standalone test-case

As opposed to HUDSON12, this configuration uses prescribed lateral boundary conditions (BDYs) for both NANUQ (sea-ice) and OCE (NEMO 3D ocean component) !

### The `nanuq.exe` to use

The NANUQ executable has been compiled in generic standalone mode, it should therefore be there:

```../../cfgs/generic/BLD/bin/nanuq.exe```


### Preparing the production directory

You have downloaded and extracted the `INPUT_NANUQ_DISTRIB.tgz` archive `<SOMEWHERE>`.

You have created a directory `ERA5_Arctic` under `<SOMEWHERE>` in which you have downloaded all the netCDF files, for the years of interest, found [here](https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/SASIP/model-forcings/atmo_forcing/ERA5_Arctic/catalog.html): 
It is an Arctic extraction of the 1-hourly ERA5 surface atmospheric state, 1 netCDF file per atmospheric variable.

Make a copy or create symbolic links, in here, of all the netCDF files found under the 3 following directories:
* `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/EGL12/`
* `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/EGL12/standalone/`
* `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/EGL12/BDY/`
* `<SOMEWHERE>/ERA5_Arctic`


Copy or create a symbolic link of the `nanuq.exe` executable compiled into `../../cfgs/generic/BLD/bin/nanuq.exe`

The `prepare_prod_dir.sh` can do all this for you, all you have to do is set the `DIR_NC_IN` variable to the full path to `INPUT_NANUQ_DISTRIB` directory, and the `FATM_DIR` variable to the full path to `ERA5_Arctic`.

Now you can launch the simulation on `N` cores.


### Launching a standalone sea-ice simulation with XIOS in server mode

Here on `N=32` cores:

```mpirun -n 2 ./xios_server.exe : -n 30 ./nanuq.exe```

Follow the progression of the simulation:

```tail -f nanuq.output```

or 

```tail -f nanuq.stat```


### Cleaning the current directory

Use the `clean.sh` with care!
