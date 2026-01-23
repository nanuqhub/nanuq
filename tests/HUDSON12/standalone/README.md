## PLaying with the HUDSON12 (standalone) test-case

The NANUQ executable has been compiled in generic standalone mode, it should therefore be there:

```../../cfgs/generic/BLD/bin/nanuq.exe```

You have downloaded and extracted the `INPUT_NANUQ_DISTRIB.tgz` archive `<SOMEWHERE>`.

You have created a directory `FATM_1997_HUDSON12` under `<SOMEWHERE>` in which you have downloaded all the netCDF files found [here](https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/MEOM/NANUQ/atmo_forcings/FATM_1997_HUDSON12/catalog.html) (16 GiB!): 
It is the 1-hourly ERA5 surface atmospheric forcing for 1997 interpolated on the HUDSON12 domain, 1 netCDF file per atmospheric variable.

Make a copy or create symbolic links, in here, of all the netCDF files found under the 3 following directories:
* `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/HUDSON12/`
* `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/HUDSON12/standalone`
* `<SOMEWHERE>/FATM_1997_HUDSON12`

Copy or create a symbolic link of the `nanuq.exe` executable compiled into `../../cfgs/generic/BLD/bin/nanuq.exe`

The `prepare_prod_dir.sh` can do all this for you, all you have to do is set the `DIR_NC_IN` variable to the full path to `INPUT_NANUQ_DISTRIB` directory, and the `FATM_DIR` variable to the full path to `FATM_1997_HUDSON12`.

Now you can launch the simulation on `N` cores.

### Using XIOS in library mode

Here on `N=30` cores:

```mpirun -n 30 ./nanuq.exe```


Follow the progression of the simulation:

```tail -f nanuq.output```

or 

```tail -f nanuq.stat```


### Same, but using XIOS in "detached" mode (server)

Copy or create a symbolic link of the `xios_server.exe` in the current directory.

Edit the `iodef.xml` to have `true` instead of `false` for the `"using_server"` variable:

```<variable id="using_server" type="bool">true</variable>```

You can launch the simulation on `N` cores for NANUQ and `M` cores for the XIOS server:

```mpirun -n <N> ./nanuq.exe: -n <M> ./xios_server.exe```


### Cleaning the current directory

Use the `clean.sh` with care!
