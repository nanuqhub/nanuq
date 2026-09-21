## PLaying with the CYCLONE test-case

The NANUQ executable has been compiled in generic standalone mode, it should therefore be there:

```../../cfgs/generic/BLD/bin/nanuq.exe```

You have downloaded and extracted the `INPUT_NANUQ_DISTRIB.tgz` archive `<SOMEWHERE>`.

Decide what resolution `<RES>` (in km) you want to use (available: 10km, 4km and 2km).

Copy or create symbolic links of all the netCDF files found under `<SOMEWHERE>/INPUT_NANUQ_DISTRIB/CYCLONE/<RES>km` in here.

Adjust the host-specific value of paths such as `DIR_NC_IN` into the file located under `<NANUQ_REPO>/tests/paths_nanuq_data.bash`; `INPUT_NANUQ_DISTRIB` should provide the full path to the `INPUT_NANUQ_DISTRIB` you have downloaded.

Create symbolic links `namelist_dom_cfg` and `namelist_ice_cfg` pointing to the appropriate namelist for the chosen resolution, _i.e._ `namelist_dom_cfg.<RES>km` and `namelist_ice_cfg.<RES>km`.

Copy or create a symbolic link of the `nanuq.exe` executable compiled into `../../cfgs/generic/BLD/bin/nanuq.exe`

The `prepare_prod_dir.sh` can do all this for you, example:

`./prepare_prod_dir.sh 10` prepares the current directory for the resolution at 10 km...

Now you can launch the simulation on `N` cores.

### Using XIOS in library mode

Here on `N=4` cores:

```mpirun -n 4 ./nanuq.exe```


Follow the progression of the simulation:

```tail -f nanuq.output```

or 

```tail -f nanuq.stat```


### Same, but using XIOS in "detached" mode (server)

Copy or create a symbolic link of the `xios_server.exe` in the current directory.

Edit the `iodef.xml` to have `true` instead of `false` for the `"using_server"` variable:

```<variable id="using_server" type="bool">true</variable>```

You can launch the simulation on `N` cores for NANUQ and `M` cores for the XIOS server:

```mpirun -n N ./nanuq.exe: -n M ./xios_server.exe```


### Cleaning the current directory

Use the `clean.sh` with care!
