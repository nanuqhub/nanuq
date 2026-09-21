## Idealized test-cases:

(In increasing order of complexity...)

### `CHANNELD17`

The channel test-case of Dansereau et al. 2017.


### `CYCLONE

The "cyclone" test-case of Mehlmann et al. 2021.


## Realistic test-cases

### `HUDSON12/standalone`

Hudson Bay at 12th of a degree (extraction from ORCA12), standalone mode, landlocked configuration.

Ocean "SSX" surface forcing (6-hourly) and ERA5 atmo forcing (1-hourly) are provided directly on the HUDSON12 horizontal domain.

### `HUDSON12/cpl_oce`

Same but couple to OCE of NEMO (`nemo.exe`) via OASIS, landlocked configuration.


### `EGL12/standalone`

East-Greenland sea at 12th of a degree (extraction from ORCA12), standalone mode, _BDY_ lateral boundary condition setup.

Ocean "SSX" surface forcing (6-hourly) and ERA5 atmo forcing (1-hourly) are provided directly on the EGL12 horizontal domain.

### `EGL12/cpl_oce`

Same but couple to OCE of NEMO (`nemo.exe`) via OASIS, _BDY_ lateral boundary condition setup for both NANUQ and NEMO/OCE.


