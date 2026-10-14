

### Adapt these paths for your machine before testing the TEST-CASES!

hostnm="$( hostname )"

case ${hostnm} in
    #
    "shaolin")
        export DIR_NC_IN="/home/laurent/data/INPUT_NANUQ_DISTRIB"

        export ERA5_HUDSON12_DIR="xxx"

        export XIOS_EXE="/opt/xios2-trunk/bin/xios_server.exe"

        export XIOS_OA3_EXE="xxx"  ; # same as above but compiled with OASIS support!

        export FATM_ERA5_DIR="xxx"
        ;;
    #   
    "frazilo")
        export DIR_NC_IN="/data/laurent/INPUT_NANUQ_DISTRIB"

        export ERA5_HUDSON12_DIR="/SUMMER/SASIP/model-forcings/nanuq/ERA5-HUDSON12"

        export XIOS_EXE="/opt/xios-trunk/bin/xios_server.exe"

        export XIOS_OA3_EXE="/home/laurent/src/xios2-trunk_oa3/bin/xios_server.exe"  ; # same as above but compiled with OASIS support!

        export FATM_ERA5_DIR="/SUMMER/SASIP/model-forcings/atmo_forcing/ERA5_Arctic"
        ;;
    *)
        echo "ERROR: unknown host: ${hostnm} !"
        exit
        ;;
esac    
