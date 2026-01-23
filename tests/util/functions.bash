#!/bin/bash

function check_on_file()
{
    if [ ! -f ${1} ]; then
        echo "PROBLEM: file ${1} does not exist!!!"
        exit
    fi
}

function n_cores_xio_nodes_tot()
{
    # Input:
    #  #1: NBCPN:      the number of cores per node
    #  #2: NCORES_NEM: the number of cores for NEMO
    #  #3: N_XIOS_PN:  the number of xios processes you want to be running on each node
    #  #4: strict (1) or not (0)
    #
    # Returns (actually UPDATES!):
    # NNODES_NEM  => the number of cores for NEMO
    # NCORES_XIO  => the number of cores for xios
    # NCORES_TOT  => the total number of cores needed for the job (NEMO + XIOS)
    # NNODES_TOT  => the total number of nodes needed for the job (NEMO + XIOS)
    # NCORES_BOOK => the total number of cores to be booked (based on NNODES_TOT)
    #
    istrict=${4}
    #
    if [ ${2} -lt 2000 ]; then
        echo " PROBLEM: * n_cores_xio_nodes_tot * => Not enough core for NEMO..."; exit
    fi
    #
    # 1 Number of nodes for NEMO:
    #NNODES_XIO=$((${NNODES_XIO} > 1 ? ${NNODES_XIO}:1))
    export NNODES_NEM=$((${2}/${1}))
    if [ $((${2}%${1})) -gt 0 ]; then export NNODES_NEM=$((${NNODES_NEM}+1)); fi
    #echo "Number of nodes for NEMO:"
    #echo ${NNODES_NEM}
    LEFTO=$((${NNODES_NEM}*${1} - ${2}))
    #echo "Free leftover cores:"
    #echo ${LEFTO}
    NCORES_XIO=$((${3}*${NNODES_NEM}))
    #echo "Initial guess of # cores XIOS:"
    #echo ${NCORES_XIO}
    if [ ${NCORES_XIO} -ge $((${LEFTO}+${1}/4)) ]; then
        nextra_c=$((${NCORES_XIO}-${LEFTO}))
        #echo "Then, we need that amount of extra cores:"
        #echo ${nextra_c}
        nextra_n=$((${nextra_c}/${1}))
        if [ $((${nextra_c}%${1})) -ge $((${1}/2)) ]; then nextra_n=$((${nextra_n}+1)); fi
        #echo " ==> nodes for xios => ${nextra_n}"
        export NCORES_XIO=$((${LEFTO}+${nextra_n}*${1}))
    else
        export NCORES_XIO=${LEFTO}
    fi
    #echo "Final # cores XIOS:"
    #echo ${NCORES_XIO}
    #echo
    #
    export NCORES_TOT=$((${2}+${NCORES_XIO}))
    if [ $((${NCORES_TOT}%${1})) -ne 0 ]; then echo "PROBLEM!!! (# of xios processes! => n_cores_xio_nodes_tot)"; exit; fi

    export NNODES_TOT=$((${NCORES_TOT}/${1}))
    if [ $((${NCORES_TOT}%${1})) -gt 0 ]; then export NNODES_TOT=$((${NNODES_TOT}+1)); fi
    export NCORES_BOOK=$((${NNODES_TOT}*${1}))


    l_special=false
    if [ ${istrict} -eq 1 ]; then

        if [ ${NCORES_XIO} -lt $((${3}*${NNODES_TOT})) ]; then
            #echo "PROBLEM: NCORES_XIO < ${3}*NNODES_TOT !!! ${NCORES_XIO} $((${3}*${NNODES_TOT}))"
            l_special=true
            #nc_to_add=$(($((${3}*${NNODES_TOT}))-${NCORES_XIO}))
            ##echo " => need to add ${nc_to_add} cores!"
            ##echo " => adding an extra node for xios!"
            #NCORES_XIO=$((${NCORES_XIO}+${1}-${nc_to_add}))
            ##echo " => NCORES_XIO=${NCORES_XIO}"
            #export NCORES_TOT=$((${2}+${NCORES_XIO}))
            #export NNODES_TOT=$((${NCORES_TOT}/${1}))
            #if [ $((${NCORES_TOT}%${1})) -gt 0 ]; then export NNODES_TOT=$((${NNODES_TOT}+1)); fi
            #export NCORES_BOOK=$((${NNODES_TOT}*${1}))

        else

            # => there is more XIOS processes than " total number of nodes booked x
            # N_XIOS_PN " and we don't want that

            #echo
            #echo " initial NCORES_XIO=${NCORES_XIO}"
            #echo " initial NCORES_TOT=${NCORES_TOT}"
            #echo " initial NCORES_BOOK=${NCORES_BOOK}"
            nx_to_remove=$((${NCORES_XIO}-${3}*${NNODES_TOT}))
            #echo "   => will remove $nx_to_remove xios procs!"
            export NCORES_XIO=$((${NCORES_XIO}-${nx_to_remove}))
            #echo " new NCORES_XIO=${NCORES_XIO}"
            export NCORES_TOT=$((${2}+${NCORES_XIO}))
            #echo " new NCORES_TOT=${NCORES_TOT}"
            nidle=$((${NNODES_TOT}*${1}-${NCORES_TOT}))
            #echo " initial nidle=${nidle}"
            while [ ${nidle} -ge ${1} ]; do
                export NNODES_TOT=$((${NNODES_TOT}-1))     ; #echo "   updated NNODES_TOT=${NNODES_TOT}"
                export NCORES_XIO=$((${3}*${NNODES_TOT}))  ; #echo "   updated NCORES_XIO=${NCORES_XIO}"
                nidle=$((${nidle}-${1}))                   ; #echo "   updated nidle=${nidle}"
            done
            if [ ${nidle} -gt 0 ]; then
                echo " *** WARNING [n_cores_xio_nodes_tot()]: we will have ${nidle} idle processors! (istrict==1)"
            fi

            #if [ $((${NCORES_TOT}%${1})) -ne $((${1}-${nx_to_remove})) ] ; then echo "PROBLEM 0!!! (# of xios processes! => n_cores_xio_nodes_tot) $((${NCORES_TOT}%${1})) $((${1}-${nx_to_remove})) "; exit; fi
        fi
    fi

    # Verification:
    if [ ${NCORES_TOT} -gt ${NCORES_BOOK} ]         ; then echo "PROBLEM 1!!! (# of xios processes! => n_cores_xio_nodes_tot) ${NCORES_TOT} ${NCORES_BOOK}"; exit; fi
    if [ $((${NCORES_BOOK}%${1})) -ne 0 ]           ; then echo "PROBLEM 2!!! (# of xios processes! => n_cores_xio_nodes_tot)"; exit; fi
    if [ $((${NCORES_BOOK}/${1})) -ne ${NNODES_TOT} ]; then echo "PROBLEM 3!!! (# of xios processes! => n_cores_xio_nodes_tot) $((${NCORES_TOT}/${1})) ${NNODES_TOT} "; exit; fi
    if [ $((${2}+${NCORES_XIO})) -ne ${NCORES_TOT} ]; then echo "PROBLEM 4!!! (# of xios processes! => n_cores_xio_nodes_tot) $((${NCORES_NEM}+${NCORES_XIO})) ${NCORES_TOT}"; exit; fi

    NEW_N_XIOS_PN=$((${NCORES_XIO}/${NNODES_TOT}))
    if [ ${NEW_N_XIOS_PN} -ne ${3} ] && [ ! $l_special ];  then echo "PROBLEM 5!!! (n_cores_xio_nodes_tot) ${NEW_N_XIOS_PN} ${3}"; exit; fi

    echo
    echo " *** As found by function n_cores_xio_nodes_tot() :"
    echo "    - total # of Nodes to book         ---> ${NNODES_TOT}"
    echo "    - of which for NEMO                ---> ${NNODES_NEM}"
    echo "    - # of xios processes              ---> ${NCORES_XIO}"
    n_xios_os=$((${NCORES_XIO}%${NNODES_TOT}))
    echo "    - How many XIOS processes per node ---> $((${NCORES_XIO}/${NNODES_TOT})) (${n_xios_os})"
    echo "    - Overshoot of XIOS processes      ---> ${n_xios_os}"
    if [ ${n_xios_os} -gt 0 ]; then
        echo "     => number of nodes with ${3} xios processes  ---> ${n_xios_os}"
        echo "     => number of nodes with only $((${3}-1)) xios processes  ---> $((${NNODES_TOT}-${n_xios_os}))"
    fi
    echo "    - How many idle processors         ---> $((${NNODES_TOT}*${1}-${NCORES_NEM}-${NCORES_XIO}))"
    echo
}


lb_is_leap()
{
    if [ "$1" = "" ]; then echo "USAGE: lb_is_leap <YEAR>"; exit; fi
    #
    i_mod_400=`expr ${1} % 400`
    i_mod_100=`expr ${1} % 100`
    i_mod_4=`expr ${1} % 4`
    #
    if [ ${i_mod_400} -eq 0 -o ${i_mod_4} -eq 0 -a ! ${i_mod_100} -eq 0 ]; then
        echo "1"
    else
        echo "0"
    fi
}



function nb_days_in_month()
{
    # Input is an integer between 1 and 366
    nm=${1}
    yyyy=${2}
    #
    if [ "${2}" = "" ]; then
        echo "ERROR in nb_days_in_month of lb_functions.sh !!!"
        echo "  ==> specify the year as 2nd argument!"; exit
    fi
    if [ ${nm} -lt 1 ] || [ ${nm} -gt 12 ]; then
        echo "ERROR in nb_days_in_month of lb_functions.sh !!!"
        echo "  ==> month # ${nm} means nothing!!!" ; exit
    fi
    #
    VMN=( 31 28 31 30 31 30 31 31 30 31 30 31 )
    if [ `lb_is_leap ${yyyy}` -eq 1 ]; then
        VMN=( 31 29 31 30 31 30 31 31 30 31 30 31 )
    fi

    jm=`expr ${nm} - 1`
    echo ${VMN[${jm}]}
    #
}



function day_year_to_MMDD()
{
    # Input is an integer between 1 and 366
    nid=${1}
    yyyy=${2}

    if [ "${2}" = "" ]; then
        echo "ERROR in day_year_to_MMDD of lb_functions.sh !!!"
        echo "  ==> specify the year as 2nd argument!"; exit
    fi

    leap=false
    ndays=365
    VMN=( 31 28 31 30 31 30 31 31 30 31 30 31 )
    if [ `lb_is_leap ${yyyy}` -eq 1 ]; then
        leap=true
        ndays=366
        VMN=( 31 29 31 30 31 30 31 31 30 31 30 31 )
    fi

    if [ ${nid} -lt 1 ] || [ ${nid} -gt ${ndays} ]; then
        echo "ERROR in day_year_to_MMDD of lb_functions.sh !!!"
        echo "  ==> day of year # ${nid} means nothing!!!" ; exit
    fi

    id=0
    jd=0
    jm=0
    while [ ${id} -lt ${nid} ]; do
        id=$((${id}+1))
        jd=$((${jd}+1))
        if [ $((${jd}-1)) -eq ${VMN[${jm}]} ]; then
            jd=1 ; jm=$((${jm}+1))
        fi
    done
    jm=$((${jm}+1))
    cm=`printf "%02d" ${jm}`
    cd=`printf "%02d" ${jd}`
    echo "${cm}${cd}"
}



function MMDD_to_day_year()
{
    # Input is 4-digit integer: MMDD
    cmmdd=${1}
    yyyy=${2}

    if [ "${2}" = "" ]; then
        echo "ERROR in MMDD_to_day_year of lb_functions.sh !!!"
        echo "  ==> specify the year as 2nd argument!"; exit
    fi

    leap=false
    ndays=365
    VMN=( 31 28 31 30 31 30 31 31 30 31 30 31 )
    VMC=( 0 31 59 90 120 151 181 212 243 273 304 334 )
    if [ `lb_is_leap ${yyyy}` -eq 1 ]; then
        leap=true
        ndays=366
        VMN=( 31 29 31 30 31 30 31 31 30 31 30 31 )
        VMC=( 0 31 60 91 121 152 182 213 244 274 305 335 )
    fi

    nbd=`echo ${cmmdd} | wc -c`
    if [ ${cmmdd} -lt 101 ] || [ ${cmmdd} -gt 1231 ] || [ ! ${nbd} -eq 5 ]; then
        echo "ERROR in MMDD_to_day_year of lb_functions.sh !!!"
        echo "  ==> day of year # ${cmmdd} means nothing!!!" ; echo "x${cmmdd}x"; exit
    fi

    im=`echo ${cmmdd} | cut -c1-2` ; im=`expr ${im} + 0`
    id=`echo ${cmmdd} | cut -c3-4` ; id=`expr ${id} + 0`

    iday=`expr ${VMC[$((${im}-1))]} + ${id}`

    echo "${iday}"
}












function lb_wait()
{
    if [ "$1" = "" ]; then
        echo "USAGE: lb_wait <run name (CASE)>"
        exit
    fi

    #cuser='bsc3232'
    cuser="${USER}"

    # Cutting to 8 characters
    run_name=`echo $1 | cut -c1-8`
    
    if [ "`which squeue`" = "" ]; then
        clrun=`ps ax | grep -i nemo.exe | grep mpirun` ; sleep 1
        clrun2=`ps ax | grep -i nemo.exe | grep mpirun`
    else        
        clrun=`squeue -u ${cuser} | grep "\ ${run_name}\ " | grep -v "\ CG\ "`; sleep 1
        clrun2=`squeue -u ${cuser} | grep "\ ${run_name}\ " | grep -v "\ CG\ "`
    fi
    #
    cmsg0="No simulation called ${run_name} is running!"
    cmsg1="A simulation called ${run_name} is running!"
    #
    if [ "${clrun}" = "" -a "${clrun2}" = "" ]; then
        echo ${cmsg0}
    else
        while [ "${clrun}" != "" ]; do
            echo "${cmsg1}! Waiting...   `date`"
            sleep 60
            if [ "`which squeue`" = "" ]; then
                clrun=`ps ax | grep -i nemo.exe | grep mpirun`
            else
                clrun=`squeue -u ${cuser} | grep "\ ${run_name}\ " | grep -v "\ CG\ "`
            fi
        done
        echo ${cmsg0}
    fi
    echo
}





function lb_wait_dtrmr()
{
    if [ "$1" = "" ]; then echo "USAGE: lb_wait_dtrmr <run name (CASE)>"; exit; fi
    cuser="${USER}"

    # Cutting to 8 characters
    run_name=`echo $1 | cut -c1-8`
    clrun1=`qstat -u ${cuser} | grep "\ ${run_name}\ "` ; # | grep -v "\ CG\ "`
    sleep 1
    clrun2=`qstat -u ${cuser} | grep "\ ${run_name}\ "` ; # | grep -v "\ CG\ "`
    cmsg0="No simulation called ${run_name} is running!"
    cmsg1="A simulation called ${run_name} is running!"
    if [ "${clrun1}" = "" -a "${clrun2}" = "" ]; then
        echo ${cmsg0}
    else
        while [ "${clrun1}" != "" ]; do
            echo "${cmsg1}! Waiting...   `date`"
            sleep 60
            clrun1=`qstat -u ${cuser} | grep "\ ${run_name}\ "` ; # | grep -v "\ CG\ "`
        done
        echo ${cmsg0}
    fi
    echo
}



#
#
lb_num_leap()
{
    #
    # number of leap years comprised in <YEAR1> (included) and <YEAR2> (excluded)...
    #
    if [ "$2" = "" ]; then echo "USAGE: num_leap <YEAR1> <YEAR2>"; exit; fi
    icpt=0
    jy=${1}
    while [ ${jy} -lt ${2} ]; do
        if [ `lb_is_leap ${jy}` -eq 1 ]; then icpt=`expr ${icpt} + 1`; fi
        jy=`expr ${jy} + 1`
    done
    echo ${icpt}
}
#
#
function lb_leap_day()
{
    # We need 2 different methods to know the current date:
    # The input argument is the file date.file
    #
    if [ "$1" = "" ]; then echo "USAGE: lb_check_leap <date.file>"; exit; fi
    DATF=$1
    #
    if [ ! -f ${DATF} ]; then echo "There should be a ${DATF} !"; exit; fi
    #
    dtg1=`paste ${DATF} | cut -c 11-18`
    dtg2=`paste ${DATF} | cut -c 20-27`
    #
    yyyy=`echo ${dtg2} | cut -c1-4`
    mmdd1=`echo ${dtg1} | cut -c5-8`
    mmdd2=`echo ${dtg2} | cut -c5-8`
    #
    ileap=`expr ${yyyy} - 2000`; ileap=`expr ${ileap} % 4`
    #
    # Last DTG done:
    export CLDTG=${dtg2}
    #
    #if [ ! ${ileap} -eq 0 ]; then echo "Not a leap year !"; fi
    #
    #
    # We were at a DTBR=
    # a regular year would make a 19900225_19900229 file
    #  (should be 19900225_19900301 actually but NEMO seems to keep the same month for a file name!)
    # so when we restart a run for a leap year and the last dtg was "19900220_19900224", we force
    # DTBR=6 !
    #
    export i_leap_day=0
    #
    if [ ${ileap} -eq 0 -a "${mmdd2}" = "0224" ]; then
        #
        # Checking if the frequency was 5 days:
        if [ ! "${mmdd1}" = "0220" ]; then echo "Problem (lb_leap_day), expecting different DTG"; exit; fi
        #
        export i_leap_day=1
    fi
    #
}
#
#
function lb_ec_nodes()
{
    # 1: tot numb of procs for nemo
    # 2: tot numb of procs for ifs
    # 3: how many procs to use on one node
    #
    if [ "$3" = "" ]; then echo "USAGE: lb_ec_proc_node <Nb Pr. NEMO> <Nb Pr. IFS> <Pr. per node>"; exit; fi
    #
    # Number of nodes for NEMO:
    if [ ! `expr ${1} % ${3}` -eq 0 ]; then
        NNN=`expr ${1} / ${3} + 1`
    else
        NNN=`expr ${1} / ${3}`
    fi
    #
    if [ ! `expr ${2} % ${3}` -eq 0 ]; then
        NNI=`expr ${2} / ${3} + 1`
    else
        NNI=`expr ${2} / ${3}`
    fi
    #
    echo
    echo " Total number of procs for NEMO = ${1}"
    echo "   => number of nodes for NEMO  = ${NNN}"
    #
    echo " Total number of procs for IFS  = ${2}"
    echo "   => number of nodes for IFS   = ${NNI}"
    echo
    #
    # Total number of nodes requested by NEMO and IFS:
    TNN=`expr ${NNN} + ${NNI}`
    #
}
#
#
function lb_loop_year_old()
{
    # 1: Y1
    # 2: Y2
    # 3: True year:
    #
    # Return fake year!
    #
    if [ "$3" = "" ]; then echo "USAGE: lb_loop_year <Y1> <Y2> <jy>"; exit; fi
    #
    Y1=${1}; Y2=${2}; jy=${3}
    #
    nby=`expr ${Y2} - ${Y1} + 1`
    jyf=${jy}
    #
    jdiff=`expr ${jy} - ${Y2}`
    #
    if [ ${jdiff} -ge 1 ]; then
        na=`expr \( ${jdiff} - 1 \) % ${nby}`
        ia=`expr \( ${jdiff} - 1 \) / ${nby}`
        jyf=`expr ${Y1} + ${jdiff} - 1 - \( ${ia} \* ${nby} \)`
    fi
    #
    echo ${jyf}
    #
}
#
#
#
function clock2hour()
{
    # convert time like hh:mm:ss
    # to decimal hour of the day
    # ex: 13:30:00  => 13.5
    #
    # Replacing ":" by " ":
    ca=`echo $1 | sed -e s/:/' '/g`
    #
    hdec=`echo ${ca} | awk '{print $1+($2*60+$3)/3600}'`
    echo ${hdec}
}
#


function lb_loop_year()
{
    Y0=0
    Y1f=1958
    Y2f=1983
    #
    nbyf=`expr ${Y2f} - ${Y1f} + 1`
    #
    if [ "$1" == "" ]; then
        echo "USAGE: $0 <model year>"; exit
    fi
    #
    YE=${1}
    #
    # Direct method:
    ir=`expr \( ${YE} - ${Y0} \) % ${nbyf}`
    #
    echo `expr ${Y1f} + ${ir}`
}
#

function nb_day_in_year()
{
    yyyy=${1}
    if [ `lb_is_leap ${yyyy}` -eq 1 ]; then
        echo "366"
    else
        echo "365"
    fi
}



function dtg_itt()
{
    DTG0=$1
    #
    VMN=( 31 28 31 30 31 30 31 31 30 31 30 31 )
    VML=( 31 29 31 30 31 30 31 31 30 31 30 31 )
    #

    YY=`echo ${DTG0}|cut -c1-4`; MM=`echo ${DTG0}|cut -c5-6`; DD=`echo ${DTG0}|cut -c7-8`
    jm=`expr ${MM} - 1` ; nbd=${VMN[${jm}]}
    if [ `lb_is_leap ${YY}` -eq 1 ]; then nbd=${VML[${jm}]}; fi
    if [ ${DD} -eq ${nbd} ]; then
        DD=1
        if [ ${MM} -eq 12 ]; then
            MM=1; jm=0; YY=`expr ${YY} + 1`
        else
            MM=`expr ${MM} + 1`; jm=`expr ${jm} + 1`
        fi
    else
        DD=`expr ${DD} + 1`
    fi
    MM=`expr ${jm} + 1`; YY=`printf "%02d" "${YY}"`; MM=`printf "%02d" "${MM}"`; DD=`printf "%02d" "${DD}"`
    echo "${YY}${MM}${DD}"
}


function fnd2m()
{
    sdate="$1"
    #
    cmmdd=`echo ${sdate} | cut -c5-8`
    jsdy=`MMDD_to_day_year ${cmmdd} ${JY}` ;# echo "LOLO: jsdy=${jsdy}"
    fnd2m=61 ; # TERRIBLE :( !!!! #fixme
    if [ "${cmmdd}" = "0101" ] || [ "${cmmdd}" = "0201" ]; then
        ndf=`nb_days_in_month 2 ${JY}`
        fnd2m=$((59-28+ndf)) ; # Normal => 59 days, Leap => 60 days
    fi
    if [ ${cmmdd} -gt 0931 ]; then
        echo "lb_functions.sh->fnd2m: FIX ME: we want 3-month chunks but case not ready!!!"
        fnd2m=-1
    fi
    echo "${fnd2m}"
}



function fnd3m()
{
    sdate="$1"
    #
    cmmdd=`echo ${sdate} | cut -c5-8`
    jsdy=`MMDD_to_day_year ${cmmdd} ${JY}` ;# echo "LOLO: jsdy=${jsdy}"
    fnd3m=89 ; # TERRIBLE :( !!!! #fixme
    if [ "${cmmdd}" = "0101" ] || [ "${cmmdd}" = "0201" ]; then
        ndf=`nb_days_in_month 2 ${JY}`
        fnd3m=$((90-28+ndf)) ; # Normal => 90 days, Leap => 91 days
    fi
    if [ ${cmmdd} -gt 0931 ]; then
        echo "lb_functions.sh->fnd3m: FIX ME: we want 3-month chunks but case not ready!!!"
        fnd3m=-1
    fi
    echo "${fnd3m}"
}

function fnd6m()
{
    sdate="$1"
    yyyy=`echo ${sdate} | cut -c1-4`
    mm=`echo ${sdate} | cut -c5-6`
    dd=`echo ${sdate} | cut -c7-8`
    #
    if   [ "${mm}" = "01" ]; then
        ndf=`nb_days_in_month 2 ${yyyy}`
        fnd6m=$((181-28+ndf)) ; # Normal => 181 days, Leap => 182 days
    elif [ "${mm}" = "07" ]; then
        fnd6m=184 ; # 1st of July to 31st of Dec., always valid!
    else
        echo "lb_functions.sh->fnd6m: ERROR: we were expected current month to start from to be '1' or '7'!"; echo " ===> ${mm} !"; exit
        fnd6m=-1
    fi
    echo "${fnd6m}"
}


function fnb_days()
{
    slfrq=$1 ; # string for length (ex: 1d, 1m, 3m, ..., 1y)
    sdate=$2
    #
    yyyy=`echo ${sdate} | cut -c1-4`
    mm=`echo ${sdate} | cut -c5-6`
    mm=`expr ${mm} + 0`
    #
    ct=${slfrq}
    if [ "`echo ${slfrq} | rev | cut -c1-1`" = "d" ]; then ct="Xd"; nXd=`echo ${slfrq} | cut -d d -f1`; fi
    #
    case ${ct} in
        #
        "1y") nnd=`nb_day_in_year ${yyyy}` ; # 1 year of simulation
              ;;
        "6m") nnd=`fnd6m ${sdate}`
              ;;
        "3m") nnd=`fnd3m ${sdate}`
              ;;
        "2m") nnd=`fnd2m ${sdate}`
              ;;        
        "1m") nnd=`nb_days_in_month ${mm} ${yyyy}` ; # 1 calendar month
              ;;
        "Xd") nnd=${nXd} ; # X days of simulation
              ;;
        *) echo "Unknow NDAYS_EXP or nbd_brs! => ${slfrq}" ; exit
           nnd=-1
           ;;
    esac
    echo ${nnd}
}
