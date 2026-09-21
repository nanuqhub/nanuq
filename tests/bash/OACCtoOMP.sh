#!/usr/bin/env bash

set +H


if [ "${VERSION_OMP}" != "4.x" ] && [ "${VERSION_OMP}" != "5.1" ] && [ "${VERSION_OMP}" != "cray" ]; then
    echo; echo "Please set the environment variable VERSION_OMP (4.x or 5.1 or cray) !"; echo
    exit
fi

if [ "${3}" = "" ]; then
    echo " $0 <file.F90> <TARGET_DIR> <file_out.F90>"
    exit
fi

FIN="${1}"
DIR_OUT="${2}"
FOUT="${3}"

dirS=`dirname ${FIN}`
finS=`basename ${FOUT}`


if [ "${DIR_OUT}" = "${dirS}" ]; then
    echo
    echo " <TARGET_DIR> must be different that the directory of your input field!"
    exit
fi

mkdir -p ${DIR_OUT}

OUTFILE="${DIR_OUT}/${finS}"

rm -f ${OUTFILE}


function ConvertClauses()
{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    # Converts a suite of OpenACC clauses to their OpenMP counterparts
    #
    #   Take 2 arguments:
    #
    #   * $1 => type of ACC directive
    #        => can be 'LOOP', 'DATA' or 'EnterDATA'
    #
    #   * $2 =>
    #
    #   Returns:
    #   * the appropriate OpenMP conterpart directive
    #    or
    #   * an error message
    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ACCcmd="${1}" ; #

    VCacc=( ${2} ) ; # making it a vector
    VComp=( ${2} ) ; # initialisation of OMP counterpart clauses array...

    if [ "${ACCcmd}" != "LOOP" ] && [ "${ACCcmd}" != "DATA" ] && [ "${ACCcmd}" != "EnterDATA" ]; then
        errMSG=" unknown ACC command passed to 'ConvertClauses()' => ${ACCcmd}"
        iexit=2
        return ${iexit}
    fi

    ncls=`echo ${VCacc[*]} | wc -w` ; # number of clauses

    listSupported="present create copy pcopy copyin pcopyin copyout pcopyout collapse private reduction"

    iexit=0  ; # exit code
    nfound=0 ; # how many clauses to be identified ?
    ks=-1    ; # counter

    while [ ${nfound} -lt ${ncls} ]; do
        #
        ks=$((ks+1))
        clause=${VCacc[${ks}]} ;                       # the clause with its arguments...
        clause_func=`echo ${clause} | cut -d'(' -f1` ; # name of the clause
        #
        if [ "`echo ${listSupported} | grep ${clause_func}`" = "" ]; then
            iexit=1
            errMSG="do not know what to do with ACC clause: ${clause} !"
            break
        fi
        #
        ifnd=0
        #
        case "${clause_func}" in
            #
            #
            "present")
                crhs=`echo ${clause} | cut -d'(' -f2-`
                #
                # Default => nothing !
                VComp[${ks}]=""
                #
                case "${VERSION_OMP}" in
                    "5.1")
                        VComp[${ks}]="map(present:${crhs}"
                        ;;
                    "cray")
                        if [ "${ACCcmd}" = "DATA" ]; then
                            VComp[${ks}]="use_device_addr(${crhs}"
                        fi
                        ;;
                esac
                #
                ifnd=1
                nfound=$((nfound+1))
                ;;
            #
            "create")
                crhs=`echo ${clause} | cut -d'(' -f2-`
                VComp[${ks}]="map(alloc:${crhs}"
                ifnd=1
                nfound=$((nfound+1))
                ;;
            #
            "copy"|"pcopy")
                crhs=`echo ${clause} | cut -d'(' -f2-`
                VComp[${ks}]="map(tofrom:${crhs}"
                ifnd=1
                nfound=$((nfound+1))
                ;;
            #
            "copyin"|"pcopyin")
                crhs=`echo ${clause} | cut -d'(' -f2-`
                VComp[${ks}]="map(to:${crhs}"
                ifnd=1
                nfound=$((nfound+1))
                ;;
            #
            "copyout"|"pcopyout")
                crhs=`echo ${clause} | cut -d'(' -f2-`
                VComp[${ks}]="map(from:${crhs}"
                ifnd=1
                nfound=$((nfound+1))
                ;;
            #
            "collapse"|"private"|"reduction")
                #  ==> value ${VComp[${ks}]} is good then!
                ifnd=1
                nfound=$((nfound+1))
                ;;
        esac
        #
    done

    if [ ${iexit} -eq 0 ]; then
        echo ${VComp[*]}
    else
        echo ${errMSG}
    fi
    return ${iexit}
}



#####################################################
# 1/ Test presence of potential problematic strings #
#####################################################

#ca=`grep '0x0' ${FIN}`
#if [ "$ca" != "" ]; then echo "PROBLEM: '0x0' exists in file!"; exit; fi
#ca=`grep '|' ${FIN}`
#if [ "$ca" != "" ]; then echo "PROBLEM: '|' exists in file!"; exit; fi




#####################################################
# 2/ Line by line... Painfull but trustworthy...
#####################################################

nbl=`cat ${FIN} | wc -l`

echo
echo " *** ${FIN} has ${nbl} lines ! ***"
echo
echo "   => launching scan!"
echo; echo

ignore_next_data_end=0
is_real_data_clause=0

jl=0

while [ ${jl} -le ${nbl} ]; do

    jl=`expr ${jl} + 1`

    rwline=`cat ${FIN} | sed -n ${jl}p`
    trline="${rwline#"${rwline%%[![:space:]]*}"}" ; # our line trimmed from any preceeding white space...
    cspace=""


    if [ "$(echo ${trline} | cut -c 1-5)" = '!$acc' ]; then


        echo; echo
        echo " *** Found an acc directive at line # ${jl} of ${FIN} ***"
        echo "      ==> raw line = >>${rwline}<<"
        echo

        # 0/ Keep in memory the length of "blank space" that was right before `!$acc`:
        cspace="${rwline%%[![:space:]]*}" ;    #    echo " ==> cspace = >>${cspace}<<"

        # A/ Delete RHS comments on an ACC directive line:
        trline=`echo ${trline} | cut -d'!' -f1-2` # second occurence because already a `!` in `!$acc` !!!

        # B/ Remove potential trailing white spaces:
        accline="${trline%"${trline##*[![:space:]]}"}"

        # C/ remove the unnecessary spaces...
        nline=`echo ${accline} | sed -e s/'  '/''/g  -e s/'( '/'('/g  -e s/' )'/')'/g  \
                                     -e s/', '/','/g  -e s/' , '/','/g  -e s/' ,'/','/g \
                                     -e s/'device ('/'device('/g -e s/'self ('/'self('/g \
                                     -e s/'present ('/'present('/g`

        accline="${nline}"

        echo " >>>${accline}<<<"

        ompline=""


        ###################################
        #  `!$acc routine`                #
        ###################################
        ctst=`echo ${accline} | cut -c 1-13`
        if [ "${ctst}" = '!$acc routine' ]; then
            ompline="!\$omp declare target"
            # ==> same for '!$acc routine seq'
        fi

        ############################################
        #  `!$acc declare create()`                #
        #  -> copy a module scalar value to device #
        ############################################
        ctst=`echo ${accline} | cut -c 1-21`
        if [ "${ctst}" = '!$acc declare create(' ]; then
            crhs=`echo ${accline} | cut -d'(' -f2-`
            ompline="!\$omp declare target(${crhs}"
        fi


        ############################
        #  `!$acc enter data ...`  #
        ############################
        ctst=`echo ${accline} | cut -c 1-17`
        if [ "${ctst}" = '!$acc enter data ' ]; then
            #
            crest=`echo ${accline} | sed -e s/'!$acc enter data '/''/g` ;# what comes after `!$acc enter data`...
            crest="${crest#"${crest%%[![:space:]]*}"}"
            #
            #DEBUG:
            #echo "LOLO: calling ConvertClauses EnterDATA ${crest}"
            #VC=( "${crest}" )
            #for cc in ${VC[*]}; do
            #    echo " *lolo cc = >>${cc}<<"
            #    #clause=${VCacc[${ks}]} ;                       # the clause with its arguments...
            #    clause_func=`echo ${cc} | cut -d'(' -f1` ; # name of the clause
            #    echo " clause_func=>>${clause_func}<<"
            #    crhs=`echo ${cc} | cut -d'(' -f2-`
            #    echo " crhs=>>${crhs}<<"
            #done
            #echo
            #DEBUG.
            comp=`ConvertClauses EnterDATA "${crest}"`
            if [ $? -ge 1 ]; then echo " ERROR: ${comp}"; echo; exit; fi

            ompline="!\$omp target enter data ${comp}"

            echo
        fi; #if [ "${ctst}" = '!$acc enter data' ]



        ######################
        #  `!$acc data ...`  #
        ######################
        ctst=`echo ${accline} | cut -c 1-11`
        if [ "${ctst}" = '!$acc data ' ]; then
            #
            crest=`echo ${accline} | sed -e s/'!$acc data '/''/g` ;# what comes after `!$acc data`...
            crest="${crest#"${crest%%[![:space:]]*}"}"
            #
            # Special case here in case `!$acc data ...` is only with a `present` clause !
            VC=( ${crest} ) ; cls=""
            ncls=`echo ${VC[*]} | wc -w` ; # number of clauses
            if [ ${ncls} -eq 1 ]; then cls=` echo "${VC[0]}" | cut -d'(' -f1`; fi
            #if [ "${VERSION_OMP}" != "5.1" ] && [ ${ncls} -eq 1 ] && [ "${cls}" = "present" ]; then
            #    ompline="!Nothing for OMP & NVIDIA compilo! (data present)"
            #    ignore_next_data_end=`expr ${ignore_next_data_end} + 1`
            #else
            comp=`ConvertClauses DATA "${crest}"`
            if [ $? -ge 1 ]; then echo " ERROR: ${comp}"; echo; exit; fi
            #
            ompline="!\$omp target data ${comp}"
            #
            #is_real_data_clause=`expr ${is_real_data_clause} + 1`
            #fi

            #
            echo
        fi; #if [ "${ctst}" = '!$acc data' ]


        ######################
        #  `!$acc end data`  #
        ######################
        ctst=`echo ${accline} | cut -c 1-14`
        if [ "${ctst}" = '!$acc end data' ]; then
            #if [ ${ignore_next_data_end} -eq 0 ] || [ ${is_real_data_clause} -gt 0 ]; then
            ompline="!\$omp end target data"
            #    is_real_data_clause=`expr ${is_real_data_clause} - 1`
            #elif [ ${ignore_next_data_end} -gt 0 ]; then
            #    ompline="!Nothing for OMP & NVIDIA compilo! (end data present)"
            #    ignore_next_data_end=`expr ${ignore_next_data_end} - 1`
            #else
            #    echo " * What the fuck are we doing here???"
            #    exit
            #fi
        fi
        if [ "${ctst}" = '!$acc loop seq' ]; then
            ompline="!Nothing for OpenMP (loop seq)"
        fi

        #################################
        #  `!$acc parallel loop ...`    #
        #  -> 2D loop                   #
        #################################
        ctst=`echo ${accline} | cut -c 1-19`
        if [ "${ctst}" = '!$acc parallel loop' ]; then
            #
            # Does this loop directive contain more than that
            crest=`echo ${accline} | sed -e s/'!$acc parallel loop'/''/g` ;# what comes after `!$acc parallel loop`...

            if [ "${crest}" = "" ]; then
                echo " --> easy :)"
                ompline="!\$omp target teams distribute parallel do"
                #
            else
                #
                comp=`ConvertClauses LOOP "${crest}"`
                if [ $? -ge 1 ]; then echo " ERROR: ${comp}"; echo; exit; fi
                #
                ompline="!\$omp target teams distribute parallel do ${comp}"
                #
            fi
            echo
        fi; #if [ "${ctst}" = '!$acc parallel loop' ]



        ###############################
        #  `!$acc end parallel loop`  #
        ###############################
        ctst=`echo ${accline} | cut -c 1-23`
        if [ "${ctst}" = '!$acc end parallel loop' ]; then
            ompline="!\$omp end target teams distribute parallel do"
        fi



        ############################
        #  `!$acc update device()`  #
        ############################
        ctst=`echo ${accline} | cut -c 1-20`
        if [ "${ctst}" = '!$acc update device(' ]; then
            crhs=`echo ${accline} | cut -d'(' -f2-`
            ompline="!\$omp target update to(${crhs}"
        fi


        ############################
        #  `!$acc update self()`  #
        ############################
        ctst=`echo ${accline} | cut -c 1-18`
        if [ "${ctst}" = '!$acc update self(' ]; then
            crhs=`echo ${accline} | cut -d'(' -f2-`
            ompline="!\$omp target update from(${crhs}"
        fi


        ######################################################################################


        if [ "${ompline}" != "" ]; then
            echo "   ===> going to be translated with:"
            echo " >>>${ompline}<<<"
        else
            echo " * ERROR  ===> DID NOT FIND AN OPENMP COUNTERPART FOR THAT!"

            exit
        fi

        echo; echo; echo

        printf "%s\n" "${cspace}${ompline}" >> ${OUTFILE}

    else

        printf "%s\n" "${rwline}"           >> ${OUTFILE}

    fi ; # if [ "${cbeg}" = '!$acc' ]


done ; # while [ ${jl} -le ${nbl} ]


echo ; echo
echo " *** File: ${OUTFILE} created !!!"
echo ; echo ; echo ; echo
