#!/usr/bin/env bash

#set +H

i_do_sc=1 ; # semicolon
i_do_1d=1
i_do_2d=1
i_do_3d=1
#
i_do_a1d=1 ; # A1Di(?) & A1Dj(?)
i_do_a2d=1 ; # A2D(?)
i_do_t2d=1 ; # T2D(?)


#(custom-set-variables
# ;; custom-set-variables was added by Custom.
# ;; If you edit it by hand, you could mess it up, so be careful.
# ;; Your init file should contain only one such instance.
# ;; If there is more than one, they won't work right.
# '(f90-auto-keyword-case (quote upcase-word))
# '(f90-continuation-indent 3)
# '(f90-do-indent 3)
# '(f90-if-indent 3)
# '(f90-program-indent 3)
# '(f90-smart-end (quote blink))
# '(f90-type-indent 3)


function indent_clean()
{
    emacs -batch ${1} --eval '(delete-trailing-whitespace)'                -f save-buffer
    emacs -batch ${1} --eval "(progn (custom-set-variables '(f90-continuation-indent 3) '(f90-do-indent 3) '(f90-if-indent 3) '(f90-program-indent 3) '(f90-type-indent 3)) (indent-region (point-min) (point-max) nil))" -f save-buffer
    emacs -batch ${1} --eval '(untabify (point-min) (point-max))'          -f save-buffer
    #
}


function extract_func_args()
{
    #------------------------------------------------------------------------------------------
    # Argument # 1 :  a string containing something like A1D(whatever) or A2D(whatever)
    # Argument # 2 :  what we are intested in => "A2D", "A1D", "A1Di", "A1Dj", etc...
    # RETURNS      :  exactly what's inside the parenthesis of the A*D(...)
    #------------------------------------------------------------------------------------------
    local s="$1"
    local func="$2"

    local i start depth=0 result=""
    local pattern="${func}("

    # locate "func("
    start=${s%%"$pattern"*}

    # function not found
    [[ "$start" == "$s" ]] && return 1

    i=${#start}

    # move past function name
    i=$(( i + ${#func} ))

    # verify and skip '('
    [[ "${s:i:1}" == "(" ]] || return 1
    ((i++))

    depth=1

    while (( i < ${#s} && depth > 0 )); do
        c="${s:i:1}"

        if [[ "$c" == "(" ]]; then
            ((depth++))
        elif [[ "$c" == ")" ]]; then
            ((depth--))
        fi

        # collect only inside outer func(...)
        if (( depth > 0 )); then
            result+="$c"
        fi

        ((i++))
    done

    echo "$result"
}


################################################
################################################





if [ "${1}" = "" ]; then
    echo " $0 <file.F90>"
    exit
fi

FIN="${1}"

OUTFILE=CLEANED_${1}

rm -rf ${OUTFILE}
mkdir -p tmp
rsync -avP ${FIN} tmp/
cd tmp


#####################################################
# 1/ Test presence of potential problematic strings #
#####################################################

# the '*' character is going to be replaced with '0x0' when working,
# so ensuring that '0x0' does not already exist in the file...
ca=$( grep '0x0' ${FIN} )
if [ "$ca" != "" ]; then echo "PROBLEM: '0x0' exists in file!"; exit; fi


# Will replace '|' with '`/`' when working, so same shit...
ca=$( grep "\`/\`" ${FIN} )
if [ "$ca" != "" ]; then echo "PROBLEM: '`/`' exists in file!"; exit; fi



i_pipe_fnd=0
ca=$( grep '|' ${FIN} )
if [ "$ca" != "" ]; then
    echo
    i_pipe_fnd=1
    rm -f ${FIN}.tmp
    echo " * INFO: '|' exists in file!"
    echo "  => will momentarily replace it with '\`/\`' !"
    sed -e "s|\||\`/\`|g" ${FIN} > ${FIN}.tmp
    mv -f ${FIN}.tmp ${FIN} ; # safe because we working into the ./tmp sub-dir...
    echo
fi



#####################################################
# 2/ General cleaning:
#####################################################

NFIN="NEW__${FIN}"

sed -e s/'IF ('/'IF('/g  -e s/'WHERE ('/'WHERE('/g  -e s/'END IF'/'ENDIF'/g \
    -e s/'\.eq\.'/'=='/gI -e s/'\.ne\.'/'\/='/gI \
    -e s/'\.le\.'/'<='/gI -e s/'\.ge\.'/'>='/gI -e s/'\.lt\.'/'<'/gI -e s/'\.gt\.'/'>'/gI \
    -e s/'# if'/'#if'/g -e s/'# endif'/'#endif'/g \
    -e s/'*'/'0x0'/g \
    ${FIN} > ${NFIN}





#####################################################
# 3/ Line by line...
#####################################################

nbl=$( cat ${NFIN} | wc -l)

echo
echo " *** ${NFIN} has ${nbl} lines ! ***"
echo

jl=0

while [ ${jl} -le ${nbl} ]; do

    jl=`expr ${jl} + 1`

    # Extract line # jl:
    line_orig=$( cat ${NFIN} | sed -n ${jl}p ) ; # extract line

    # Remove leading empty space characters of the line:
    line="${line_orig#"${line_orig%%[![:space:]]*}"}"

    # Replace any sequence of spaces/tabs/newlines with a single space:
    line="$(echo "$line" | tr -s '[:space:]' ' ')"

    # Replace ' ' with '|':
    line="$( echo ${line} | sed -e s/' '/'|'/g )"



    # Ignore lines starting with `!`:
    if [ "`echo ${line} | cut -c1-1`" != "!" ] && [ "${line}" != "" ]; then

        iOK=0
        ipass=1

        ct_old=${line}
        ct=${line}


        #echo ""
        #echo " ===> line :${ct}"

        while [ ${iOK} -ne 1 ]; do

            #if [ ${ipass} -eq 2 ]; then echo " SECOND PASS !!!"; fi

            #echo; echo
            #echo " Pass # ${ipass}, line => ${ct}"
            #echo

            idone=0 ; # the line has not been modified during current pass !



            # DO_1D
            # ~~~~~
            if [ ${i_do_1d} -eq 1 ]; then

                ca=`echo "${ct}"  | grep '^DO_1D'`
                if [ "${ca}" != "" ]; then
                    echo " 'DO_1D' found at line #${jl} ... (pass #${ipass})"
                    cb=`echo "${ca}" | sed -e s/'|'/''/g | cut -d'!' -f1`
                    c1=`echo "${cb}" | grep 'DO_1Di(0,0)'`
                    c2=`echo "${cb}" | grep 'DO_1Dj(0,0)'`
                    #echo "LOLO: cb => ${cb}"
                    #echo "LOLO: c1 => ${c1}"
                    #echo "LOLO: c2 => ${c2}"
                    #exit

                    if   [ "${c1}" != "" ]; then
                        echo "   ==> 'DO_1Di(0,0)'"
                        ct=`echo "${cb}" | sed -e s/'DO_1Di(0,0)'/'DO ji=Nis0, Nie0'/g`
                        #
                    elif [ "${c2}" != "" ]; then
                        echo "   ==> 'DO_1Dj(0,0)'"
                        ct=`echo "${cb}" | sed -e s/'DO_1Dj(0,0)'/'DO jj=Njs0, Nje0'/g`
                        #
                    else
                        echo "   ==> '${cb}' UNKNOWN!!!"
                        exit
                        # Need to find whats inside the parenthesis!
                        #echo "cb = |${cb}|"
                        cx=`echo ${cb} | cut -d'(' -f2 | cut -d')' -f1 | sed -e s/','/' '/g`
                        #echo "cx = |${cx}|"
                        vx=(${cx})
                        #echo ${vx[1]}
                        k1=${vx[0]}; k2=${vx[1]}
                        ct=`echo "${cb}" | sed -e s/"${cb}"/"DO jj=Njs0-${k3}, Nje0+${k4}\n   DO ji=Nis0-${k1}, Nie0+${k2}"/g`

                    fi
                    idone=1
                    #
                else
                    ca=`echo "${ct}"  | grep '^END_1D'`
                    if [ "${ca}" != "" ]; then
                        echo " 'END_1D' found at line #${jl} ... (pass #${ipass})"
                        cb=`echo "${ca}" | sed -e s/'|'/''/g`
                        ct=`echo "${cb}" | sed -e s/'END_1D'/'    END DO'/g`
                        idone=1
                    fi

                fi ;#if [ "${ca}" != "" ]


            fi ;#if [ ${i_do_1d} -eq 1 ]




            # DO_2D
            # ~~~~~
            if [ ${i_do_2d} -eq 1 ]; then

                ca=`echo "${ct}"  | grep '^DO_2D'`
                if [ "${ca}" != "" ]; then
                    echo " 'DO_2D' found at line #${jl} ... (pass #${ipass})"
                    cb=`echo "${ca}" | sed -e s/'|'/''/g | cut -d'!' -f1`

                    c1=`echo "${cb}" | grep 'DO_2D(0,0,0,0)'`
                    c2=`echo "${cb}" | grep 'DO_2D(nn_hls,nn_hls,nn_hls,nn_hls)'`
                    c3=`echo "${cb}" | grep 'DO_2D(ihls,ihls,ihls,ihls)'`
                    c4=`echo "${cb}" | grep 'DO_2D(1,1,1,1)'`
                    c5=`echo "${cb}" | grep 'DO_2D(ihls+1,ihls+1,ihls+1,ihls+1)'`
                    c6=`echo "${cb}" | grep 'DO_1Dj(0,0)'`
                    c7=`echo "${cb}" | grep 'DO_1Di(0,0)'`

                    if    [ "${c1}" != "" ]; then
                        echo "   ==> 'DO_2D(0,0,0,0)'"
                        ct=`echo "${cb}" | sed -e s/'DO_2D(0,0,0,0)'/'DO jj=Njs0, Nje0\n   DO ji=Nis0, Nie0'/g`
                        #
                    elif [ "${c2}" != "" ]; then
                        echo "   ==> 'DO_2D(nn_hls,nn_hls,nn_hls,nn_hls)'"
                        ct=`echo "${cb}" | sed -e s/'DO_2D(nn_hls,nn_hls,nn_hls,nn_hls)'/'DO jj=Njs0-nn_hls, Nje0+nn_hls\n   DO ji=Nis0-nn_hls, Nie0+nn_hls'/g`
                        #
                    elif [ "${c3}" != "" ]; then
                        echo "   ==> 'DO_2D(ihls,ihls,ihls,ihls)'"
                        ct=`echo "${cb}" | sed -e s/'DO_2D(ihls,ihls,ihls,ihls)'/'DO jj=Njs0-ihls, Nje0+ihls\n   DO ji=Nis0-ihls, Nie0+ihls'/g`
                        #
                    elif [ "${c4}" != "" ]; then
                        echo "   ==> 'DO_2D(1,1,1,1)'"
                        ct=`echo "${cb}" | sed -e s/'DO_2D(1,1,1,1)'/'DO jj=Njs0-1, Nje0+1\n   DO ji=Nis0-1, Nie0+1'/g`
                        #
                    elif [ "${c5}" != "" ]; then
                        echo "   ==> 'DO_2D(ihls+1,ihls+1,ihls+1,ihls+1)'"
                        ct=`echo "${cb}" | sed -e s/'DO_2D(ihls+1,ihls+1,ihls+1,ihls+1)'/'DO jj=Njs0-ihls+1, Nje0+ihls+1\n   DO ji=Nis0-ihls+1, Nie0+ihls+1'/g`
                        #
                    else
                        echo "   ==> '${cb}'"
                        # Need to find whats inside the parenthesis!
                        #echo "cb = |${cb}|"
                        cx=`echo ${cb} | cut -d'(' -f2 | cut -d')' -f1 | sed -e s/','/' '/g`
                        #echo "cx = |${cx}|"
                        vx=(${cx})
                        #echo ${vx[1]}
                        k1=${vx[0]}; k2=${vx[1]}; k3=${vx[2]}; k4=${vx[3]}
                        ct=`echo "${cb}" | sed -e s/"${cb}"/"DO jj=Njs0-${k3}, Nje0+${k4}\n   DO ji=Nis0-${k1}, Nie0+${k2}"/g`
                        #echo "${ct}"
                        #exit
                        #
                        #else
                        #    echo; echo "PROBLEM: unknow DO_2D pattern !!!"
                        #    exit
                    fi
                    idone=1
                    #
                else
                    ca=`echo "${ct}"  | grep '^END_2D'`
                    if [ "${ca}" != "" ]; then
                        echo " 'END_2D' found at line #${jl} ... (pass #${ipass})"
                        cb=`echo "${ca}" | sed -e s/'|'/''/g`
                        ct=`echo "${cb}" | sed -e s/'END_2D'/'    END DO\n END DO'/g`
                        idone=1
                    fi

                fi ;#if [ "${ca}" != "" ]


            fi ;#if [ ${i_do_2d} -eq 1 ]




            # DO_3D
            # ~~~~~
            if [ ${i_do_3d} -eq 1 ]; then

                ca=`echo "${ct}"  | grep '^DO_3D'`
                if [ "${ca}" != "" ]; then
                    cb=`echo "${ca}" | sed -e s/'|'/''/g | cut -d'!' -f1`
                    #
                    for CLT in "nlay_i" "nlay_s" "jpka" "ipk" "jpk" "jpkm1"; do
                        cc=`echo "${cb}" | grep ",${CLT})"`
                        if [ "${cc}" != "" ]; then
                            echo " 'DO_3D'/${CLT} found at line #${jl} ... (pass #${ipass})"

                            c1=`echo "${cb}" | grep "DO_3D(0,0,0,0,1,${CLT})"`
                            c2=`echo "${cb}" | grep "DO_3D(nn_hls,nn_hls,nn_hls,nn_hls,1,${CLT})"`
                            c3=`echo "${cb}" | grep "DO_3D(ihls,ihls,ihls,ihls,1,${CLT})"`
                            c4=`echo "${cb}" | grep "DO_3D(1,1,1,1,1,${CLT})"`
                            c5=`echo "${cb}" | grep "DO_3D(ihls+1,ihls+1,ihls+1,ihls+1,1,${CLT})"`

                            if    [ "${c1}" != "" ]; then
                                echo "   ==> DO_3D(0,0,0,0,1,${CLT})"
                                ct=`echo "${cb}" | sed -e s/"DO_3D(0,0,0,0,1,${CLT})"/"DO jj=Njs0, Nje0\n   DO ji=Nis0, Nie0\n      DO jk=1, ${CLT}"/g`
                                #
                            elif [ "${c2}" != "" ]; then
                                echo "   ==> DO_3D(nn_hls,nn_hls,nn_hls,nn_hls,1,${CLT})"
                                ct=`echo "${cb}" | sed -e s/"DO_3D(nn_hls,nn_hls,nn_hls,nn_hls,1,${CLT})"/"DO jj=Njs0-nn_hls, Nje0+nn_hls\n   DO ji=Nis0-nn_hls, Nie0+nn_hls\n      DO jk=1, ${CLT}"/g`
                                #
                            elif [ "${c3}" != "" ]; then
                                echo "   ==> DO_3D(ihls,ihls,ihls,ihls,1,${CLT})"
                                ct=`echo "${cb}" | sed -e s/"DO_3D(ihls,ihls,ihls,ihls,1,${CLT})"/"DO jj=Njs0-ihls, Nje0+ihls\n   DO ji=Nis0-ihls, Nie0+ihls\n      DO jk=1, ${CLT}"/g`
                                #
                            elif [ "${c4}" != "" ]; then
                                echo "   ==> DO_3D(1,1,1,1,1,${CLT})"
                                ct=`echo "${cb}" | sed -e s/"DO_3D(1,1,1,1,1,${CLT})"/"DO jj=Njs0-1, Nje0+1\n   DO ji=Nis0-1, Nie0+1\n      DO jk=1, ${CLT}"/g`
                                #
                            elif [ "${c5}" != "" ]; then
                                echo "   ==> DO_3D(ihls+1,ihls+1,ihls+1,ihls+1,1,${CLT})"
                                ct=`echo "${cb}" | sed -e s/"DO_3D(ihls+1,ihls+1,ihls+1,ihls+1,1,${CLT})"/"DO jj=Njs0-ihls+1, Nje0+ihls+1\n   DO ji=Nis0-ihls+1, Nie0+ihls+1\n      DO jk=1, ${CLT}"/g`
                                #
                            else
                                echo; echo "PROBLEM: unknow DO_3D pattern !!! ==> ${cb}"
                                exit
                            fi

                        fi
                    done
                    idone=1
                    #
                else
                    #
                    ca=`echo "${ct}"  | grep '^END_3D'`
                    if [ "${ca}" != "" ]; then
                        echo " 'END_3D' found at line #${jl} ... (pass #${ipass})"
                        cb=`echo "${ca}" | sed -e s/'|'/''/g`
                        ct=`echo "${cb}" | sed -e s/'END_3D'/'       END DO\n    END DO\n END DO'/g`
                        idone=1
                    fi

                fi ;#if [ "${ca}" != "" ]

            fi ;#if [ ${i_do_3d} -eq 1 ]




            # A1D*(*)
            # ~~~~~
            if [ ${i_do_a1d} -eq 1 ]; then

                ca=`echo "${ct}" | grep 'A1Di(' | cut -d'!' -f1`
                if [ "${ca}" != "" ]; then
                    # What's inside the parenthesis:
                    #cin=$( echo "${ca}" | grep -oP 'A1Di\(\K[^)]*\)[^)]*(?=\))' ) ; # 2026 :D
                    #cin=$( extract_AXD ${ca} A1Di )
                    cin=$( extract_func_args ${ca} A1Di )
                    echo " A1Di(${cin}) found at line #${jl} ... (pass #${ipass})"
                    #
                    cb4=`echo ${ca} | sed -e s/'|'/''/g | cut -c1-4`
                    #
                    #if [ "${cb4}" = "ALLO" ] || [ "${cb4}" = "LOGI" ] || [ "${cb4}" = "REAL" ] || [ "${cb4}" = "INTE" ] || [ "${cb4}" = "CHAR" ]; then
                    #    ct=`echo "${ca}" | sed -e s/"A1Di(${cin})"/"jpi"/g`
                    #else
                    if   [ "${cin}" = "0" ]; then
                        ct=`echo "${ca}" | sed -e s/"A1Di(${cin})"/"Nis0:Nie0"/g`
                        echo ; echo " *** WARNING: A1Di(${cin}) to be changed to Nis0:Nie0 ***"; echo
                    else
                        ct=`echo "${ca}" | sed -e s/"A1Di(${cin})"/"Nis0-${cin}:Nie0+${cin}"/g`
                        echo ; echo " *** WARNING: A1Di(${cin}) to be changed to Nis0-${cin}:Nie0+${cin} ***"; echo
                    fi
                    idone=1
                    #fi
                fi ;#if [ "${ca}" != "" ]


                ca=`echo "${ct}" | grep 'A1Dj(' | cut -d'!' -f1`
                if [ "${ca}" != "" ]; then
                    # What's inside the parenthesis:
                    #cin=$( echo "${ca}" | grep -oP 'A1Dj\(\K[^)]*\)[^)]*(?=\))' ) ; # 2026 :D
                    #cin=$( extract_AXD ${ca} A1Dj )
                    cin=$( extract_func_args ${ca} A1Dj )
                    echo " A1Dj(${cin}) found at line #${jl} ... (pass #${ipass})"
                    #
                    cb4=`echo ${ca} | sed -e s/'|'/''/g | cut -c1-4`
                    #
                    if [ "${cb4}" = "ALLO" ] || [ "${cb4}" = "LOGI" ] || [ "${cb4}" = "REAL" ] || [ "${cb4}" = "INTE" ] || [ "${cb4}" = "CHAR" ]; then
                        ct=`echo "${ca}" | sed -e s/"A1Dj(${cin})"/"jpj"/g`
                    else
                        if   [ "${cin}" = "0" ]; then
                            ct=`echo "${ca}" | sed -e s/"A1Dj(${cin})"/"Njs0:Nje0"/g`
                            echo ; echo " *** WARNING: A1Dj(${cin}) to be changed to Njs0:Nje0 ***"; echo
                        else
                            ct=`echo "${ca}" | sed -e s/"A1Dj(${cin})"/"Njs0-${cin}:Nje0+${cin}"/g`
                            echo ; echo " *** WARNING: A1Dj(${cin}) to be changed to Njs0-${cin}:Nje0+${cin} ***"; echo
                        fi
                    fi
                    idone=1
                fi ;#if [ "${ca}" != "" ]

            fi ;#if [ ${i_do_2d} -eq 1 ]




            # A2D()
            # ~~~~~
            if [ ${i_do_a2d} -eq 1 ]; then

                ca=`echo "${ct}" | grep 'A2D('`
                if [ "${ca}" != "" ]; then
                    # What's inside the parenthesis:
                    cin=$( extract_func_args ${ca} A2D )
                    echo " A2D(${cin}) found at line #${jl} ... (pass #${ipass}) [ ${ct} ] "
                    cb4=`echo ${ca} | sed -e s/'|'/''/g | cut -c1-4`
                    #
                    if   [ "${cin}" = "0" ]; then
                        ct=`echo "${ca}" | sed -e s/"A2D(${cin})"/"Nis0:Nie0,Njs0:Nje0"/g`
                        echo ; echo " *** WARNING: A2D(${cin}) to be changed to Nis0:Nie0,Njs0:Nje0 ***"; echo
                    else
                        ct=`echo "${ca}" | sed -e s/"A2D(${cin})"/"Nis0-${cin}:Nie0+${cin},Njs0-${cin}:Nje0+${cin}"/g`
                        echo ; echo " *** WARNING: A2D(${cin}) to be changed to Nis0-${cin}:Nie0+${cin},Njs0-${cin}:Nje0+${cin} ***"; echo
                    fi
                    #fi
                    idone=1

                fi ;#if [ "${ca}" != "" ]

            fi ;#if [ ${i_do_2d} -eq 1 ]


            # T2D()
            # ~~~~~
            if [ ${i_do_t2d} -eq 1 ]; then

                ca=`echo "${ct}" | grep 'T2D('`
                if [ "${ca}" != "" ]; then
                    # What's inside the parenthesis:
                    cin=$( extract_func_args ${ca} T2D )
                    echo " T2D(${cin}) found at line #${jl} ... (pass #${ipass}) [ ${ct} ] "
                    cb4=`echo ${ca} | sed -e s/'|'/''/g | cut -c1-4`
                    #
                    if   [ "${cin}" = "0" ]; then
                        ct=`echo "${ca}" | sed -e s/"T2D(${cin})"/"ntsi:ntei,ntsj:ntej"/g`
                        echo ; echo " *** WARNING: T2D(${cin}) to be changed to ntsi:ntei,ntsj:ntej ***"; echo
                    else
                        ct=`echo "${ca}" | sed -e s/"T2D(${cin})"/"ntsi-${cin}:ntei+${cin},ntsj-${cin}:ntej+${cin}"/g`
                        echo ; echo " *** WARNING: T2D(${cin}) to be changed to ntsi-${cin}:ntei+${cin},ntsj-${cin}:ntej+${cin} ***"; echo
                    fi
                    #
                    idone=1

                fi ;#if [ "${ca}" != "" ]

            fi ;#if [ ${i_do_2d} -eq 1 ]










            # Go to next pass (while loop) if something has been changed...
            if [ ${idone} -ne 1 ]; then

                # Semicolon to line break
                # ~~~~~~~~~~~~~~~~~~~~~~~
                #  ==> must be in last position because creates 2 lines !!!!

                if [ ${i_do_sc} -eq 1 ]; then

                    # Stupid `IF` statement starting after something and a `;` :
                    ca=`echo "${ct}"  | grep -e ';|IF' -e ';IF'`
                    if [ "${ca}" != "" ]; then
                        echo "Stupid 'IF-statement' starting after a ';'' found at line #${jl} ... (pass #${ipass}) [ ${ct} ]"
                        ct=`echo ${ca} | sed -e s/';'/'\n'/g`
                        idone=1
                    fi

                    ca=`echo "${ct}"  | grep '^IF' | grep -e '|THEN' | grep -e '\;'`
                    if [ "${ca}" != "" ]; then
                        echo " 'THEN ;' found at line #${jl}... (pass #${ipass}) [ ${ct} ]"
                        #cb=`echo "${ca}" | sed -e s/'|'/''/g | grep 'THEN;'`
                        ct=`echo ${ca} | sed -e s/';'/'\n'/g`
                        idone=1
                    fi

                    ca=`echo "${ct}"  | grep '^ELSE' | grep -e '\;'`
                    if [ "${ca}" != "" ]; then
                        echo " 'ELSE ;' found at line #${jl}... (pass #${ipass}) [ ${ct} ]"
                        #cb=`echo "${ca}" | sed -e s/'|'/''/g | grep 'THEN;'`
                        ct=`echo ${ca} | sed -e s/';'/'\n'/g`
                        idone=1
                    fi

                    ca=`echo "${ct}"  | grep '^WHERE' | grep -e '\;'`
                    if [ "${ca}" != "" ]; then
                        echo " 'WHERE ;' found at line #${jl}... (pass #${ipass}) [ ${ct} ]"
                        #cb=`echo "${ca}" | sed -e s/'|'/''/g | grep 'THEN;'`
                        ct=`echo ${ca} | sed -e s/';'/'\n'/g`
                        idone=1
                    fi


                fi


            fi ; #if [ ${idone} -ne 1 ]



            #echo "Line after pass ${ipass} => ${ct}"



            if [ "${ct}" != "${ct_old}" ] && [ ${ipass} -eq 1 ]; then #
                # A modif has been don
                iOK=0
                ipass=2
            else
                iOK=1
            fi


        done


        if [ "${ct}" != "${ct_old}" ]; then
            # Line has been modified, use the new line:
            echo "${ct}" | sed -e s/'|'/' '/g   -e s/'0x0'/'*'/g >> ../${OUTFILE}
        else
            # Nothing has been modified, keep the original line:
            echo "${line_orig}" | sed -e s/'|'/' '/g   -e s/'0x0'/'*'/g >> ../${OUTFILE}
        fi

    else
        # Line that has been ignored for the search (comment or blank line):
        echo "${line_orig}" | sed -e s/'|'/' '/g   -e s/'0x0'/'*'/g >> ../${OUTFILE}


        #
    fi #if [ "`echo ${line} | cut -c1-1`" != "!" ] && [ "${line}" != "" ]



done


cd ../


if [ ${i_pipe_fnd} -eq 1 ]; then
    rm -f ${OUTFILE}.tmp
    sed -e "s|\`/\`|\||g" ${OUTFILE} > ${OUTFILE}.tmp
    mv -f ${OUTFILE}.tmp ${OUTFILE}
fi


indent_clean ${OUTFILE}


echo; echo
echo " Adapated file '${OUTFILE}' generated !!!"
echo
