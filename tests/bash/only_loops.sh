#!/bin/bash

i_keep_hls=1


if [ "${1}" = "" ]; then
    echo " $0 <file.F90>"
    exit
fi

FIN="${1}"

# 1st, backing up file:
if [ -f ./${FIN}.orig ]; then
    echo "PROBLEM: there is already a ${FIN}.orig in here!"; exit
fi

rsync -avP ${FIN} ${FIN}.orig


cscrpt=do.tmp


cat > ${cscrpt} <<EOF
#!/bin/bash

cat ${FIN} | sed \\
EOF




### DO_2D substitutions:
#
# * DO_2D( nn_hls, nn_hls, nn_hls, nn_hls ) => DO jj = ntsj-( nn_hls), ntej+( nn_hls ) ; DO ji = ntsi-( nn_hls), ntei+( nn_hls)
#      ==> DO jj = ntsj-nn_hls, ntej+nn_hls
#             DO ji = ntsi-nn_hls, ntei+nn_hls
#
# * DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 ) => DO jj = ntsj-( nn_hls-1), ntej+( nn_hls-1 ) ; DO ji = ntsi-( nn_hls-1), ntei+( nn_hls-1)
#      ==> DO jj = ntsj-nn_hls+1, ntej+nn_hls-1
#             DO ji = ntsi-nn_hls+1, ntei+nn_hls-1
#
# * DO_2D( 1, 1, 1, 1 ) => DO jj = ntsj-( 1), ntej+( 1 ) ; DO ji = ntsi-( 1), ntei+( 1)
#
# * DO_2D( 1, 0, 1, 0 ) => DO jj = ntsj-( 1), ntej+( 0 ) ; DO ji = ntsi-( 1), ntei+( 0)
#
# * DO_2D( 0, 1, 0, 0 ) => DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 1)    ! CONFIRMS that in `DO_2D` the 2 arg 0 & 1 are for ji and 2 & 3 for jj !!

vva=`cat ${FIN} | grep 'DO_2D(' | sed -e s/' '/'.'/g`
vva=`echo ${vva} | cut -d '!' -f1`
echo ${vva}


noc=`echo ${vva} | wc -w`
echo "  ==> ${noc} occurences!"
echo

vva=( ${vva} )

k=0
cc=""
while [ ${k} -lt ${noc} ]; do

    str=${vva[${k}]}
    echo ${str}
    pspace=`echo ${str}| cut -d'D' -f1`
    echo " ==> ${pspace}"
    ppar="DO_2D`echo ${str}| cut -d'D' -f3`"
    echo " ppar ==> ${ppar}"
    csp=`echo ${pspace} | sed -e s/'.'/' '/g`
    echo " * preceeding spaces ==>x${csp}x"
    pnums=`echo ${str}| cut -d'(' -f2 | sed -e s/'\.'/''/g  -e s/')'/''/g`
    echo " * pnums ==> ${pnums}"
    vnums=(`echo ${pnums} | sed -e s/','/' '/g`)
    echo " * vnums ==> ${vnums[*]}"

    # Important:
    # ntsi=2 ; ntei=jpi-1
    # ntsj=2 ; ntej=jpj-1

    # *** ji ***
    case "${vnums[0]}" in
        #
        "0")        i1="Nis0";; #"ntsi";;
        #
        "1")        i1="Nis0-1";; #"ntsi-1";;
        #
        "nn_hls-1") i1="Nis0";; #""ntsi-nn_hls+1";;
        #
        "nn_hls")   i1="Nis0-1";; #"ntsi-nn_hls";;
        #
        *)          i1="Nis0-(${vnums[0]})" #"ntsi-${vnums[0]}"
                    echo "WARNING: i1-substituting ${vnums[0]} with ${i1}!"
                    ;;
    esac

    case "${vnums[1]}" in
        #
        "0")        i2="Nie0";; #"ntei";;
        #
        "1")        i2="Nie0+1";; #"ntei+1";;
        #
        "nn_hls-1") i2="Nie0";; #"ntei+nn_hls-1";;
        #
        "nn_hls")   i2="Nie0+1";; #"ntei+nn_hls";;
        #
        *)          i2="Nie0+(${vnums[1]})"  #"ntei+${vnums[1]}"
                    echo "WARNING: i2-substituting ${vnums[1]} with ${i2}!"
                    ;;
    esac

    # *** jj ***
    case "${vnums[2]}" in
        #
        "0")        j1="Njs0";; #"ntsj";;                              2
        #
        "1")        j1="Njs0-1";; #"ntsj-1";;                            1
        #
        "nn_hls-1") j1="Njs0";; #"ntsj-nn_hls+1";;            2
        #
        "nn_hls")   j1="Njs0-1";; #"ntsj-nn_hls";;                1
        #
        *)          j1="Njs0-(${vnums[2]})" #j1="ntsj-${vnums[2]}"
                    echo "WARNING: j1-substituting ${vnums[2]} with ${j1}!"
                    ;;
    esac

    case "${vnums[3]}" in
        #
        "0")        j2="Nje0";; #"ntej";;            
        #
        "1")        j2="Nje0+1";; #"ntej+1";; 
        #
        "nn_hls-1") j2="Nje0";; #"ntej+nn_hls-1";;
        #
        "nn_hls")   j2="Nje0+1";; #"ntej+nn_hls";; 
        #
        *)          j2="Nje0+(${vnums[3]})" #"ntej+${vnums[3]}"
                    echo "WARNING: j2-substituting ${vnums[3]} with ${j2}!"
                    ;;
    esac

    #exit;#lolo

    #if [ "${vnums[0]}" = "1" ]; then i1="1";   fi
    #if [ "${vnums[1]}" = "1" ]; then i2="jpi"; fi
    #if [ "${vnums[2]}" = "1" ]; then j1="1";   fi
    #if [ "${vnums[3]}" = "1" ]; then j2="jpj"; fi


    #if [ ${i_keep_hls} -eq 1 ]; then

    #if [ "${vnums[0]}" = "nn_hls" ]; then i1="ntsi-nn_hls"; fi
    #if [ "${vnums[1]}" = "nn_hls" ]; then i2="ntei+nn_hls"; fi
    #if [ "${vnums[2]}" = "nn_hls" ]; then j1="ntsj-nn_hls"; fi
    #if [ "${vnums[3]}" = "nn_hls" ]; then j2="ntej+nn_hls"; fi

    #if [ "${vnums[0]}" = "nn_hls-1" ]; then i1="ntsi-nn_hls+1"; fi
    #if [ "${vnums[1]}" = "nn_hls-1" ]; then i2="ntei+nn_hls-1"; fi
    #if [ "${vnums[2]}" = "nn_hls-1" ]; then j1="ntsj-nn_hls+1"; fi
    #if [ "${vnums[3]}" = "nn_hls-1" ]; then j2="ntej+nn_hls-1"; fi

    #else

    #if [ "${vnums[0]}" = "nn_hls" ]; then i1="1";   fi
    #if [ "${vnums[1]}" = "nn_hls" ]; then i2="jpi"; fi
    #if [ "${vnums[2]}" = "nn_hls" ]; then j1="1";   fi
    #if [ "${vnums[3]}" = "nn_hls" ]; then j2="jpj"; fi

    #if [ "${vnums[0]}" = "nn_hls-1" ]; then i1="2";     fi
    #if [ "${vnums[1]}" = "nn_hls-1" ]; then i2="jpi-1"; fi
    #if [ "${vnums[2]}" = "nn_hls-1" ]; then j1="2";     fi
    #if [ "${vnums[3]}" = "nn_hls-1" ]; then j2="jpj-1"; fi

    #fi

    #exit


    cat >> ${cscrpt} <<EOF
    -e s/"${csp}${ppar}"/"${csp}DO jj=${j1}, ${j2}\n${csp}   DO ji=${i1}, ${i2}"/g \\
    -e s/"${csp}END_2D"/"   ${csp}END DO\n${csp}END DO"/g \\
EOF


    k=`expr ${k} + 1`
    echo
done

cat >> ${cscrpt} <<EOF
        > OUT.F90
EOF

exit
chmod +x ${cscrpt}


./${cscrpt}

#exit

indent_clean OUT.F90

rm OUT.F90~ ${cscrpt}

mv OUT.F90 ${FIN}
