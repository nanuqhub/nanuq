#!/bin/bash
set -u
#
vname=timing_step
doplot=1
dopng=0
t1st=1    ; tend=0
vminup=-1 ; vmaxup=-1
vmindn=-1 ; vmaxdn=-1
dinput="not defined"
finput="not defined"
rnk=0
#
while [ $# -gt 0 ]
do
    case $( echo $1 | tr '[:upper:]' '[:lower:]' ) in
     -h|--help)
         echo "   Purpose:"
         echo "       Print statistics and plot of NEMO stats time series in micro-second"
         echo "   Usage:"
         echo "       ./timing_nanuq_gnuplot.sh"
         echo "   Options:"
         echo "       --png               # save plots in png"
         echo "       --noplot            # print only the statistics (no plot)"
         echo "       --vname vname       # NetCDF variable name (timing_step by default)"
         echo "       --dir dir_name      # dir_name: directory where is located (./ by default)"
         echo "       --finput file_name  # file_name: NetCDF input file (timing_step.nc by default)"
         echo "       --rnk               # MPI rank who wrote the NetCDF input file (0 by default)"
         echo "       --t1st time1        # time1: first time step to select ($t1st by default)"
         echo "       --tend time2        # time2: last time step to select (last by default)"
         echo "       --rgup vmin vmax    # [vmin, vmax]: upper plot range ([0, max] by default)"
         echo "       --rgdn vmin vmax    # [vmin, vmax]: lower plot range"
         exit 0 ;;
     --png)     dopng=1 ;;
     --noplot)  doplot=0 ;;
     --vname)   vname=${2}  ; shift ;;
     --dir)     dinput=${2} ; shift ;;
     --finput)  finput=${2} ; shift ;;
     --rnk)     rnk=${2}    ; shift ;;
     --t1st)    t1st=${2}   ; shift ;;
     --tend)    tend=${2}   ; shift ;;
     --rgup)    vminup=${2} ; shift ; vmaxup=${2} ; shift ;;
     --rgdn)    vmindn=${2} ; shift ; vmaxdn=${2} ; shift ;;
    esac
    shift
done
#
[ "$finput" = "not defined" ] && finput=${vname}.nc
if [ "$dinput" != "not defined" ] && finput=$dinput/$finput
then
   [ ! -d $dinput ] && echo "ERROR: $dinput is not a directory" && exit 1
   finput=$dinput/$finput
fi
if [ $rnk -gt 0 ]
then
   finput=$( ls -1 ${finput/.nc/_*.nc} | grep -e "_0*${rnk}.nc$" )
   rnk=$( echo $finput | sed -e "s/.*_/_/" -e "s/\.nc//" )
else
   rnk=""
fi
[ ! -f $finput ] && echo "ERROR: $finput not found" && exit 1
#
ok=$( ncdump -h $finput | grep -c "double *${vname}(" )
[ $ok -eq 0 ] && echo "ERROR: $vname not found in $finput" && exit 1
#
[ $tend -eq 0 ] && tend=$( ncdump -h $finput | grep UNLIMITED | sed -e "s/[^0-9]*//g" )
#
which gnuplot &> /dev/null
[ $? -ne 0 ] && echo "$( basename $0 ) uses gnuplot which was not found" && exit 1
#
ncdump -f f -v $vname $finput | sed -n -e "/\/\/ *${vname}(${t1st})/,/\/\/ *${vname}(${tend})/p" | sed -e "s/[;,] .*//" -e "s/.* /scale=6 ; 1000000. * /" | bc > timing_gnuplot.$$
#
gnuplot -persist << EOF
dopng = $dopng
doplot = $doplot
vminup = $vminup ; vmaxup = $vmaxup
vmindn = $vmindn ; vmaxdn = $vmaxdn
if ( dopng == 0 ) {
   set terminal x11 size 800,800
} else {
   set terminal png size 800,800
   set output "${finput/${rnk}.nc/_t${t1st}_t${tend}${rnk}.png}"
}
stats "timing_gnuplot.$$" name "ST"
if ( doplot == 1 ) {
   mn = ST_mean ; md = ST_median ; std = ST_stddev
   if ( vminup == -1 ) { ; vminup = ST_min ; }
   if ( vmaxup == -1 ) { ; vmaxup = ST_max ; }
   lo = ST_lo_quartile ; up = ST_up_quartile ; iqr = up-lo
   set xrange [0:ST_records]
   set xlabel "${vname//_/ } number"
   set ylabel " ${vname//_/ } elapsed time (microsecond)"
   set multiplot layout 2,1
   set yrange [vminup:vmaxup]
   set title sprintf("FULL RANGE: ${vname//_/ } (microsecond), mean = %f, median = %f", mn, md)
   plot "timing_gnuplot.$$" notitle, mn title "Mean" lw 2, md title "Median" lw 2, mn-std title "Mean-StdDev" lw 2, mn+std title "Mean+StdDev" lw 2
   if ( vmindn == -1 ) {
      vmindn = lo - 1.5*iqr
      if ( vmindn < 0      ) { ; vmindn = 0 ; }
   }
   if ( vmaxdn == -1 ) {
      vmaxdn = up + 1.5*iqr
      if ( vmaxdn > ST_max ) { ; vmaxdn = ST_max ; }
   }
   set yrange [vmindn:vmaxdn]
   set title sprintf("ZOOM: ${vname//_/ } (microsecond), mean = %f, median = %f", mn, md)
   plot "timing_gnuplot.$$" notitle, mn title "Mean" lw 2, md title "Median" lw 2, lo title "1st Quatile" lw 2, up title "3rd Quartile" lw 2
}
EOF
#
rm -f timing_gnuplot.$$ &>/dev/null
