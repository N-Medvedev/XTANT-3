#!/bin/bash

NAME=OUTPUT_diffraction_peaks_Ta_COMPARE.png
LW=3.0
M=3.34641e-26
T_D=745
q=2.95229
LABL="Diffraction peak"
TICSIZ=100.00
echo " 
set terminal pngcairo dashed font \"arial,14\" 
set output \"$NAME\"
set xlabel \"Time (ps)\" font \"arial,18\" 
set ylabel \"Atomic tmperature (K)\" font \"arial,18\" 
set key left top
set xtics \"$TICSIZ\" 
p [  -500.0:500][0:] \"OUTPUT_temperatures.dat\"u (\$1/1):3 w l lw \"$LW\" title \"XTANT-3 actual T_a\" ,\
\"OUTPUT_diffraction_peaks.dat\" u (\$1/1):(300-log(\$2) * $M * 1.38064852e-23 * $T_D**2 / (2*3*(1.05457162853e-34)**2 * ($q*1e10)**2) ) w l lw \"$LW\" dashtype \"_\" title \"DW analysis of (200)\"
" | gnuplot 
