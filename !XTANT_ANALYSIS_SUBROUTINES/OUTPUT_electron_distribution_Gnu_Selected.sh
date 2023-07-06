#!/bin/bash

NAME=OUTPUT_electron_distribution_Selected.gif
LW=1.0
LABL="Distribution"
TICSIZ=5.00
echo " 
set terminal gif font \"arial,14\" 
set output \"$NAME\"
set xlabel \"Energy (eV) \" 
set ylabel \"Electron distribution (a.u.) \" 
set key right top 
set xtics \"$TICSIZ\" 
stats \"OUTPUT_electron_distribution.dat\" nooutput
#do for [i=1:int(STATS_blocks)] {
p [:25.0000000000000000][0:2]  \"OUTPUT_electron_distribution.dat\" index 50 u 1:3 w l lw 2 lt rgb \"green\" notitle  ,\
 \"OUTPUT_electron_distribution.dat\" index 50 u 1:2 pt 7 ps 1 title \"-150 fs\" ,\
 \"OUTPUT_electron_distribution.dat\" index 200 u 1:3 w l lw 2 lt rgb \"yellow\" notitle  ,\
 \"OUTPUT_electron_distribution.dat\" index 200 u 1:2 pt 6 ps 1 title \"0 fs\"  ,\
 \"OUTPUT_electron_distribution.dat\" index 400 u 1:3 w l lw 2 lt rgb \"blue\" notitle  ,\
 \"OUTPUT_electron_distribution.dat\" index 400 u 1:2 pt 5 ps 0.5 title \"200 fs\" ,\
 \"OUTPUT_electron_distribution.dat\" index 1200 u 1:3 w l lw 2 lt rgb \"grey\" notitle  ,\
 \"OUTPUT_electron_distribution.dat\" index 1200 u 1:2 pt 4 ps 0.5 title \"1000 fs\"
#}
reset
" | gnuplot 


