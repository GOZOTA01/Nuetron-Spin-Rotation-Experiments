# plot.plt
set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*"
set title "OUTCOIL-thetax"
set nokey
set grid
set xlabel "position(cm)"
set ylabel "angle"
m="mydata.txt"
plot m using 1:2 with linespoints
