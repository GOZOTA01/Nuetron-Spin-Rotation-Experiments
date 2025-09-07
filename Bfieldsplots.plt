#
# Set overall margins for the combined set of plots and size them
# to generate a requested inter-plot spacing
#
if (!exists("MP_LEFT"))   MP_LEFT = .1
if (!exists("MP_RIGHT"))  MP_RIGHT = .95
if (!exists("MP_BOTTOM")) MP_BOTTOM = .14
if (!exists("MP_TOP"))    MP_TOP = .9
if (!exists("MP_xGAP"))   MP_xGAP = 0.05
if (!exists("MP_yGAP"))   MP_yGAP = 0.1

set multiplot layout 2,2 columnsfirst title "{/:Bold=15 Multiplot with explicit page margins}" \
              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_xGAP, MP_yGAP

set format y "%.1f"
set key box opaque
set ylabel 'Bfield(G)'
set xrange [0:200]

Bx="Bfieldsx.txt"
plot Bx lt 4


set xlabel 'Position(cm)'
By="Bfieldsy.txt"
plot By lt 3

unset ylabel
unset xlabel

set xrange [-5:5]
set yrange [-5:5]
quadrant="quadrants.txt"
plot quadrant 

unset yrange

set xlabel 'Position(cm)'
set xrange [0:200]
Bz="Bfieldsz.txt"
plot Bz lt 2






unset multiplot
