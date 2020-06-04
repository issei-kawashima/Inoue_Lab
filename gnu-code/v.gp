reset
set term x11
#se pm3d map
se view 90,270,1,1
set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"

set xl "x"
set yl "y"
set zl "v"
set grid

#set xrange [0:1]
set yrange [-1:1]
#set zrange [-1.5:1.5]

set term x11 1
splot '/Users/isseyshome/Documents/Lab/thesis/result_DT/v00200.d' w pm3d
#se term jpeg
#se out "/Users/isseyshome/Documents/Lab/thesis/result_DT/v00000.jpg"
#rep
#set term x11
#splot '/Users/isseyshome/Documents/Lab/thesis/result_DT/v01000.d' w pm3d
