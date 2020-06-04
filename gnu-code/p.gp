reset
set term x11 1
se pm3d map

set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"
#set colorbox horiz user origin .25,.07 size .65,.04
set palette defined (0.34"blue",0.37"white",0.4"red")
set cbrange [0.34:0.4]


set xl "x"
set yl "y"
set zl "parameter"
set grid
se contour#等高線の表示
unset clabel
se cntrparam levels incremental 0.365,0.0022,0.374#等高線の初期値、間隔、最終値の設定

set xrange [0:24]
set yrange [-8:8]
#set zrange [0:1.1]

set term x11
splot '/Users/isseyshome/Documents/Lab/2dimensions/result_super/parameter030000.d' u 1:2:5 w pm3d t ""
se term jpeg
se out "/Users/isseyshome/Documents/Lab/2dimensions/p030000.jpg"
rep
#set term x11
#splot '/Users/isseyshome/Documents/Lab/ex_20/result_file/p01000.d' w pm3d
