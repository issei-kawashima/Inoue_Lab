reset
set term x11
#se pm3d map
#se view 0,90,1,1
set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"
#set colorbox horiz user origin .25,.07 size .65,.04
set palette defined (-0.5"blue",0.5"white",1.5"red")
set cbrange [-0.5:1.5]

set xl "x"
set yl "y"
set zl "u"
set grid

set xrange [0:24]
set yrange [-8:8]
#set zrange [-1.5:1.5]
se contour#等高線の表示
unset clabel
se cntrparam levels incremental -0.2,0.2,1.2#等高線の初期値、間隔、最終値の設定

set term x11 1
splot '/Users/isseyshome/Documents/Lab/2dimensions/result_super/parameter000000.d' u 1:2:6 w pm3d t ""

set term x11 2
splot '/Users/isseyshome/Documents/Lab/2dimensions/result_super/parameter000100.d' u 1:2:6 w pm3d t ""

set term x11 3
splot '/Users/isseyshome/Documents/Lab/2dimensions/result_super/parameter000500.d' u 1:2:6 w pm3d t ""

set term x11 4
splot '/Users/isseyshome/Documents/Lab/2dimensions/result_super/parameter000600.d' u 1:2:6 w pm3d t ""

#se term jpeg
#se out "/Users/isseyshome/Documents/Lab/2dimensions/u000600.jpg"
#rep
