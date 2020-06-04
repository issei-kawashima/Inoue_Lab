reset
set term x11 1
se pm3d map
se contour#等高線の表示
unset clabel#等高線の凡例の消去
se cntrparam levels incremental -0.001,0.001,0.0005#等高線の初期値、間隔、最終値の設定
set colorbox horiz user origin .25,.07 size .65,.04
set palette defined (-1"blue",0"white",1"red")
set cbrange [-1:1]
#se view 90,270,1,1
set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"
set title "dp/dt"
set title font "Arial,12"

set xl "x"
set yl "y"
set zl "dp/dt"
set grid

set xrange [0:24]
set yrange [-8:8]
#set zrange [-1.5:1.5]

set term x11
splot '/Users/isseyshome/Documents/Lab/2dimensions/result_super/parameter002000.d' u 1:2:5 w pm3d t ""
#se term jpeg
#se out "/Users/isseyshome/Documents/Lab/thesis/result_DT/dp000000.jpg"
#rep
#set term x11
#splot '/Users/isseyshome/Documents/Lab/thesis/result_DT/dp010000.d' w pm3d
