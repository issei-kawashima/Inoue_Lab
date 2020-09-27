reset
set term x11
se pm3d map
#se size square
#se view 90,270,1,1
se contour#等高線の表示
se cntrparam levels incremental -0.2,0.03,0.3#等高線の初期値、間隔、最終値の設定
unset clabel#等高線の凡例の消去
set colorbox horiz user origin .25,.07 size .65,.04
set cbrange [-1:1]
set palette defined (-0.5"blue",0"white",0.5"red")

set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"

set xl "x"
set yl "y"
set zl "2nd Invaritant"
set grid

set xrange [0:24]
set yrange [-8:8]
#set zrange [0:1]

set term x11 1
#3次元　通常 (2020/08/26 Invariant_2第二不変量を投入する前のデータのみ)
splot '/Users/isseyshome/Downloads/result_omp_1.6/parameter015000_10.txt' u 1:2:6  w pm3d t ""

#2次元
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/2dimensions/result_super_2d-3Ma1.6/parameter001500.txt' u 1:2:4  w pm3d t ""

se term jpeg
se out "/Users/isseyshome/Documents/Lab/#3 修士論文/Visual_result/Ma1.6_random/Q_15000_10.jpg"
rep
