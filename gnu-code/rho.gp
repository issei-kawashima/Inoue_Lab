reset
set term x11
#se pm3d map
#se size square
#se view 56,333,1,1
#se view 90,270,1,1

#se contour  surface#等高線の表示
#se cntrparam levels incremental 0.7,0.02,1#等高線の初期値、間隔、最終値の設定
#se cntrparam levels incremental 0.5,0.05,1.5#等高線の初期値、間隔、最終値の設定
#unset clabel#等高線の凡例の消去

set colorbox horiz user origin .25,.03 size .65,.04
set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"


#set palette rgbformulae 22,13,-31	#!コントラストが強い色
#set palette defined (0.5"blue",1"white",1.5"red")
#set cbrange [0:1.5]



set xl "x"
set yl "y"
set zl "rho"
set grid

#set xrange [0:24]
#set yrange [-8:8]
set xrange [0:36]
set yrange [-10:10]
#set zrange [0:0.2]

#set dgrid3d
set term x11 2
splot '/Users/isseyshome/Downloads/result_square/parameter000041_44.txt'u 1:2:4 w pm3d t ""

#se term jpeg
#unset colorbox
#se out "/Users/isseyshome/Documents/Lab/#3 修士論文/Visual_result/Ma1.4_random_10%/rho015356.jpg"
#rep
