reset
set term x11
se pm3d map
#se size square
#se view 56,333,1,1
#se view 90,270,1,1

#se contour  surface#等高線の表示
#se cntrparam levels incremental 0.7,0.02,1#等高線の初期値、間隔、最終値の設定
#se cntrparam levels incremental 0.5,0.04,1.2#等高線の初期値、間隔、最終値の設定
#unset clabel#等高線の凡例の消去

set colorbox horiz user origin .25,.03 size .65,.04
#set cbrange [0:2]
set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"


#set palette rgbformulae 22,13,-31	#!コントラストが強い色
#set palette defined (-0.6"blue",0"white",0.6."red")
#set cbrange [0.4:1.2]


set xl "x"
set yl "y"
set zl "rho"
set grid

set xrange [0:24]
set yrange [-8:8]
#set xrange [-18:18]
#set yrange [-18:18]#表示範囲を正方形にしたいため?そうすればx,y方向にきれいに広がっているかわかるはず
#set zrange [0:0.2]

#set dgrid3d
set term x11 1
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_3D/parameter_initial000000.d'u 1:2:4 w pm3d t ""
splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_3D/parameter000000_00.txt'u 1:2:4 w pm3d t ""
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/comparison_result.d'u 1:2:3 w pm3d t ""

#set term x11 2
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_3D/parameter000005.d'u 1:2:4 w pm3d t ""
#splot '/Users/isseyshome/Documents/Lab/2dimensions/result_super/parameter000400.d'u 1:2:3 w pm3d t ""

#set term x11 3
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_3D/parameter00009.d'u 1:2:4 w pm3d t ""

#set term x11 4
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_3D/parameter000019.d'u 1:2:4 w pm3d t ""

#se term jpeg
#unset colorbox
#se out "/Users/isseyshome/Documents/Lab/2dimensions/Ma=2.4/rho008000.jpg"
#rep

#set term x11 1
#splot '/Users/isseyshome/Documents/Lab/thesis/result_DT/rho05000.d' w pm3d
