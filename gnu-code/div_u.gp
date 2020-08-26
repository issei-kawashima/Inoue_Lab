reset
set term x11
se pm3d map
#se size square
#se view 90,270,1,1
se contour#等高線の表示
#se cntrparam levels incremental -0.2,0.05,0.2#等高線の初期値、間隔、最終値の設定
se cntrparam levels incremental -2,0.3,2#等高線の初期値、間隔、最終値の設定
unset clabel#等高線の凡例の消去
set colorbox horiz user origin .25,.07 size .65,.04
#set cbrange [-2.5:5]
#set palette defined (-2.5"blue",0"white",5"red")

set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"

set xl "x"
set yl "y"
set zl "Sound element div_u"
set grid

set xrange [0:10]
set yrange [-3:3]
#set zrange [0:1]

set term x11 1
#3次元
splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_omp/parameter001000_18.txt' u 1:2:5  w pm3d t ""

####dpの出力位置を変更する前のデータ用(8.25~26)#######################
#3次元　通常
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_omp/parameter001000_18.txt' u 1:2:7  w pm3d t ""
#3次元 Nan直前
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_omp/parameter001000_18.txt' u 1:2:6  w pm3d t ""
#################################################################


set term x11 2
#2次元
splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/2dimensions/result_super_2d-3Ma2/parameter001000.txt' u 1:2:5  w pm3d t ""

####dpの出力位置を変更する前のデータ用(8.25~26)#######################
#NaNの場合は上のコードで大丈夫
#2次元　通常
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/2dimensions/result_super_2d-3Ma2/parameter001000.txt' u 1:2:6  w pm3d t ""
#################################################################
#se term jpeg
#se out "/Users/isseyshome/Documents/Lab/3-dimensions/result_2D/vorticity000000.jpg"
#rep
#set term x11
#splot '/Users/isseyshome/Documents/Lab/3-dimensions/dp000100.d' w pm3d
