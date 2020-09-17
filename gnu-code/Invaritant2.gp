reset

#For Mac
set term x11

#For Windows
#se term wxt

se pm3d map
#se size square
#se view 90,270,1,1
se contour#等高線の表示
se cntrparam levels incremental -0.2,0.05,0.2#等高線の初期値、間隔、最終値の設定
#se cntrparam levels incremental -2,0.1,2#等高線の初期値、間隔、最終値の設定
unset clabel#等高線の凡例の消去
set colorbox horiz user origin .25,.07 size .65,.04
set cbrange [-0.1:0.1]
set palette defined (-0.1"blue",0"white",0.1"red")

set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"

set xl "x"
set yl "y"
set zl "Invaritant 2"
set grid

set xrange [0:24]
set yrange [-8:8]
#set zrange [0:1]

#For Mac
set term x11 1

#For Windows
#se term wxt 1


#(通常&NaN)3次元
#For Mac
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_omp/parameter001500_18.txt' u 1:2:6  w pm3d t ""
splot "/Users/isseyshome/Downloads/parameter012900_07.txt" u 1:2:6  w pm3d t ""


#2次元
#For Mac
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/2dimensions/result_super_2d-3Ma1.5/parameter075000.txt' u 1:2:7  w pm3d t ""

#For Windows
#splot 'C:\Users\Kawashima\Desktop\2-dimensions\result_super_2d-3Ma1.5\parameter005000.txt' u 1:2:7  w pm3d t ""

#############JPEGで出力#############################################
#se term jpeg
#For Mac
#se out '/Users/isseyshome/Documents/Lab/#3 修士論文/Visual_result/Invaritant2_005000.jpg'

#For Windows
#se out 'C:\Users\Kawashima\Desktop\2-dimensions\visual-result\Invaritant2_005000.jpg'
#rep
#############JPEGで出力#############################################

#卒論のデータ用
#set term x11
#splot '/Users/isseyshome/Documents/Lab/3-dimensions/dp000100.d' w pm3d
