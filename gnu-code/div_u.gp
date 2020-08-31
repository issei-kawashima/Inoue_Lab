reset
#For Mac
set term x11

#For Windows
#set term wxt

se pm3d map
#se size square
#se view 90,270,1,1
se contour#等高線の表示
se cntrparam levels incremental -4,0.5,4#等高線の初期値、間隔、最終値の設定
#se cntrparam levels incremental -1,0.03,0.5#等高線の初期値、間隔、最終値の設定
unset clabel#等高線の凡例の消去
set colorbox horiz user origin .25,.07 size .65,.04
set cbrange [-2:2]
set palette defined (-2"blue",0"white", 2"red")

set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"

set xl "x"
set yl "y"
set zl "Sound element div_u"
set grid

set xrange [0:24]
set yrange [-10:10]
#set zrange [0:1]

#For Mac
set term x11 1

#For Windows
#set term wxt 1

#3次元
#For Mac
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_omp/parameter001000_18.txt' u 1:2:5  w pm3d t ""

####dpの出力位置を変更する前のデータ用(8.25~26)#######################
#3次元　通常
#For Mac
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_omp/parameter005000_18.txt' u 1:2:7  w pm3d t ""
splot '/Users/isseyshome/Downloads/result_omp/parameter006000_18.txt' u 1:2:7  w pm3d t ""
#3次元 Nan直前
#For Mac
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_omp/parameter001000_18.txt' u 1:2:6  w pm3d t ""
#################################################################


#set term x11 2
#2次元
#For Mac
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/2dimensions/result_super_2d-3Ma2/parameter001000.txt' u 1:2:5  w pm3d t ""

####dpの出力位置を変更する前のデータ用(8.25~26)#######################
#NaNの場合は上のコードで大丈夫
#2次元　通常
#For Mac
#splot '/Users/isseyshome/Documents/GitHub/Inoue_Lab/2dimensions/result_super_2d-3Ma2/parameter001000.txt' u 1:2:6  w pm3d t ""

#For Windows
#splot 'C:\Users\Kawashima\Desktop\2-dimensions\result_super_2d-3Ma1.5\parameter025000.txt' u 1:2:6  w pm3d t ""
#################################################################

#############JPEGで出力#############################################
se term jpeg
#For Mac
se out '/Users/isseyshome/Documents/Lab/#3 修士論文/Visual_result/div_u_005000.jpg'


#se out 'C:\Users\Kawashima\Desktop\2-dimensions\visual-result\div_u025000.jpg' #For Windows
rep
#############JPEGで出力#############################################

#卒論のデータ用
#set term x11
#For Mac
#splot '/Users/isseyshome/Documents/Lab/3-dimensions/dp000100.d' w pm3d
