reset
#set term x11
se term wxt
se pm3d map
#se size square
#se view 90,270,1,1
se contour#等高線の表示
se cntrparam levels incremental -1,0.05,1#等高線の初期値、間隔、最終値の設定
unset clabel#等高線の凡例の消去
set colorbox horiz user origin .25,.07 size .65,.04
set palette defined (-0.5'blue',0'white',0.5'red')
set cbrange [-0.5:0.5]

set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set tics font "Arial,12"
set key below right box font "Arial,12"

set xl "x"
set yl "y"
set zl "Invaritant 2"
set grid

set xrange [0:36]
set yrange [-10:10]
#set zrange [0:1]


#3次元
set term wxt 1
set title "vorticity t=160 z=0" offset -30,0
splot '\Users\Kawashima\Desktop\3-dimensions\result_square\parameter016000_00.txt' using 1:2:6 w pm3d t ""
set term wxt 2
set title "vorticity t=160 z=10" offset -30,0
splot '\Users\Kawashima\Desktop\3-dimensions\result_square\parameter016000_10.txt' using 1:2:6 w pm3d t ""
set term wxt 3
set title "vorticity t=160 z=17" offset -30,0
splot '\Users\Kawashima\Desktop\3-dimensions\result_square\parameter016000_17.txt' using 1:2:6 w pm3d t ""
set term wxt 4
set title "vorticity t=160 z=33" offset -30,0
splot '\Users\Kawashima\Desktop\3-dimensions\result_square\parameter016000_33.txt' using 1:2:6 w pm3d t ""
set term wxt 5
set title "vorticity t=160 z=40" offset -30,0
splot '\Users\Kawashima\Desktop\3-dimensions\result_square\parameter016000_40.txt' using 1:2:6 w pm3d t ""
set term wxt 6
set title "vorticity t=160 z=50" offset -30,0
splot '\Users\Kawashima\Desktop\3-dimensions\result_square\parameter016000_50.txt' using 1:2:6 w pm3d t ""

#2次元
#splot 'C:\Users\Kawashima\Desktop\2-dimensions\result_super_2d-3Ma1.5\parameter005000.txt' u 1:2:7  w pm3d t ""

#se term jpeg
#se out 'C:\Users\Kawashima\Desktop\2-dimensions\visual-result\Invaritant2_005000.jpg'
#rep
