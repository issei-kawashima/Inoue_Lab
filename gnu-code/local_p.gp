x#色を自動で設定させた
#今までの設定をリセット
reset
set term x11 1
#aquaの後に数字をつけているのはグラフを別ウィンドウで表示させるため
#グラフタイトルので設定
#se tit "Local Pressure (x,y)=(2.5,0)"
#se tit font "Arial 20"
#軸設定
se xl font "Arial 12"
se yl font "Arial 12"
set tics font "Arial 12"
set key below right box font "Arial 10"
set xrange [-10:10]
set yrange [-0.1:1.1]
se xl "y"
se yl "u(y)"
unset key

p '/Users/isseyshome/Documents/Lab/thesis/result_SC/u.d'  w l t "" lc "black" dt (7,3)

#JPEG形式で出力
se term jpeg
se out "/Users/isseyshome/Documents/Lab/thesis/result_jpeg/u(y).jpg"
rep
