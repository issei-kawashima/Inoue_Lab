####!計算結果を出力してアニメを作製するコマンド
####!カレントディレクトリは適宜変更
####!参考：gnuplotで連番ファイルをプロットする↓
####!<http://qiita.com/NatsukiLab/items/f08445731bdb473c52ee>

reset


se term gif size 840, 480	#! 7:4 の大きさで出力
se term gif animate	#!アニメの生成
#!出力先ファイル名の設定
#se output "/Users/isseyshome/Documents/Lab/2dimensions/result_analysis/div_u.gif"
se output '\Users\Kawashima\Desktop\2-dimensions\visual-result\div_u.gif'

se pm3d map
#!表示領域の設定
se zr [-8:8]
se xr [0:24]
se yr [-8:8]
set cbrange [-5:5]
set xlabel 'x'
set ylabel 'y'
set zlabel 'div u'
set ticslevel 0
set size ratio -1

##!フォントの設定--------------
se xlabel font "Arial,12"
se ylabel font "Arial,12"
se zlabel font "Arial,12"
se tics font "Arial,12"
se key font "Arial,12"

set key outside top left #right top##!時間の出力場所
unset clabel #!等高線の凡例表示が消えた。やったぜ！

#!//////////出力するカラーのセット，気分で好きなやつ使って////////////
#set palette defined (0"blue",0.25"white",0.4"red")
#set palette rgbformulae 33,13,10	#!パステルっぽい色
#set palette rgbformulae 22,13,-31	#!コントラストが強い色
#set palette defined(0"#00008b",1"#2ca9e1",2"#38b48b",3.5"#ffff00",5"#eb6101",5.3"#c9171e") #!研究室でよく使ってるやつ
#set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)
#! ↑2017/6/9ネットで見つけたよさげなやつ
set palette define (-5"blue",-3"purple",-1"skyblue",0"white",1"gold",3"orange",5"dark-orange") #!0を白にしてる

set contour surface			#!等高線の出力
set cntrparam levels auto 5 	#!間隔
#se cntrparam levels incremental -5,0.06,5#等高線の初期値、間隔、最終値の設定
#set palette maxcolors 1
#!出力のループ========================================================
st = 0 #常に0
to = 75000 #計算ステップ数
dt = 0.0020 #時間ステップ刻み幅
offset = 1/dt #使用するコマの間隔(1秒毎に設定)
do for[i = st : to : offset]{
	#input1 = sprintf("/Users/isseyshome/Documents/Lab/ex_20/result_gradient/rho%05d.d", i)
	input1 = sprintf('\Users\Kawashima\Desktop\2-dimensions\result_super_2d-3Ma1.5\parameter%06d.txt', i)
	time =sprintf("t = %.1f",i/500.0)
	splot input1 u 1:2:6 w pm3d ti time
	pause 10
}
