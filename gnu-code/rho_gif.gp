####!計算結果を出力してアニメを作製するコマンド
####!カレントディレクトリは適宜変更
####!参考：gnuplotで連番ファイルをプロットする↓
####!<http://qiita.com/NatsukiLab/items/f08445731bdb473c52ee>

reset


se term gif size 840, 480	#! 7:4 の大きさで出力
se term gif animate	#!アニメの生成
se output "/Users/isseyshome/Documents/Lab/#3 修士論文/Visual_result/subsonic_rho.gif"	#!出力先ファイル名の設定

se pm3d map
#!表示領域の設定
se zr [0.65:1.4]
set cbrange [0.65:1.4]
set xlabel 'x'
set ylabel 'y'
set zlabel 'rho'
set ticslevel 0
set size ratio -1

##!フォントの設定--------------
se xlabel font "Arial,12"
se ylabel font "Arial,12"
se zlabel font "Arial,12"
se tics font "Arial,12"
se key font "Arial,12"

set key at 1,11	#right top##!時間の出力場所
unset clabel #!等高線の凡例表示が消えた

#!//////////出力するカラーのセット////////////
set palette define (0.8"blue",1"white",1.2"orange") #1を白にしてる

set contour surface			#!等高線の出力
set cntrparam levels auto 6	#!間隔
#!出力のループ========================================================
st = 0 #開始したい計算ステップ数
to = 15000 #計算ステップ数
dt = 0.01 #時間ステップ刻み幅
offset = 1/dt #使用するコマの間隔(1秒毎に設定)
do for[i = st : to : offset]{
	input1 = sprintf("/Users/isseyshome/Downloads/result_sub/parameter%06d_07.txt", i)
	#input1 = sprintf('\Users\Kawashima\Desktop\2-dimensions\result_super_2d-3Ma1.5\parameter%06d.txt', i)
	time =sprintf("t = %.1f",i/100.0)
	splot input1 u 1:2:4 w pm3d ti time
	pause 10
}
