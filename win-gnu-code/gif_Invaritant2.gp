####!計算結果を出力してアニメを作製するコマンド
####!カレントディレクトリは適宜変更
####!参考：gnuplotで連番ファイルをプロットする↓
####!<http://qiita.com/NatsukiLab/items/f08445731bdb473c52ee>

reset

se term gif size 840, 480	#! 7:4 の大きさで出力
se term gif animate	#!アニメの生成
#!出力先ファイル名の設定
se output '\Users\Kawashima\Desktop\3-dimensions\result_square\Invaritant2.gif'
#set view 70, 10, 1, 1

se pm3d map
#!表示領域の設定
se zr [0.65:1.4]
set cbrange [0.65:1.4]
set xl 'x'
set ylabel 'y'
set zl "Invaritant 2"
set ticslevel 0
set size ratio -1

se xr[0:36]
se yr[-10:10]

##!フォントの設定--------------
se xlabel font 'Arial,12'
se ylabel font 'Arial,12'
se zlabel font 'Arial,12'
se tics font 'Arial,12'
se key font 'Arial,12'

set colorbox horiz user origin .40,.04 size .50,.05
set key at 1,11	#right top##!時間の出力場所
unset clabel #!等高線の凡例表示が消えた

#!//////////出力するカラーのセット////////////
set palette defined (-5'blue',0'white',5'red')

set contour surface
se cntrparam levels incremental -2,0.2,2 #手動で等高線の最大値・最小値を設定して、間隔も決める
#set cntrparam levels auto 5 #自動で等高線を設定する(間隔は自分で設定)

#!出力のループ========================================================
st = 0 #開始したい計算ステップ数
to = 15000 #計算ステップ数
dt = 0.01 #時間ステップ刻み幅
offset = 1/dt #使用するコマの間隔(1秒毎に設定)
do for[i = st : to : offset]{
	input1 = sprintf('\Users\Kawashima\Desktop\3-dimensions\result_square\parameter%06d_10.txt', i)
	time =sprintf("t = %.1f",i/100.0)#dtの逆数でiを割る
	splot input1 u 1:2:6 w pm3d ti time
	pause 10
}
