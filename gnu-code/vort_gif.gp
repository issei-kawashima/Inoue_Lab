####!計算結果を出力してアニメを作製するコマンド
####!カレントディレクトリは適宜変更
####!参考：gnuplotで連番ファイルをプロットする↓
####!<http://qiita.com/NatsukiLab/items/f08445731bdb473c52ee>

reset


se term gif size 840, 480	#! 7:4 の大きさで出力
se term gif animate	#!アニメの生成
se output "/Users/isseyshome/Documents/Lab/thesis/vorticity-100.gif"	#!出力先ファイル名の設定
#set view 70, 10, 1, 1 #視点の変更
#!表示領域の設定
se pm3d map
se zr [-10:10]
se xr[0:10]
se yr[-4:4]
set cbrange [-4:4]#色分けする範囲の設定
set xlabel 'x'
set ylabel 'y'
set zlabel 'vorticity'
set ticslevel 0

##!フォントの設定--------------
se xlabel font "Arial,12"
se ylabel font "Arial,12"
se zlabel font "Arial,12"
se tics font "Arial,12"
se key font "Arial,12"

set key outside top left #!時間の出力場所
set colorbox horiz user origin .25,.07 size .65,.04
unset clabel #!等高線の凡例表示が消えた。やったぜ！
set contour surface			#!等高線の出力
se cntrparam levels incremental -4,0.3,4#等高線の初期値、間隔、最終値の設定
#set cntrparam levels auto 5 	#!間隔
#set palette maxcolors 1

#!//////////出力するカラーのセット，気分で好きなやつ使って////////////
set palette defined (-2.5"blue",0"white",2.5"red")
#set palette rgbformulae 33,13,10	#!パステルっぽい色
#set palette rgbformulae 22,13,-31	#!コントラストが強い色
#set palette defined(0"#00008b",1"#2ca9e1",2"#38b48b",3.5"#ffff00",5"#eb6101",5.3"#c9171e") #!研究室でよく使ってるやつ
#set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)
#! ↑2017/6/9ネットで見つけたよさげなやつ
#set palette define (-0.3"blue",0"white",0.3"yellow",0.7"orange") #!0を白にしてる


#!出力のループ========================================================
st = 0
to = 15000
offset = 300
do for[i = st : to : offset]{
	input1 = sprintf("/Users/isseyshome/Documents/Lab/thesis/result_KORe100-vortex/parameter%06d.d", i)
	time =sprintf("t = %.1f",i/100.0)
	splot input1 u 1:2:4 w pm3d ti time #!'' w pm3d,¥
	#pause 10
}
