!paraviewでの可視化のために.dファイルを.txtファイルへ変換するプログラム
!2020.06.03 今後の計算プログラムでの出力方式を.dから.txtに変更したので過去データをpara viewで可視化したい時にのみ使用する
program main
implicit none
integer,parameter :: Nx = 181
integer,parameter :: Ny = 100
integer,parameter :: Nz = 20
integer i, j,k,x,y!変数の宣言
double precision,dimension(0:Nx,0:Ny) :: compare_dim
double precision,dimension(0:Nx,0:Ny,3) :: dim2
double precision,dimension(0:Nx,0:Ny,4) :: dim3,dim310
double precision,dimension(0:Nx,0:Ny,6,0:Nz-1) :: dim_kukei
open(41, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_00.d")
open(42, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_01.d")
open(43, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_02.d")
open(44, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_03.d")
open(45, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_04.d")
open(46, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_05.d")
open(47, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_06.d")
open(48, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_07.d")
open(49, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_08.d")
open(50, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_09.d")
open(51, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_10.d")
open(52, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_11.d")
open(53, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_12.d")
open(54, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_13.d")
open(55, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_14.d")
open(56, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_15.d")
open(57, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_16.d")
open(58, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_17.d")
open(59, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_18.d")
open(60, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_kukei/parameter050500_19.d")

open(71, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_00.txt")
open(72, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_01.txt")
open(73, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_02.txt")
open(74, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_03.txt")
open(75, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_04.txt")
open(76, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_05.txt")
open(77, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_06.txt")
open(78, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_07.txt")
open(79, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_08.txt")
open(80, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_09.txt")
open(81, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_10.txt")
open(82, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_11.txt")
open(83, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_12.txt")
open(84, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_13.txt")
open(85, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_14.txt")
open(86, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_15.txt")
open(87, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_16.txt")
open(88, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_17.txt")
open(89, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_18.txt")
open(90, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_analysis/parameter050500_19.txt")

!入力データの読み込みと確認
do i=0,Ny
  do j = 0,Nx
  read(41,'(6f24.16)') dim_kukei(j,i,1:6,0)
  read(42,'(6f24.16)') dim_kukei(j,i,1:6,1)
  read(43,'(6f24.16)') dim_kukei(j,i,1:6,2)
  read(44,'(6f24.16)') dim_kukei(j,i,1:6,3)
  read(45,'(6f24.16)') dim_kukei(j,i,1:6,4)
  read(46,'(6f24.16)') dim_kukei(j,i,1:6,5)
  read(47,'(6f24.16)') dim_kukei(j,i,1:6,6)
  read(48,'(6f24.16)') dim_kukei(j,i,1:6,7)
  read(49,'(6f24.16)') dim_kukei(j,i,1:6,8)
  read(50,'(6f24.16)') dim_kukei(j,i,1:6,9)
  read(51,'(6f24.16)') dim_kukei(j,i,1:6,10)
  read(52,'(6f24.16)') dim_kukei(j,i,1:6,11)
  read(53,'(6f24.16)') dim_kukei(j,i,1:6,12)
  read(54,'(6f24.16)') dim_kukei(j,i,1:6,13)
  read(55,'(6f24.16)') dim_kukei(j,i,1:6,14)
  read(56,'(6f24.16)') dim_kukei(j,i,1:6,15)
  read(57,'(6f24.16)') dim_kukei(j,i,1:6,16)
  read(58,'(6f24.16)') dim_kukei(j,i,1:6,17)
  read(59,'(6f24.16)') dim_kukei(j,i,1:6,18)
  read(60,'(6f24.16)') dim_kukei(j,i,1:6,19)
  enddo
enddo
close(41)
close(42)
close(43)
close(44)
close(45)
close(46)
close(47)
close(48)
close(49)
close(50)
close(51)
close(52)
close(53)
close(54)
close(55)
close(56)
close(57)
close(58)
close(59)
close(60)

do i =0,Ny
  do j =0,Nx
    if((dim_kukei(j,i,3,0) == 0.d0) .and. (dim_kukei(j,i,1,0) == 0.d0) .and.&
       (dim_kukei(j,i,2,0) == 0.d0) .and. (dim_kukei(j,i,4,0) == 0.d0)) then
    else
      write(71,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,0),dim_kukei(j,i,2,0),dim_kukei(j,i,3,0)&
      ,dim_kukei(j,i,4,0),dim_kukei(j,i,5,0),dim_kukei(j,i,6,0)
      write(72,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,1),dim_kukei(j,i,2,1),dim_kukei(j,i,3,1)&
      ,dim_kukei(j,i,4,1),dim_kukei(j,i,5,1),dim_kukei(j,i,6,1)
      write(73,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,2),dim_kukei(j,i,2,2),dim_kukei(j,i,3,2)&
      ,dim_kukei(j,i,4,2),dim_kukei(j,i,5,2),dim_kukei(j,i,6,2)
      write(74,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,3),dim_kukei(j,i,2,3),dim_kukei(j,i,3,3)&
      ,dim_kukei(j,i,4,3),dim_kukei(j,i,5,3),dim_kukei(j,i,6,3)
      write(75,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,4),dim_kukei(j,i,2,4),dim_kukei(j,i,3,4)&
      ,dim_kukei(j,i,4,4),dim_kukei(j,i,5,4),dim_kukei(j,i,6,4)
      write(76,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,5),dim_kukei(j,i,2,5),dim_kukei(j,i,3,5)&
      ,dim_kukei(j,i,4,5),dim_kukei(j,i,5,5),dim_kukei(j,i,6,5)
      write(77,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,6),dim_kukei(j,i,2,6),dim_kukei(j,i,3,6)&
      ,dim_kukei(j,i,4,6),dim_kukei(j,i,5,6),dim_kukei(j,i,6,6)
      write(78,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,7),dim_kukei(j,i,2,7),dim_kukei(j,i,3,7)&
      ,dim_kukei(j,i,4,7),dim_kukei(j,i,5,7),dim_kukei(j,i,6,7)
      write(79,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,8),dim_kukei(j,i,2,8),dim_kukei(j,i,3,8)&
      ,dim_kukei(j,i,4,8),dim_kukei(j,i,5,8),dim_kukei(j,i,6,8)
      write(80,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,9),dim_kukei(j,i,2,9),dim_kukei(j,i,3,9)&
      ,dim_kukei(j,i,4,9),dim_kukei(j,i,5,9),dim_kukei(j,i,6,9)
      write(81,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,10),dim_kukei(j,i,2,10),dim_kukei(j,i,3,10)&
      ,dim_kukei(j,i,4,10),dim_kukei(j,i,5,10),dim_kukei(j,i,6,10)
      write(82,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,11),dim_kukei(j,i,2,11),dim_kukei(j,i,3,11)&
      ,dim_kukei(j,i,4,11),dim_kukei(j,i,5,11),dim_kukei(j,i,6,11)
      write(83,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,12),dim_kukei(j,i,2,12),dim_kukei(j,i,3,12)&
      ,dim_kukei(j,i,4,12),dim_kukei(j,i,5,12),dim_kukei(j,i,6,12)
      write(84,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,13),dim_kukei(j,i,2,13),dim_kukei(j,i,3,13)&
      ,dim_kukei(j,i,4,13),dim_kukei(j,i,5,13),dim_kukei(j,i,6,13)
      write(85,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,14),dim_kukei(j,i,2,14),dim_kukei(j,i,3,14)&
      ,dim_kukei(j,i,4,14),dim_kukei(j,i,5,14),dim_kukei(j,i,6,14)
      write(86,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,15),dim_kukei(j,i,2,15),dim_kukei(j,i,3,15)&
      ,dim_kukei(j,i,4,15),dim_kukei(j,i,5,15),dim_kukei(j,i,6,15)
      write(87,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,16),dim_kukei(j,i,2,16),dim_kukei(j,i,3,16)&
      ,dim_kukei(j,i,4,16),dim_kukei(j,i,5,16),dim_kukei(j,i,6,16)
      write(88,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,17),dim_kukei(j,i,2,17),dim_kukei(j,i,3,17)&
      ,dim_kukei(j,i,4,17),dim_kukei(j,i,5,17),dim_kukei(j,i,6,17)
      write(89,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,18),dim_kukei(j,i,2,18),dim_kukei(j,i,3,18)&
      ,dim_kukei(j,i,4,18),dim_kukei(j,i,5,18),dim_kukei(j,i,6,18)
      write(90,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') dim_kukei(j,i,1,19),dim_kukei(j,i,2,19),dim_kukei(j,i,3,19)&
      ,dim_kukei(j,i,4,19),dim_kukei(j,i,5,19),dim_kukei(j,i,6,19)
    endif
  end do
end do
close(71)
close(72)
close(73)
close(74)
close(75)
close(76)
close(77)
close(78)
close(79)
close(80)
close(81)
close(82)
close(83)
close(84)
close(85)
close(86)
close(87)
close(88)
close(89)
close(90)


end program main
