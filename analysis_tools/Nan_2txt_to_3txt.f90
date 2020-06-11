!2次元と3次元の計算結果を比較するためのプログラム
!2020.06.03 3次元コードのデバックが完了したので2次元と3次元の比較は行なわない。
!今後は_0~_19(Nz-1)までの複数の2次元txtファイルを1つのtxtファイルにまとめるプログラム。
!それにより、oara viewで3次元での可視化を可能にするファイルが生成される
program main
implicit none
integer,parameter :: Nx = 181
integer,parameter :: Ny = 100
integer,parameter :: Nz = 20
integer i, j,k,x,y!変数の宣言
double precision,dimension(0:Nx,0:Ny) :: compare_dim
double precision,dimension(0:Nx,0:Ny,3) :: dim2
double precision,dimension(0:Nx,0:Ny,4) :: dim3,dim310
double precision,dimension(0:Nx,0:Ny,10,0:Nz-1) :: dim_kukei
!条件に応じて、フォルダ名とCNTを変更する
integer,parameter :: CNT = 2199
character(len = 16) folder_name
character(len = 16) M_number
write(folder_name, '(a)') "super_3"
write(M_number, '(i6.6)') CNT
open(41, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_00.txt")
open(42, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_01.txt")
open(43, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_02.txt")
open(44, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_03.txt")
open(45, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_04.txt")
open(46, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_05.txt")
open(47, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_06.txt")
open(48, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_07.txt")
open(49, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_08.txt")
open(50, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_09.txt")
open(51, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_10.txt")
open(52, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_11.txt")
open(53, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_12.txt")
open(54, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_13.txt")
open(55, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_14.txt")
open(56, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_15.txt")
open(57, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_16.txt")
open(58, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_17.txt")
open(59, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_18.txt")
open(60, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                    //trim(folder_name)//"/parameter"//trim(M_number)//"_19.txt")


open(99, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_analysis/parameter"&
                //trim(M_number)//".txt")
!入力データの読み込みと確認
do i=0,Ny
  do j = 0,Nx
  read(41,'(10f24.16)')dim_kukei(j,i,1:10,0)
  read(42,'(10f24.16)')dim_kukei(j,i,1:10,1)
  read(43,'(10f24.16)')dim_kukei(j,i,1:10,2)
  read(44,'(10f24.16)')dim_kukei(j,i,1:10,3)
  read(45,'(10f24.16)')dim_kukei(j,i,1:10,4)
  read(46,'(10f24.16)')dim_kukei(j,i,1:10,5)
  read(47,'(10f24.16)')dim_kukei(j,i,1:10,6)
  read(48,'(10f24.16)')dim_kukei(j,i,1:10,7)
  read(49,'(10f24.16)')dim_kukei(j,i,1:10,8)
  read(50,'(10f24.16)')dim_kukei(j,i,1:10,9)
  read(51,'(10f24.16)')dim_kukei(j,i,1:10,10)
  read(52,'(10f24.16)')dim_kukei(j,i,1:10,11)
  read(53,'(10f24.16)')dim_kukei(j,i,1:10,12)
  read(54,'(10f24.16)')dim_kukei(j,i,1:10,13)
  read(55,'(10f24.16)')dim_kukei(j,i,1:10,14)
  read(56,'(10f24.16)')dim_kukei(j,i,1:10,15)
  read(57,'(10f24.16)')dim_kukei(j,i,1:10,16)
  read(58,'(10f24.16)')dim_kukei(j,i,1:10,17)
  read(59,'(10f24.16)')dim_kukei(j,i,1:10,18)
  read(60,'(10f24.16)')dim_kukei(j,i,1:10,19)
  enddo
enddo
! close(10)
! close(20)
! close(30)
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

do k = 0,Nz-1
  do i =0,Ny
    do j =0,Nx
      write(*,*) dim_kukei(j,i,1,0),dim_kukei(j,i,3,0),dim_kukei(j,i,5,0),dim_kukei(j,i,7,0),dim_kukei(j,i,9,0)
      if((dim_kukei(j,i,1,k) == 0.d0) .and. (dim_kukei(j,i,3,k) == 0.d0) .and.&
       (dim_kukei(j,i,5,k) == 0.d0) .and. (dim_kukei(j,i,7,k) == 0.d0)) then
      else
        write(99,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16)') &
        dim_kukei(j,i,1,k),dim_kukei(j,i,3,k),dim_kukei(j,i,5,k)&
        ,dim_kukei(j,i,7,k),dim_kukei(j,i,9,k)
      endif
    end do
  end do
enddo
close(99)

end program main
