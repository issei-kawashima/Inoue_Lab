!2次元と3次元の計算結果を比較するためのプログラム
!2020.06.03 3次元コードのデバックが完了したので2次元と3次元の比較は行なわない。
!今後は_0~_19(Nz-1)までの複数の2次元txtファイルを1つのtxtファイルにまとめるプログラム。
!それにより、oara viewで3次元での可視化を可能にするファイルが生成される
program main
implicit none
!integer,parameter :: Nx = 180 !#Mを出力結果から削除した最新版ではこちら
integer,parameter :: Nx = 181 !#Mが入っている旧バージョンでの計算結果はこちら。Nx+1にして#Mのせい読み込みがでずれないようにする
integer,parameter :: Ny = 100
integer,parameter :: Nz = 20
integer i, j,k,x,y!変数の宣言
double precision,dimension(0:Nx,0:Ny) :: compare_dim
double precision,dimension(0:Nx,0:Ny,3) :: dim2
double precision,dimension(0:Nx,0:Ny,4) :: dim3,dim310
double precision,dimension(0:Nx,0:Ny,12,0:Nz-1) :: dim_kukei
!条件に応じて、フォルダ名とCNTを変更する
integer,parameter :: CNT = 2000
character(len = 16) folder_name
character(len = 16) M_number
write(folder_name, '(a)') "super_3"
write(M_number, '(i6.6)') CNT
! open(10, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_2D/parameter000300.d")!入力ファイルを選択
! open(20, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_3D/parameter000300.d")!1/3角処理
! open(30, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_3D_1/parameter000000.d")!1/2角処理
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
!上はdim3のz方向同士の比較もできるように書いた
! open(3, file = "comparison_result.d")!出力ファイルの選択
  !入力データの読み込みと確認
do i=0,Ny
  do j = 0,Nx
  ! read(10,'(6f24.16)') dim2(j,i,1:3)!
  !write(*,*) dim2(j,i,1:3)!確認
  ! read(20,'(6f24.16)') dim3(j,i,1:4)
  ! read(30,'(6f24.16)') dim310(j,i,1:4)
  read(41,'(12f24.16)')dim_kukei(j,i,1:12,0)
  read(42,'(12f24.16)')dim_kukei(j,i,1:12,1)
  read(43,'(12f24.16)')dim_kukei(j,i,1:12,2)
  read(44,'(12f24.16)')dim_kukei(j,i,1:12,3)
  read(45,'(12f24.16)')dim_kukei(j,i,1:12,4)
  read(46,'(12f24.16)')dim_kukei(j,i,1:12,5)
  read(47,'(12f24.16)')dim_kukei(j,i,1:12,6)
  read(48,'(12f24.16)')dim_kukei(j,i,1:12,7)
  read(49,'(12f24.16)')dim_kukei(j,i,1:12,8)
  read(50,'(12f24.16)')dim_kukei(j,i,1:12,9)
  read(51,'(12f24.16)')dim_kukei(j,i,1:12,10)
  read(52,'(12f24.16)')dim_kukei(j,i,1:12,11)
  read(53,'(12f24.16)')dim_kukei(j,i,1:12,12)
  read(54,'(12f24.16)')dim_kukei(j,i,1:12,13)
  read(55,'(12f24.16)')dim_kukei(j,i,1:12,14)
  read(56,'(12f24.16)')dim_kukei(j,i,1:12,15)
  read(57,'(12f24.16)')dim_kukei(j,i,1:12,16)
  read(58,'(12f24.16)')dim_kukei(j,i,1:12,17)
  read(59,'(12f24.16)')dim_kukei(j,i,1:12,18)
  read(60,'(12f24.16)')dim_kukei(j,i,1:12,19)

  ! write(*,*) dim3(j,i,1:4)!確認
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
      if((dim_kukei(j,i,1,k) == 0.d0) .and. (dim_kukei(j,i,3,k) == 0.d0) .and.&
       (dim_kukei(j,i,5,k) == 0.d0) .and. (dim_kukei(j,i,7,k) == 0.d0)) then
      else
        write(99,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
        f24.16)') &
        dim_kukei(j,i,1,k),dim_kukei(j,i,3,k),dim_kukei(j,i,5,k)&
        ,dim_kukei(j,i,7,k),dim_kukei(j,i,9,k),dim_kukei(j,i,11,k)
      endif
    end do
  end do
enddo
close(99)



! compare_dim(:,:) = dim2(:,:,3) - dim3(:,:,4)
! compare_dim(:,:) = dim310(:,:,4) - dim3(:,:,4)

!特に差が大きいもののみを出力する
! do i =0,Ny
!   do j =0,Nx
!     if(abs(compare_dim(j,i)) >= 5.d-3) then
!       write(*,*) dim3(j,i,1),dim3(j,i,2),compare_dim(j,i)
!     endif
!   end do
! end do

!全ての結果を出力する
! do i =0,Ny
!   do j =0,Nx
!     if((compare_dim(j,i) == 0.d0) .and. (dim3(j,i,1) == 0.d0) .and. (dim3(j,i,2) == 0.d0)) then
!       else
!         write(3,'(6f24.16)') dim3(j,i,1),dim3(j,i,2),compare_dim(j,i)
!     endif
!   end do
!   write(3,*)
! end do
! close(3)

end program main
