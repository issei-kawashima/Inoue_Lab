!2次元と3次元の計算結果を比較するためのプログラム
program main
implicit none
integer,parameter :: Nx = 181
integer,parameter :: Ny = 100
integer i, j,x,y!変数の宣言
double precision,dimension(0:Nx,0:Ny) :: compare_dim
double precision,dimension(0:Nx,0:Ny,3) :: dim2
double precision,dimension(0:Nx,0:Ny,4) :: dim3,dim310
open(10, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_2D/parameter000300.d")!入力ファイルを選択
open(20, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_3D/parameter000300.d")!1/3角処理
! open(30, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/result_3D_1/parameter000000.d")!1/2角処理
!上はdim3のz方向同士の比較もできるように書いた
! open(3, file = "comparison_result.d")!出力ファイルの選択
  !入力データの読み込みと確認
do i=0,Ny
  do j = 0,Nx
    read(10,'(6f24.16)') dim2(j,i,1:3)!
    !write(*,*) dim2(j,i,1:3)!確認
    read(20,'(6f24.16)') dim3(j,i,1:4)
    ! read(30,'(6f24.16)') dim310(j,i,1:4)
    ! write(*,*) dim3(j,i,1:4)!確認
  enddo
enddo
close(10)
close(20)
! close(30)

compare_dim(:,:) = dim2(:,:,3) - dim3(:,:,4)
! compare_dim(:,:) = dim310(:,:,4) - dim3(:,:,4)

!特に差が大きいもののみを出力する
do i =0,Ny
  do j =0,Nx
    if(abs(compare_dim(j,i)) >= 5.d-3) then
      write(*,*) dim3(j,i,1),dim3(j,i,2),compare_dim(j,i)
    endif
  end do
end do

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
