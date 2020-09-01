!3次元同士の任意のz面での結果を比較するコード
!3次元コードを高速化などの改修を加えた際に問題なく計算できているかの比較をするコード
program main
implicit none
!!ここで各計算結果の格子点数を設定!!!!!!!!!!
integer,parameter :: Nx1 = 181
integer,parameter :: Ny1 = 100
integer,parameter :: Nx2 = 181
integer,parameter :: Ny2 = 100
integer i,ii,j,jj,k
double precision,dimension(0:Nx1,0:Ny1,12) :: in_dim1
double precision,dimension(0:Nx2,0:Ny2,12) :: in_dim2
double precision,dimension(0:Nx1,0:Ny1,12) :: compare_dim
double precision result_1,result_3,result_7
result_1=0.d0;result_3=0.d0;result_7=0.d0
!!!比較するファイルのz面は同じものを選ぶこと！！！
!!現状ではそうでないとcompare_dimを計算してくれない！！！
open(41, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_super_3_2.4/parameter002000_02.txt")
open(42, file = "/Users/isseyshome/Downloads/parameter002000_02.txt")

open(99, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_analysis/rapid_dif.txt")
write(*,*) 1
do i=0,Ny1
  do j = 0,Nx1
  read(41,'(12f24.16)') in_dim1(j,i,1:12)
  enddo
enddo
close(41)
write(*,*) 2
do i=0,Ny2
  do j = 0,Nx2
  read(42,'(12f24.16)') in_dim2(j,i,1:12)
  enddo
enddo
close(42)
write(*,*) 3
do ii = 0,Ny2
  do i = 0,Ny1
    do jj = 0,Nx2
      do j=0,Nx1
        if ((in_dim1(j,i,1) == in_dim2(jj,ii,1)).and.&
            (in_dim1(j,i,3) == in_dim2(jj,ii,3)).and.&
            (in_dim1(j,i,5) == in_dim2(jj,ii,5))) then
          do k=1,6
            compare_dim(j,i,k) = in_dim1(j,i,k)
          end do
          do k =7,12
            compare_dim(j,i,k) = in_dim1(j,i,k) - in_dim2(jj,ii,k)
          end do
        endif
      end do
    end do
  end do
end do
write(*,*) 4

do i =0,Ny1
  do j =0,Nx1
    if((compare_dim(j,i,1) == 0.d0) .and. (compare_dim(j,i,3) == 0.d0) .and.&
     (compare_dim(j,i,5) == 0.d0)) then
    else
      write(99,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
      f24.16)') &
      compare_dim(j,i,1),compare_dim(j,i,3),compare_dim(j,i,5)&
      ,compare_dim(j,i,7),compare_dim(j,i,9),compare_dim(j,i,11)
    endif
  end do
  write(99,*) ''
end do
close(99)

do i =0,Ny1
  do j =0,Nx1
    ! if(abs(compare_dim(j,i,7)) >= 0.2d0) then
    !   write(*,*)compare_dim(j,i,1),compare_dim(j,i,3),compare_dim(j,i,7)
    ! endif
    result_1 = result_1 + compare_dim(j,i,7)
    result_3 = result_3 + compare_dim(j,i,9)
    result_7 = result_7 + compare_dim(j,i,11)
  end do
end do
write(*,*) "2つの計算結果の誤差の合計を出力"
write(*,'(4A)') "rho,","omega3,","dp,","(dpはNaNの計算結果では必ず0になる)"
write(*,*) result_1,result_3,result_7
end program main
