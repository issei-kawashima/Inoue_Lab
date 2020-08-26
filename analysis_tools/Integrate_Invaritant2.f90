!3次元結果の第二不変量Qのみを複数z面から1つのtxtファイルにまとめるコード
!2次元txtを1つのtxtにまとめる昨日も兼ねている
!2020.08.26 Q=xxxの特定面のみの情報を切り出す処理は未実装
program main
implicit none
integer i,j,k,l,time
integer,parameter :: Nx = 360+1 !一行改行分を含めてる
integer,parameter :: Ny = 200
integer,parameter :: Nz = 20
integer,parameter :: CNT = 1500
double precision,allocatable,dimension(:,:,:,:) :: dim_kukei
character(len = 16) z_name
character(len = 16) folder_name
character(len = 16) M_number

allocate(dim_kukei(0:Nx,0:Ny,0:Nz-1,14))
write(folder_name, '(a)') "ompMa2.0"
write(M_number, '(i6.6)') CNT

do k=0,Nz-1
  write(z_name, '(i2.2)') k
  open(11, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
                      //trim(folder_name)//"/parameter"//trim(M_number)//"_"//trim(z_name)//".txt")
  !入力データの読み込みと確認
  do i=0,Ny
    do j = 0,Nx
      read(11,'(14f24.16)') dim_kukei(j,i,k,1:14)
    end do
  end do
  close(11)
enddo

!複数の2次元txtファイルから第二不変量と座標のみを1つのtxtファイルにまとめる時はこちら
open(99, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_analysis/"&
          //trim(folder_name)//"_2ndInvaritant"//trim(M_number)//".txt")
!複数の2次元txtファイルを1つのtxtファイルにまとめる時はこちら
! open(99, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_analysis/"&
!           //trim(folder_name)//"_All"//trim(M_number)//".txt")
do k = 0,Nz-1
  do i =0,Ny
    do j =0,Nx
      if((dim_kukei(j,i,k,1) == 0.d0) .and. (dim_kukei(j,i,k,3) == 0.d0) .and.&
       (dim_kukei(j,i,k,5) == 0.d0) .and. (dim_kukei(j,i,k,7) == 0.d0)) then
      else
        write(99,'(f24.16,",",f24.16,",",f24.16,",",f24.16)') &
        dim_kukei(j,i,k,1),dim_kukei(j,i,k,3),dim_kukei(j,i,k,5),&
        dim_kukei(j,i,k,13)!x,y,z,第二不変量Q(~2020.08.26)
        ! dim_kukei(j,i,k,11)!x,y,z,第二不変量Q(2020.08.27~出力順を変更した)

!全ての変数を1つのtxtファイルにまとめる時はこちら
        ! write(99,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,&
        ! ",",f24.16,",",f24.16)') &
        ! dim_kukei(j,i,k,1),dim_kukei(j,i,k,3),dim_kukei(j,i,k,5),&
        ! dim_kukei(j,i,k,7),dim_kukei(j,i,k,9),dim_kukei(j,i,k,11),&
        ! dim_kukei(j,i,k,13)
      endif
    end do
  end do
enddo
close(99)
deallocate(dim_kukei)
end program main