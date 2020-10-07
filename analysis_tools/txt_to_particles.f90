!3次元結果の第二不変量Qのみを複数z面から1つのtxtファイルにまとめるコード
!2次元txtを1つのtxtにまとめる機能も兼ねている
!NaNの時でも対応できるように、読み込む変数の数を調整したコードも書いてある(コメントアウト済み)
!等値面をparaviewで可視化したいときは、particlesファイルで出力する。
!正直等値面表示以外もparticlesでできるので、txtファイルへの出力する機能はいらないかもしない

program main
implicit none
integer i,j,k,l,time
integer,parameter :: Nx = 360+1 !一行改行分を含めてる
integer,parameter :: Ny = 200
integer,parameter :: Nz = 20
integer,parameter :: CNT = 15356
double precision,allocatable,dimension(:,:,:,:) :: dim_kukei
character(len = 16) z_name
character(len = 16) folder_name
character(len = 16) M_number

allocate(dim_kukei(0:Nx,0:Ny,0:Nz-1,14))
write(folder_name, '(a)') "1.4_random_10"!ここは長すぎるとエラーになる
write(M_number, '(i6.6)') CNT
do k=0,Nz-1
  write(z_name, '(i2.2)') k
  ! open(11, file = "/Users/isseyshome/Documents/GitHub/Inoue_Lab/3dimensions/result_"&
  !                     //trim(folder_name)//"/parameter"//trim(M_number)//"_"//trim(z_name)//".txt")
  ! open(11, file = "/Users/isseyshome/Dropbox/Lab/3-dimensions/result_"&
  !                     //trim(folder_name)//"/parameter"//trim(M_number)//"_"//trim(z_name)//".txt")
  open(11, file = "/Users/isseyshome/Downloads/result_"&
                      //trim(folder_name)//"/parameter"//trim(M_number)//"_"//trim(z_name)//".txt")
  !入力データの読み込みと確認
  do i=0,Ny
    do j = 0,Nx
      !通常はこちら
      ! read(11,'(14f24.16)') dim_kukei(j,i,k,1:14)
      !NaN直前ファイルはこちら
      read(11,'(12f24.16)') dim_kukei(j,i,k,1:12)
    end do
  end do
  close(11)
enddo

!paraviewでの等値面可視化用ファイル
!複数の2次元txtファイルから第二不変量/音響成分/密度と座標のみを1つのparticlesファイルにまとめる時はこちら
open(97, file = "/Users/isseyshome/Paper/result_analysis/"&
          //trim(folder_name)//"_2ndInvaritant"//trim(M_number)//".particles")
open(98, file = "/Users/isseyshome/Paper/result_analysis/"&
          //trim(folder_name)//"_div_u"//trim(M_number)//".particles")
open(99, file = "/Users/isseyshome/Paper/result_analysis/"&
          //trim(folder_name)//"_rho"//trim(M_number)//".particles")
!複数の2次元txtファイルを1つのtxtファイルにまとめる時はこちら
! open(99, file = "/Users/isseyshome/Paper/result_analysis/"&
!           //trim(folder_name)//"_All"//trim(M_number)//".txt")
do k = 0,Nz-1
  do i =0,Ny
    do j =0,Nx
      if((dim_kukei(j,i,k,1) == 0.d0) .and. (dim_kukei(j,i,k,3) == 0.d0) .and.&
       (dim_kukei(j,i,k,5) == 0.d0) .and. (dim_kukei(j,i,k,7) == 0.d0)) then
      else
!(通常&NaN直前)第二不変量のみの時のparticlesファイル出力はこちら
        write(97,'(4f24.16)') &
        dim_kukei(j,i,k,1),dim_kukei(j,i,k,3),dim_kukei(j,i,k,5),&
        dim_kukei(j,i,k,11)!x,y,z,第二不変量Q(2020.08.27~出力順を変更した)
!(通常&NaN直前)div_uのみの時のparticlesファイル出力はこちら
        write(98,'(4f24.16)') &
        dim_kukei(j,i,k,1),dim_kukei(j,i,k,3),dim_kukei(j,i,k,5),&
        dim_kukei(j,i,k,9)!x,y,z,div_u
!(通常&NaN直前)rhoのみの時のparticlesファイル出力はこちら
        write(99,'(4f24.16)') &
        dim_kukei(j,i,k,1),dim_kukei(j,i,k,3),dim_kukei(j,i,k,5),&
        dim_kukei(j,i,k,7)!x,y,z,rho


!(通常)全ての変数を1つのtxtファイルにまとめる時はこちら
        ! write(99,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,&
        ! ",",f24.16,",",f24.16)') &
        ! dim_kukei(j,i,k,1),dim_kukei(j,i,k,3),dim_kukei(j,i,k,5),&
        ! dim_kukei(j,i,k,7),dim_kukei(j,i,k,9),dim_kukei(j,i,k,11),&
        ! dim_kukei(j,i,k,13)

!(NaN直前)全ての変数を1つのtxtファイルにまとめる時はこちら
        ! write(99,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,&
        ! ",",f24.16)') &
        ! dim_kukei(j,i,k,1),dim_kukei(j,i,k,3),dim_kukei(j,i,k,5),&
        ! dim_kukei(j,i,k,7),dim_kukei(j,i,k,9),dim_kukei(j,i,k,11)
      endif
    end do
  end do
enddo
close(99)
deallocate(dim_kukei)
end program main
