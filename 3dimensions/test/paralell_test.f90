program main
  implicit none

  integer,parameter :: Ny = 2000000
  integer,parameter :: Nz = 20
  double precision,parameter :: b = 1.d0!Jet半径は1で固定してしまう
  double precision,parameter :: Cy = 8.d0*b !y軸の幅の設定
  double precision,parameter :: Wry = 2.d0*b!Buffer領域y方向右側の幅
  double precision,parameter :: Wly = Wry!Buffer領域y方向左側の幅
  double precision,parameter :: Ly = 2.d0*Cy+Wry+Wly!y方向の長さを定義 計算領域がy軸対称なのでCyは*2にしている
  double precision,parameter :: Temp = 1.d0
  double precision,parameter :: Tjet = 1.4d0*Temp
  double precision,parameter :: ujet = 1.d0
  double precision,parameter :: Ma = 2.4d0
  double precision,parameter :: gamma = 1.4d0
  double precision,allocatable,dimension(:) :: ur,zeta_fy,Tu,y,com
  double precision,allocatable,dimension(:,:,:) :: in_G
  double precision dy,Ymin,width,a1,a2,b1,tc0,tc1
  integer i,t0,t1,tr,k
  call system_clock(t0)
  call cpu_time(tc0)
  width=3.d0;a1=1d0/14d0;a2=7d0;b1=1.d0/1.4d0
  allocate(ur(0:Ny),zeta_fy(0:Ny),Tu(0:Ny),y(0:Ny),in_G(0:3,0:Ny,0:Nz-1),com(0:Ny))
  ur=0.d0;Tu=0.d0;in_G=0.d0;com=0.d0

  dy = Ly /dble(Ny)
  ur(Ny/2) = ujet
  Tu(Ny/2) = Tjet

  Ymin = -(Ly/2.d0)
  !$omp parallel do
    do i= 0,Ny
      zeta_fy(i) = (b1 * ((1.7d0*(Ymin + dy*dble(i))) - a1 * &
    (-dlog(dcosh(a2*((Ymin + dy*dble(i)) - width))) + dlog(dcosh(a2*((Ymin + dy*dble(i)) + width))))))
    enddo
  !$omp end parallel do

  !====urとTuは同時に並列で求めることができなく、Doループを別にしないと並列化できないので、しない====
  !$omp parallel do
    do i = (Ny/2)+1,Ny
      !Top-hat型の分布になるような式を設定
      ur(i) = ujet/2.d0*(1.d0 - dtanh((12.5d0/4.d0)*((zeta_fy(i)/b)- (b/zeta_fy(i)))))
    enddo
  !$omp end parallel do
  !$omp parallel do
    do i = (Ny/2)+1,Ny
      !Crocco-Busemanの関係式と主流速度分布を用いて温度分布Tuを設定
      Tu(i) = Ma**2.d0*(gamma-1.d0)/2.d0*(ur(i)*ujet-ur(i)**2.d0)/ujet+&
              Tjet*ur(i)/ujet+Temp*(ujet-ur(i))/ujet
    enddo
  !$omp end parallel do
  !====並列化しない================================================================

  !初期分布をx軸対象になるようにする。
  !そのためにy軸正の範囲の値を負の範囲に軸対象になるようにコピーする
  !$omp parallel do
    do i = 0,(Ny/2)-1
      ur(i) = ur(Ny-i)
      Tu(i) = Tu(Ny-i)
    enddo
  !$omp end parallel do
  !$omp parallel do
    do k =0,Nz-1
      do i=0,Ny
        in_G(0,i,k) = 1.d0/Tu(i)!密度ρは理想気体状態方程式に従うから
        in_G(1,i,k) = ur(i)! NSCBC流入条件で使う流入値のみを保存する配列
        in_G(2,i,k) = 0.d0! NSCBC流入条件で使う流入値のみを保存する配列
        in_G(3,i,k) = 0.d0! NSCBC流入条件で使う流入値のみを保存する配列
      end do
    enddo
  !$omp end parallel do
  !$omp parallel do
    do i=0,Ny
      com(i) = in_G(1,i,0) -ur(i)
    enddo
  !$omp end parallel do
  ! open(10, file="tu_o_do.csv")
  !   do i=0,Ny
  !     write(10,*) ur(i),",",Tu(i)
  !   end do
  ! close(10)
call cpu_time(tc1)
call system_clock(t1,tr)
write(*,'(f10.3,A)')(t1-t0)/dble(tr),"[s],",tc1-tc0,"[s]"
end program main
