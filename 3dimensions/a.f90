program main
  implicit none

  integer,parameter :: Ny = 4
  integer,parameter :: Nz = 2
  double precision,parameter :: b = 1.d0!Jet半径は1で固定してしまう
  double precision,parameter :: Cy = 8.d0*b !y軸の幅の設定
  double precision,parameter :: Wry = 2.d0*b!Buffer領域y方向右側の幅
  double precision,parameter :: Wly = Wry!Buffer領域y方向左側の幅
  double precision,parameter :: Ly = 2.d0*Cy+Wry+Wly!y方向の長さを定義 計算領域がy軸対称なのでCyは*2にしている
  double precision,allocatable,dimension(:) :: ur,yu
  double precision,allocatable,dimension(:,:) :: Tu
  double precision dy,aaa,aa1a
  integer i,k
  allocate(ur(0:Ny),Tu(0:Ny,0:Nz),yu(0:Nz))
  ur=0.d0;Tu=0.d0;yu=0.d0

  yu(0) = 0.d0
  yu(1) = 1.d0
  yu(2) = 3.d0
  Tu(0,0) = 2.d0
  TU(1,0) = 3.d0
  TU(2,0) = 5.d0
  TU(3,0) = 7.d0
  TU(4,0) = 9.d0
  Tu(0,1) = 11.d0
  TU(1,1) = 13.d0
  TU(2,1) = 17.d0
  TU(3,1) = 19.d0
  TU(4,1) = 23.d0
  Tu(0,2) = 29.d0
  TU(1,2) = 31.d0
  TU(2,2) = 37.d0
  TU(3,2) = 41.d0
  TU(4,2) = 43.d0
  aaa =0.d0
  !$omp parallel do
  do k =0,Nz
    write(*,*) "1",aaa
    aaa = aaa + Tu(1,k) * yu(k)
    write(*,*) "2",aaa
  enddo
  !$omp end parallel do
  write(*,*) "last",aaa
end program main
