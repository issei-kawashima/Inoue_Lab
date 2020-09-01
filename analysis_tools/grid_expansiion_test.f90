program main
  integer Nx,i
  double precision,allocatable,dimension(:) :: dzeta,dzeta_in
  double precision,allocatable,dimension(:) :: zeta,ygs,dygs,dygs_in
  double precision dx,x,Lx
  double precision,parameter :: b1 = 10.8d0
  double precision,parameter :: b2 = 10.8d0
  real,parameter :: a1 = 7.d0;a2=14.d0
  Nx=180;Lx=36.d0;dx=Lx/dble(Nx)
  allocate(dzeta(0:Nx),dzeta_in(0:Nx),zeta(0:Nx),ygs(0:Nx),dygs(0:Nx),dygs_in(0:Nx))
  dzeta=0.d0;dzeta_in=0.d0;zeta=0.d0;ygs=0.d0;dygs_in=0.d0;dygs=0.d0
  do i= 0,Nx
    x =dx*dble(i)
    zeta(i) = (1.d0/1.4d0) * ((1.7d0*x) - (1.d0/a2) * &
    (-dlog(dcosh(a1*(x - b1))) + dlog(dcosh(a1*(x + b2)))))
    dzeta(i) = (1.d0/1.4d0) * (1.7d0 - (a1/a2) * &
    (-dtanh(a1*(x - b1)) + dtanh(a1*(x + b2))))
  enddo
  dzeta_in(:) = 1.d0/dzeta(:)
  open(1,file ="xzeta.csv")
  open(2,file ="xdzeta.csv")
  do i =0,Nx
    x =dx*dble(i)
    write(1,*) x,",",zeta(i)
    write(2,*) x,",",dzeta(i)
  enddo
  close(1);close(2)
  write(*,*) a1,a2,zeta(Nx),dzeta(0),dzeta(Nx/2)
endprogram
