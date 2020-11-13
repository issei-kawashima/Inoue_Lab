!格子伸長の効き具合を調整するコード
!a1,a2を調整する事で、元のLzなどの長さになるように調整する。
!したがってLzなどを変更したらこのコードでパラメータの調整をしなくてはいけない
!widthが狭ければ格子幅を狭める倍率は0.1くらいに小さくなるが、widthが広ければ0.5倍などがせいぜいになる。
program main
  integer Nx,i
  double precision,allocatable,dimension(:) :: dzeta,dzeta_in
  double precision,allocatable,dimension(:) :: zeta
  double precision dx,x,Lx
  double precision,parameter :: width = 5.d-1
  ! double precision p1,p2,p3,bai1,bai2
  ! bai1=1.4d0;bai2=1.7/bai1;p1=-0.5d0;p2=5.d0
  double precision,parameter :: a1 = 5d0
  double precision,parameter :: a2= 20.d0/3.d0
  Nx=50;Lx=5.d0;dx=Lx/dble(Nx)
  allocate(dzeta(0:Nx),dzeta_in(0:Nx),zeta(0:Nx))
  dzeta=0.d0;dzeta_in=0.d0;zeta=0.d0
  do i= 0,Nx
    x =-Lx/2.d0+dx*dble(i)
    zeta(i) = (1.d0/1.4d0) * ((1.7d0*x) - (1.d0/a2) * &
    (-dlog(dcosh(a1*(x - width))) + dlog(dcosh(a1*(x + width)))))
    dzeta(i) = (1.d0/1.4d0) * (1.7d0 - (a1/a2) * &
    (-dtanh(a1*(x - width)) + dtanh(a1*(x + width))))
    ! zeta(i) = bai2*x+(1.d0/bai1)*(p1/p2)*(-dlog(dcosh(p2*(x-width)))&
    ! +dlog(dcosh(p2*(x+width))))
    ! dzeta(i) = bai2+p1/bai1*(-dtanh(p2*(x-width))+dtanh(p2*(x+width)))
  enddo
  dzeta_in(:) = 1.d0/dzeta(:)
  open(1,file ="z_zeta.csv")
  open(2,file ="z_dzeta.csv")
  do i =0,Nx
    x =-Lx/2.d0+dx*dble(i)
    write(1,*) x,",",zeta(i)
    write(2,*) x,",",dzeta(i)
  enddo
  close(1);close(2)
  write(*,*) a1,a2,zeta(0),zeta(Nx),dzeta(0),dzeta(Nx/2)
endprogram
