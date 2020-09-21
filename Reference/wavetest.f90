program iryuu6
implicit none
integer i,M,j,k,m1,im1,im2,ip1,ip2
double precision dz,dt,pi,Z_LU,s_1,s_2,s_3,s_4,a,a3,ad,ad3,dzinv
double precision alpha,alpha3,alphaD,b,bd,sigma,cd,dxinv,alpha5,Dalpha,usum,lsum
integer,parameter::Nz=100

double  precision,parameter :: ra = 14.d0/9.d0, rb = 1.d0/9.d0&
&,da = 4.d0 / 9.d0,db = 2.d0 / 9.d0 !5次精度のDCSとなるための係数設定

double precision,allocatable::L(:,:),U(:,:),AA(:,:),LU(:,:),f(:)
double precision,allocatable::y(:),df(:),Q(:),z(:),QR(:),RHS_z(:),x(:)
double precision,allocatable::D2(:),D4(:),D6(:),D8(:)

allocate(QR(0:NZ-1),f(0:NZ-1),y(0:NZ-1),Q(0:NZ-1),z(0:NZ-1),df(0:NZ-1),RHS_z(0:NZ-1))
allocate(AA(0:NZ-1,0:NZ-1),L(0:NZ-1,0:NZ-1),U(0:NZ-1,0:NZ-1))
allocate(LU(-2:2,0:Nz-1))
allocate(D2(0:NZ-1),D4(0:NZ-1),D6(0:NZ-1),D8(0:NZ-1))
allocate(x(0:NZ-1))

pi=dacos(-1d0)
dz=1d0/dble(NZ)

AA=0d0
L=0d0
U=0d0
LU=0d0
Q=0d0;f=0d0;y=0d0;x=0d0
alpha5 = 1.d0/3.d0;Dalpha = 1.d0;sigma=025.d0
D2=0.d0;D4=0.d0;D6=0.d0;D8=0.d0
!!!--------------------------------------A行列

!周期条件なので5次精度DCSを使う(sigmaを0にすれば6次CCSとなる)
!0行目
AA(0,0) = 1.d0
AA(0,1) = alpha5*(1.d0+sigma*Dalpha)
AA(0,NZ-1) = alpha5*(1.d0-sigma*Dalpha)
!1からN-2行目まで
do i = 1,NZ-2
  AA(i,i-1) = (1.d0 - Dalpha * sigma) * alpha5
  AA(i,i) = 1.d0
  AA(i,i+1) = (1.d0 + Dalpha * sigma) * alpha5
enddo
!N行目
AA(NZ-1,0) = alpha5 * (1.d0 + Dalpha * sigma)
AA(NZ-1,NZ-2) = alpha5*(1.d0-sigma*Dalpha)
AA(NZ-1,NZ-1) = 1.d0

!!!---------------------------------------LU分解
do i = 0,NZ-1
  U(0,i) = AA(0,i)
  L(i,0) = AA(i,0) / U(0,0)
  L(i,i) = 1.d0
enddo
!Uは行ごとに、Lは列ごとに求めていく。
!ただしUの２行目、 Lの２列目、Uの3行目、Lの３列目といった順番
do i = 1,NZ-1 !初期条件の結果を利用してUのi列、Lのi行の順に求めていく
  do j = i,NZ-1 !Uのi行の列要素を求めていく。i<=jを考慮する
    Usum = 0.d0
    do k = 0, i-1
      Usum = Usum + L(i,k) * U(k,j)
    enddo
    U(i,j) = AA(i,j) - Usum
  enddo

  do j = i+1,NZ-1 !Lのi列の行要素を求めていく。j<iを考慮
    Lsum = 0.d0
    do k = 0,i-1
      Lsum = Lsum + L(j,k) * U(k,i)
    enddo
    L(j,i) = (AA(j,i) - Lsum) / U(i,i)
  enddo
enddo

!圧縮=========================================================
LU(-2,NZ-2:NZ-1) = 0.d0
LU(-1,0) = 0.d0
LU(0,NZ-1) = U(NZ-1,NZ-1)
LU(1,NZ-1) = 0.d0
LU(2,NZ-2:NZ-1) = 0.d0
do i= 0,NZ-3
  LU(-2,i) = U(i,Nz-1)
  LU(2,i)  = L(Nz-1,i)
enddo
do i = 1,NZ-1
  LU(-1,i) = L(i,i-1)
enddo
do j = 0,1
  do i = 0,NZ-2
    LU(j,i) = U(i,i+j)
  enddo
enddo


do m=0,Nz-1
	z(m)=dz*dble(m)
	f(m)=sin(2.d0*pi*z(m))
end do

!do M=0,500

dzinv = 1.d0/dz
  !片側DCSの右辺設定
  do i = 0,Nz-1
    im1 = mod(Nz+i-1, Nz)
    im2 = mod(Nz+i-2, Nz)
    ip1 = mod(Nz+i+1, Nz)
    ip2 = mod(Nz+i+2, Nz)
    D2(i) = (-F(im1)+F(ip1)) * (0.5d0*dzinv)
    D4(i) = (-F(im2)+F(ip2)) * (0.25d0*dzinv)

    D6(i) = (F(im1)+F(ip1)- 2.d0* F(i)) * dzinv
    D8(i) = (F(im2)+F(ip2)- 2.d0* F(i)) * (0.25d0*dzinv)
  RHS_z(i)=ra*D2(i)+rb*D4(i)+sigma*(da*D6(i)+db*D8(i))
  enddo
     y(0) = RHS_z(0)
     do i = 1,Nz-2
       y(i) = RHS_z(i) - LU(-1,i)*y(i-1)
     enddo


     Lsum = 0.d0
     do i = 0,Nz-3
       Lsum = Lsum +LU(2,i)*y(i)
     enddo
       y(Nz-1) = RHS_z(Nz-1) - LU(-1,Nz-1)*y(Nz-2)-Lsum

     x(Nz-1) = y(Nz-1) / LU(0,Nz-1)
     x(Nz-2) = (y(Nz-2) -LU(1,Nz-2)*x(Nz-1))/ LU(0,Nz-2)
     do i = Nz-3, 0, -1
       x(i) = (y(i) - LU(1,i)*x(i+1)-LU(-2,i)*x(Nz-1)) / LU(0,i)
     enddo

open(10,file ='testz.csv')
do m=0,Nz-1
write(10,*) z(m),',',f(m),',',x(m)/(2.d0*pi)
end do
close(10)
end program iryuu6
