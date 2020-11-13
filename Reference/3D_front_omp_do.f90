!!!!==================================================================
!!!3次元化およびNSCBC法
!!!森山史孝　2017/11/06
!!!===================================================================
!!!x方向格子伸長未実装
!!!平均化および計算手法を省略化

module sub_mod
  implicit none
  integer,parameter::NX=500+1,NY=200+1,NZ=40
  integer,parameter::NT=300,MMM=650*NT
  integer,parameter::MM=450*NT,M1=250*NT
  integer,parameter::Mmod=NT*10,Mmod1=NT
  integer,parameter::Mlarge=0,n0=0,NTT=120
  integer,parameter::Ndt=61,loop=NT*2!18
  double precision,parameter::pi=dacos(-1d0)
  double precision,parameter::RE=1250d0,Pr=0.71d0
  double precision,parameter::zeta=1d0,sig=0.27d0,gamma=1.4d0
  double precision,parameter::S=120d0,Ma=0.3d0,Tc=S/(291.15d0)
  integer,parameter::bx2=34,by=5,bx1=8
  integer,parameter::Mk=20*NT
!!!Maごとにsigmaを設定すること
 !double precision,parameter::sigma=0.025d0
  double precision,parameter::sigma_x=0.25d0,sigma_y=0.25d0
  double precision,parameter::Lx=150d0,Ly=70d0,Lyb=70d0,Lya=0d0,Lz=4.5d0
  double precision,parameter::Lxr=75d0,Lxl=-75d0,B1=5d0
  double precision,parameter::dx=Lx/dble(NX-1),dy=Ly/dble(NY-1)
  double precision,parameter::dz=Lz/dble(NZ),dt=1d0/dble(NT)
  double precision,parameter::gdx=1d0/dx
  double precision,parameter::gdy=1d0/dy
  double precision,parameter::gdz=1d0/dz
  double precision,parameter::rho0=1d0,T0=1d0
  double precision,parameter::p0=T0*rho0/(ma*ma*gamma),p_inf=p0
  double precision,parameter::c_inf=dsqrt(1d0/(Ma**2))
  double precision,parameter::a1=0.15d0
  integer,parameter::kizonsiyou=1
  integer,parameter::data2dim=1,data3dim=1,randim=0,zydanmen=1,xydanmen=1
  double precision,parameter::a111=-17d0/6d0,a112=1.5d0,a113=1.5d0,a114=-1d0/6d0
  double precision,parameter::a121=-21d0/6d0,a122=0.5d0,a123=3.5d0,a124=-0.5d0
  double precision,parameter::a51=7d0/9d0,a52=1d0/36d0,b51=4d0/9d0,b52=1d0/18d0
  double precision,parameter::a31=0.75d0,b31=0.5d0
  double precision,parameter::tau1=2d0/3d0,gRE=1d0/RE,Kq=1d0/((gamma-1d0)*RE*Pr*Ma**2)
  double precision,parameter::tvd1=1d0/3d0,tvd2=2d0/3d0,g6=1d0/6d0
  integer,parameter::xk=210,yk=56
  integer,parameter::hx=6,hz=6
  double precision,parameter::parav=0.50d1
  integer,parameter::Nz2=NZ/2,hz2=hz/2
  double precision,parameter::gmamaga=1d0/(gamma*ma*ma)

contains
!!!!====================================================================================

  subroutine zahyou_henkan(x_k,z_k,y,y_k,dx_k,dy_k) !!!初期座標および格子伸長
    implicit none
    double precision,dimension(:)::y_k,dy_k,y !!格子点関数
    double precision,dimension(:)::x_k,dx_k
    double precision,allocatable,dimension(:)::x
    double precision,dimension(:)::z_k
    integer::i,j,k
	allocate(x(NX))
    do i=1,NX
       x(i)=dx*(dble(i-1))+Lxl
    end do

    do i=1,NZ
       z_k(i)=dz*(dble(i-1))
    end do
!dx_k=1d0
    do i=1,NY
       y(i)=dy*dble(i-1)+Lya
    end do
    y_k = (Lyb)*dexp(-a1*(Lyb-y)) - (Lyb - y)*dexp(-a1*(Lyb))+y/20d0!tanh((y/10d0)**3)*5d0
    y_k = Y_k-Y_k(1)
  ! y_k(1)=0d0
    dy_k =1d0/(a1*(Lyb)*dexp(-a1*(Lyb-y)) +dexp(-a1*(Lyb))+1d0/20d0 )! +(1-(tanh((y/10d0)**3))**2)*5d0*3d0*(y/10d0)**2 )
    !Y = Y_a + ( y_b - y_a )*dexp(-a1*( y_b - y_a - Ys )) - ( y_b - y_a - Ys )*dexp(-a1*( y_b - y_a ))
    !Y(1) = 0.0d0
!    y_k=y
!    dy_k=1d0
	x_k=(X/Lxl)**3*Lxl+(X/B1)
	dx_k=1d0/(3d0*(1d0/Lxl)**3d0*X**2*Lxl+1d0/B1)

    !    dy_k=1d0
    open(10,file='x_grid.txt')
    open(11,file='y_grid.txt')
    open(12,file='z_grid.txt')
    do i=1,NX
 	write(10,*)i,x_k(i)
    end do
    do j=1,NY
        write(11,*)j,y_k(j)
    end do
    do k=1,NZ
        write(12,*)k,z_k(k)*1d0
    end do
    close(10)
    close(11)
    close(12)


  end subroutine zahyou_henkan

!!!=====Qからrho,u,v,pを求める（亜音速流入用）=======
  subroutine Q_g_in(Q_1,Q_2,Q_3,Q_4,Q_5,g_1,g_2,g_3,g_4,g_5,T_k,mu,mugRE,mukq,t,y_k,g_in,vin)
    implicit none
    double precision,dimension(:,:,:)::Q_1,Q_2,Q_3,Q_4,Q_5
    double precision,dimension(:,:,:)::g_1,g_2,g_3,g_4,g_5
    double precision,dimension(:,:,:)::T_k,mu,mukq,mugRE
    double precision,dimension(:)::y_k
    double precision,dimension(:,:)::g_in
    double precision,dimension(1:hx,1:hz)::vin
    integer::i,j,k,t,time,karit,karis

  time=dble(t)/dble(NT*10)
  if(time>=1d0)time=1d0
   time=1d0
   !$omp parallel do
    do k=1,NZ
    	do i=1,NY
	Q_2(i,1,k)=g_in(2,i)*Q_1(i,1,k)
   	Q_3(i,1,k)=g_in(3,i)*Q_1(i,1,k)
    	Q_4(i,1,k)=0d0
    	T_k(i,1,k)=g_in(5,i)
    	g_5(i,1,k)=Q_1(i,1,k)*g_in(5,i)*gmamaga

    	end do
    	do j=1,xk
    	Q_2(1,j,k) =0d0
    	Q_3(1,j,k) =0d0
    	Q_4(1,j,k) =0d0
        T_k(1,j,k) =T_k(2,j,k)
        Q_1(1,j,k) =Q_1(2,j,k)
        g_5(1,j,k) =T_k(1,j,k)*Q_1(1,j,k)*gmamaga
        Q_5(1,j,k) =g_5(1,j,k) /(gamma-1d0)
    	end do
    	do j=xk,NX
   	Q_2(yk,j,k) =0d0
	Q_3(yk,j,k) =0d0
   	Q_4(yk,j,k) =0d0
   	T_k(yk,j,k) =T_k(yk+1,j,k)
        Q_1(yk,j,k) =Q_1(yk+1,j,k)
        g_5(yk,j,k) =T_k(yk,j,k)*Q_1(yk,j,k)*gmamaga
        Q_5(yk,j,k) =g_5(yk,j,k) /(gamma-1d0)
   	end do
	do i=1,yk
        Q_2(i,xk,k) =0d0
        Q_3(i,xk,k) =0d0
        Q_4(i,xk,k) =0d0
        T_k(i,xk,k) =T_k(i,xk-1,k)
        Q_1(i,xk,k) =Q_1(i,xk-1,k)
        g_5(i,xk,k) =T_k(i,xk,k)*Q_1(i,xk,k)*gmamaga
        Q_5(i,xk,k) =g_5(i,xk,k) /(gamma-1d0)
	end do

    Q_1(yk,xk,k)=(Q_1(yk,xk-1,k)*0.4d0+Q_1(yk+1,xk,k)*0.6d0)
    T_k(yk,xk,k)=(T_k(yk,xk-1,k)*0.4d0+T_k(yk+1,xk,k)*0.6d0)
    g_5(yk,xk,k) =T_k(yk,xk,k)*Q_1(yk,xk,k)*gmamaga
    Q_5(yk,xk,k) =g_5(yk,xk,k) /(gamma-1d0)

   end do
	!$omp end parallel do
!!!----------------------------------------------------------------------
!!!hunryujouken
   Q_3(1,bx1+1:bx1+1+hx-1,5+1:hz+5)  =0.95d0*vin(1:hx,1:hz)*Q_1(1,bx1+1:bx1+1+hx-1,5+1:hz+5) *time
!   Q_3(1,bx1+1:bx1+1+hx-1,NZ2-hz2+1:NZ2+hz2)=0.3d0*vin(1:hx,1:hz)*Q_1(1,bx1+1:bx1+1+hx-1,NZ2-hz2+1:NZ2+hz2) *time
   Q_3(1,bx1+1:bx1+1+hx-1,NZ-hz-4:NZ-5)=1.05d0*vin(1:hx,1:hz)*Q_1(1,bx1+1:bx1+1+hx-1,NZ-hz-4:NZ-5) *time

   Q_4(1,bx1+1:bx1+1+hx-1,5+1:hz+5)  =-0.1d0*vin(1:hx,1:hz)*Q_1(1,bx1+1:bx1+1+hx-1,5+1:hz+5) *time
!   Q_4(1,bx1+1:bx1+1+hx-1,NZ2-hz2+1:NZ2+hz2)=-0.d0*vin(:,:)*Q_1(1,bx1+1:bx1+1+hx-1,NZ2-hz2+1:NZ2+hz2) *time
   Q_4(1,bx1+1:bx1+1+hx-1,NZ-hz-4:NZ-5)=0.05d0*vin(1:hx,1:hz)*Q_1(1,bx1+1:bx1+1+hx-1,NZ-hz-4:NZ-5) *time

    T_k(1,bx1+1:bx1+1+hx-1,5+1:hz+5)  =1d0
!    T_k(1,bx1+1:bx1+1+hx-1,NZ2-hz2+1:NZ2+hz2)=1d0
    T_k(1,bx1+1:bx1+1+hx-1,NZ-hz-4:NZ-5)=1d0

!!!-----------------------------------------------------------------
    !$omp parallel do
    do k=1,NZ
    	do j=1,NX
         Q_5(1,j,k)=g_5(1,j,k)/(gamma-1d0)&
         +0.5d0/Q_1(1,j,k)*(Q_2(1,j,k)**2+Q_3(1,j,k)**2+Q_4(1,j,k)**2)
        end do
    	do i=1,NY
    	Q_5(i,1,k)=g_5(i,1,k)/(gamma-1d0)&
        	 +(1d0/(Q_1(i,1,k))*(Q_2(i,1,k)**2+Q_3(i,1,k)**2+Q_4(i,1,k)**2))*0.5d0
        end do

    	do j=1,NX
    		do i=1,NY
    g_1(i,j,k)=Q_1(i,j,k)
    g_2(i,j,k)=Q_2(i,j,k)/Q_1(i,j,k)
    g_3(i,j,k)=Q_3(i,j,k)/Q_1(i,j,k)
    g_4(i,j,k)=Q_4(i,j,k)/Q_1(i,j,k)
    g_5(i,j,k)=(gamma-1d0)*&
         (Q_5(i,j,k)-(Q_2(i,j,k)**2+Q_3(i,j,k)**2+Q_4(i,j,k)**2)/(2d0*Q_1(i,j,k)))
    T_k(i,j,k)=gamma/g_1(i,j,k)*(g_5(i,j,k))*Ma*Ma
    mu(i,j,k)=T_k(i,j,k)**(tvd2)
    mugRE(i,j,k)=mu(i,j,k)*gRE
    mukq(i,j,k)=-mu(i,j,k)*Kq
    		end do
    	end do
    end do
    !$omp end parallel do
  end subroutine Q_g_in

  function fBDCS(n,NN,gdh,F,LUcom,sig), result(fBDCS)!!!nからNNまで計算
    !,resultが大事！
    !返す変数にはintent(in)をつけず、仮引数(function()で読み込む引数)にはintent(in)をつける
    implicit none
    integer::i,j,k,NN,n
    double precision,dimension(n:NN),intent(in)::F
    double precision,dimension(n:NN,1:3),intent(in)::LUcom
    double precision,dimension(n:NN)::fBDCS
    double precision::sig,gdh

    do j=n+2,NN-2
       fBDCS(j)=(a51*(F(j+1)-F(j-1))&
            +a52*(F(j+2)-F(j-2))&
            +sig*(b51*(F(j+1)-2.d0*F(j)+F(j-1))&
            +b52*(F(j+2)-2.d0*F(j)+F(j-2))))*gdh
    end do
    if (sig>0d0)then
       fBDCS(n)=(-17d0/6d0*F(n)+1.5d0*F(n+1)&
            +1.5d0*F(n+2)-1.d0/6d0*F(n+3))*gdh
       fBDCS(NN)=-((-21d0/6d0*F(NN)+0.5d0*F(NN-1)&
            +3.5d0*F(NN-2)-0.5d0*F(NN-3)))*gdh
    else
       fBDCS(n)=(-21d0/6d0*F(n)+0.5d0*F(n+1)&
            +3.5d0*F(n+2)-0.5d0*F(n+3))*gdh
       fBDCS(NN)=-((-17d0/6d0*F(NN)+1.5d0*F(NN-1)&
            +1.5d0*F(NN-2)-1.d0/6d0*F(NN-3)))*gdh

    end if
   !DCS_3
    fBDCS(n+1)=(a31*(F(n+2)-F(n))&
         +sig*b31*(F(n)-2d0*F(n+1)+F(n+2)))*gdh
    fBDCS(NN-1)=(a31*(F(NN)-F(NN-2))&
         +sig*b31*(F(NN-2)-2d0*F(NN-1)+F(NN)))*gdh

    do j=n+1,NN
       fBDCS(j)=fBDCS(j)-LUcom(j-1,1)*fBDCS(j-1)
    end do
    fBDCS(NN)=fBDCS(NN)/LUcom(NN,2)
    do j=NN-1,n,-1
       fBDCS(j)=(fBDCS(j)-LUcom(j,3)*fBDCS(j+1))/LUcom(j,2)
    end do
  end function fBDCS
!!!===========================================================
!!!CCS関数の定義
!!!===========================================================
  function CCS(n,NN,gdh,F,LUc)
    implicit none
    integer,intent(in)::NN,n
    integer::i,j
    double precision,dimension(n:NN),intent(in)::F
    double precision,dimension(n:NN,1:3),intent(in)::LUc
    double precision,dimension(n:NN)::CCS
    double precision::gdh
    !CCS_6次
    do j=n+2,NN-2
       CCS(j)=(a51*(F(j+1)-F(j-1))&
            +a52*(F(j+2)-F(j-2)))*gdh
    end do

!!!片側CCS
    CCS(n)=(a111*F(n)+a112*F(n+1)&
         +a113*F(n+2)+a114*F(n+3))*gdh
    CCS(NN)=-((a111*F(NN)+a112*F(NN-1)&
         +a113*F(NN-2)+a114*F(NN-3)))*gdh
!!!CCS_4
    CCS(n+1)=a31*(F(n+2)-F(n))*gdh
    CCS(NN-1)=a31*(F(NN)-F(NN-2))*gdh

    do j=n+1,NN
       CCS(j)=(CCS(j)-LUc(j-1,1)*CCS(j-1))
    end do

    CCS(NN)=CCS(NN)/LUc(NN,2)
    do j=NN-1,n,-1
       CCS(j)=(CCS(j)-LUc(j,3)*CCS(j+1))/LUc(j,2)
    end do
  end function CCS

  !!!-------------------------------------------------------------
  function DCSz(F,L,U,sig)!!!nからNNまで計算
    implicit none
    integer::i,j,k
    double precision,dimension(:,:)::L,U
    double precision,dimension(:)::F
    double precision,allocatable,dimension(:)::DCSz
    double precision::sig,sum
    allocate(DCSz(NZ))
    sum=0d0

    do j=3,NZ-2
       DCSz(j)=(a51*(F(j+1)-F(j-1))&
            +a52*(F(j+2)-F(j-2))&
            +sig*(b51*(F(j+1)-2.d0*F(j)+F(j-1))&
            +b52*(F(j+2)-2.d0*F(j)+F(j-2))))*gdz
    end do

       DCSz(1)=(a51*(F(2)-F(NZ))&
            +a52*(F(3)-F(NZ-1))&
            +sig*(b51*(F(2)-2.d0*F(1)+F(NZ))&
            +b52*(F(3)-2.d0*F(1)+F(NZ-1))))*gdz

       DCSz(2)=(a51*(F(3)-F(1))&
            +a52*(F(4)-F(NZ))&
            +sig*(b51*(F(3)-2.d0*F(2)+F(1))&
            +b52*(F(4)-2.d0*F(2)+F(NZ))))*gdz

       DCSz(NZ-1)=(a51*(F(NZ)-F(NZ-2))&
            +a52*(F(1)-F(NZ-3))&
            +sig*(b51*(F(NZ)-2.d0*F(NZ-1)+F(NZ-2))&
            +b52*(F(1)-2.d0*F(NZ-1)+F(NZ-3))))*gdz

       DCSz(NZ)=(a51*(F(1)-F(NZ-1))&
            +a52*(F(2)-F(NZ-2))&
            +sig*(b51*(F(1)-2.d0*F(NZ)+F(NZ-1))&
            +b52*(F(2)-2.d0*F(NZ)+F(NZ-2))))*gdz

    do i=2,NZ
       DCSz(i)=-DCSz(i-1)*L(i,i-1)+DCSz(i)
    end do

    do i=1,NZ-2
       DCSz(NZ)=DCSz(NZ)-DCSz(i)*L(NZ,i)
    end do

    DCSz(NZ)=DCSz(NZ)/U(NZ,NZ)
    DCSz(NZ-1)=(DCSz(NZ-1)-DCSz(NZ)*U(NZ-1,NZ))/U(NZ-1,NZ-1)
    do i=NZ-2,1,-1
       DCSz(i)=(-DCSz(NZ)*U(i,NZ)-DCSz(i+1)*U(i,i+1)+DCSz(i))/(U(i,i))
    end do
  end function DCSz
!!!===========================================================
!!!CCS関数の定義
!!!===========================================================
  function CCSz(F,L,U)
    implicit none
    integer::i,j,k,n,n1
    double precision,dimension(:),intent(in)::F
    double precision,allocatable,dimension(:)::CCSz
    double precision,dimension(:,:),intent(in)::L,U
    double precision::sum
    allocate(CCSz(NZ))
    sum=0d0
    !CCS_6次
    do j=3,NZ-2
       CCSz(j)=(a51*(F(j+1)-F(j-1))&
            +a52*(F(j+2)-F(j-2)))*gdz
    end do

    CCSz(1)=(a51*(F(2)-F(NZ))&
         +a52*(F(3)-F(NZ-1)))*gdz

    CCSz(2)=(a51*(F(3)-F(1))&
         +a52*(F(4)-F(NZ)))*gdz

    CCSz(NZ-1)=(a51*(F(NZ)-F(NZ-2))&
         +a52*(F(1)-F(NZ-3)))*gdz

    CCSz(NZ)=(a51*(F(1)-F(NZ-1))&
         +a52*(F(2)-F(NZ-2)))*gdz

    CCSz(1)=CCSz(1)
    do i=2,NZ
       CCSz(i)=-CCSz(i-1)*L(i,i-1)+CCSz(i)
    end do

    do i=1,NZ-2
       CCSz(NZ)=CCSz(NZ)-CCSz(i)*L(NZ,i)
    end do

    CCSz(NZ)=CCSz(NZ)/U(NZ,NZ)
    CCSz(NZ-1)=(CCSz(NZ-1)-CCSz(NZ)*U(NZ-1,NZ))/U(NZ-1,NZ-1)
    do i=NZ-2,1,-1
       CCSz(i)=(-CCSz(NZ)*U(i,NZ)-CCSz(i+1)*U(i,i+1)+CCSz(i))/(U(i,i))
    end do
   end function CCSz

  subroutine NSCBC_xhoukou(g_1,g_2,g_3,g_4,g_5,dxg_1,dxg_2,dxg_3,dxg_4,dxg_5,dF1,dF2,dF3,dF4,dF5,c)
    implicit none
    integer ::i,j,k,t
    double precision,dimension(:,:,:),intent(inout)::g_1,g_2,g_3,g_4,g_5
    double precision,dimension(:,:,:),intent(inout)::dxg_1,dxg_2,dxg_3,dxg_4,dxg_5
    double precision,dimension(:,:,:),intent(inout)::dF1,dF2,dF3,dF4,dF5
    double precision,allocatable,dimension(:,:,:,:)::L,d  !L＿番号＿NX＿左側と右側

    double precision,dimension(:,:,:),intent(in)::c
    double precision,allocatable,dimension(:,:)::Ma1
    allocate(L(1:5,NY,1,NZ),d(1:5,NY,1,NZ),Ma1(NY,NZ))

!!!亜音速流入条件
 !$omp parallel do
	do k=1,NZ
	    do i=1,NY
    Ma1(i,k)=g_2(i,NX,k)/c(i,NX,k)
    L(1,i,1,k)=(g_2(i,1,k)-c(i,1,k))&
         *(-g_1(i,1,k)*c(i,1,k)*dxg_2(i,1,k)+dxg_5(i,1,k))
    L(5,i,1,k)=L(1,i,1,k)
    L(2,i,1,k)=0.5d0*(gamma-1d0)*(L(1,i,1,k)+L(5,i,1,k))
    !            &rho(:,0)*c(:,0)**2d0/T(:,0)*dT/d0=0d0
!!!流入d行列の計算i=0
    dF1(i,1,k)=1d0/(c(i,1,k)**2d0)*(L(1,i,1,k)+L(2,i,1,k))

!!!=========================================
    L(:,i,1,k)=0d0
end do
!!!流入d行列の計算
	do i=yk+1,NY
    L(1,i,1,k)=sigma_x*(1d0-Ma1(i,k)**2d0)*c(i,NX,k)*(g_5(i,NX,k)-p_inf)/Lx
 !   L(1,i,1,k)=sigma*(1d0-Ma**2d0)*c(i,NX,k)*(g_5(i,NX,k)-p_inf)  !/Lx
    L(2,i,1,k)=g_2(i,NX,k)*(c(i,NX,k)**2d0*dxg_1(i,NX,k)-dxg_5(i,NX,k))
    L(3,i,1,k)=g_2(i,NX,k)*dxg_3(i,NX,k)
    L(4,i,1,k)=g_2(i,NX,k)*dxg_4(i,NX,k)
    L(5,i,1,k)=(g_2(i,NX,k)+c(i,NX,k))*(g_1(i,NX,k)*c(i,NX,k)&
                  *dxg_2(i,NX,k)+dxg_5(i,NX,k))

    d(1,i,1,k)=1d0/(c(i,NX,k)**2d0)*(0.5d0*(L(1,i,1,k)+L(5,i,1,k))+L(2,i,1,k))
    d(2,i,1,k)=0.5d0*(L(1,i,1,k)+L(5,i,1,k))
    d(3,i,1,k)=0.5d0/(g_1(i,NX,k)*c(i,NX,k))*(-L(1,i,1,k)+L(5,i,1,k))
    d(4,i,1,k)=L(3,i,1,k)
    d(5,i,1,k)=L(4,i,1,k)
!!!i=NX
    dF1(i,NX,k)=d(1,i,1,k)
    dF2(i,NX,k)=g_2(i,NX,k)*d(1,i,1,k)+g_1(i,NX,k)*d(3,i,1,k)
    dF3(i,NX,k)=g_3(i,NX,k)*d(1,i,1,k)+g_1(i,NX,k)*d(4,i,1,k)
    dF4(i,NX,k)=g_4(i,NX,k)*d(1,i,1,k)+g_1(i,NX,k)*d(5,i,1,k)
    dF5(i,NX,k)=0.5d0*(g_2(i,NX,k)**2+g_3(i,NX,k)**2+g_4(i,NX,k)**2)&
         *d(1,i,1,k)+d(2,i,1,k)/(gamma-1d0)+&
             g_1(i,NX,k)*(g_2(i,NX,k)*d(3,i,1,k)+&
             g_3(i,NX,k)*d(4,i,1,k)+g_4(i,NX,k)*d(5,i,1,k))
	end do
	end do
!$omp end parallel do
  end subroutine NSCBC_XHOUKOU

  subroutine NSCBC_yhoukou(g_1,g_2,g_3,g_4,g_5,dyg_1,dyg_2,dyg_3,dyg_4,dyg_5,dF1,dF2,dF3,dF4,dF5,c)
    implicit none
    integer ::i,j,k
    double precision,allocatable,dimension(:,:,:)::g_1,g_2,g_3,g_4,g_5,dyg_1,dyg_2,dyg_3,dyg_4,dyg_5,dF1,dF2,dF3,dF4,dF5
    double precision,allocatable,dimension(:,:,:,:)::L,d
    double precision,dimension(:,:,:),intent(in)::c
    double precision,allocatable,dimension(:,:,:)::Ma1
    allocate(L(1:5,1,NX,NZ),d(1:5,1,NX,NZ),Ma1(NY,NX,NZ))
    Ma1=g_3(:,:,:)/c(:,:,:)
!!!等温壁面条件=================================================
!!!NSCBCの適用i=0

!    dF1(1,1:xk-1,:)=1d0/(c(1,1:xk-1,:)**2)*(g_3(1,1:xk-1,:) - c(1,1:xk-1,:))&
!         *(-g_1(1,1:xk-1,:)*c(1,1:xk-1,:)*dyg_3(1,1:xk-1,:)+dyg_5(1,1:xk-1,:))
!    dF1(yk,xk:NX,:)=1d0/(c(yk,xk:NX,:)**2)*(g_3(yk,xk:NX,:) - c(yk,xk:NX,:))&
!         *(-g_1(yk,xk:NX,:)*c(yk,xk:NX,:)*dyg_3(yk,xk:NX,:)+dyg_5(yk,xk:NX,:))
!    dF1(yk,xk,:)=dF1(yk,xk,:)

!    dF1(1,bx1+1:bx1+1+hx-1,5+1:hz+5)=dF1(1,bx1+1:bx1+1+hx-1,5+1:hz+5)*1.4d0
!    dF1(1,bx1+1:bx1+1+hx-1,NZ2-hz2+1:NZ2+hz2)=dF1(1,bx1+1:bx1+1+hx-1,NZ2-hz2+1:NZ2+hz2)*1.4d0
!    dF1(1,bx1+1:bx1+1+hx-1,NZ-hz-4:NZ-5)=dF1(1,bx1+1:bx1+1+hx-1,NZ-hz-4:NZ-5)*1.4d0

!!!    i=NY　流出条件
!if(MA<1d0)then
!$omp parallel do
	do k=1,NZ
	do j=1,NX
    L(1,1,j,k)=sigma_y*(1d0-Ma1(NY,j,k)**2d0)*c(NY,j,k)*(g_5(NY,j,k)-p_inf)/Ly
    L(5,1,j,k)=(g_3(NY,j,k) + c(NY,j,k))&
         *(g_1(NY,j,k)*c(NY,j,k)*dyg_3(NY,j,k)+dyg_5(NY,j,k))
    L(2,1,j,k)=g_3(NY,j,k)*(c(NY,j,k)**2d0*dyg_1(NY,j,k)-dyg_5(NY,j,k))
    L(3,1,j,k)=g_3(NY,j,k)*dyg_2(NY,j,k)
    L(4,1,j,k)=g_3(NY,j,k)*dyg_4(NY,j,k)
!!!流入d行列の計算i=NY
    d(1,1,j,k)=1d0/(c(NY,j,k)**2)*(0.5d0*(L(1,1,j,k)+L(5,1,j,k))+L(2,1,j,k))
    d(2,1,j,k)=0.5d0*(L(1,1,j,k)+L(5,1,j,k))
    d(3,1,j,k)=L(3,1,j,k)
    d(4,1,j,k)=0.5d0/(g_1(NY,j,k)*c(NY,j,k))*(-L(1,1,j,k)+L(5,1,j,k))
    d(5,1,j,k)=L(4,1,j,k)
!!!i=NY
    dF1(NY,j,k)=d(1,1,j,k)
    dF2(NY,j,k)=g_2(NY,j,k)*d(1,1,j,k)+g_1(NY,j,k)*d(3,1,j,k)
    dF3(NY,j,k)=g_3(NY,j,k)*d(1,1,j,k)+g_1(NY,j,k)*d(4,1,j,k)
    dF4(NY,j,k)=g_4(NY,j,k)*d(1,1,j,k)+g_1(NY,j,k)*d(5,1,j,k)
    dF5(NY,j,k)=0.5d0*(g_2(NY,j,k)**2d0+g_3(NY,j,k)**2d0+g_4(NY,j,k)**2d0)*d(1,1,j,k)&
         +d(2,1,j,k)/(gamma-1d0)&
         +g_1(NY,j,k)*(g_2(NY,j,k)*d(3,1,j,k)+g_3(NY,j,k)*d(4,1,j,k)+g_4(NY,j,k)*d(5,1,j,k))
   ! dF(5,NY,j,k)=L(5,1,j,k)/(gamma-1d0)
	end do
end do
!$omp end parallel do

  end subroutine NSCBC_YHOUKOU

     subroutine result_2(t,T_k,g_old,g_1,g_2,g_3,g_4,g_5,&
     			dxg_1,dxg_2,dxg_3,dxg_4,dxg_5,&
     			dyg_1,dyg_2,dyg_3,dyg_4,dyg_5,&
     			dzg_1,dzg_2,dzg_3,dzg_4,dzg_5,LU_x,LU_y,LU_xk,LU_yk,x_k,dx_k,y_k,dy_k,avg_p)
    implicit none
    integer,intent(in)::t
    double precision,intent(in),dimension(:,:,:)::T_k
    double precision,dimension(:,:)::g_old
    double precision,dimension(:,:,:),allocatable::uzu_w
    double precision,dimension(:,:),allocatable::avg_u,q_y,q_x,q_z,uv,uu,avg_p
    double precision,dimension(:,:,:),intent(inout)::g_1,g_2,g_3,g_4,g_5
    double precision,dimension(:,:,:),intent(inout)::dxg_1,dxg_2,dxg_3,dxg_4,dxg_5
    double precision,dimension(:,:,:),intent(inout)::dyg_1,dyg_2,dyg_3,dyg_4,dyg_5
    double precision,dimension(:,:,:),intent(inout)::dzg_1,dzg_2,dzg_3,dzg_4,dzg_5
    double precision,dimension(:),intent(in)::y_k,dy_k
    double precision,dimension(:),intent(in)::dx_k,x_k
    double precision,dimension(:,:),intent(in)::LU_x,LU_y
    double precision,dimension(1:xk,1:3),intent(in)::LU_xk
    double precision,dimension(yk:NY,1:3),intent(in)::LU_yk
    integer::n7,n14,n21,n28,n35
    double precision::p
    integer::i,j,k,handan=1
    character(len=13)name
    character(len=13)name1
!!!xydanmen
    character(len=17)xy1_name
    character(len=17)xy2_name
    character(len=17)xz1_name
    allocate(avg_u(NY,NX),q_y(NY,NX),q_z(ny,nx),q_x(ny,nx),uv(NY,NX),uu(NY,NX))
    allocate(uzu_w(NY,NX,NZ))

    avg_u=0d0
     do k=1,NZ
     avg_u(:,:)=avg_u+g_2(:,:,k)
     end do
     avg_u=avg_u/dble(NZ)

    !$omp parallel do
    do k=1,NZ
    do i=1,yk-1
    	do j=xk+1,NX
    	g_1(i,j,k)=1d0
    	g_2(i,j,k)=0d0
    	g_3(i,j,k)=0d0
    	g_4(i,j,k)=0d0
    	g_5(i,j,:)=p_inf
    	dxg_1(i,j,k)=0d0
    	dyg_1(i,j,k)=0d0
    	dzg_1(i,j,k)=0d0
    	dxg_2(i,j,k)=0d0
    	dyg_2(i,j,k)=0d0
    	dzg_2(i,j,k)=0d0
    	dxg_3(i,j,k)=0d0
    	dyg_3(i,j,k)=0d0
    	dzg_3(i,j,k)=0d0
    	dxg_4(i,j,k)=0d0
    	dyg_4(i,j,k)=0d0
    	dzg_4(i,j,k)=0d0
    	dxg_5(i,j,k)=0d0
    	dyg_5(i,j,k)=0d0
    	dzg_5(i,j,k)=0d0
    	    	avg_u(i,j)=0d0
	end do
    	end do
   end do
   !$omp end parallel do

       uzu_w=dxg_2(:,:,:)*dyg_3(:,:,:)+dyg_3(:,:,:)*dzg_4(:,:,:)+dzg_4(:,:,:)*dxg_2(:,:,:)&
            -(dxg_3(:,:,:)*dyg_2(:,:,:)+dzg_2(:,:,:)*dxg_4(:,:,:)+dyg_4(:,:,:)*dzg_3(:,:,:))
	do j=1,xk
	          q_y(:,j)=CCS(1,NY,gdy,T_k(:,j,5),LU_y)*dy_k
	          uv(:,j)=CCS(1,NY,gdy,g_2(:,j,5)*g_3(:,j,5),LU_y)*dy_k
	end do
	do j=xk+1,NX
	          q_y(yk:NY,j)=CCS(yk,NY,gdy,T_k(yk:NY,j,5),LU_yk)*dy_k(yk:NY)
	          uv(yk:NY,j)=CCS(yk,NY,gdy,g_2(yk:NY,j,5)*g_3(yk:NY,j,5),LU_yk)*dy_k(yk:NY)
	end do
	do i=1,yk
		uu(i,1:xk)=CCS(1,xk,gdx,g_2(i,1:xk,5)**2d0,LU_xk)*dx_k(1:xk)
	        q_x(i,1:xk)=CCS(1,xk,gdx,T_k(i,1:xk,5),LU_xk)*dx_k(1:xk)
	end do
	do i=yk,NY
		uu(i,:)=CCS(1,NX,gdx,g_2(i,:,5)**2d0,LU_x)*dx_k
	        q_x(i,:)=CCS(1,Nx,gdx,T_k(i,:,5),LU_x)*dx_k
	end do

    if (xydanmen==1)then

       k=NZ/2
       xy1_name='xy1_000000_nn.txt'
       xy2_name='xy2_000000_nn.txt'
       xz1_name='xz1_000000_nn.txt'
       write(xy1_name(5:10),"(i6.6)")t
       write(xy1_name(12:13),"(i2.2)")NZ/2
       write(xz1_name(5:10),"(i6.6)")t
       write(xz1_name(12:13),"(i2.2)")NZ/2
       write(xy2_name(5:10),"(i6.6)")t
       write(xy2_name(12:13),"(i2.2)")Nz/2
       open(111,file=xy1_name)
  !     open(112,file=xy2_name)
       open(113,file=xz1_name)

       !!!ライトヒルテンソル音源項
       do j=1,NX
          do i=1,NY
             write(111,'(19E25.15e3)') x_k(j),y_k(i),&
                  g_2(i,j,k),g_2(i,j,k)-avg_u(i,j),dxg_2(i,j,k)+dyg_3(i,j,k)+dzg_4(i,j,k),&
		 g_1(i,j,k),g_3(i,j,k),g_4(i,j,k),g_5(i,j,k)-g_old(i,j)&
		 ,g_5(i,j,k),T_k(i,j,k),dxg_2(i,j,k),dyg_3(i,j,k),dzg_4(i,j,k),g_5(i,j,k)-avg_p(i,j),&
     dyg_2(i,j,k),dxg_3(i,j,k)-dyg_2(i,j,k),uzu_w(i,j,k)
          end do
          write(111,'(14E23.13e3)')
       end do
       close(111)
        i=yk+15
       do j=1,NX
          do k=1,NZ
             write(113,'(9E20.10e3)') x_k(j),dble(k-1)*dz,g_2(i,j,k),&
             g_2(i,j,k)-avg_u(i,j),dxg_2(i,j,k)+dyg_3(i,j,k)+dzg_4(i,j,k),&
             g_5(i,j,k)-avg_p(i,j),dxg_3(i,j,k)-dyg_2(i,j,k),uzu_w(i,j,k),dyg_2(i,j,k)
          end do
         write(113,'(7E20.10e3)')
     end do
    close(113)

    end if

     if(mod(t,Mmod*3)==0)then
    open(31,file='data_rho3D1.txt',status='replace')
    open(32,file='data_u3D1.txt',status='replace')
    open(33,file='data_v3D1.txt',status='replace')
    open(34,file='data_w3D1.txt',status='replace')
    open(35,file='data_p3D1.txt',status='replace')
    open(36,file='data_T3D1.txt',status='replace')
do i=1,NY
 do j=1,NX
   do k=1,NZ
    write(31,*) g_1(i,j,k)
    write(32,*) g_2(i,j,k)
    write(33,*) g_3(i,j,k)
    write(34,*) g_4(i,j,k)
    write(35,*) g_5(i,j,k)
    write(36,*) T_k(i,j,k)
end do
end do
end do
    close(31)
    close(32)
    close(33)
    close(34)
    close(35)
    close(36)
    end if
end subroutine result_2
!!!====================================

subroutine result_rms(g1,g2,g3,g4,g5,T_k,t,g_h1,g_h2,g_h3,g_h4,g_h5,g_h6,x_k,y_k,rmsg1,rmsg2,rmsg3,rmsg4,rmsg5,rmsg6)
  implicit none
  integer::i,j,k,t
  double precision,dimension(:),intent(in)::y_k,x_k
  double precision,dimension(:,:,:),intent(in)::g1,g2,g3,g4,g5
    double precision,dimension(:,:,:)::T_k
    double precision,dimension(:,:)::rmsg1,rmsg2,rmsg3,rmsg4,rmsg5,rmsg6
    double precision,dimension(:,:)::g_h1,g_h2,g_h3,g_h4,g_h5,g_h6
  character(len=8)name1
!!!!ファブル平均
!$omp parallel do
  do k=1,NZ
  	do j=1,NX
  	do i=1,NY
     rmsg1(i,j) =rmsg1(i,j)+abs(g1(i,j,k)-g_h1(i,j))
     rmsg2(i,j) =rmsg2(i,j)+abs(g2(i,j,k)-g_h2(i,j))
     rmsg3(i,j) =rmsg3(i,j)+abs(g3(i,j,k)-g_h3(i,j))
     rmsg4(i,j) =rmsg4(i,j)+abs(g4(i,j,k)-g_h4(i,j))
     rmsg5(i,j) =rmsg5(i,j)+abs(g5(i,j,k)-g_h5(i,j))
     rmsg6(i,j) =rmsg6(i,j)+abs(T_k(i,j,k)-g_h6(i,j))
  	end do
  	end do
  end do
!$omp end parallel do

if(t==MMM)then
  name1='rmsg.txt'
 rmsg1=rmsg1/dble(NZ)/dble(MMM-MM+1)
 rmsg2=rmsg2/dble(NZ)/dble(MMM-MM+1)
 rmsg3=rmsg3/dble(NZ)/dble(MMM-MM+1)
 rmsg4=rmsg4/dble(NZ)/dble(MMM-MM+1)
 rmsg5=rmsg5/dble(NZ)/dble(MMM-MM+1)
 rmsg6=rmsg6/dble(NZ)/dble(MMM-MM+1)

  open(21,file=name1)
  do j=1,NX
     do i=1,NY
        write(21,'(17E21.10e3)')x_k(j),y_k(i),rmsg1(i,j),rmsg2(i,j),rmsg3(i,j)&
        ,rmsg4(i,j),rmsg5(i,j),rmsg1(i,j)/g_h1(i,j),rmsg5(i,j)/g_h5(i,j),&
        -rmsg2(i,j)*rmsg3(i,j),-rmsg3(i,j)*rmsg4(i,j),-rmsg2(i,j)*rmsg4(i,j),&
       -g_h1(i,j)*rmsg2(i,j)*rmsg3(i,j),-g_h1(i,j)*rmsg3(i,j)*rmsg4(i,j),-g_h1(i,j)*rmsg2(i,j)*rmsg4(i,j),&
        rmsg6(i,j),rmsg6(i,j)/g_h6(i,j)
     end do
     write(21,'(17E21.10e3)')
  end do
  close(21)
end if
end subroutine result_rms
!!!=============================================
  subroutine result_final_2(g_1,g_2,g_3,g_4,g_5,T_k,t,g_h1,g_h2,g_h3,g_h4,g_h5,g_h6,x_k,y_k,avg_p)
    implicit none
	double precision,dimension(:),intent(in)::y_k,x_k
	double precision,dimension(:,:,:),intent(in)::g_1,g_2,g_3,g_4,g_5
	double precision,dimension(:,:)::g_h1,g_h2,g_h3,g_h4,g_h5,g_h6
	double precision,dimension(:,:,:)::T_k
	double precision,dimension(:,:)::avg_p
    integer::i,j,k,t
    !$omp parallel do
  do k=1,NZ
  	do j=1,NX
  		do i=1,NY
     g_h1(i,j) =g_h1(i,j)+g_1(i,j,k)
     g_h2(i,j) =g_h2(i,j)+g_2(i,j,k)
     g_h3(i,j) =g_h3(i,j)+g_3(i,j,k)
     g_h4(i,j) =g_h4(i,j)+g_4(i,j,k)
     g_h5(i,j) =g_h5(i,j)+g_5(i,j,k)
     g_h6(i,j) =g_h6(i,j)+T_k(i,j,k)
  end do
  end do
  end do
  !$omp end parallel do
    if(t==MM)then
	g_h1=g_h1/dble(NZ)/dble(MM-M1+1)
	g_h2=g_h2/dble(NZ)/dble(MM-M1+1)
	g_h3=g_h3/dble(NZ)/dble(MM-M1+1)
	g_h4=g_h4/dble(NZ)/dble(MM-M1+1)
	g_h5=g_h5/dble(NZ)/dble(MM-M1+1)
	g_h6=g_h6/dble(NZ)/dble(MM-M1+1)

       open(31,file='average_g.txt')
       do j=1,NX
          do i=1,NY
             write(31,'(8E30.15e3)')x_k(j),y_k(i),g_h1(i,j),g_h2(i,j),&
             	g_h3(i,j),g_h4(i,j),g_h5(i,j),g_h6(i,j)
          end do
          write(31,'(8E30.15e3)')

       end do
       close(31)
       avg_p=g_h5(:,:)
    end if
  end subroutine result_final_2

  subroutine result_bunseki(G_H1,G_H2,G_H3,G_H4,G_H5,G_h6,LUy,LUyk,x_k,y_k,dy_k,rmsg1,rmsg2,rmsg3,rmsg4,rmsg5,rmsg6)
    implicit none
    integer::i,j,k,n
    double precision,dimension(:,:)::G_H1,G_H2,G_H3,G_H4,G_H5,G_h6,rmsg1,rmsg2,rmsg3,rmsg4,rmsg5,rmsg6
    double precision,dimension(:,:),allocatable::dyUD
    double precision::mu1
    double precision,dimension(:),allocatable::tauw
    double precision,dimension(:),allocatable::U_t,U_t1,x_kRE
    double precision,dimension(1:NY,1:3)::LUy
    double precision,dimension(yk:NY,1:3)::LUyk
    double precision,dimension(:,:),allocatable::U_p,y_p
    double precision,dimension(:)::dy_k,y_k,x_k
    double precision,dimension(:,:),allocatable::omomi
    double precision,dimension(:,:,:),allocatable::yobi1,yobi2
    double precision,dimension(:),allocatable::theta1,theta2
    allocate(yobi1(6,NY,NX),yobi2(6,NY,NX),omomi(NY,NX),y_p(NY,NX))
    allocate(U_t(nx),U_p(NY,NX),tauw(NX),dyUD(NY,NX),x_kRE(NX),U_t1(NX))
    allocate(theta1(Nx),theta2(NX))
    mu1=1d0/RE*1.1d0
	U_p=0d0
	U_t=0d0
	dyUD=0d0

    do j=1,xk-1
    	do i=1,NY
       		omomi(i,j)=dsqrt(G_h1(i,j)/G_h1(1,j))
	end do
    end do
    do j=xk,NY
    	do i=1,NY
       		omomi(i,j)=dsqrt(G_h1(i,j)/G_h1(yk,j))
	end do
    end do

    do i=1,xk-1
       dyUD(:,i)=(dy_k)*CCS(1,NY,gdy,G_h2(:,i),LUy)
    end do
    do i=xk,NX
       dyUD(yk:NY,i)=(dy_k(yk:NY))*CCS(yk,NY,gdy,G_h2(yk:NY,i),LUyk)
    end do
    do i=1,xk-1
    tauw(i)=mu1*dyUD(1,i)
    end do
    do i=xk,NX
    tauw(i)=mu1*dyUD(yk,i)
    end do

    U_t=dsqrt(tauw)

    do j=1,xk-1
	    do i=1,NY
       		U_p(i,j)=G_H2(i,j)/U_t(j)
	    end do
    end do
    do j=xk,NX
	    do i=yk,NY
       		U_p(i,j)=G_H2(i,j)/U_t(j)
	    end do
    end do

    do i=1,xk-1
       x_kRE(i)=(x_k(i)-x_k(1))*RE+(RE/1.7208d0)**2d0
       y_p(:,i)=y_k(:)*U_t(i)/mu1
    end do
    do i=xk,NX
       x_kRE(i)=(x_k(i)-x_k(1))*RE+(RE/1.7208d0)**2d0
       y_p(yk:NY,i)=(y_k(yk:NY)-y_k(yk))*U_t(i)/mu1
    end do

    open(200,file='bunsekiy=-1.csv')
    do j=1,xk-1
       write(200,*)'x=',dble(j)*dx
       do i=1,NY
          write(200,*)y_p(i,j),',',U_p(i,j),',',omomi(i,j)*U_p(i,j)
       end do
       write(200,*)
    end do
    do j=xk,NX
       write(200,*)'x=',dble(j)*dx
       do i=yk,NY
          write(200,*)y_p(i,j),',',U_p(i,j),',',omomi(i,j)*U_p(i,j)
       end do
       write(200,*)
    end do
    close(200)

    open(201,file='uvw_bunseki_1.dat')
    do j=1,xk-1
       do i=1,NY
          write(201,'(6E20.10e3)')x_k(j),y_k(i),y_p(i,j),omomi(i,j)*rmsg2(i,j)/U_t(j),omomi(i,j)*rmsg3(i,j)/U_t(j)&
          ,omomi(i,j)*rmsg4(i,j)/U_t(j)
       end do

       write(201,*)
    end do
    close(201)

    open(201,file='uvw_bunseki_2.dat')
    do j=xk,NX
       do i=yk,NY
          write(201,'(6E20.10e3)')x_k(j),y_k(i),y_p(i,j),omomi(i,j)*rmsg2(i,j)/U_t(j),omomi(i,j)*rmsg3(i,j)/U_t(j)&
          ,omomi(i,j)*rmsg4(i,j)/U_t(j)
       end do
       write(201,*)
    end do
    close(201)

theta1=0d0
theta2=0d0
	do i=1,190
		theta1(:)=theta1(:)+(y_k(i+1)-y_k(i))*(G_h2(190,:)-g_h2(i,:))
		theta2(:)=theta2(:)+(y_k(i+1)-y_k(i))*(G_h2(190,:)-g_h2(i,:))*g_h2(i,:)
	end do
	theta1=theta1/g_h2(190,:)
	theta2=theta2/(g_h2(190,:)**2d0)


    open(100,file='REtau_Ut_tauw_Cf.csv')
    write(100,*)'!Retau,Ut,tw,Cf'
    do i=1,NX
       write(100,'(10E20.10e3)')x_kre(i),x_k(i),theta2(i)/mu1,theta1(i)/mu1/1.4d0,U_t(i)&
       ,tauw(i)*2d0,0.455d0/(LOG10(x_kre(i))**(2.58d0)),theta1(i),theta2(i),theta1(i)/theta2(i)
    end do
    close (100)
  end subroutine result_bunseki

  !!!=====================================================================
!!!Subroutine LU行列の計算
!!!=====================================================================
  subroutine LU(fn,NN,LUcom,sig)
    implicit none
    integer,intent(in)::NN,fn
    double precision,intent(in)::sig
    double precision,intent(inout),dimension(fn:NN,1:3)::LUcom
    double precision,dimension(fn:NN,fn:NN)::A,L,U
    integer::i,j,k,n
    double precision::alpha_5,alphaD_5 !DCS5次精度alpha
    double precision::alpha_3,alphaD_3 !DCS3次精度alpha
    double precision::sum,alphaOS !片側差分
!!!-----------------------------
!!!初期値入力
!!!-----------------------------
    LUcom(:,:)=0d0
    A(:,:)=0d0
    L(:,:)=0d0
    U(:,:)=0d0
    sum=0d0
    alpha_5=1d0/3d0
    alphaD_5=1d0
    alpha_3=0.25d0
    alphaD_3=1d0
    alphaOS=3d0
!!!======================================================
!!!LUの値の計算
!!!======================================================
!!!A行列およびLUの初期設定
    do i=fn,NN
       A(i,i)=1d0
    end do
    do i=fn+2,NN-2
       A(i,i-1)=alpha_5*(1d0-sig*alphaD_5)
       A(i,i+1)=alpha_5*(1d0+sig*alphaD_5)
    end do
!!!--------------------
!!!片側差分
!!!--------------------!片側CCS
       A(fn,fn+1)=alphaOS
       A(NN,NN-1)=A(fn,fn+1)
    if(sig==0d0)then
       A(fn,fn+1)=alphaOS
       A(NN,NN-1)=A(fn,fn+1)
    else if(sig<-0.01d0)then
       A(fn,fn+1)=5d0
       A(NN,NN-1)=3d0
    else if(sig>0.01d0)then
       A(fn,fn+1)=3d0
       A(NN,NN-1)=5d0
   end if
 !!---------------------------------
    A(fn+1,fn)=alpha_3*(1d0-sig*alphaD_3)
    A(fn+1,fn+2)=alpha_3*(1d0+sig*alphaD_3)
    A(NN-1,NN-2)=alpha_3*(1d0-sig*alphaD_3)
    A(NN-1,NN)=alpha_3*(1d0+sig*alphaD_3)
!!!---------------------------------------------
!!!LU行列の計算開始
    do j=fn,NN
       U(fn,j)=A(fn,j)
    end do
    do i=fn,NN
       L(i,fn)=A(i,fn)/U(fn,fn)
       L(i,i)=1.d0
    end do
!!!----------------------
    do i=fn+1,NN
       do j=fn+1,NN
          do k=fn,i-1
             sum=sum+L(i,k)*U(k,j)
          end do
          U(i,j)=A(i,j)-sum
          sum=0d0
       end do
       do n=i+1,NN
          do k=fn,i-1
             sum=sum+L(n,k)*U(k,i)
          end do
          L(n,i)=(A(n,i)-sum)/U(i,i)
          sum=0d0
       end do
    end do


    LUcom(NN,2)=U(NN,NN)
    do i=fn,NN-1
       LUcom(i,1)=L(i+1,i)
    end do
    do i=fn,NN-1
       do k=2,3
          LUcom(i,k)=U(i,i+k-2)
       end do
    end do
!!!!!!!!!!!!!!!!!LUの確認
  end subroutine LU
  !!!==========================================================
!!!周期境界条件
!!!==========================================================
    subroutine LUsyuuki(L,U,sig)
    implicit none
    double precision,intent(in)::sig
    double precision,dimension(:,:)::L,U
    double precision,allocatable,dimension(:,:)::A
    integer::i,j,k,n,im1,im2,ip1,ip2
    double precision::alpha_5,alphaD_5,sum !DCS5次精度alpha
    allocate(A(NZ,NZ))
!!!-----------------------------
!!!初期値入力
!!!-----------------------------
    A(:,:)=0d0
    L(:,:)=0d0
    U(:,:)=0d0
    sum=0d0
    alpha_5=1d0/3d0
    alphaD_5=1d0
!!!======================================================
!!!LUの値の計算
!!!======================================================
!!!A行列およびLUの初期設定
    do i=1,NZ
       A(i,i)=1d0
    end do
    do i=1,NZ
       im1=mod((NZ-1)+i-1,NZ)+1
       ip1=mod((NZ-1)+i+1,NZ)+1
       A(i,im1)=alpha_5*(1d0-sig*alphaD_5)
       A(i,ip1)=alpha_5*(1d0+sig*alphaD_5)
    end do
!!!LU行列の計算開始
    do j=1,NZ
       U(1,j)=A(1,j)
    end do
    do i=1,NZ
       L(i,1)=A(i,1)/U(1,1)
       L(i,i)=1.d0
    end do
!!!----------------------
    do i=2,NZ
       do j=2,NZ
          do k=1,i-1
             sum=sum+L(i,k)*U(k,j)
          end do
          U(i,j)=A(i,j)-sum
          sum=0d0
       end do
       do n=i+1,NZ
          do k=1,i-1
             sum=sum+L(n,k)*U(k,i)
          end do
          L(n,i)=(A(n,i)-sum)/U(i,i)
          sum=0d0
       end do
    end do

  end subroutine LUsyuuki
!!-------------------------------------------------------------
!!-------------------------------------------------------------
  subroutine buffer(Ux,Uy,sigxy,x,y,tt)
	implicit none
    double precision,dimension(:,:,:),intent(out)::Ux,Uy
    double precision,dimension(:,:,:),intent(out)::sigxy
    double precision,dimension(:),intent(in)::x,y
    double precision::alU,als,berl,tt
    integer::i,j,k
    alU=1.0d0
    als=1.05d0
    berl=0.01d0
    !!!x
    do i=1,NX/2
       Ux(:,i,:)=alU*c_inf*(dtanh(datanh(berl/alU-1d0)*(x(i)-x(NX))/(x(NX-bx1)-x(NX)))&
            -dtanh(datanh(berl/alU-1d0)*(x(i)-x(1))/(x(bx1)-x(1))))
    end do
    do i=NX/2,NX
       Ux(:,i,:)=alU*c_inf*(dtanh(datanh(berl/alU-1d0)*(x(i)-x(NX))/(x(NX-bx2)-x(NX)))&
            -dtanh(datanh(berl/alU-1d0)*(x(i)-x(1))/(x(bx2)-x(1))))
    end do
    Ux=-Ux

    do i=1,NX/2
  	Ux(:,i,:)=0.5d0*UX(:,i,:)
    end do

     do i=NX/2,NX
 	Ux(:,i,:)=0.7d0*UX(:,i,:)
    end do

   do i=1,bx1
       sigxy(:,i,:)=-als*c_inf*((x(bx1)-x(i))/(x(bx1)-x(1)))**3d0*0.5d0
   end do

    do i=NX-bx2,NX
       sigxy(:,i,:)=-als*c_inf*((x(i)-x(NX-bx2))/(x(NX)-x(NX-bx2)))**3d0*0.7d0
    end do

!!!y
   !
    do i=1,NY
       Uy(i,:,:)=alU*c_inf*(dtanh(datanh(berl/alU-1d0)*(y(i)-y(Ny))/(y(Ny-by)-y(NY)))&
            -dtanh(datanh(berl/alU-1d0)*(y(i)-y(1))/(y(by)-y(1))))
    end do

    do i=1,NY/2
       Uy(i,:,:)=0d0
    end do
    Uy=-Uy!*0.d0

    do i=Ny-by,Ny
       sigxy(i,:,:)=sigxy(i,:,:)-als*c_inf*((y(i)-y(Ny-by))/(y(Ny)-y(Ny-by)))**3d0!*0.d0
    end do

    open(10,file='tasikame1.csv')
    do i=1,NX
       write(10,*)dble(i)*dx,Ux(1,i,1),sigxy(1,i,1)
    end do
    open(11,file='tasikame2.csv')
    do i=1,NY
       write(11,*)y(i),Uy(i,1,1),sigxy(i,1,1)
    end do
    close(10)
    close(11)
 ! if(buffer_h==0)then
 !   Ux=0d0*Ux
!Ux=Ux/3d0
  ! Uy=0d0*Uy

  !sigxy=0d0*sigxy

 !sigxy=sigxy/3d0
 !end if
  end subroutine buffer
end module sub_mod

!!!===========================================================================
!!!===========================================================================
!!!===========================================================================
program main
  !$use omp_lib
  use sub_mod
 implicit none
  integer::i,j,k,t,n,Mmod2
  double precision,allocatable,dimension(:,:,:)::Q_1,Q1_1,Q2_1,Q01
  double precision,allocatable,dimension(:,:,:)::Q_2,Q1_2,Q2_2,Q02
  double precision,allocatable,dimension(:,:,:)::Q_3,Q1_3,Q2_3,Q03
  double precision,allocatable,dimension(:,:,:)::Q_4,Q1_4,Q2_4,Q04
  double precision,allocatable,dimension(:,:,:)::Q_5,Q1_5,Q2_5,Q05

  double precision,allocatable,dimension(:,:,:)::g_1,dxg_1,dyg_1,dzg_1
  double precision,allocatable,dimension(:,:,:)::g_2,dxg_2,dyg_2,dzg_2
  double precision,allocatable,dimension(:,:,:)::g_3,dxg_3,dyg_3,dzg_3
  double precision,allocatable,dimension(:,:,:)::g_4,dxg_4,dyg_4,dzg_4
  double precision,allocatable,dimension(:,:,:)::g_5,dxg_5,dyg_5,dzg_5
  double precision,allocatable,dimension(:,:,:)::Fxp1,Fxm1,Fx1
  double precision,allocatable,dimension(:,:,:)::Fxp2,Fxm2,Fx2
  double precision,allocatable,dimension(:,:,:)::Fxp3,Fxm3,Fx3
  double precision,allocatable,dimension(:,:,:)::Fxp4,Fxm4,Fx4
  double precision,allocatable,dimension(:,:,:)::Fxp5,Fxm5,Fx5
  double precision,allocatable,dimension(:,:,:)::Fyp1,Fym1,Fy1
  double precision,allocatable,dimension(:,:,:)::Fyp2,Fym2,Fy2
  double precision,allocatable,dimension(:,:,:)::Fyp3,Fym3,Fy3
  double precision,allocatable,dimension(:,:,:)::Fyp4,Fym4,Fy4
  double precision,allocatable,dimension(:,:,:)::Fyp5,Fym5,Fy5

  double precision,allocatable,dimension(:,:,:)::Fzp1,Fzm1,Fz1
  double precision,allocatable,dimension(:,:,:)::Fzp2,Fzm2,Fz2
  double precision,allocatable,dimension(:,:,:)::Fzp3,Fzm3,Fz3
  double precision,allocatable,dimension(:,:,:)::Fzp4,Fzm4,Fz4
  double precision,allocatable,dimension(:,:,:)::Fzp5,Fzm5,Fz5

  double precision,allocatable,dimension(:,:)::LU_D_p_nx,LU_D_m_nx
  double precision,allocatable,dimension(:,:)::LU_D_p_ny,LU_D_m_ny
  double precision,allocatable,dimension(:,:):: LU_C_nx,LU_C_ny
  double precision,dimension(1:xk,1:3)::LU_D_p_nxk,LU_D_m_nxk
  double precision,dimension(yk:NY,1:3)::LU_D_p_nyk,LU_D_m_nyk
  double precision,dimension(1:xk,1:3):: LU_C_nxk
  double precision,dimension(yk:NY,1:3):: LU_C_nyk
  double precision,allocatable,dimension(:,:)::Lz_Dp,Uz_Dp,Lz_Dm,Uz_Dm,Lz_C,Uz_C

  double precision,allocatable,dimension(:,:,:)::T_k,c,dxT,dyT,dzT
  double precision,allocatable,dimension(:,:,:)::iryu_x_1,iryu_y_1,iryu_z_1
  double precision,allocatable,dimension(:,:,:)::iryu_x_2,iryu_y_2,iryu_z_2
  double precision,allocatable,dimension(:,:,:)::iryu_x_3,iryu_y_3,iryu_z_3
  double precision,allocatable,dimension(:,:,:)::iryu_x_4,iryu_y_4,iryu_z_4
  double precision,allocatable,dimension(:,:,:)::iryu_x_5,iryu_y_5,iryu_z_5

  double precision,allocatable,dimension(:,:,:) :: mu,mukq,mugRE,q_x,q_y,q_z
  double precision,allocatable,dimension(:,:,:) :: txx,txy,txz
  double precision,allocatable,dimension(:,:,:) :: tyy,tyz
  double precision,allocatable,dimension(:,:,:) :: tzz,dxQ5,dyQ5
  double precision,allocatable,dimension(:,:,:) :: dxtxx,dxtxy,dxtxz
  double precision,allocatable,dimension(:,:,:) :: dytxy,dytyy,dytyz
  double precision,allocatable,dimension(:,:,:) :: dztxz,dztyz,dztzz
  double precision,allocatable,dimension(:,:,:) :: dxq_x,dyq_y,dzq_z

  double precision,allocatable,dimension(:,:,:) ::dFx1,dFy1,dFz1
  double precision,allocatable,dimension(:,:,:) ::dFx2,dFy2,dFz2
  double precision,allocatable,dimension(:,:,:) ::dFx3,dFy3,dFz3
  double precision,allocatable,dimension(:,:,:) ::dFx4,dFy4,dFz4
  double precision,allocatable,dimension(:,:,:) ::dFx5,dFy5,dFz5

  double precision,allocatable,dimension(:,:,:) ::sigxy
  double precision::time
  double precision::t1,t2
  double precision,allocatable,dimension(:)::y_k,dy_k,dydy_k,y_moto
  double precision,allocatable,dimension(:)::x_k,dx_k
  double precision,allocatable,dimension(:)::z_k
  double precision,allocatable,dimension(:,:)::D_rho,D_u,D_v,D_p,D_T,D_w
  double precision,allocatable,dimension(:,:,:)::ran_u,ran_v,ran_w,data1,data2,data3,data4,data5,datak
  double precision,allocatable,dimension(:,:,:)::avg_u,avg_uzu,avg_u2
  double precision,allocatable,dimension(:,:)::rho_h,rho_h2,u_h,v_h,w_h,p_h,T_h,zenT
  double precision,allocatable,dimension(:,:)::D2_rho,D2_u,D2_v,D2_p,D2_T

  double precision,allocatable,dimension(:,:,:)::Ux,Uy
  double precision,allocatable,dimension(:,:)::rmsg1,rmsg2,rmsg3,rmsg4,rmsg5,rmsg6
  double precision,allocatable,dimension(:,:)::g_h1,g_h2,g_h3,g_h4,g_h5,g_h6
  double precision,allocatable,dimension(:,:)::g_in,g_old
  double precision,allocatable,dimension(:,:)::Vin,avg_p
  integer::loopMK


  allocate(Q_1(NY,NX,NZ),Q1_1(NY,NX,NZ),Q2_1(NY,NX,NZ),Q01(NY,NX,NZ))
  allocate(Q_2(NY,NX,NZ),Q1_2(NY,NX,NZ),Q2_2(NY,NX,NZ),Q02(NY,NX,NZ))
  allocate(Q_3(NY,NX,NZ),Q1_3(NY,NX,NZ),Q2_3(NY,NX,NZ),Q03(NY,NX,NZ))
  allocate(Q_4(NY,NX,NZ),Q1_4(NY,NX,NZ),Q2_4(NY,NX,NZ),Q04(NY,NX,NZ))
  allocate(Q_5(NY,NX,NZ),Q1_5(NY,NX,NZ),Q2_5(NY,NX,NZ),Q05(NY,NX,NZ))
  allocate(g_1(NY,NX,NZ),dxg_1(NY,NX,NZ),dyg_1(NY,NX,NZ),dzg_1(NY,NX,NZ))
  allocate(g_2(NY,NX,NZ),dxg_2(NY,NX,NZ),dyg_2(NY,NX,NZ),dzg_2(NY,NX,NZ))
  allocate(g_3(NY,NX,NZ),dxg_3(NY,NX,NZ),dyg_3(NY,NX,NZ),dzg_3(NY,NX,NZ))
  allocate(g_4(NY,NX,NZ),dxg_4(NY,NX,NZ),dyg_4(NY,NX,NZ),dzg_4(NY,NX,NZ))
  allocate(g_5(NY,NX,NZ),dxg_5(NY,NX,NZ),dyg_5(NY,NX,NZ),dzg_5(NY,NX,NZ))
  allocate(Fxp1(NY,NX,NZ),Fxm1(NY,NX,NZ),Fx1(NY,NX,NZ))
  allocate(Fxp2(NY,NX,NZ),Fxm2(NY,NX,NZ),Fx2(NY,NX,NZ))
  allocate(Fxp3(NY,NX,NZ),Fxm3(NY,NX,NZ),Fx3(NY,NX,NZ))
  allocate(Fxp4(NY,NX,NZ),Fxm4(NY,NX,NZ),Fx4(NY,NX,NZ))
  allocate(Fxp5(NY,NX,NZ),Fxm5(NY,NX,NZ),Fx5(NY,NX,NZ))
  allocate(Fyp1(NY,NX,NZ),Fym1(NY,NX,NZ),Fy1(NY,NX,NZ))
  allocate(Fyp2(NY,NX,NZ),Fym2(NY,NX,NZ),Fy2(NY,NX,NZ))
  allocate(Fyp3(NY,NX,NZ),Fym3(NY,NX,NZ),Fy3(NY,NX,NZ))
  allocate(Fyp4(NY,NX,NZ),Fym4(NY,NX,NZ),Fy4(NY,NX,NZ))
  allocate(Fyp5(NY,NX,NZ),Fym5(NY,NX,NZ),Fy5(NY,NX,NZ))
  allocate(Fzp1(NY,NX,NZ),Fzm1(NY,NX,NZ),Fz1(NY,NX,NZ))
  allocate(Fzp2(NY,NX,NZ),Fzm2(NY,NX,NZ),Fz2(NY,NX,NZ))
  allocate(Fzp3(NY,NX,NZ),Fzm3(NY,NX,NZ),Fz3(NY,NX,NZ))
  allocate(Fzp4(NY,NX,NZ),Fzm4(NY,NX,NZ),Fz4(NY,NX,NZ))
  allocate(Fzp5(NY,NX,NZ),Fzm5(NY,NX,NZ),Fz5(NY,NX,NZ))
  allocate(iryu_x_1(NY,NX,NZ),iryu_y_1(NY,NX,NZ),iryu_z_1(NY,NX,NZ))
  allocate(iryu_x_2(NY,NX,NZ),iryu_y_2(NY,NX,NZ),iryu_z_2(NY,NX,NZ))
  allocate(iryu_x_3(NY,NX,NZ),iryu_y_3(NY,NX,NZ),iryu_z_3(NY,NX,NZ))
  allocate(iryu_x_4(NY,NX,NZ),iryu_y_4(NY,NX,NZ),iryu_z_4(NY,NX,NZ))
  allocate(iryu_x_5(NY,NX,NZ),iryu_y_5(NY,NX,NZ),iryu_z_5(NY,NX,NZ))
  allocate(dFx1(NY,NX,NZ),dFy1(NY,NX,NZ),dFz1(NY,NX,NZ))
  allocate(dFx2(NY,NX,NZ),dFy2(NY,NX,NZ),dFz2(NY,NX,NZ))
  allocate(dFx3(NY,NX,NZ),dFy3(NY,NX,NZ),dFz3(NY,NX,NZ))
  allocate(dFx4(NY,NX,NZ),dFy4(NY,NX,NZ),dFz4(NY,NX,NZ))
  allocate(dFx5(NY,NX,NZ),dFy5(NY,NX,NZ),dFz5(NY,NX,NZ))

  allocate(g_old(NY,NX))
  allocate(Ux(NY,NX,NZ),UY(NY,NX,NZ),sigxy(NY,NX,NZ),dxQ5(NY,NX,NZ),dyQ5(Ny,NX,NZ))
  allocate(LU_D_p_nx(NX,3),LU_D_m_nx(NX,3),LU_C_nx(NX,3))
  allocate(LU_D_p_ny(NY,3),LU_D_m_ny(NY,3),LU_C_ny(NY,3))
  allocate(Lz_Dp(NZ,NZ),Uz_Dp(NZ,NZ),Lz_Dm(NZ,NZ),Uz_Dm(NZ,NZ)&
       ,Lz_C(NZ,NZ),Uz_C(NZ,NZ))
  allocate(T_k(NY,NX,NZ),c(NY,NX,NZ),dxT(NY,NX,NZ),dyT(NY,NX,NZ),dzT(NY,NX,NZ))
  allocate(mu(NY,NX,NZ),mukq(NY,NX,NZ),mugRE(NY,NX,NZ),q_x(NY,NX,NZ),q_y(NY,NX,NZ),q_z(NY,NX,NZ))
  allocate(txx(NY,NX,NZ),txy(NY,NX,NZ),txz(NY,NX,NZ))
  allocate(tyy(NY,NX,NZ),tyz(NY,NX,NZ))
  allocate(tzz(NY,NX,NZ))
  allocate(y_k(NY),dy_k(NY),dydy_k(NY),y_moto(NY))
  allocate(x_k(NX),dx_k(NX))
  allocate(z_k(NZ))
  allocate(D2_rho(NY,NX),D2_u(NY,NX),D2_v(NY,NX),D2_p(NY,NX),D2_T(NY,NX))
  allocate(D_rho(NY,NX),D_u(NY,NX),D_v(NY,NX),D_p(NY,NX),D_T(NY,NX),D_w(NY,NX))
  allocate(ran_u(Ndt,NY,NZ),ran_v(Ndt,NY,NZ),ran_w(Ndt,NY,NZ))
  allocate(avg_u(NY,NX,NZ),avg_uzu(NY,NX,NZ),avg_p(NY,NX),avg_u2(NY,NX,NZ))
  allocate(rho_h(NY,NX),rho_h2(NY,NX),u_h(NY,NX),v_h(NY,NX),w_h(NY,NX)&
       ,p_h(NY,NX),T_h(NY,NX),zenT(NY,NX))
  allocate(g_in(5,NY))
  allocate(dxtxx(NY,NX,NZ),dxtxy(NY,NX,NZ),dxtxz(NY,NX,NZ),dxq_x(NY,NX,NZ))
  allocate(dytxy(NY,NX,NZ),dytyy(NY,NX,NZ),dytyz(NY,NX,NZ),dyq_y(NY,NX,NZ))
  allocate(dztzz(NY,NX,NZ),dztxz(NY,NX,NZ),dztyz(NY,NX,NZ),dzq_z(NY,NX,NZ))
  allocate(Vin(hx,hz))
  allocate(rmsg1(NY,NX),rmsg2(NY,NX),rmsg3(NY,NX),rmsg4(NY,NX),rmsg5(NY,NX),rmsg6(NY,NX))
  allocate(g_h1(NY,NX),g_h2(NY,NX),g_h3(NY,NX),g_h4(NY,NX),g_h5(NY,NX),g_h6(NY,NX))
open(222,file='p_keisoku.dat')
open(223,file='p_keisoku_1.dat')
open(224,file='p_keisoku_2.dat')
!!!====================================================
!!!初期値設定
!!!====================================================
  call zahyou_henkan(x_k,z_k,y_moto,y_k,dx_k,dy_k)  !!!初期座標および格子伸長
  t=0
  g_h1=0d0;g_h2=0d0;g_h3=0d0;g_h4=0d0;g_h5=0d0;g_h6=0d0
  rmsg1=0d0;rmsg2=0d0;rmsg3=0d0;rmsg4=0d0;rmsg5=0d0;rmsg6=0d0
sigxy=0d0
Ux=0d0
Uy=0d0
  call buffer(Ux,Uy,sigxy,x_k,y_k,0d0)

  open(200,file='hun.csv')
   do i=1,hx
	do j=1,hz
    vin(i,j)=parav*(((dble(i-1)/dble(hx-1)-0.5d0)**2-0.5d0**2)*&
    			((dble(j-1)/dble(hz-1)-0.5d0)**2-0.5d0**2))
  write(200,*)i,j,vin(i,j)
   end do
   end do
	close(200)
!!!-----------------------------
!!!LU行列の計算
!!!------------------------------
!  open(44,file='kakuniin_LU.csv')
!  open(45,file='kakunin_LUsyuuki.csv')
  call LU(1,NX,LU_D_p_nx,-sig)  !!DCS_plus_x
  call LU(1,NY,LU_D_p_ny,-sig)  !!DCS_plus_y
  call LU(1,NX,LU_D_m_nx,sig)   !!DCS_minus_x
  call LU(1,NY,LU_D_m_ny,sig)   !!DCS_minus_y
  call LU(1,NX,LU_C_nx,0d0)   !!CCS_x
  call LU(1,NY,LU_C_ny,0d0)   !!CCS_y
  call LU(1,xk,LU_D_p_nxk,-sig)  !!DCS_plus_x
  call LU(yk,NY,LU_D_p_nyk,-sig)  !!DCS_plus_y
  call LU(1,xk,LU_D_m_nxk,sig)   !!DCS_minus_x
  call LU(yk,NY,LU_D_m_nyk,sig)   !!DCS_minus_y
  call LU(1,xk,LU_C_nxk,0d0)   !!CCS_x
  call LU(yk,NY,LU_C_nyk,0d0)   !!CCS_y
  call LUsyuuki(Lz_Dp,Uz_Dp,-sig)
  call LUsyuuki(Lz_Dm,Uz_Dm,sig)
  call LUsyuuki(Lz_C,Uz_C,0d0)

    if (data2dim==1)then
     open(101,file='data_rho1.txt',status='old')
     open(102,file='data_u1.txt',status='old')
     open(103,file='data_v1.txt',status='old')
     open(104,file='data_p1.txt',status='old')
     open(105,file='data_T1.txt',status='old')
     do i=1,NY
        do j=1,NX
     read(101,*)D2_rho(i,j)
     read(102,*)D2_u(i,j)
     read(103,*)D2_v(i,j)
     read(104,*)D2_p(i,j)
     read(105,*)D2_T(i,j)
end do
end do
     close(101)
     close(102)
     close(103)
     close(104)
     close(105)
     do k=1,NZ
        g_1(:,:,k)=D2_rho
        g_2(:,:,k)=D2_u
        g_3(:,:,k)=D2_v
        g_5(:,:,k)=D2_p
        T_k(:,:,k)=D2_T
     end do
     g_4(:,:,:)=0d0
     g_in(1,:)=D2_rho(:,1)
     g_in(2,:)=D2_u(:,1)
     g_in(3,:)=D2_v(:,1)
     g_in(4,:)=D2_p(:,1)
     g_in(5,:)=D2_T(:,1)

     !$omp parallel do
     do k=1,NZ
     	do j=1,NX
     		do i=1,NY
    Q_1(i,j,k)=g_1(i,j,k)
    Q_2(i,j,k)=g_1(i,j,k)*g_2(i,j,k)
    Q_3(i,j,k)=g_1(i,j,k)*g_3(i,j,k)
    Q_4(i,j,k)=g_1(i,j,k)*g_4(i,j,k)
    Q_5(i,j,k)=g_5(i,j,k)/(gamma-1d0)&
         +(g_1(i,j,k)*(g_2(i,j,k)**2+g_3(i,j,k)**2+g_4(i,j,k)**2))*0.5d0
         	end do
         end do
     end do
     !$omp end parallel do
  call Q_g_in(Q_1,Q_2,Q_3,Q_4,Q_5,g_1,g_2,g_3,g_4,g_5,T_k,mu,mugRE,mukq,t,y_k,g_in,vin)
  c=dsqrt(T_k/(Ma**2))

  Q01=Q_1*(sigxy)
  Q02=Q_2*(sigxy)
  Q03=Q_3*(sigxy)
  Q04=Q_4*(sigxy)
  Q05=Q_5*(sigxy)
  end if

  if (data3dim==1)then
     open(101,file='data_rho3D1.txt',status='old')
     open(102,file='data_u3D1.txt',status='old')
     open(103,file='data_v3D1.txt',status='old')
     open(104,file='data_w3D1.txt',status='old')
     open(105,file='data_p3D1.txt',status='old')
     open(106,file='data_T3D1.txt',status='old')

     do i=1,NY
        do j=1,NX
           do k=1,NZ
              read(101,*)g_1(i,j,k)
              read(102,*)g_2(i,j,k)
              read(103,*)g_3(i,j,k)
              read(104,*)g_4(i,j,k)
              read(105,*)g_5(i,j,k)
              read(106,*)T_k(i,j,k)
           end do
        end do
     end do
     g_in(1,:)=g_1(:,1,1)
     g_in(2,:)=g_2(:,1,1)
     g_in(3,:)=g_3(:,1,1)
     g_in(4,:)=g_5(:,1,1)
     g_in(5,:)=T_k(:,1,1)

     close(101)
     close(102)
     close(103)
     close(104)
     close(105)
     close(106)

          !$omp parallel do
     do k=1,NZ
     	do j=1,NX
     		do i=1,NY
    Q_1(i,j,k)=g_1(i,j,k)
    Q_2(i,j,k)=g_1(i,j,k)*g_2(i,j,k)
    Q_3(i,j,k)=g_1(i,j,k)*g_3(i,j,k)
    Q_4(i,j,k)=g_1(i,j,k)*g_4(i,j,k)
    Q_5(i,j,k)=g_5(i,j,k)/(gamma-1d0)&
         +(g_1(i,j,k)*(g_2(i,j,k)**2+g_3(i,j,k)**2+g_4(i,j,k)**2))*0.5d0
         	end do
         end do
     end do
     !$omp end parallel do
  call Q_g_in(Q_1,Q_2,Q_3,Q_4,Q_5,g_1,g_2,g_3,g_4,g_5,T_k,mu,mugRE,mukq,t,y_k,g_in,vin)
  c=dsqrt(T_k/(Ma**2))
  end if

write(*,*)'a'

!!!=============================================================
  do t=1,MMM
  do j=1,NX
	do i=1,NY
  g_old(i,j)=g_5(i,j,NZ/2)
	end do
  end do
  time=dble(t)*dt

!!!===================================
!!!TVD method
!!!===================================

	!$omp parallel do
    do k=1,NZ
    	do j=1,NX
    		do i=1,NY
!!!x
 	 c(i,j,k)=dsqrt(T_k(i,j,k)/(Ma**2))
       Fx1(i,j,k)=Q_2(i,j,k)
       Fx2(i,j,k)=Q_2(i,j,k)*g_2(i,j,k)+g_5(i,j,k)
       Fx3(i,j,k)=Q_2(i,j,k)*g_3(i,j,k)
       Fx4(i,j,k)=Q_2(i,j,k)*g_4(i,j,k)
       Fx5(i,j,k)=(Q_5(i,j,k)+g_5(i,j,k))*g_2(i,j,k)

       Fy1(i,j,k)=Q_3(i,j,k)
       Fy2(i,j,k)=Q_3(i,j,k)*g_2(i,j,k)
       Fy3(i,j,k)=Q_3(i,j,k)*g_3(i,j,k)+g_5(i,j,k)
       Fy4(i,j,k)=Q_3(i,j,k)*g_4(i,j,k)
       Fy5(i,j,k)=(Q_5(i,j,k)+g_5(i,j,k))*g_3(i,j,k)
    !!!Lax F
       Fz1(i,j,k)=Q_4(i,j,k)
       Fz2(i,j,k)=Q_4(i,j,k)*g_2(i,j,k)
       Fz3(i,j,k)=Q_4(i,j,k)*g_3(i,j,k)
       Fz4(i,j,k)=Q_4(i,j,k)*g_4(i,j,k)+g_5(i,j,k)
       Fz5(i,j,k)=(Q_5(i,j,k)+g_5(i,j,k))*g_4(i,j,k)

       Fxp1(i,j,k)=(Fx1(i,j,k)+Q_1(i,j,k))*0.5d0
       Fxp2(i,j,k)=(Fx2(i,j,k)+Q_2(i,j,k))*0.5d0
       Fxp3(i,j,k)=(Fx3(i,j,k)+Q_3(i,j,k))*0.5d0
       Fxp4(i,j,k)=(Fx4(i,j,k)+Q_4(i,j,k))*0.5d0
       Fxp5(i,j,k)=(Fx5(i,j,k)+Q_5(i,j,k))*0.5d0

       Fxm1(i,j,k)=(Fx1(i,j,k)-Q_1(i,j,k))*0.5d0
       Fxm2(i,j,k)=(Fx2(i,j,k)-Q_2(i,j,k))*0.5d0
       Fxm3(i,j,k)=(Fx3(i,j,k)-Q_3(i,j,k))*0.5d0
       Fxm4(i,j,k)=(Fx4(i,j,k)-Q_4(i,j,k))*0.5d0
       Fxm5(i,j,k)=(Fx5(i,j,k)-Q_5(i,j,k))*0.5d0


       Fyp1(i,j,k)=(Fy1(i,j,k)+Q_1(i,j,k))*0.5d0
       Fyp2(i,j,k)=(Fy2(i,j,k)+Q_2(i,j,k))*0.5d0
       Fyp3(i,j,k)=(Fy3(i,j,k)+Q_3(i,j,k))*0.5d0
       Fyp4(i,j,k)=(Fy4(i,j,k)+Q_4(i,j,k))*0.5d0
       Fyp5(i,j,k)=(Fy5(i,j,k)+Q_5(i,j,k))*0.5d0

       Fym1(i,j,k)=(Fy1(i,j,k)-Q_1(i,j,k))*0.5d0
       Fym2(i,j,k)=(Fy2(i,j,k)-Q_2(i,j,k))*0.5d0
       Fym3(i,j,k)=(Fy3(i,j,k)-Q_3(i,j,k))*0.5d0
       Fym4(i,j,k)=(Fy4(i,j,k)-Q_4(i,j,k))*0.5d0
       Fym5(i,j,k)=(Fy5(i,j,k)-Q_5(i,j,k))*0.5d0

!!!!z
       Fzp1(i,j,k)=(Fz1(i,j,k)+Q_1(i,j,k))*0.5d0
       Fzp2(i,j,k)=(Fz2(i,j,k)+Q_2(i,j,k))*0.5d0
       Fzp3(i,j,k)=(Fz3(i,j,k)+Q_3(i,j,k))*0.5d0
       Fzp4(i,j,k)=(Fz4(i,j,k)+Q_4(i,j,k))*0.5d0
       Fzp5(i,j,k)=(Fz5(i,j,k)+Q_5(i,j,k))*0.5d0

       Fzm1(i,j,k)=(Fz1(i,j,k)-Q_1(i,j,k))*0.5d0
       Fzm2(i,j,k)=(Fz2(i,j,k)-Q_2(i,j,k))*0.5d0
       Fzm3(i,j,k)=(Fz3(i,j,k)-Q_3(i,j,k))*0.5d0
       Fzm4(i,j,k)=(Fz4(i,j,k)-Q_4(i,j,k))*0.5d0
       Fzm5(i,j,k)=(Fz5(i,j,k)-Q_5(i,j,k))*0.5d0
       			end do
       		end do
       end do
   !$omp end parallel do
   !$omp parallel do
    do k=1,NZ
        do i=1,yk
           dxg_1(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_1(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_2(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_2(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_3(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_3(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_4(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_4(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_5(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_5(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           q_x(i,1:xk,k)=dx_k(1:xk)*mukq(i,1:xk,k)*CCS(1,xk,gdx,T_k(i,1:xk,k),LU_C_nxk(1:xk,1:3))
            q_x(i,xk,k)=0d0
	 end do
        do i=yk+1,NY
           dxg_1(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_1(i,:,k),LU_C_nx)
           dxg_2(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_2(i,:,k),LU_C_nx)
           dxg_3(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_3(i,:,k),LU_C_nx)
           dxg_4(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_4(i,:,k),LU_C_nx)
           dxg_5(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_5(i,:,k),LU_C_nx)
           q_x(i,:,k)=dx_k(:)*mukq(i,:,k)*CCS(1,NX,gdx,T_k(i,:,k),LU_C_nx)
        end do
    end do
       !$omp end parallel do
    !$omp parallel do
    do k=1,NZ
	    do i=1,xk-1
              dyg_1(:,i,k)=dy_k*CCS(1,NY,gdy,g_1(:,i,k),LU_C_ny)
              dyg_2(:,i,k)=dy_k*CCS(1,NY,gdy,g_2(:,i,k),LU_C_ny)
              dyg_3(:,i,k)=dy_k*CCS(1,NY,gdy,g_3(:,i,k),LU_C_ny)
              dyg_4(:,i,k)=dy_k*CCS(1,NY,gdy,g_4(:,i,k),LU_C_ny)
              dyg_5(:,i,k)=dy_k*CCS(1,NY,gdy,g_5(:,i,k),LU_C_ny)
              q_y(:,i,k)=mukq(:,i,k)*dy_k*CCS(1,NY,gdy,T_k(:,i,k),LU_C_ny)
     	      q_y(1,i,k)=0d0
	end do

        do i=xk,NX
              dyg_1(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_1(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_2(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_2(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_3(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_3(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_4(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_4(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_5(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_5(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              q_y(yk:NY,i,k)=mukq(yk:NY,i,k)*dy_k(yk:NY)*CCS(yk,NY,gdy,T_k(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
     	      q_y(yk,i,k)=0d0
        end do
     end do
        !$omp end parallel do
     !$omp parallel do
     do i=1,NX
        do j=1,NY
              dzg_1(j,i,:)=CCSz(g_1(j,i,:),Lz_C,Uz_C)
              dzg_2(j,i,:)=CCSz(g_2(j,i,:),Lz_C,Uz_C)
              dzg_3(j,i,:)=CCSz(g_3(j,i,:),Lz_C,Uz_C)
              dzg_4(j,i,:)=CCSz(g_4(j,i,:),Lz_C,Uz_C)
              dzg_5(j,i,:)=CCSz(g_5(j,i,:),Lz_C,Uz_C)
           q_z(j,i,:)=mukq(j,i,:)*CCSz(T_k(j,i,:),Lz_C,Uz_C)
        end do
	end do
    !$omp end parallel do

    !$omp parallel do
    do k=1,NZ
    	do j=1,NX
    		do i=1,NY
	        txy(i,j,k)=(dyg_2(i,j,k)+dxg_3(i,j,k))*mugRE(i,j,k)
         	txz(i,j,k)=(dzg_2(i,j,k)+dxg_4(i,j,k))*mugRE(i,j,k)
    	    	tyz(i,j,k)=(dzg_3(i,j,k)+dyg_4(i,j,k))*mugRE(i,j,k)
    		txx(i,j,k)=tvd2*(2d0*dxg_2(i,j,k)-dyg_3(i,j,k)-dzg_4(i,j,k))*mugRE(i,j,k)
  		tyy(i,j,k)=tvd2*(2d0*dyg_3(i,j,k)-dxg_2(i,j,k)-dzg_4(i,j,k))*mugRE(i,j,k)
		tzz(i,j,k)=tvd2*(2d0*dzg_4(i,j,k)-dyg_3(i,j,k)-dxg_2(i,j,k))*mugRe(i,j,k)
			end do
		end do
	end do
	!$omp end parallel do
	!$omp parallel do
    do k=1,NZ
       do i=1,yk      !!!xhoukouw
          iryu_x_1(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp1(i,1:xk,k),LU_D_p_nxk(1:xk,:),-sig)&
          						+fBDCS(1,xk,gdx,Fxm1(i,1:xk,k),LU_D_m_nxk(1:xk,:),sig))
          iryu_x_2(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp2(i,1:xk,k),LU_D_p_nxk(1:xk,:),-sig)&
          						+fBDCS(1,xk,gdx,Fxm2(i,1:xk,k),LU_D_m_nxk(1:xk,:),sig))
          iryu_x_3(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp3(i,1:xk,k),LU_D_p_nxk(1:xk,:),-sig)&
          						+fBDCS(1,xk,gdx,Fxm3(i,1:xk,k),LU_D_m_nxk(1:xk,:),sig))
          iryu_x_4(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp4(i,1:xk,k),LU_D_p_nxk(1:xk,:),-sig)&
          						+fBDCS(1,xk,gdx,Fxm4(i,1:xk,k),LU_D_m_nxk(1:xk,:),sig))
          iryu_x_5(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp5(i,1:xk,k),LU_D_p_nxk(1:xk,:),-sig)&
          						+fBDCS(1,xk,gdx,Fxm5(i,1:xk,k),LU_D_m_nxk(1:xk,:),sig))
       end do
       do i=yk+1,NY       !!!xhoukouw
          iryu_x_1(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp1(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm1(i,:,k),LU_D_m_nx,sig))
          iryu_x_2(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp2(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm2(i,:,k),LU_D_m_nx,sig))
          iryu_x_3(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp3(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm3(i,:,k),LU_D_m_nx,sig))
          iryu_x_4(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp4(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm4(i,:,k),LU_D_m_nx,sig))
          iryu_x_5(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp5(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm5(i,:,k),LU_D_m_nx,sig))
       end do
      end do
!$omp end parallel do
!$omp parallel do
    do k=1,NZ
       do j=1,xk-1       !!!yhoukou
          iryu_y_1(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp1(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym1(:,j,k),LU_D_m_ny,sig))
          iryu_y_2(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp2(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym2(:,j,k),LU_D_m_ny,sig))
          iryu_y_3(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp3(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym3(:,j,k),LU_D_m_ny,sig))
          iryu_y_4(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp4(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym4(:,j,k),LU_D_m_ny,sig))
          iryu_y_5(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp5(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym5(:,j,k),LU_D_m_ny,sig))
       end do
       do j=xk,NX       !!!yhoukou
          iryu_y_1(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp1(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym1(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_2(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp2(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym2(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_3(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp3(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym3(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_4(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp4(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym4(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_5(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp5(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym5(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
       end do
    end do
    !$omp end parallel do
    !$omp parallel do
    do j=1,NX
       do i=1,NY          !!!z方向移流項
          iryu_z_1(i,j,:)=DCSz(Fzp1(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm1(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_2(i,j,:)=DCSz(Fzp2(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm2(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_3(i,j,:)=DCSz(Fzp3(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm3(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_4(i,j,:)=DCSz(Fzp4(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm4(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_5(i,j,:)=DCSz(Fzp5(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm5(i,j,:),Lz_Dm,Uz_Dm,sig)
       end do
    end do
    !$omp end parallel do
	!$omp parallel do
    do k=1,NZ
       do i=1,yk
          dxtxx(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,txx(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxtxy(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,txy(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxtxz(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,txz(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxq_x(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,q_x(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxQ5(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,Q_5(i,1:xk,k),LU_C_NXk(1:xk,1:3))
       end do
       do i=yk+1,NY
          dxtxx(i,:,k)=dx_k*CCS(1,NX,gdx,txx(i,:,k),LU_C_nx)
          dxtxy(i,:,k)=dx_k*CCS(1,NX,gdx,txy(i,:,k),LU_C_nx)
          dxtxz(i,:,k)=dx_k*CCS(1,NX,gdx,txz(i,:,k),LU_C_nx)
          dxq_x(i,:,k)=dx_k*CCS(1,NX,gdx,q_x(i,:,k),LU_C_nx)
          dxQ5(i,:,k)=dx_k*CCS(1,NX,gdx,Q_5(i,:,k),LU_C_NX)
       end do
    end do
    !$omp end parallel do
    !$omp parallel do
   do k=1,NZ
       do j=1,xk-1
          dytyy(:,j,k)=dy_k*CCS(1,NY,gdy,tyy(:,j,k),LU_C_ny)
          dytxy(:,j,k)=dy_k*CCS(1,NY,gdy,txy(:,j,k),LU_C_ny)
          dytyz(:,j,k)=dy_k*CCS(1,NY,gdy,tyz(:,j,k),LU_C_ny)
          dyq_y(:,j,k)=dy_k*CCS(1,NY,gdy,q_y(:,j,k),LU_C_ny)
          dyQ5(:,j,k)=dy_k*CCS(1,NY,gdy,Q_5(:,j,k),LU_C_ny)
       end do
       do j=xk,NX
          dytyy(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,tyy(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dytxy(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,txy(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dytyz(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,tyz(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dyq_y(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,q_y(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dyQ5(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,Q_5(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
       end do
    end do

	!$omp end parallel do
	!$omp parallel do
    do j=1,NX
       do i=1,NY
          dztxz(i,j,:)=CCSz(txz(i,j,:),Lz_C,Uz_C)
          dztyz(i,j,:)=CCSz(tyz(i,j,:),Lz_C,Uz_C)
          dztzz(i,j,:)=CCSz(tzz(i,j,:),Lz_C,Uz_C)
          dzq_z(i,j,:)=CCSz(q_z(i,j,:),Lz_C,Uz_C)
       end do
    end do
   !$omp end parallel do
   !$omp parallel do
    do k=1,NZ
    	do j=1,NX
    		do i=1,NY
	    dFx2(i,j,k)=dxtxx(i,j,k)
    	dFx3(i,j,k)=dxtxy(i,j,k)
    	dFx4(i,j,k)=dxtxz(i,j,k)
    	dFx5(i,j,k)=dxg_2(i,j,k)*txx(i,j,k)+g_2(i,j,k)*dxtxx(i,j,k)&
         +dxg_3(i,j,k)*txy(i,j,k)+g_3(i,j,k)*dxtxy(i,j,k)&
         +dxg_4(i,j,k)*txz(i,j,k)+g_4(i,j,k)*dxtxz(i,j,k) -dxq_x(i,j,k)
    !!ryuusyutu jouken
    dFx3(i,NX,k)=0d0
    dFx4(i,NX,k)=0d0
    dFx5(i,NX,k)=dFx5(i,NX,k)+dxq_x(i,NX,k)-g_3(i,NX,k)*dxtxy(i,NX,k)-g_4(i,NX,k)*dxtxz(i,NX,k)

    dFx1(i,j,k)=dxg_1(i,j,k)*Ux(i,j,k)
    dFx2(i,j,k)=dFx2(i,j,k)+(dxg_1(i,j,k)*g_2(i,j,k)+dxg_2(i,j,k)*g_1(i,j,k))*Ux(i,j,k)
    dFx3(i,j,k)=dFx3(i,j,k)+(dxg_1(i,j,k)*g_3(i,j,k)+dxg_3(i,j,k)*g_1(i,j,k))*Ux(i,j,k)
    dFx4(i,j,k)=dFx4(i,j,k)+(dxg_1(i,j,k)*g_4(i,j,k)+dxg_4(i,j,k)*g_1(i,j,k))*Ux(i,j,k)
    dFx5(i,j,k)=dFx5(i,j,k)+dxQ5(i,j,k)*Ux(i,j,k)

    dFy2(i,j,k)=dytxy(i,j,k)
    dFy3(i,j,k)=dytyy(i,j,k)
    dFy4(i,j,k)=dytyz(i,j,k)
    dFy5(i,j,k)=dyg_2(i,j,k)*txy(i,j,k)+g_2(i,j,k)*dytxy(i,j,k)&
         +dyg_3(i,j,k)*tyy(i,j,k) + g_3(i,j,k)*dytyy(i,j,k)&
         +dyg_4(i,j,k)*tyz(i,j,k) + g_4(i,j,k)*dytyz(i,j,k) -dyq_y(i,j,k)
	!!ryusyutu jouken
	dFy2(NY,j,k)=0d0
	dFy4(NY,j,k)=0d0
	dFy5(NY,j,k)=dFy5(NY,j,k)+dyq_y(NY,j,k)-g_2(NY,j,k)*dytxy(NY,j,k)-g_4(NY,j,k)*dytyz(NY,j,k)

    dFy1(i,j,k)=dyg_1(i,j,k)*Uy(i,j,k)
    dFy2(i,j,k)=dFy2(i,j,k)+(dyg_1(i,j,k)*g_2(i,j,k)+dyg_2(i,j,k)*g_1(i,j,k))*Uy(i,j,k)
    dFy3(i,j,k)=dFy3(i,j,k)+(dyg_1(i,j,k)*g_3(i,j,k)+dyg_3(i,j,k)*g_1(i,j,k))*Uy(i,j,k)
    dFy4(i,j,k)=dFy4(i,j,k)+(dyg_1(i,j,k)*g_4(i,j,k)+dyg_4(i,j,k)*g_1(i,j,k))*Uy(i,j,k)
    dFy5(i,j,k)=dFy5(i,j,k)+(dyQ5(i,j,k))*Uy(i,j,k)


    dFz1(i,j,k)=0d0
    dFz2(i,j,k)=dztxz(i,j,k)
    dFz3(i,j,k)=dztyz(i,j,k)
    dFz4(i,j,k)=dztzz(i,j,k)
    dFz5(i,j,k)=dzg_2(i,j,k)*txz(i,j,k)+g_2(i,j,k)*dztxz(i,j,k)&
         +dzg_3(i,j,k)*tyz(i,j,k)+g_3(i,j,k)*dztyz(i,j,k)&
         +dzg_4(i,j,k)*tzz(i,j,k)+g_4(i,j,k)*dztzz(i,j,k)-dzq_z(i,j,k)
    		end do
    	end do
    end do

!$omp end parallel do
    call NSCBC_xhoukou(g_1,g_2,g_3,g_4,g_5,dxg_1,dxg_2,dxg_3,dxg_4,dxg_5,iryu_x_1,iryu_x_2,iryu_x_3,iryu_x_4,iryu_x_5,c)
    call NSCBC_yhoukou(g_1,g_2,g_3,g_4,g_5,dyg_1,dyg_2,dyg_3,dyg_4,dyg_5,iryu_y_1,iryu_y_2,iryu_y_3,iryu_y_4,iryu_y_5,c)

!!!TVD1---------------------------------------------------------------------------------------------------
    !Q1=Q-dt*(iryu_x+iryu_y+iryu_z-(dFx+dFy+dFz)-sigxy*Q+Q0)
!$omp parallel do
  do k=1,NZ
  	do j=1,NX
  		do i=1,NY
	  	Q1_1(i,j,k)=Q_1(i,j,k)-dt*(iryu_x_1(i,j,k)+iryu_y_1(i,j,k)+iryu_z_1(i,j,k)&
	  				-dFx1(i,j,k)-dFy1(i,j,k)-dFz1(i,j,k)-sigxy(i,j,k)*Q_1(i,j,k)+Q01(i,j,k))
	  	Q1_2(i,j,k)=Q_2(i,j,k)-dt*(iryu_x_2(i,j,k)+iryu_y_2(i,j,k)+iryu_z_2(i,j,k)&
	  				-dFx2(i,j,k)-dFy2(i,j,k)-dFz2(i,j,k)-sigxy(i,j,k)*Q_2(i,j,k)+Q02(i,j,k))
	  	Q1_3(i,j,k)=Q_3(i,j,k)-dt*(iryu_x_3(i,j,k)+iryu_y_3(i,j,k)+iryu_z_3(i,j,k)&
	  				-dFx3(i,j,k)-dFy3(i,j,k)-dFz3(i,j,k)-sigxy(i,j,k)*Q_3(i,j,k)+Q03(i,j,k))
	  	Q1_4(i,j,k)=Q_4(i,j,k)-dt*(iryu_x_4(i,j,k)+iryu_y_4(i,j,k)+iryu_z_4(i,j,k)&
	  				-dFx4(i,j,k)-dFy4(i,j,k)-dFz4(i,j,k)-sigxy(i,j,k)*Q_4(i,j,k)+Q04(i,j,k))
	  	Q1_5(i,j,k)=Q_5(i,j,k)-dt*(iryu_x_5(i,j,k)+iryu_y_5(i,j,k)+iryu_z_5(i,j,k)&
	  				-dFx5(i,j,k)-dFy5(i,j,k)-dFz5(i,j,k)-sigxy(i,j,k)*Q_5(i,j,k)+Q05(i,j,k))
		end do
	end do
   end do
  !$omp end parallel do

  call Q_g_in(Q1_1,Q1_2,Q1_3,Q1_4,Q1_5,g_1,g_2,g_3,g_4,g_5,T_k,mu,mugRE,mukq,t,y_k,g_in,vin)

	!$omp parallel do
    do k=1,NZ
    	do j=1,NX
    		do i=1,NY
    		  c(i,j,k)=dsqrt(T_k(i,j,k)/(0.09d0))
!!!x
       Fx1(i,j,k)=Q1_2(i,j,k)
       Fx2(i,j,k)=Q1_2(i,j,k)*g_2(i,j,k)+g_5(i,j,k)
       Fx3(i,j,k)=Q1_2(i,j,k)*g_3(i,j,k)
       Fx4(i,j,k)=Q1_2(i,j,k)*g_4(i,j,k)
       Fx5(i,j,k)=(Q1_5(i,j,k)+g_5(i,j,k))*g_2(i,j,k)

       Fy1(i,j,k)=Q1_3(i,j,k)
       Fy2(i,j,k)=Q1_3(i,j,k)*g_2(i,j,k)
       Fy3(i,j,k)=Q1_3(i,j,k)*g_3(i,j,k)+g_5(i,j,k)
       Fy4(i,j,k)=Q1_3(i,j,k)*g_4(i,j,k)
       Fy5(i,j,k)=(Q1_5(i,j,k)+g_5(i,j,k))*g_3(i,j,k)
    !!!Lax F
       Fz1(i,j,k)=Q1_4(i,j,k)
       Fz2(i,j,k)=g_2(i,j,k)*Q1_4(i,j,k)
       Fz3(i,j,k)=g_3(i,j,k)*Q1_4(i,j,k)
       Fz4(i,j,k)=Q1_4(i,j,k)*g_4(i,j,k)+g_5(i,j,k)
       Fz5(i,j,k)=(Q1_5(i,j,k)+g_5(i,j,k))*g_4(i,j,k)

       Fxp1(i,j,k)=(Fx1(i,j,k)+Q1_1(i,j,k))*0.5d0
       Fxp2(i,j,k)=(Fx2(i,j,k)+Q1_2(i,j,k))*0.5d0
       Fxp3(i,j,k)=(Fx3(i,j,k)+Q1_3(i,j,k))*0.5d0
       Fxp4(i,j,k)=(Fx4(i,j,k)+Q1_4(i,j,k))*0.5d0
       Fxp5(i,j,k)=(Fx5(i,j,k)+Q1_5(i,j,k))*0.5d0

       Fxm1(i,j,k)=(Fx1(i,j,k)-Q1_1(i,j,k))*0.5d0
       Fxm2(i,j,k)=(Fx2(i,j,k)-Q1_2(i,j,k))*0.5d0
       Fxm3(i,j,k)=(Fx3(i,j,k)-Q1_3(i,j,k))*0.5d0
       Fxm4(i,j,k)=(Fx4(i,j,k)-Q1_4(i,j,k))*0.5d0
       Fxm5(i,j,k)=(Fx5(i,j,k)-Q1_5(i,j,k))*0.5d0


       Fyp1(i,j,k)=(Fy1(i,j,k)+Q1_1(i,j,k))*0.5d0
       Fyp2(i,j,k)=(Fy2(i,j,k)+Q1_2(i,j,k))*0.5d0
       Fyp3(i,j,k)=(Fy3(i,j,k)+Q1_3(i,j,k))*0.5d0
       Fyp4(i,j,k)=(Fy4(i,j,k)+Q1_4(i,j,k))*0.5d0
       Fyp5(i,j,k)=(Fy5(i,j,k)+Q1_5(i,j,k))*0.5d0

       Fym1(i,j,k)=(Fy1(i,j,k)-Q1_1(i,j,k))*0.5d0
       Fym2(i,j,k)=(Fy2(i,j,k)-Q1_2(i,j,k))*0.5d0
       Fym3(i,j,k)=(Fy3(i,j,k)-Q1_3(i,j,k))*0.5d0
       Fym4(i,j,k)=(Fy4(i,j,k)-Q1_4(i,j,k))*0.5d0
       Fym5(i,j,k)=(Fy5(i,j,k)-Q1_5(i,j,k))*0.5d0

!!!!z
       Fzp1(i,j,k)=(Fz1(i,j,k)+Q1_1(i,j,k))*0.5d0
       Fzp2(i,j,k)=(Fz2(i,j,k)+Q1_2(i,j,k))*0.5d0
       Fzp3(i,j,k)=(Fz3(i,j,k)+Q1_3(i,j,k))*0.5d0
       Fzp4(i,j,k)=(Fz4(i,j,k)+Q1_4(i,j,k))*0.5d0
       Fzp5(i,j,k)=(Fz5(i,j,k)+Q1_5(i,j,k))*0.5d0

       Fzm1(i,j,k)=(Fz1(i,j,k)-Q1_1(i,j,k))*0.5d0
       Fzm2(i,j,k)=(Fz2(i,j,k)-Q1_2(i,j,k))*0.5d0
       Fzm3(i,j,k)=(Fz3(i,j,k)-Q1_3(i,j,k))*0.5d0
       Fzm4(i,j,k)=(Fz4(i,j,k)-Q1_4(i,j,k))*0.5d0
       Fzm5(i,j,k)=(Fz5(i,j,k)-Q1_5(i,j,k))*0.5d0
       			end do
       		end do
       end do
	!$omp end parallel do
	!$omp parallel do

    do k=1,NZ
        do i=1,yk
           dxg_1(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_1(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_2(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_2(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_3(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_3(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_4(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_4(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_5(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_5(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           q_x(i,1:xk,k)=dx_k(1:xk)*mukq(i,1:xk,k)*CCS(1,xk,gdx,T_k(i,1:xk,k),LU_C_nxk(1:xk,1:3))
      	   q_x(i,xk,k)=0d0
	  end do
        do i=yk+1,NY
           dxg_1(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_1(i,:,k),LU_C_nx)
           dxg_2(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_2(i,:,k),LU_C_nx)
           dxg_3(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_3(i,:,k),LU_C_nx)
           dxg_4(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_4(i,:,k),LU_C_nx)
           dxg_5(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_5(i,:,k),LU_C_nx)
           q_x(i,:,k)=dx_k(:)*mukq(i,:,k)*CCS(1,NX,gdx,T_k(i,:,k),LU_C_nx)
        end do
    end do
    !$omp end parallel do
    !$omp parallel do
    do k=1,NZ
	    do i=1,xk-1
              dyg_1(:,i,k)=dy_k*CCS(1,NY,gdy,g_1(:,i,k),LU_C_ny)
              dyg_2(:,i,k)=dy_k*CCS(1,NY,gdy,g_2(:,i,k),LU_C_ny)
              dyg_3(:,i,k)=dy_k*CCS(1,NY,gdy,g_3(:,i,k),LU_C_ny)
              dyg_4(:,i,k)=dy_k*CCS(1,NY,gdy,g_4(:,i,k),LU_C_ny)
              dyg_5(:,i,k)=dy_k*CCS(1,NY,gdy,g_5(:,i,k),LU_C_ny)
           	  q_y(:,i,k)=mukq(:,i,k)*dy_k*CCS(1,NY,gdy,T_k(:,i,k),LU_C_ny)
		q_y(1,i,k)=0d0
        end do
        do i=xk,NX
              dyg_1(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_1(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_2(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_2(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_3(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_3(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_4(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_4(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_5(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_5(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              q_y(yk:NY,i,k)=mukq(yk:NY,i,k)*dy_k(yk:NY)*CCS(yk,NY,gdy,T_k(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
		q_y(yk,i,k)=0d0
    end do
     end do

	!$omp end parallel do
	!$omp parallel do

     do i=1,NX
        do j=1,NY
              dzg_1(j,i,:)=CCSz(g_1(j,i,:),Lz_C,Uz_C)
              dzg_2(j,i,:)=CCSz(g_2(j,i,:),Lz_C,Uz_C)
              dzg_3(j,i,:)=CCSz(g_3(j,i,:),Lz_C,Uz_C)
              dzg_4(j,i,:)=CCSz(g_4(j,i,:),Lz_C,Uz_C)
              dzg_5(j,i,:)=CCSz(g_5(j,i,:),Lz_C,Uz_C)
           q_z(j,i,:)=mukq(j,i,:)*CCSz(T_k(j,i,:),Lz_C,Uz_C)
        end do
     end do
    !$omp end parallel do

    !$omp parallel do
    do k=1,NZ
    	do j=1,NX
    		do i=1,NY
	        txy(i,j,k)=(dyg_2(i,j,k)+dxg_3(i,j,k))*mugRE(i,j,k)
         	txz(i,j,k)=(dzg_2(i,j,k)+dxg_4(i,j,k))*mugRE(i,j,k)
    	  	tyz(i,j,k)=(dzg_3(i,j,k)+dyg_4(i,j,k))*mugRE(i,j,k)
    		txx(i,j,k)=tvd2*(2d0*dxg_2(i,j,k)-dyg_3(i,j,k)-dzg_4(i,j,k))*mugRE(i,j,k)
  		tyy(i,j,k)=tvd2*(2d0*dyg_3(i,j,k)-dxg_2(i,j,k)-dzg_4(i,j,k))*mugRE(i,j,k)
		tzz(i,j,k)=tvd2*(2d0*dzg_4(i,j,k)-dyg_3(i,j,k)-dxg_2(i,j,k))*mugRe(i,j,k)
			end do
		end do
	end do
	!$omp end parallel do
	!$omp parallel do
    do k=1,NZ
       do i=1,yk      !!!xhoukouw
          iryu_x_1(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp1(i,1:xk,k),LU_D_p_nxk(1:xk,:),-sig)&
          						+fBDCS(1,xk,gdx,Fxm1(i,1:xk,k),LU_D_m_nxk(1:xk,:),sig))
          iryu_x_2(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp2(i,1:xk,k),LU_D_p_nxk(1:xk,:),-sig)&
          						+fBDCS(1,xk,gdx,Fxm2(i,1:xk,k),LU_D_m_nxk(1:xk,:),sig))
          iryu_x_3(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp3(i,1:xk,k),LU_D_p_nxk(1:xk,:),-sig)&
          						+fBDCS(1,xk,gdx,Fxm3(i,1:xk,k),LU_D_m_nxk(1:xk,:),sig))
          iryu_x_4(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp4(i,1:xk,k),LU_D_p_nxk(1:xk,:),-sig)&
          						+fBDCS(1,xk,gdx,Fxm4(i,1:xk,k),LU_D_m_nxk(1:xk,:),sig))
          iryu_x_5(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp5(i,1:xk,k),LU_D_p_nxk(1:xk,:),-sig)&
          						+fBDCS(1,xk,gdx,Fxm5(i,1:xk,k),LU_D_m_nxk(1:xk,:),sig))
       end do
       do i=yk+1,NY       !!!xhoukouw
          iryu_x_1(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp1(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm1(i,:,k),LU_D_m_nx,sig))
          iryu_x_2(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp2(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm2(i,:,k),LU_D_m_nx,sig))
          iryu_x_3(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp3(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm3(i,:,k),LU_D_m_nx,sig))
          iryu_x_4(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp4(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm4(i,:,k),LU_D_m_nx,sig))
          iryu_x_5(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp5(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm5(i,:,k),LU_D_m_nx,sig))
       end do
      end do
	!$omp end parallel do
     !$omp parallel do
    do k=1,NZ
       do j=1,xk-1       !!!yhoukou
          iryu_y_1(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp1(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym1(:,j,k),LU_D_m_ny,sig))
          iryu_y_2(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp2(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym2(:,j,k),LU_D_m_ny,sig))
          iryu_y_3(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp3(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym3(:,j,k),LU_D_m_ny,sig))
          iryu_y_4(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp4(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym4(:,j,k),LU_D_m_ny,sig))
          iryu_y_5(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp5(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym5(:,j,k),LU_D_m_ny,sig))
       end do
       do j=xk,NX       !!!yhoukou
          iryu_y_1(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp1(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym1(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_2(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp2(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym2(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_3(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp3(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym3(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_4(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp4(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym4(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_5(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp5(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym5(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
       end do
    end do
	!$omp end parallel do
	!$omp parallel do

    do j=1,NX
       do i=1,NY          !!!z方向移流項
          iryu_z_1(i,j,:)=DCSz(Fzp1(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm1(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_2(i,j,:)=DCSz(Fzp2(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm2(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_3(i,j,:)=DCSz(Fzp3(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm3(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_4(i,j,:)=DCSz(Fzp4(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm4(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_5(i,j,:)=DCSz(Fzp5(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm5(i,j,:),Lz_Dm,Uz_Dm,sig)
       end do
    end do
	!$omp end parallel do
	!$omp parallel do

    do k=1,NZ
       do i=1,yk
          dxtxx(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,txx(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxtxy(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,txy(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxtxz(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,txz(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxq_x(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,q_x(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxQ5(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,Q1_5(i,1:xk,k),LU_C_NXk(1:xk,1:3))
       end do
       do i=yk+1,NY
          dxtxx(i,:,k)=dx_k*CCS(1,NX,gdx,txx(i,:,k),LU_C_nx)
          dxtxy(i,:,k)=dx_k*CCS(1,NX,gdx,txy(i,:,k),LU_C_nx)
          dxtxz(i,:,k)=dx_k*CCS(1,NX,gdx,txz(i,:,k),LU_C_nx)
          dxq_x(i,:,k)=dx_k*CCS(1,NX,gdx,q_x(i,:,k),LU_C_nx)
          dxQ5(i,:,k)=dx_k*CCS(1,NX,gdx,Q1_5(i,:,k),LU_C_NX)
       end do
    end do

    !$omp end parallel do
    !$omp parallel do
    do k=1,NZ
       do j=1,xk-1
          dytyy(:,j,k)=dy_k*CCS(1,NY,gdy,tyy(:,j,k),LU_C_ny)
          dytxy(:,j,k)=dy_k*CCS(1,NY,gdy,txy(:,j,k),LU_C_ny)
          dytyz(:,j,k)=dy_k*CCS(1,NY,gdy,tyz(:,j,k),LU_C_ny)
          dyq_y(:,j,k)=dy_k*CCS(1,NY,gdy,q_y(:,j,k),LU_C_ny)
          dyQ5(:,j,k)=dy_k*CCS(1,NY,gdy,Q1_5(:,j,k),LU_C_ny)
       end do
       do j=xk,NX
          dytyy(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,tyy(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dytxy(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,txy(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dytyz(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,tyz(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dyq_y(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,q_y(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dyQ5(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,Q1_5(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
       end do
    end do
	!$omp end parallel do
	!$omp parallel do

    do j=1,NX
       do i=1,NY
          dztxz(i,j,:)=CCSz(txz(i,j,:),Lz_C,Uz_C)
          dztyz(i,j,:)=CCSz(tyz(i,j,:),Lz_C,Uz_C)
          dztzz(i,j,:)=CCSz(tzz(i,j,:),Lz_C,Uz_C)
          dzq_z(i,j,:)=CCSz(q_z(i,j,:),Lz_C,Uz_C)
       end do
    end do
	!$omp end parallel do
	!$omp parallel do


    do k=1,NZ
    	do j=1,NX
    		do i=1,NY
	    dFx2(i,j,k)=dxtxx(i,j,k)
    	dFx3(i,j,k)=dxtxy(i,j,k)
    	dFx4(i,j,k)=dxtxz(i,j,k)
    	dFx5(i,j,k)=dxg_2(i,j,k)*txx(i,j,k)+g_2(i,j,k)*dxtxx(i,j,k)&
         +dxg_3(i,j,k)*txy(i,j,k)+g_3(i,j,k)*dxtxy(i,j,k)&
         +dxg_4(i,j,k)*txz(i,j,k)+g_4(i,j,k)*dxtxz(i,j,k) -dxq_x(i,j,k)
    !!ryuusyutu jouken
    dFx3(i,NX,k)=0d0
    dFx4(i,NX,k)=0d0
    dFx5(i,NX,k)=dFx5(i,NX,k)+dxq_x(i,NX,k)-g_3(i,NX,k)*dxtxy(i,NX,k)-g_4(i,NX,k)*dxtxz(i,NX,k)

    dFx1(i,j,k)=dxg_1(i,j,k)*Ux(i,j,k)
    dFx2(i,j,k)=dFx2(i,j,k)+(dxg_1(i,j,k)*g_2(i,j,k)+dxg_2(i,j,k)*g_1(i,j,k))*Ux(i,j,k)
    dFx3(i,j,k)=dFx3(i,j,k)+(dxg_1(i,j,k)*g_3(i,j,k)+dxg_3(i,j,k)*g_1(i,j,k))*Ux(i,j,k)
    dFx4(i,j,k)=dFx4(i,j,k)+(dxg_1(i,j,k)*g_4(i,j,k)+dxg_4(i,j,k)*g_1(i,j,k))*Ux(i,j,k)
    dFx5(i,j,k)=dFx5(i,j,k)+dxQ5(i,j,k)*Ux(i,j,k)

    dFy2(i,j,k)=dytxy(i,j,k)
    dFy3(i,j,k)=dytyy(i,j,k)
    dFy4(i,j,k)=dytyz(i,j,k)
    dFy5(i,j,k)=dyg_2(i,j,k)*txy(i,j,k)+g_2(i,j,k)*dytxy(i,j,k)&
         +dyg_3(i,j,k)*tyy(i,j,k) + g_3(i,j,k)*dytyy(i,j,k)&
         +dyg_4(i,j,k)*tyz(i,j,k) + g_4(i,j,k)*dytyz(i,j,k) -dyq_y(i,j,k)
	!!ryusyutu jouken
	dFy2(NY,j,k)=0d0
	dFy4(NY,j,k)=0d0
	dFy5(NY,j,k)=dFy5(NY,j,k)+dyq_y(NY,j,k)-g_2(NY,j,k)*dytxy(NY,j,k)-g_4(NY,j,k)*dytyz(NY,j,k)

    dFy1(i,j,k)=dyg_1(i,j,k)*Uy(i,j,k)
    dFy2(i,j,k)=dFy2(i,j,k)+(dyg_1(i,j,k)*g_2(i,j,k)+dyg_2(i,j,k)*g_1(i,j,k))*Uy(i,j,k)
    dFy3(i,j,k)=dFy3(i,j,k)+(dyg_1(i,j,k)*g_3(i,j,k)+dyg_3(i,j,k)*g_1(i,j,k))*Uy(i,j,k)
    dFy4(i,j,k)=dFy4(i,j,k)+(dyg_1(i,j,k)*g_4(i,j,k)+dyg_4(i,j,k)*g_1(i,j,k))*Uy(i,j,k)
    dFy5(i,j,k)=dFy5(i,j,k)+(dyQ5(i,j,k))*Uy(i,j,k)


	dFz1(i,j,k)=0d0
    dFz2(i,j,k)=dztxz(i,j,k)
    dFz3(i,j,k)=dztyz(i,j,k)
    dFz4(i,j,k)=dztzz(i,j,k)
    dFz5(i,j,k)=dzg_2(i,j,k)*txz(i,j,k)+g_2(i,j,k)*dztxz(i,j,k)&
         +dzg_3(i,j,k)*tyz(i,j,k)+g_3(i,j,k)*dztyz(i,j,k)&
         +dzg_4(i,j,k)*tzz(i,j,k)+g_4(i,j,k)*dztzz(i,j,k)-dzq_z(i,j,k)
    		end do
    	end do
    end do

!$omp end parallel do
    call NSCBC_xhoukou(g_1,g_2,g_3,g_4,g_5,dxg_1,dxg_2,dxg_3,dxg_4,dxg_5,iryu_x_1,iryu_x_2,iryu_x_3,iryu_x_4,iryu_x_5,c)
    call NSCBC_yhoukou(g_1,g_2,g_3,g_4,g_5,dyg_1,dyg_2,dyg_3,dyg_4,dyg_5,iryu_y_1,iryu_y_2,iryu_y_3,iryu_y_4,iryu_y_5,c)

!    Q2=3d0/4d0*Q+0.25d0*(Q1-dt*(iryu_x+iryu_y+iryu_z-(dFx+dFy+dFz)-Q1*sigxy+Q0))

 !$omp parallel do
  do k=1,NZ
  	do j=1,NX
  		do i=1,NY
	  	Q2_1(i,j,k)=0.75d0*Q_1(i,j,k)+0.25d0*Q1_1(i,j,k)-0.25d0*dt*(iryu_x_1(i,j,k)+iryu_y_1(i,j,k)+iryu_z_1(i,j,k)&
	  				-dFx1(i,j,k)-dFy1(i,j,k)-dFz1(i,j,k)-sigxy(i,j,k)*Q1_1(i,j,k)+Q01(i,j,k))
	  	Q2_2(i,j,k)=0.75d0*Q_2(i,j,k)+0.25d0*Q1_2(i,j,k)-0.25d0*dt*(iryu_x_2(i,j,k)+iryu_y_2(i,j,k)+iryu_z_2(i,j,k)&
	  				-dFx2(i,j,k)-dFy2(i,j,k)-dFz2(i,j,k)-sigxy(i,j,k)*Q1_2(i,j,k)+Q02(i,j,k))
	  	Q2_3(i,j,k)=0.75d0*Q_3(i,j,k)+0.25d0*Q1_3(i,j,k)-0.25d0*dt*(iryu_x_3(i,j,k)+iryu_y_3(i,j,k)+iryu_z_3(i,j,k)&
	  				-dFx3(i,j,k)-dFy3(i,j,k)-dFz3(i,j,k)-sigxy(i,j,k)*Q1_3(i,j,k)+Q03(i,j,k))
	  	Q2_4(i,j,k)=0.75d0*Q_4(i,j,k)+0.25d0*Q1_4(i,j,k)-0.25d0*dt*(iryu_x_4(i,j,k)+iryu_y_4(i,j,k)+iryu_z_4(i,j,k)&
	  				-dFx4(i,j,k)-dFy4(i,j,k)-dFz4(i,j,k)-sigxy(i,j,k)*Q1_4(i,j,k)+Q04(i,j,k))
	  	Q2_5(i,j,k)=0.75d0*Q_5(i,j,k)+0.25d0*Q1_5(i,j,k)-0.25d0*dt*(iryu_x_5(i,j,k)+iryu_y_5(i,j,k)+iryu_z_5(i,j,k)&
	  				-dFx5(i,j,k)-dFy5(i,j,k)-dFz5(i,j,k)-sigxy(i,j,k)*Q1_5(i,j,k)+Q05(i,j,k))
		end do
	end do
   end do
  !$omp end parallel do


  call Q_g_in(Q2_1,Q2_2,Q2_3,Q2_4,Q2_5,g_1,g_2,g_3,g_4,g_5,T_k,mu,mugRE,mukq,t,y_k,g_in,vin)

	!$omp parallel do
    do k=1,NZ
    	do j=1,NX
    		do i=1,NY
    		  c(i,j,k)=dsqrt(T_k(i,j,k)/(0.09d0))
!!!x
       Fx1(i,j,k)=Q2_2(i,j,k)
       Fx2(i,j,k)=Q2_2(i,j,k)*g_2(i,j,k)+g_5(i,j,k)
       Fx3(i,j,k)=Q2_2(i,j,k)*g_3(i,j,k)
       Fx4(i,j,k)=Q2_2(i,j,k)*g_4(i,j,k)
       Fx5(i,j,k)=(Q2_5(i,j,k)+g_5(i,j,k))*g_2(i,j,k)

       Fy1(i,j,k)=Q2_3(i,j,k)
       Fy2(i,j,k)=Q2_3(i,j,k)*g_2(i,j,k)
       Fy3(i,j,k)=Q2_3(i,j,k)*g_3(i,j,k)+g_5(i,j,k)
       Fy4(i,j,k)=Q2_3(i,j,k)*g_4(i,j,k)
       Fy5(i,j,k)=(Q2_5(i,j,k)+g_5(i,j,k))*g_3(i,j,k)
    !!!Lax F
       Fz1(i,j,k)=Q2_4(i,j,k)
       Fz2(i,j,k)=g_2(i,j,k)*Q2_4(i,j,k)
       Fz3(i,j,k)=g_3(i,j,k)*Q2_4(i,j,k)
       Fz4(i,j,k)=Q2_4(i,j,k)*g_4(i,j,k)+g_5(i,j,k)
       Fz5(i,j,k)=(Q2_5(i,j,k)+g_5(i,j,k))*g_4(i,j,k)

       Fxp1(i,j,k)=(Fx1(i,j,k)+Q2_1(i,j,k))*0.5d0
       Fxp2(i,j,k)=(Fx2(i,j,k)+Q2_2(i,j,k))*0.5d0
       Fxp3(i,j,k)=(Fx3(i,j,k)+Q2_3(i,j,k))*0.5d0
       Fxp4(i,j,k)=(Fx4(i,j,k)+Q2_4(i,j,k))*0.5d0
       Fxp5(i,j,k)=(Fx5(i,j,k)+Q2_5(i,j,k))*0.5d0

       Fxm1(i,j,k)=(Fx1(i,j,k)-Q2_1(i,j,k))*0.5d0
       Fxm2(i,j,k)=(Fx2(i,j,k)-Q2_2(i,j,k))*0.5d0
       Fxm3(i,j,k)=(Fx3(i,j,k)-Q2_3(i,j,k))*0.5d0
       Fxm4(i,j,k)=(Fx4(i,j,k)-Q2_4(i,j,k))*0.5d0
       Fxm5(i,j,k)=(Fx5(i,j,k)-Q2_5(i,j,k))*0.5d0


       Fyp1(i,j,k)=(Fy1(i,j,k)+Q2_1(i,j,k))*0.5d0
       Fyp2(i,j,k)=(Fy2(i,j,k)+Q2_2(i,j,k))*0.5d0
       Fyp3(i,j,k)=(Fy3(i,j,k)+Q2_3(i,j,k))*0.5d0
       Fyp4(i,j,k)=(Fy4(i,j,k)+Q2_4(i,j,k))*0.5d0
       Fyp5(i,j,k)=(Fy5(i,j,k)+Q2_5(i,j,k))*0.5d0

       Fym1(i,j,k)=(Fy1(i,j,k)-Q2_1(i,j,k))*0.5d0
       Fym2(i,j,k)=(Fy2(i,j,k)-Q2_2(i,j,k))*0.5d0
       Fym3(i,j,k)=(Fy3(i,j,k)-Q2_3(i,j,k))*0.5d0
       Fym4(i,j,k)=(Fy4(i,j,k)-Q2_4(i,j,k))*0.5d0
       Fym5(i,j,k)=(Fy5(i,j,k)-Q2_5(i,j,k))*0.5d0

!!!!z
       Fzp1(i,j,k)=(Fz1(i,j,k)+Q2_1(i,j,k))*0.5d0
       Fzp2(i,j,k)=(Fz2(i,j,k)+Q2_2(i,j,k))*0.5d0
       Fzp3(i,j,k)=(Fz3(i,j,k)+Q2_3(i,j,k))*0.5d0
       Fzp4(i,j,k)=(Fz4(i,j,k)+Q2_4(i,j,k))*0.5d0
       Fzp5(i,j,k)=(Fz5(i,j,k)+Q2_5(i,j,k))*0.5d0

       Fzm1(i,j,k)=(Fz1(i,j,k)-Q2_1(i,j,k))*0.5d0
       Fzm2(i,j,k)=(Fz2(i,j,k)-Q2_2(i,j,k))*0.5d0
       Fzm3(i,j,k)=(Fz3(i,j,k)-Q2_3(i,j,k))*0.5d0
       Fzm4(i,j,k)=(Fz4(i,j,k)-Q2_4(i,j,k))*0.5d0
       Fzm5(i,j,k)=(Fz5(i,j,k)-Q2_5(i,j,k))*0.5d0
       			end do
       		end do
       end do
	!$omp end parallel do
	!$omp parallel do

    do k=1,NZ
        do i=1,yk
           dxg_1(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_1(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_2(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_2(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_3(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_3(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_4(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_4(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_5(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_5(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           q_x(i,1:xk,k)=dx_k(1:xk)*mukq(i,1:xk,k)*CCS(1,xk,gdx,T_k(i,1:xk,k),LU_C_nxk(1:xk,1:3))
     q_x(i,xk,k)=0d0

        end do
        do i=yk+1,NY
           dxg_1(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_1(i,:,k),LU_C_nx)
           dxg_2(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_2(i,:,k),LU_C_nx)
           dxg_3(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_3(i,:,k),LU_C_nx)
           dxg_4(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_4(i,:,k),LU_C_nx)
           dxg_5(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_5(i,:,k),LU_C_nx)
           q_x(i,:,k)=dx_k(:)*mukq(i,:,k)*CCS(1,NX,gdx,T_k(i,:,k),LU_C_nx)
        end do
    end do

    !$omp end parallel do
    !$omp parallel do
    do k=1,NZ
	    do i=1,xk-1
              dyg_1(:,i,k)=dy_k*CCS(1,NY,gdy,g_1(:,i,k),LU_C_ny)
              dyg_2(:,i,k)=dy_k*CCS(1,NY,gdy,g_2(:,i,k),LU_C_ny)
              dyg_3(:,i,k)=dy_k*CCS(1,NY,gdy,g_3(:,i,k),LU_C_ny)
              dyg_4(:,i,k)=dy_k*CCS(1,NY,gdy,g_4(:,i,k),LU_C_ny)
              dyg_5(:,i,k)=dy_k*CCS(1,NY,gdy,g_5(:,i,k),LU_C_ny)
           	  q_y(:,i,k)=mukq(:,i,k)*dy_k*CCS(1,NY,gdy,T_k(:,i,k),LU_C_ny)
     q_y(1,i,k)=0d0
        end do
        do i=xk,NX
              dyg_1(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_1(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_2(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_2(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_3(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_3(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_4(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_4(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_5(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_5(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              q_y(yk:NY,i,k)=mukq(yk:NY,i,k)*dy_k(yk:NY)*CCS(yk,NY,gdy,T_k(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
     q_y(yk,i,k)=0d0
        end do
     end do

	!$omp end parallel do
	!$omp parallel do

     do i=1,NX
        do j=1,NY
              dzg_1(j,i,:)=CCSz(g_1(j,i,:),Lz_C,Uz_C)
              dzg_2(j,i,:)=CCSz(g_2(j,i,:),Lz_C,Uz_C)
              dzg_3(j,i,:)=CCSz(g_3(j,i,:),Lz_C,Uz_C)
              dzg_4(j,i,:)=CCSz(g_4(j,i,:),Lz_C,Uz_C)
              dzg_5(j,i,:)=CCSz(g_5(j,i,:),Lz_C,Uz_C)
           q_z(j,i,:)=mukq(j,i,:)*CCSz(T_k(j,i,:),Lz_C,Uz_C)
        end do
     end do
    !$omp end parallel do
	!$omp parallel do

    do k=1,NZ
    	do j=1,NX
    		do i=1,NY
	        txy(i,j,k)=(dyg_2(i,j,k)+dxg_3(i,j,k))*mugRE(i,j,k)
         	txz(i,j,k)=(dzg_2(i,j,k)+dxg_4(i,j,k))*mugRE(i,j,k)
    	 	tyz(i,j,k)=(dzg_3(i,j,k)+dyg_4(i,j,k))*mugRE(i,j,k)
    		txx(i,j,k)=tvd2*(2d0*dxg_2(i,j,k)-dyg_3(i,j,k)-dzg_4(i,j,k))*mugRE(i,j,k)
  		  	tyy(i,j,k)=tvd2*(2d0*dyg_3(i,j,k)-dxg_2(i,j,k)-dzg_4(i,j,k))*mugRE(i,j,k)
		    tzz(i,j,k)=tvd2*(2d0*dzg_4(i,j,k)-dyg_3(i,j,k)-dxg_2(i,j,k))*mugRe(i,j,k)
			end do
		end do
	end do
	!$omp end parallel do
	!$omp parallel do
    do k=1,NZ
       do i=1,yk      !!!xhoukouw
          iryu_x_1(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp1(i,1:xk,k),LU_D_p_nxk,-sig)&
          						+fBDCS(1,xk,gdx,Fxm1(i,1:xk,k),LU_D_m_nxk,sig))
          iryu_x_2(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp2(i,1:xk,k),LU_D_p_nxk,-sig)&
          						+fBDCS(1,xk,gdx,Fxm2(i,1:xk,k),LU_D_m_nxk,sig))
          iryu_x_3(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp3(i,1:xk,k),LU_D_p_nxk,-sig)&
          						+fBDCS(1,xk,gdx,Fxm3(i,1:xk,k),LU_D_m_nxk,sig))
          iryu_x_4(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp4(i,1:xk,k),LU_D_p_nxk,-sig)&
          						+fBDCS(1,xk,gdx,Fxm4(i,1:xk,k),LU_D_m_nxk,sig))
          iryu_x_5(i,1:xk,k)=dx_k(1:xk)*(fBDCS(1,xk,gdx,Fxp5(i,1:xk,k),LU_D_p_nxk,-sig)&
          						+fBDCS(1,xk,gdx,Fxm5(i,1:xk,k),LU_D_m_nxk,sig))
       end do
       do i=yk+1,NY       !!!xhoukouw
          iryu_x_1(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp1(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm1(i,:,k),LU_D_m_nx,sig))
          iryu_x_2(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp2(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm2(i,:,k),LU_D_m_nx,sig))
          iryu_x_3(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp3(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm3(i,:,k),LU_D_m_nx,sig))
          iryu_x_4(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp4(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm4(i,:,k),LU_D_m_nx,sig))
          iryu_x_5(i,:,k)=dx_k(:)*(fBDCS(1,NX,gdx,Fxp5(i,:,k),LU_D_p_nx,-sig)+fBDCS(1,NX,gdx,Fxm5(i,:,k),LU_D_m_nx,sig))
       end do
      end do

    !$omp end parallel do
    !$omp parallel do
    do k=1,NZ
       do j=1,xk-1       !!!yhoukou
          iryu_y_1(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp1(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym1(:,j,k),LU_D_m_ny,sig))
          iryu_y_2(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp2(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym2(:,j,k),LU_D_m_ny,sig))
          iryu_y_3(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp3(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym3(:,j,k),LU_D_m_ny,sig))
          iryu_y_4(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp4(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym4(:,j,k),LU_D_m_ny,sig))
          iryu_y_5(:,j,k)=dy_k(:)*(fBDCS(1,NY,gdy,Fyp5(:,j,k),LU_D_p_ny,-sig)+fBDCS(1,NY,gdy,Fym5(:,j,k),LU_D_m_ny,sig))
       end do
       do j=xk,NX       !!!yhoukou
          iryu_y_1(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp1(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym1(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_2(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp2(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym2(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_3(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp3(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym3(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_4(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp4(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym4(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
          iryu_y_5(yk:NY,j,k)=dy_k(yk:NY)*(fBDCS(yk,NY,gdy,Fyp5(yk:NY,j,k),LU_D_p_nyk(yk:NY,1:3),-sig)&
          						+fBDCS(yk,NY,gdy,Fym5(yk:NY,j,k),LU_D_m_nyk(yk:NY,1:3),sig))
       end do
    end do
	!$omp end parallel do
	!$omp parallel do

    do j=1,NX
       do i=1,NY          !!!z方向移流項
          iryu_z_1(i,j,:)=DCSz(Fzp1(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm1(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_2(i,j,:)=DCSz(Fzp2(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm2(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_3(i,j,:)=DCSz(Fzp3(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm3(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_4(i,j,:)=DCSz(Fzp4(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm4(i,j,:),Lz_Dm,Uz_Dm,sig)
          iryu_z_5(i,j,:)=DCSz(Fzp5(i,j,:),Lz_Dp,Uz_Dp,-sig)+DCSz(Fzm5(i,j,:),Lz_Dm,Uz_Dm,sig)
       end do
    end do
    !$omp end parallel do
    !$omp parallel do
    do k=1,NZ
       do i=1,yk
          dxtxx(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,txx(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxtxy(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,txy(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxtxz(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,txz(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxq_x(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,q_x(i,1:xk,k),LU_C_nxk(1:xk,1:3))
          dxQ5(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,Q2_5(i,1:xk,k),LU_C_NXk(1:xk,1:3))
       end do
       do i=yk+1,NY
          dxtxx(i,:,k)=dx_k*CCS(1,NX,gdx,txx(i,:,k),LU_C_nx)
          dxtxy(i,:,k)=dx_k*CCS(1,NX,gdx,txy(i,:,k),LU_C_nx)
          dxtxz(i,:,k)=dx_k*CCS(1,NX,gdx,txz(i,:,k),LU_C_nx)
          dxq_x(i,:,k)=dx_k*CCS(1,NX,gdx,q_x(i,:,k),LU_C_nx)
          dxQ5(i,:,k)=dx_k*CCS(1,NX,gdx,Q2_5(i,:,k),LU_C_NX)
       end do
    end do

    !$omp end parallel do
    !$omp parallel do
    do k=1,NZ
       do j=1,xk-1
          dytyy(:,j,k)=dy_k*CCS(1,NY,gdy,tyy(:,j,k),LU_C_ny)
          dytxy(:,j,k)=dy_k*CCS(1,NY,gdy,txy(:,j,k),LU_C_ny)
          dytyz(:,j,k)=dy_k*CCS(1,NY,gdy,tyz(:,j,k),LU_C_ny)
          dyq_y(:,j,k)=dy_k*CCS(1,NY,gdy,q_y(:,j,k),LU_C_ny)
          dyQ5(:,j,k)=dy_k*CCS(1,NY,gdy,Q2_5(:,j,k),LU_C_ny)
       end do
       do j=xk,NX
          dytyy(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,tyy(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dytxy(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,txy(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dytyz(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,tyz(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dyq_y(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,q_y(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
          dyQ5(yk:NY,j,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,Q2_5(yk:NY,j,k),LU_C_nyk(yk:NY,1:3))
       end do
    end do
	!$omp end parallel do
	!$omp parallel do

    do j=1,NX
       do i=1,NY
          dztxz(i,j,:)=CCSz(txz(i,j,:),Lz_C,Uz_C)
          dztyz(i,j,:)=CCSz(tyz(i,j,:),Lz_C,Uz_C)
          dztzz(i,j,:)=CCSz(tzz(i,j,:),Lz_C,Uz_C)
          dzq_z(i,j,:)=CCSz(q_z(i,j,:),Lz_C,Uz_C)
       end do
    end do
	!$omp end parallel do
	!$omp parallel do
    do k=1,NZ
    	do j=1,NX
    		do i=1,NY
	    dFx2(i,j,k)=dxtxx(i,j,k)
    	dFx3(i,j,k)=dxtxy(i,j,k)
    	dFx4(i,j,k)=dxtxz(i,j,k)
    	dFx5(i,j,k)=dxg_2(i,j,k)*txx(i,j,k)+g_2(i,j,k)*dxtxx(i,j,k)&
         +dxg_3(i,j,k)*txy(i,j,k)+g_3(i,j,k)*dxtxy(i,j,k)&
         +dxg_4(i,j,k)*txz(i,j,k)+g_4(i,j,k)*dxtxz(i,j,k) -dxq_x(i,j,k)
    !!ryuusyutu jouken
    dFx3(i,NX,k)=0d0
    dFx4(i,NX,k)=0d0
    dFx5(i,NX,k)=dFx5(i,NX,k)+dxq_x(i,NX,k)-g_3(i,NX,k)*dxtxy(i,NX,k)-g_4(i,NX,k)*dxtxz(i,NX,k)

    dFx1(i,j,k)=dxg_1(i,j,k)*Ux(i,j,k)
    dFx2(i,j,k)=dFx2(i,j,k)+(dxg_1(i,j,k)*g_2(i,j,k)+dxg_2(i,j,k)*g_1(i,j,k))*Ux(i,j,k)
    dFx3(i,j,k)=dFx3(i,j,k)+(dxg_1(i,j,k)*g_3(i,j,k)+dxg_3(i,j,k)*g_1(i,j,k))*Ux(i,j,k)
    dFx4(i,j,k)=dFx4(i,j,k)+(dxg_1(i,j,k)*g_4(i,j,k)+dxg_4(i,j,k)*g_1(i,j,k))*Ux(i,j,k)
    dFx5(i,j,k)=dFx5(i,j,k)+dxQ5(i,j,k)*Ux(i,j,k)

    dFy2(i,j,k)=dytxy(i,j,k)
    dFy3(i,j,k)=dytyy(i,j,k)
    dFy4(i,j,k)=dytyz(i,j,k)
    dFy5(i,j,k)=dyg_2(i,j,k)*txy(i,j,k)+g_2(i,j,k)*dytxy(i,j,k)&
         +dyg_3(i,j,k)*tyy(i,j,k) + g_3(i,j,k)*dytyy(i,j,k)&
         +dyg_4(i,j,k)*tyz(i,j,k) + g_4(i,j,k)*dytyz(i,j,k) -dyq_y(i,j,k)
	!!ryusyutu jouken
	dFy2(NY,j,k)=0d0
	dFy4(NY,j,k)=0d0
	dFy5(NY,j,k)=dFy5(NY,j,k)+dyq_y(NY,j,k)-g_2(NY,j,k)*dytxy(NY,j,k)-g_4(NY,j,k)*dytyz(NY,j,k)

    dFy1(i,j,k)=dyg_1(i,j,k)*Uy(i,j,k)
    dFy2(i,j,k)=dFy2(i,j,k)+(dyg_1(i,j,k)*g_2(i,j,k)+dyg_2(i,j,k)*g_1(i,j,k))*Uy(i,j,k)
    dFy3(i,j,k)=dFy3(i,j,k)+(dyg_1(i,j,k)*g_3(i,j,k)+dyg_3(i,j,k)*g_1(i,j,k))*Uy(i,j,k)
    dFy4(i,j,k)=dFy4(i,j,k)+(dyg_1(i,j,k)*g_4(i,j,k)+dyg_4(i,j,k)*g_1(i,j,k))*Uy(i,j,k)
    dFy5(i,j,k)=dFy5(i,j,k)+(dyQ5(i,j,k))*Uy(i,j,k)


	dFz1(i,j,k)=0d0
    dFz2(i,j,k)=dztxz(i,j,k)
    dFz3(i,j,k)=dztyz(i,j,k)
    dFz4(i,j,k)=dztzz(i,j,k)
    dFz5(i,j,k)=dzg_2(i,j,k)*txz(i,j,k)+g_2(i,j,k)*dztxz(i,j,k)&
         +dzg_3(i,j,k)*tyz(i,j,k)+g_3(i,j,k)*dztyz(i,j,k)&
         +dzg_4(i,j,k)*tzz(i,j,k)+g_4(i,j,k)*dztzz(i,j,k)-dzq_z(i,j,k)
    		end do
    	end do
    end do

!$omp end parallel do
    call NSCBC_xhoukou(g_1,g_2,g_3,g_4,g_5,dxg_1,dxg_2,dxg_3,dxg_4,dxg_5,iryu_x_1,iryu_x_2,iryu_x_3,iryu_x_4,iryu_x_5,c)
    call NSCBC_yhoukou(g_1,g_2,g_3,g_4,g_5,dyg_1,dyg_2,dyg_3,dyg_4,dyg_5,iryu_y_1,iryu_y_2,iryu_y_3,iryu_y_4,iryu_y_5,c)

  !!   Q=tvd1*Q+tvd2*(Q2-dt*(iryu_x+iryu_y+iryu_z-(dFx+dFy+dFz)-sigxy*Q2+Q0))
!$omp parallel do
  do k=1,NZ
  	do j=1,NX
  		do i=1,NY
	  	Q_1(i,j,k)=tvd1*Q_1(i,j,k) +tvd2*Q2_1(i,j,k)-tvd2*dt*(iryu_x_1(i,j,k)+iryu_y_1(i,j,k)+iryu_z_1(i,j,k)&
	  				-dFx1(i,j,k)-dFy1(i,j,k)-dFz1(i,j,k)-sigxy(i,j,k)*Q2_1(i,j,k)+Q01(i,j,k))
	  	Q_2(i,j,k)=tvd1*Q_2(i,j,k)+tvd2*Q2_2(i,j,k)-tvd2*dt*(iryu_x_2(i,j,k)+iryu_y_2(i,j,k)+iryu_z_2(i,j,k)&
	  				-dFx2(i,j,k)-dFy2(i,j,k)-dFz2(i,j,k)-sigxy(i,j,k)*Q2_2(i,j,k)+Q02(i,j,k))
	  	Q_3(i,j,k)=tvd1*Q_3(i,j,k)+tvd2*Q2_3(i,j,k)-tvd2*dt*(iryu_x_3(i,j,k)+iryu_y_3(i,j,k)+iryu_z_3(i,j,k)&
	  				-dFx3(i,j,k)-dFy3(i,j,k)-dFz3(i,j,k)-sigxy(i,j,k)*Q2_3(i,j,k)+Q03(i,j,k))
	  	Q_4(i,j,k)=tvd1*Q_4(i,j,k)+tvd2*Q2_4(i,j,k)-tvd2*dt*(iryu_x_4(i,j,k)+iryu_y_4(i,j,k)+iryu_z_4(i,j,k)&
	  				-dFx4(i,j,k)-dFy4(i,j,k)-dFz4(i,j,k)-sigxy(i,j,k)*Q2_4(i,j,k)+Q04(i,j,k))
	  	Q_5(i,j,k)=tvd1*Q_5(i,j,k)+tvd2*Q2_5(i,j,k)-tvd2*dt*(iryu_x_5(i,j,k)+iryu_y_5(i,j,k)+iryu_z_5(i,j,k)&
	  				-dFx5(i,j,k)-dFy5(i,j,k)-dFz5(i,j,k)-sigxy(i,j,k)*Q2_5(i,j,k)+Q05(i,j,k))
		end do
	end do
   end do
  !$omp end parallel do

   call Q_g_in(Q_1,Q_2,Q_3,Q_4,Q_5,g_1,g_2,g_3,g_4,g_5,T_k,mu,mugRE,mukq,t,y_k,g_in,vin)


!!!===========================================================
     write(*,*)t
     write(222,'(8E26.15e3)')g_5(1,92,5),g_5(30,92,5),g_5(100,92,5),g_5(140,92,5)&
    ,g_5(170,92,5),g_5(NY-by-5,92,5),g_5(NY-by-2,92,5),g_5(NY-by,92,5)
     write(223,'(8E26.15e3)')g_5(1,410,5),g_5(30,410,5),g_5(100,410,5),g_5(140,410,5)&
    ,g_5(170,410,5),g_5(NY-by-5,410,5),g_5(NY-by-2,410,5),g_5(NY-by,410,5)
     write(224,'(8E26.15e3)')g_5(1,xk,5),g_5(30,xk,5),g_5(100,xk,5),g_5(140,xk,5)&
    ,g_5(170,xk,5),g_5(NY-by-5,xk,5),g_5(NY-by-2,xk,5),g_5(NY-by,xk,5)
     !!!===============================================
!!!計算ループ終点および結果出力
!!!===============================================!!!=====================================================
Mmod2=Mmod
if(t>MM)Mmod2=Mmod1
     if(M1<=t.and.t<=MM)then
  call result_final_2(g_1,g_2,g_3,g_4,g_5,T_k,t,g_h1,g_h2,g_h3,g_h4,g_h5,g_h6,x_k,y_k,avg_p)
      end if
     if(mod(t,Mmod2)==0)then
	!$omp parallel do
    do k=1,NZ
        do i=1,yk
           dxg_1(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_1(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_2(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_2(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_3(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_3(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_4(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_4(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_5(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_5(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           q_x(i,1:xk,k)=dx_k(1:xk)*mukq(i,1:xk,k)*CCS(1,xk,gdx,T_k(i,1:xk,k),LU_C_nxk(1:xk,1:3))
        end do
        do i=yk+1,NY
           dxg_1(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_1(i,:,k),LU_C_nx)
           dxg_2(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_2(i,:,k),LU_C_nx)
           dxg_3(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_3(i,:,k),LU_C_nx)
           dxg_4(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_4(i,:,k),LU_C_nx)
           dxg_5(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_5(i,:,k),LU_C_nx)
           q_x(i,:,k)=dx_k(:)*mukq(i,:,k)*CCS(1,NX,gdx,T_k(i,:,k),LU_C_nx)
        end do
!    end do
!!!    !?$omp do
!    do k=1,NZ
	    do i=1,xk-1
              dyg_1(:,i,k)=dy_k*CCS(1,NY,gdy,g_1(:,i,k),LU_C_ny)
              dyg_2(:,i,k)=dy_k*CCS(1,NY,gdy,g_2(:,i,k),LU_C_ny)
              dyg_3(:,i,k)=dy_k*CCS(1,NY,gdy,g_3(:,i,k),LU_C_ny)
              dyg_4(:,i,k)=dy_k*CCS(1,NY,gdy,g_4(:,i,k),LU_C_ny)
              dyg_5(:,i,k)=dy_k*CCS(1,NY,gdy,g_5(:,i,k),LU_C_ny)
           	  q_y(:,i,k)=mukq(:,i,k)*dy_k*CCS(1,NY,gdy,T_k(:,i,k),LU_C_ny)
        end do
        do i=xk,NX
              dyg_1(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_1(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_2(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_2(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_3(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_3(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_4(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_4(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_5(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_5(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              q_y(yk:NY,i,k)=mukq(yk:NY,i,k)*dy_k(yk:NY)*CCS(yk,NY,gdy,T_k(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
        end do
     end do
   	!$omp end parallel do
	!$omp parallel do
  do i=1,NX
        do j=1,NY
              dzg_1(j,i,:)=CCSz(g_1(j,i,:),Lz_C,Uz_C)
              dzg_2(j,i,:)=CCSz(g_2(j,i,:),Lz_C,Uz_C)
              dzg_3(j,i,:)=CCSz(g_3(j,i,:),Lz_C,Uz_C)
              dzg_4(j,i,:)=CCSz(g_4(j,i,:),Lz_C,Uz_C)
              dzg_5(j,i,:)=CCSz(g_5(j,i,:),Lz_C,Uz_C)
           q_z(j,i,:)=mukq(j,i,:)*CCSz(T_k(j,i,:),Lz_C,Uz_C)
        end do
     end do
    !$omp end parallel do

     call  result_2(t,T_k,g_old,g_1,g_2,g_3,g_4,g_5,&
     			dxg_1,dxg_2,dxg_3,dxg_4,dxg_5,&
     			dyg_1,dyg_2,dyg_3,dyg_4,dyg_5,&
     			dzg_1,dzg_2,dzg_3,dzg_4,dzg_5,LU_C_nx,LU_C_ny,LU_C_nxk,LU_C_nyk,x_k,dx_k,y_k,dy_k,avg_p)

  end if
  if(t==MM) call result_bunseki(G_H1,G_H2,G_H3,G_H4,G_H5,G_h6,LU_C_ny,LU_C_nyk,x_k,y_k,dy_k,rmsg1,rmsg2,rmsg3,rmsg4,rmsg5,rmsg6)
       if(t>=MM)then
        call result_rms(g_1,g_2,g_3,g_4,g_5,T_k,t,g_h1,g_h2,g_h3,g_h4,g_h5,g_h6,x_k,y_k,rmsg1,rmsg2,rmsg3,rmsg4,rmsg5,rmsg6)

     end if
     if(t==MMM)then
         do k=1,Nz
		g_1(:,:,k)=g_h1(:,:)
		g_2(:,:,k)=g_h2(:,:)
		g_3(:,:,k)=g_h3(:,:)
		g_4(:,:,k)=g_h4(:,:)
		g_5(:,:,k)=g_h5(:,:)
		T_k(:,:,k)=g_h6(:,:)
    	end do
	!$omp parallel do
  do k=1,NZ
        do i=1,yk
           dxg_1(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_1(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_2(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_2(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_3(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_3(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_4(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_4(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           dxg_5(i,1:xk,k)=dx_k(1:xk)*CCS(1,xk,gdx,g_5(i,1:xk,k),LU_C_nxk(1:xk,1:3))
           q_x(i,1:xk,k)=dx_k(1:xk)*mukq(i,1:xk,k)*CCS(1,xk,gdx,T_k(i,1:xk,k),LU_C_nxk(1:xk,1:3))
        end do
        do i=yk+1,NY
           dxg_1(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_1(i,:,k),LU_C_nx)
           dxg_2(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_2(i,:,k),LU_C_nx)
           dxg_3(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_3(i,:,k),LU_C_nx)
           dxg_4(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_4(i,:,k),LU_C_nx)
           dxg_5(i,:,k)=dx_k(:)*CCS(1,NX,gdx,g_5(i,:,k),LU_C_nx)
           q_x(i,:,k)=dx_k(:)*mukq(i,:,k)*CCS(1,NX,gdx,T_k(i,:,k),LU_C_nx)
        end do
!    end do
!!!    !?$omp do
!    do k=1,NZ
	    do i=1,xk-1
              dyg_1(:,i,k)=dy_k*CCS(1,NY,gdy,g_1(:,i,k),LU_C_ny)
              dyg_2(:,i,k)=dy_k*CCS(1,NY,gdy,g_2(:,i,k),LU_C_ny)
              dyg_3(:,i,k)=dy_k*CCS(1,NY,gdy,g_3(:,i,k),LU_C_ny)
              dyg_4(:,i,k)=dy_k*CCS(1,NY,gdy,g_4(:,i,k),LU_C_ny)
              dyg_5(:,i,k)=dy_k*CCS(1,NY,gdy,g_5(:,i,k),LU_C_ny)
           	  q_y(:,i,k)=mukq(:,i,k)*dy_k*CCS(1,NY,gdy,T_k(:,i,k),LU_C_ny)
        end do
        do i=xk,NX
              dyg_1(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_1(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_2(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_2(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_3(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_3(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_4(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_4(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              dyg_5(yk:NY,i,k)=dy_k(yk:NY)*CCS(yk,NY,gdy,g_5(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
              q_y(yk:NY,i,k)=mukq(yk:NY,i,k)*dy_k(yk:NY)*CCS(yk,NY,gdy,T_k(yk:NY,i,k),LU_C_nyk(yk:NY,1:3))
        end do
     end do

  	!$omp end parallel do
	!$omp parallel do
   do i=1,NX
        do j=1,NY
              dzg_1(j,i,:)=CCSz(g_1(j,i,:),Lz_C,Uz_C)
              dzg_2(j,i,:)=CCSz(g_2(j,i,:),Lz_C,Uz_C)
              dzg_3(j,i,:)=CCSz(g_3(j,i,:),Lz_C,Uz_C)
              dzg_4(j,i,:)=CCSz(g_4(j,i,:),Lz_C,Uz_C)
              dzg_5(j,i,:)=CCSz(g_5(j,i,:),Lz_C,Uz_C)
           q_z(j,i,:)=mukq(j,i,:)*CCSz(T_k(j,i,:),Lz_C,Uz_C)
        end do
     end do
    !$omp end parallel do



  	 call result_bunseki(G_H1,G_H2,G_H3,G_H4,G_H5,G_h6,LU_C_ny,LU_C_nyk,x_k,y_k,dy_k,rmsg1,rmsg2,rmsg3,rmsg4,rmsg5,rmsg6)
     call result_2(t,T_k,g_old,g_1,g_2,g_3,g_4,g_5,&
     			dxg_1,dxg_2,dxg_3,dxg_4,dxg_5,&
     			dyg_1,dyg_2,dyg_3,dyg_4,dyg_5,&
     			dzg_1,dzg_2,dzg_3,dzg_4,dzg_5,LU_C_nx,LU_C_ny,LU_C_nxk,LU_C_nyk,x_k,dx_k,y_k,dy_k,avg_p)

end if
!!!!NAN用　結果出力コード
     if(isnan(Q_1(1,1,1))) then
        write(*,*)'NAN'
        write(*,*)'NAN';stop
     end if

              end do
            end program main
