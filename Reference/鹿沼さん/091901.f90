module globals
implicit none
	real(8), parameter :: x0   = 0.0d0  !Xの最大値
!	real(8), parameter :: xmax = 1.0d0  !Xの最大値
	real(8), parameter :: y0   = 0.0d0  !yの最大値
!	real(8), parameter :: ymax = 1.0d0  !yの最大値
	real(8), parameter :: a1   = 1.0d0  !格子伸長（伸ばし）
	real(8), parameter :: b1   = 0.5d0  !格子伸長（細分）
!	integer, parameter  ::n1   =   200  ! 第1交点
	real(8), parameter :: Mach = 0.1d0  !最大マッハ数
	real(8), parameter :: p0   = 0.71552d0  !無限遠方の圧力
!	real(8), parameter :: T1   =1.4d0     !壁面温度
	real(8), parameter :: Temp0= 1.0d0   !入り口の無次元温度			!7/7
	real(8), parameter :: T0   = 288.0d0   !入り口の有次元温度			!7/7
	real(8), parameter :: T1   =1.0d0     !壁面温度						!7/7
	real(8), parameter :: ST   =111.0d0     !サザーランド定数			!7/7
	real(8), parameter :: gam   =1.4d0     !比熱比
	real(8), parameter :: Re   = 2.328d6   !Raynolds数
end module



module sub_compactP
use globals
!use DCS_Subroutine
implicit none
contains

subroutine CompactDCS_NP(LcNP, UcNP, f, x, vmax,sig, n)
real(8) ::   LcNP(0:n,0:1    )  !配列宣言　LcNP(n,n) 
real(8) ::   UcNP(0:n,0:1    )  !配列宣言　UcNP(n,n) 
real(8) ::     f(0:n)  !配列宣言　f(n) をもらう
real(8) ::     x(0:n)       !配列宣言　x(n) を返す
integer ::     n           !整数　nをもらう

real(8) AlNP,       arNP,brNP,crNP,drNP
real(8) AL4 , Be4,  ar4, br4 ,cr4, AL4_DCS , Be4_DCS,  ar4_DCS, br4_DCS ,cr4_DCS
real(8) Al  , Be ,  ar , br  ,cr , Al_DCS  , Be_DCS ,  ar_DCS , br_DCS  ,cr_DCS
real(8) xr, crs, brs, ars
!real(8) h, hinv, sig
real(8) h, vmax, hinv, sig

integer  i,nn

real(8) ::    a(0:n)
real(8) ::   x1(0:n)
 a(0:n) =0.0d0
x1(0:n) =0.0d0
 x(0:n) =0.0d0

!DCS1行目と最終行用(Non-Periodic)
	AlNP=3.0d0
	arNP=-2.8333333333333333d0	!=-17/6
	brNP=1.5d0					!=3/2
	crNP=1.5d0					!=3/2
	drNP=-0.1666666666666666d0	!-1/6

!DCS2行目と最終-1行用(Non-Periodic)
	Al4=0.25d0		!=1/4
	Be4=0.0d0		!=0
	ar4=1.5d0		!=3/2
	br4=0.0d0		!=0
	cr4=0.0d0		!=0
	Al4_DCS=1.0d0	!=1
	Be4_DCS=0.0d0	!=0
	ar4_DCS=0.5d0	!=1/2
	br4_DCS=0.0d0	!=0
	cr4_DCS=0.0d0	!=0
!DCS一般行用
	Al=1.0d0/3.0d0
	Be=0.0d0
	ar=14.0d0/9.0d0
	br=1.0d0/9.0d0
	cr=0.0d0
	Al_DCS=1.0d0
	Be_DCS=0.0d0
	ar_DCS=4.0d0/9.0d0
	br_DCS=2.0d0/9.0d0
	cr_DCS=0.0d0
!xmax=1.0d0 !後で修正
h=vmax/n
nn=dble(n)
hinv=nn/vmax

	i=0
			a(i)=(arNP*f(0)+brNP*f(1)+crNP*f(2)+drNP*f(3))/h
	i=1
			a(i)=ar4*(f(i+1)-f(i-1))/2.0d0/h+sig*ar4_DCS*(f(i-1)-2.0d0*f(i)+f(i+1))/h
	do i=2, n-2
			a(i)=ar*(f(i+1)-f(i-1))/2.0d0/h+br*(f(i+2)-f(i-2))/4.0d0/h	&
				+sig*(ar_DCS*(f(i-1)-2.0d0*f(i)+f(i+1))/h+br_DCS*(f(i-2)-2.0d0*f(i)+f(i+2))/4.0d0/h)
	enddo
	i=n-1
			a(i)=ar4*(f(i+1)-f(i-1))/2.0d0/h+sig*ar4_DCS*(f(i-1)-2.0d0*f(i)+f(i+1))/h
	i=n
			a(i)=-(arNP*f(n)+brNP*f(n-1)+crNP*f(n-2)+drNP*f(n-3))/h


!作成したUL行列から解を計算
!まずは下三角行列に掛かる解を出す
		x1(0)=a(0)
	do i=1,n
		x1(i)=a(i)-(LcNP(i,0)*x1(i-1))
	enddo


!上三角行列を使って本当の解を出す
		x(n)=x1(n)/UcNP(n,0)
	do i=n-1,0,-1
		x(i)=(x1(i)-UcNP(i,1)*x(i+1))/UcNP(i,0)
	enddo


end subroutine CompactDCS_NP






subroutine LU_DCS_NPcomp(n, LcNP, UcNP, sig)
integer ::    n           !整数　nをもらう
real(8) ::     LcNP(0:n, 0:1)       !配列宣言　L(n) を返す
real(8) ::     UcNP(0:n, 0:1)       !配列宣言　U(n) を返す
real(8) AlNP, AL4, AL4_DCS, AL4p, AL4m, AL, AL_DCS, ALp, ALm, sig
integer  i

!CS1行目と最終行用(Non-Periodic)
	AlNP=3.0d0
!CS2行目と最終-1行用(Non-Periodic)
	Al4=0.25d0
	Al4_DCS=1.0d0	!=1
	AL4p=AL4*(1.0d0+sig*Al4_DCS)
	AL4m=AL4*(1.0d0-sig*Al4_DCS)
!CS一般行用
	Al=1.0d0/3.0d0
	Al_DCS=1.0d0	!=1
	ALp=AL*(1.0d0+sig*Al_DCS)
	ALm=AL*(1.0d0-sig*Al_DCS)
	i=0
		LcNP(i,0)=1.0d0
		LcNP(i,1)=0.0d0

		UcNP(i,0)=1.0d0
		UcNP(i,1)=ALNP

	i=1
		LcNP(i,0)=AL4m
		LcNP(i,1)=1.0d0

		UcNP(i,0)=1.0d0-LcNP(i,0)*UcNP(i-1,1)
		UcNP(i,1)=AL4p

	do i=2, N-2
		LcNP(i,0)=ALm/UcNP(i-1,0)
		LcNP(i,1)=1.0d0

		UcNP(i,0)=1.0d0-LcNP(i,0)*UcNP(i-1,1)
		UcNP(i,1)=ALp
	enddo

	i=N-1
		LcNP(i,0)=AL4m/UcNP(i-1,0)
		LcNP(i,1)=1.0d0

		UcNP(i,0)=1.0d0-LcNP(i,0)*UcNP(i-1,1)
		UcNP(i,1)=AL4p

	i=N
		LcNP(i,0)=ALNP/UcNP(i-1,0)
		LcNP(i,1)=1.0d0

		UcNP(i,0)=1.0d0-LcNP(i,0)*UcNP(i-1,1)
		UcNP(i,1)=0.0d0

end subroutine LU_DCS_NPcomp





subroutine CompactNP(LcNP, UcNP, f, x, vmax,n)
real(8) ::   LcNP(0:n,0:1    )  !配列宣言　LcNP(n,n) 
real(8) ::   UcNP(0:n,0:1    )  !配列宣言　UcNP(n,n) 
real(8) ::     f(0:)  !配列宣言　f(n) をもらう
real(8) ::     x(0:)       !配列宣言　x(n) を返す
integer ::     n           !整数　nをもらう

real(8) AlNP,       arNP,brNP,crNP,drNP
real(8) AL4 , Be4,  ar4, br4 ,cr4
real(8) Al  , Be ,  ar , br  ,cr
real(8) xr, crs, brs, ars
real(8) h, vmax, hinv
!real(8) h, hinv

integer  i,nn

real(8) ::    a(0:n)
real(8) ::   x1(0:n)
 a(0:n) =0.0d0
x1(0:n) =0.0d0
 x(0:n) =0.0d0

!CS1行目と最終行用(Non-Periodic)
	AlNP=3.0d0
	arNP=-2.8333333333333333d0
	brNP=1.5d0
	crNP=1.5d0
	drNP=-0.1666666666666666d0

!CS2行目と最終-1行用(Non-Periodic)
	Al4=0.25d0
	Be4=0.0d0
	ar4=1.5d0
	br4=0.0d0
	cr4=0.0d0

!CS一般行用
	Al=1.0d0/3.0d0
	Be=0.0d0
	ar=14.0d0/9.0d0
	br=1.0d0/9.0d0
	cr=0.0d0
!xmax=1.0d0 !後で修正
h=vmax/n
nn=dble(n)
hinv=nn/vmax

	i=0
			a(i)=(arNP*f(0)+brNP*f(1)+crNP*f(2)+drNP*f(3))/h
	i=1
			a(i)=ar4*(f(i+1)-f(i-1))/2.0d0/h
	do i=2, n-2
			a(i)=ar*(f(i+1)-f(i-1))/2.0d0/h+br*(f(i+2)-f(i-2))/4.0d0/h
	enddo
	i=n-1
			a(i)=ar4*(f(i+1)-f(i-1))/2.0d0/h
	i=n
			a(i)=-(arNP*f(n)+brNP*f(n-1)+crNP*f(n-2)+drNP*f(n-3))/h


!作成したUL行列から解を計算
!まずは下三角行列に掛かる解を出す
		x1(0)=a(0)
	do i=1,n
		x1(i)=a(i)-(LcNP(i,0)*x1(i-1))
	enddo


!上三角行列を使って本当の解を出す
		x(n)=x1(n)/UcNP(n,0)
	do i=n-1,0,-1
		x(i)=(x1(i)-UcNP(i,1)*x(i+1))/UcNP(i,0)
	enddo


end subroutine CompactNP






subroutine LUNPcomp(n, LcNP, UcNP)
integer ::    n           !整数　nをもらう
real(8) ::     LcNP(0:n, 0:1)       !配列宣言　L(n) を返す
real(8) ::     UcNP(0:n, 0:1)       !配列宣言　U(n) を返す
real(8) AlNP, AL4, AL
integer  i

!CS1行目と最終行用(Non-Periodic)
	AlNP=3.0d0
!CS2行目と最終-1行用(Non-Periodic)
	Al4=0.25d0
!CS一般行用
	Al=1.0d0/3.0d0

	i=0
		LcNP(i,0)=1.0d0
		LcNP(i,1)=0.0d0

		UcNP(i,0)=1.0d0
		UcNP(i,1)=ALNP

	i=1
		LcNP(i,0)=AL4
		LcNP(i,1)=1.0d0

		UcNP(i,0)=1.0d0-LcNP(1,0)*ALNP
		UcNP(i,1)=AL4

	do i=2, N-2
		LcNP(i,0)=AL/UcNP(i-1,0)
		LcNP(i,1)=1.0d0

		UcNP(i,0)=1.0d0-LcNP(i,0)*UcNP(i-1,1)
		UcNP(i,1)=AL
	enddo
	i=N-1
		LcNP(i,0)=AL4/UcNP(i-1,0)
		LcNP(i,1)=1.0d0

		UcNP(i,0)=1.0d0-LcNP(i,0)*UcNP(i-1,1)
		UcNP(i,1)=AL4

	i=N
		LcNP(i,0)=ALNP/UcNP(i-1,0)
		LcNP(i,1)=1.0d0

		UcNP(i,0)=1.0d0-LcNP(i,0)*AL4
		UcNP(i,1)=0.0d0

end subroutine LUNPcomp


subroutine Newton_max(dYmax, D, JL, kappa)
integer :: i,   JL           !整数 JLをもらう
real(8) dYmax, D, kappa, Beta
integer,parameter :: itrmax = 10

kappa=1.0d0
Beta=dble(JL-1)/dble(JL)

	do i=1, itrmax
		kappa=kappa-(dYmax-D*(1.0d0-((exp(Beta*kappa)-1.0d0)/(exp(kappa)-1.0d0))))&
		/(D*(Beta*exp(Beta*kappa)/(exp(kappa)-1.0d0)+(exp(Beta*kappa)-1.0d0)*(-exp(kappa))/((exp(kappa)-1.0d0)**2)))
	enddo
!	write(*,*)  "kappa_in=", kappa
!	write(*,*)  "Beta=", beta
end subroutine Newton_max



subroutine Newton_min(dYmin, D, JL, kappa)
integer :: i,   JL           !整数 JLをもらう
real(8) dYmin, D, kappa, Alpha
integer,parameter :: itrmax = 10

kappa=1.0d0
Alpha=1.0d0/dble(JL)

	do i=1, itrmax
		kappa=kappa-(dYmin-D*((exp(Alpha*kappa)-1.0d0)/(exp(kappa)-1.0d0)))&
		/(-D*(Alpha*exp(Alpha*kappa)/(exp(kappa)-1.0d0)+(exp(Alpha*kappa)-1.0d0)*(-exp(kappa))/((exp(kappa)-1.0d0)**2)))
	enddo
!	write(*,*)  "kappa_in=", kappa
!	write(*,*)  "Alpha=", Alpha
end subroutine Newton_min


subroutine CompactDefNP(LcNP, UcNP, f, Defx, Defdf,vmax, n)
real(8) ::   LcNP(0:n,0:1    )  !配列宣言　LcNP(n,n) 
real(8) ::   UcNP(0:n,0:1    )  !配列宣言　UcNP(n,n)
real(8) ::        f(0:n)    !配列宣言　f(n) をもらう
real(8) ::     Defx(0:n)    !配列宣言　Defx(n) をもらう
real(8) ::    Defdf(0:n)    !配列宣言　Defdf(n) を返す
real(8) ::    vmax         !配列宣言　vmax をもらう
integer ::    n            !整数　nをもらう
integer ::    i            !配列宣言　整数l
real(8) ::    Defdx1(0:n)   !subroutine CompactNPで配列宣言　Defdx1(n)
real(8) ::         x(0:n)   !subroutine CompactNPで配列宣言　x(n) を返す

!			call CompactNP(LcNP, UcNP, Defx, Defdx1,       n)
			call CompactNP(LcNP, UcNP,    f,      x, vmax, n)
	do i=0,n		
			Defdf(i)=1.0d0/Defx(i)*x(i)
	enddo		
end subroutine CompactDefNP


subroutine CompactDefDCS_NP(L_d_m, U_d_m, f, Defx, Defdf,vmax,sig, n)
real(8) ::   L_d_m(0:n,0:1    )  !配列宣言　LcNP(n,n) 
real(8) ::   U_d_m(0:n,0:1    )  !配列宣言　UcNP(n,n)
real(8) ::        f(0:)    !配列宣言　f(n) をもらう
real(8) ::     Defx(0:)    !配列宣言　Defx(n) をもらう
real(8) ::    Defdf(0:)    !配列宣言　Defdf(n) を返す
real(8) ::    vmax         !vmax をもらう
integer ::    n            !整数　nをもらう
integer ::    i            !配列宣言　整数l
real(8) ::    Defdx1(0:n)   !subroutine CompactNPで配列宣言　Defdx1(n)
real(8) ::         x(0:n)   !subroutine CompactNPで配列宣言　x(n) を返す
real(8) ::    sig           !sigをもらう

!			call CompactNP(LcNP, UcNP,    f,      x, vmax, n)
			call CompactDCS_NP(L_d_m,U_d_m,f, x,vmax, sig,n)
	do i=0,n		
			Defdf(i)=1.0d0/Defx(i)*x(i)
	enddo		
end subroutine CompactDefDCS_NP

!サブステップの時間進行  ------------------------------------------------------------------------------------------------------

subroutine RK_substep(Q,&
			LcN1,  UcN1,  LcN2,  UcN2,  LcN3,  UcN3,  LcNT,  UcNT,&
			LcM1,  UcM1,  LcM2,  UcM2,  LcM3,  UcM3,  LcMT,  UcMT,&
			LcN1dp,UcN1dp,LcN2dp,UcN2dp,LcN3dp,UcN3dp,LcNTdp,UcNTdp,&
			LcM1dp,UcM1dp,LcM2dp,UcM2dp,LcM3dp,UcM3dp,LcMTdp,UcMTdp,&
			LcN1dm,UcN1dm,LcN2dm,UcN2dm,LcN3dm,UcN3dm,LcNTdm,UcNTdm,&
			LcM1dm,UcM1dm,LcM2dm,UcM2dm,LcM3dm,UcM3dm,LcMTdm,UcMTdm,&
			n1,n2,n3,nt,m1,m2,m3,mt,dt,sk,s_in,sout,dxd1,dxd11,dxd13,dyd1,dyd11,dyd13,c0,RK_Step,&
			xmax,ximax,ximax1,ximax3,ymax,ytmax,ytmax1,ytmax3)	

	real(8) :: Q(   1:4, 0:nt, 0:mt)	!もらう
	
	real(8) ::     LcN1(  0:  n1, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcN2(  0:  n2, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcN3(  0:  n3, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcNT(  0:  nt, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcM1(  0:  m1, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcM2(  0:  m2, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcM3(  0:  m3, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcMT(  0:  mt, 0:  1)       !配列宣言　LcNP(m)

	real(8) ::     UcN1(  0:  n1, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcN2(  0:  n2, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcN3(  0:  n3, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcNT(  0:  nt, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcM1(  0:  m1, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcM2(  0:  m2, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcM3(  0:  m3, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcMT(  0:  mt, 0:  1)       !配列宣言　UcNP(m)
	
	real(8) ::     LcN1dp(  0:  n1, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcN2dp(  0:  n2, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcN3dp(  0:  n3, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcNTdp(  0:  nt, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcM1dp(  0:  m1, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcM2dp(  0:  m2, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcM3dp(  0:  m3, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcMTdp(  0:  mt, 0:  1)       !配列宣言　LcNP(m)

	real(8) ::     LcN1dm(  0:  n1, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcN2dm(  0:  n2, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcN3dm(  0:  n3, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcNTdm(  0:  nt, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcM1dm(  0:  m1, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcM2dm(  0:  m2, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcM3dm(  0:  m3, 0:  1)       !配列宣言　LcNP(m)
	real(8) ::     LcMTdm(  0:  mt, 0:  1)       !配列宣言　LcNP(m)

	real(8) ::     UcN1dp(  0:  n1, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcN2dp(  0:  n2, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcN3dp(  0:  n3, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcNTdp(  0:  nt, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcM1dp(  0:  m1, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcM2dp(  0:  m2, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcM3dp(  0:  m3, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcMTdp(  0:  mt, 0:  1)       !配列宣言　UcNP(m)

	real(8) ::     UcN1dm(  0:  n1, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcN2dm(  0:  n2, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcN3dm(  0:  n3, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcNTdm(  0:  nt, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcM1dm(  0:  m1, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcM2dm(  0:  m2, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcM3dm(  0:  m3, 0:  1)       !配列宣言　UcNP(m)
	real(8) ::     UcMTdm(  0:  mt, 0:  1)       !配列宣言　UcNP(m)


real(8) :: sig
!real(8) :: dxd1(0:nt), dxd11(0:n1),dxd13(0:n3), dyd1(0:mt),dyd11(0:m1),dyd13(0:m3)
real(8) :: xmax,xmax1,xmax3,Ximax, Ximax1, Ximax2, Ximax3
real(8) :: ymax,ymax1,ymax3,Ytmax, Ytmax1, Ytmax2, Ytmax3

integer ::   n1,n2,n3,nt,m1,m2,m3,mt                 !整数　n, m をもらう------------------------------------------------- 

!RK用


real(8) :: F(   1:4, 0:nt, 0:mt)
real(8) :: G(   1:4, 0:nt, 0:mt)
real(8) :: Q23(      0:nt, 0:mt)

real(8) :: ff_F_p(   0:nt)
real(8) :: ff_F_m(   0:nt)
real(8) :: Fp(       0:nt)
real(8) :: Fm(       0:nt)
real(8) :: Fp3( 1:4, 0:nt, 0:mt)
real(8) :: Fm3( 1:4, 0:nt, 0:mt)
real(8) :: Fpm3(1:4, 0:nt, 0:mt)

real(8) :: ff_G_p(        0:mt)
real(8) :: ff_G_m(        0:mt)
real(8) :: Gp(            0:mt)
real(8) :: Gm(            0:mt)
real(8) :: Gp3( 1:4, 0:nt, 0:mt)
real(8) :: Gm3( 1:4, 0:nt, 0:mt)
real(8) :: Gpm3(1:4, 0:nt, 0:mt)

real(8) :: s_in( 1:4, 0:nt, 0:mt)
real(8) :: sout( 1:4, 0:nt, 0:mt)
real(8) ::   QS( 1:4, 0:nt, 0:mt)

!real(8)  gam, lam, sk, dt
real(8)  lam, sk, dt
real(8)   Pr, a0

!エラー検知
real(8)		daP, daRho	!圧力の平方根内を正に

integer :: RK_Step

integer  j, k, i, xi, yi, jt
integer  ix, iy
integer  ixshift, iyshift

!粘性用

real(8) ::dvx2(0:nt, 0:mt)
real(8) ::duy2(0:nt, 0:mt)
real(8) ::vx1(0:nt), uy1(0:mt), dvx1(0:nt), duy1(0:mt)

real(8) ::dux2(0:nt, 0:mt)
real(8) ::dvy2(0:nt, 0:mt)
real(8) ::ux1(0:nt), vy1(0:mt), dux1(0:nt), dvy1(0:mt)
real(8) ::dTx2(0:nt, 0:mt)
real(8) ::dTy2(0:nt, 0:mt)
real(8) ::Tx1(0:nt), Ty1(0:mt), dTx1(0:nt), dTy1(0:mt)


!粘性用
real(8) :: Fv(1:4, 0:nt, 0:mt)
real(8) :: Gv(1:4, 0:nt, 0:mt)
real(8) :: Temp(   0:nt, 0:mt)
real(8) :: dTemp(  0:nt, 0:mt)
real(8) :: myu(    0:nt, 0:mt)
real(8) Temp0,myu0

real(8) :: Fv1(0:nt),  dFv1(0:nt)
real(8) :: Gv1(0:mt),  dGv1(0:mt)
real(8) :: Fv_0(  4, 0:nt, 0:mt)
real(8) :: dFv3(1:4, 0:nt, 0:mt)
real(8) :: dFv3_1(4, 0:nt, 0:mt)
real(8) :: dFv3_2(4, 0:nt, 0:mt)
real(8) :: Gv_0(  4, 0:nt, 0:mt)
real(8) :: dGv3(1:4, 0:nt, 0:mt)
real(8) :: dGv3_1(4, 0:nt, 0:mt)
real(8) :: dGv3_2(4, 0:nt, 0:mt)

!NSCBC境界用
real(8) ::rhox1(0:nt), px1(0:nt), drhox1(0:nt), drhox2(0:nt, 0:mt), dpx1(0:nt), dpx2(0:nt, 0:mt)
real(8) ::rhoy1(0:mt), py1(0:mt), drhoy1(0:mt), drhoy2(0:nt, 0:mt), dpy1(0:mt), dpy2(0:nt, 0:mt)

real(8) ::LLx1(0:nt, 0:mt), LLx2(0:nt, 0:mt), LLx3(0:nt, 0:mt), LLx4(0:nt, 0:mt)
real(8) ::LLy1(0:nt, 0:mt), LLy2(0:nt, 0:mt), LLy3(0:nt, 0:mt), LLy4(0:nt, 0:mt)
real(8) ::Dx1(0:nt, 0:mt), Dx2(0:nt, 0:mt), Dx3(0:nt, 0:mt), Dx4(0:nt, 0:mt)
real(8) ::Dy1(0:nt, 0:mt), Dy2(0:nt, 0:mt), Dy3(0:nt, 0:mt), Dy4(0:nt, 0:mt)
real(8) :: rho(0:nt, 0:mt), u(0:nt, 0:mt), v(0:nt, 0:mt), p(0:nt, 0:mt), c0(0:nt, 0:mt)
real(8) ::LK, Sigma

!不等長格子用
real(8) :: dxd1(0:nt), dxd11(0:n1), dxd13(0:n3), dyd1(0:mt), dyd11(0:m1), dyd13(0:m3)


F(  1:4, 0:nt, 0:mt)=0.0d0
G(  1:4, 0:nt, 0:mt)=0.0d0
Q23(     0:nt, 0:mt)=0.0d0

ff_F_p(  0:nt)       =0.0d0
ff_F_m(  0:nt)       =0.0d0
Fp(      0:nt)       =0.0d0
Fm(      0:nt)       =0.0d0
Fp3(1:4, 0:nt, 0:mt)=0.0d0
Fm3(1:4, 0:nt, 0:mt)=0.0d0

ff_G_p(       0:mt)=0.0d0
ff_G_m(       0:mt)=0.0d0
Gp(           0:mt)=0.0d0
Gm(           0:mt)=0.0d0
Gp3(1:4, 0:nt, 0:mt)=0.0d0
Gm3(1:4, 0:nt, 0:mt)=0.0d0

QS( 1:4, 0:nt, 0:mt)=0.0d0

!NSCBC境界用
rho(     0:nt, 0:mt)=0.0d0
u(       0:nt, 0:mt)=0.0d0
v(       0:nt, 0:mt)=0.0d0
p(       0:nt, 0:mt)=0.0d0
c0(       0:nt, 0:mt)=0.0d0

LLx1(0:nt, 0:mt)=0.0d0
LLX2(0:nt, 0:mt)=0.0d0
LLX3(0:nt, 0:mt)=0.0d0
LLX4(0:nt, 0:mt)=0.0d0
 Dx1(0:nt, 0:mt)=0.0d0
 Dx2(0:nt, 0:mt)=0.0d0
 Dx3(0:nt, 0:mt)=0.0d0
 Dx4(0:nt, 0:mt)=0.0d0

LLy1(0:nt, 0:mt)=0.0d0
LLy2(0:nt, 0:mt)=0.0d0
LLy3(0:nt, 0:mt)=0.0d0
LLy4(0:nt, 0:mt)=0.0d0
 Dy1(0:nt, 0:mt)=0.0d0
 Dy2(0:nt, 0:mt)=0.0d0
 Dy3(0:nt, 0:mt)=0.0d0
 Dy4(0:nt, 0:mt)=0.0d0

!粘性用
Fv( 1:4, 0:nt, 0:mt)=0.0d0
Gv( 1:4, 0:nt, 0:mt)=0.0d0
Temp(    0:nt, 0:mt)=0.0d0
dTemp(   0:nt, 0:mt)=0.0d0
myu(     0:nt, 0:mt)=0.0d0

!粘性用
 vx1(    0:nt      )=0.0d0
 uy1(          0:mt)=0.0d0
 ux1(    0:nt      )=0.0d0
 vy1(          0:mt)=0.0d0
dvx2(    0:nt, 0:mt)=0.0d0
duy2(    0:nt, 0:mt)=0.0d0
dux2(    0:nt, 0:mt)=0.0d0
dvy2(    0:nt, 0:mt)=0.0d0

 Tx1(    0:nt      )=0.0d0
 Ty1(          0:mt)=0.0d0
DTx1(    0:nt      )=0.0d0
DTy1(          0:mt)=0.0d0
dTx2(    0:nt, 0:mt)=0.0d0
dTy2(    0:nt, 0:mt)=0.0d0



!gam=   1.4d0
!lam=   5.0d0
lam=   1.0d0
!Mach=  1.1d0
!Temp0= 1.4d0			!7/7
Temp0= 1.0d0			!7/7
myu0=  1.0d0
!Re=2.328d6
!Re=1.164d6
Pr=    0.7d0
a0=dsqrt(gam)

Sigma=0.25d0
!Sigma=0.0d0
LK=Sigma*(1.0d0-Mach*Mach)*a0   !--------------------注意
!p0=1.0d0


			xmax1=xmax*dble(n1)/dble(nt)
			xmax3=xmax*dble(n3)/dble(nt)
			ymax1=ymax*dble(m1)/dble(mt)
			ymax3=ymax*dble(m3)/dble(mt)


	!S1～4作成
	     !QS(k,i,j)の作成
!$omp parallel do private(iy,ix)   ! 7/7 
		do k=1,4
		   do iy=0,mt
			do ix=0,nt
				QS(k,ix,iy)=Q(k,ix,iy)+sk*dt*s_in(k,ix,iy)
			end do
		   end do
		end do


	!NSCBC用にQを一旦、元に戻す
!$omp parallel do private(iy)   ! 7/7 
	do ix=0, nt
		do iy=0, mt
			rho(ix, iy) =QS(1,ix, iy)
			  u(ix, iy) =QS(2,ix, iy)/QS(1,ix, iy)
			  v(ix, iy) =QS(3,ix, iy)/QS(1,ix, iy)
			  p(ix, iy) =(gam-1.0d0)*(QS(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
		enddo
	enddo

	!Tempの作成
		!一般部
!$omp parallel do private(i, daP, daRho)   ! 7/7
		do j=0,mt
			do i=0,nt
		   	   Temp(i,j)=gam*(gam-1.0d0)/QS(1,i,j)*(QS(4,i,j)-0.5d0*QS(2,i,j)*QS(2,i,j)/QS(1,i,j) &
		   	   -0.5d0*QS(3,i,j)*QS(3,i,j)/QS(1,i,j))
!		   	   c0(i,j)=dsqrt(gam*p(i, j)/rho(i, j))

				daP=dabs(p(i, j))
				daRho=dabs(rho(i, j))
		   	   c0(i,j)=dsqrt(gam*daP/daRho)

		if( p(i,j)<0.0d0)then
		write(*,*) i,j,p(i,j),daP,"p for Temp Error hassei!", RK_Step

		p(i,j)=daP		!要注意

		endif

		if( rho(i,j)<=0.0d0)then
		write(*,*) i,j, rho(i,j), daRho,"rho Error hassei!" , RK_Step

		rho(i,j)=daRho		!要注意
		QS(1,i,j)=daRho		!要注意

		endif

		if( isnan(c0(i,j)))then
		write(*,*)"c0 Error hassei!", RK_Step
		write(*,*) i,j,daP, p(i, j), rho(i, j);stop
		endif


			end do
		end do
!!ここから要チェック		
!!		!入力側境界
!		do iy=0,m
!			ix=0
!		   	   Temp(ix,iy)=gam    !p0もρも初期値は1のため
!!		   	   Temp(ix,iy)=gam*p(ix, iy)/rho(ix, iy)
!		end do
!!要チェックここまで

	     !dv/dx, du/dx, dT/dx, d(rho)/dx, dp/dx,d(ρ×u×v)/dx,d(ρ×u×u)/dx,d((e+p)×u)/dxの作成

!各領域毎の微分設定
	!XLW部
!$omp parallel do private(ix,vx1,ux1,Tx1,rhox1,px1,dvx1,dux1,dTx1,drhox1,dpx1) ! 7/7 
		do iy=0,m1-1
			do ix=0,nt
 			   vx1(ix)=   v(ix, iy)
 			   ux1(ix)=   u(ix, iy)
 			   Tx1(ix)=Temp(ix, iy)
 		     rhox1(ix)= rho(ix, iy) 
			   px1(ix)=   p(ix, iy)
			end do

			call CompactDefNP(LcNT, UcNT,     vx1, dxd1,    dvx1,ximax,nt)
			call CompactDefNP(LcNT, UcNT,     ux1, dxd1,    dux1,ximax,nt)
			call CompactDefNP(LcNT, UcNT,     Tx1, dxd1,    dTx1,ximax,nt)
			call CompactDefNP(LcNT, UcNT,   rhox1, dxd1,  drhox1,ximax,nt)
			call CompactDefNP(LcNT, UcNT,     px1, dxd1,    dpx1,ximax,nt)

			do ix=0,nt
			   dvx2(ix,iy)=  dvx1(ix)
			   dux2(ix,iy)=  dux1(ix)
			   dTx2(ix,iy)=  dTx1(ix)
			 drhox2(ix,iy)=drhox1(ix)
			   dpx2(ix,iy)=  dpx1(ix)
			end do
		end do

	!XML部
!$omp parallel do private(ix,vx1,ux1,Tx1,rhox1,px1,dvx1,dux1,dTx1,drhox1,dpx1) ! 7/7 
		do iy=m1,m1+m2
			do ix=0,n1
 			   vx1(ix)=   v(ix, iy)
 			   ux1(ix)=   u(ix, iy)
 			   Tx1(ix)=Temp(ix, iy)
 		     rhox1(ix)= rho(ix, iy) 
			   px1(ix)=   p(ix, iy)
			end do

			call CompactDefNP(LcN1, UcN1,     vx1, dxd11,    dvx1,ximax1,n1)
			call CompactDefNP(LcN1, UcN1,     ux1, dxd11,    dux1,ximax1,n1)
			call CompactDefNP(LcN1, UcN1,     Tx1, dxd11,    dTx1,ximax1,n1)
			call CompactDefNP(LcN1, UcN1,   rhox1, dxd11,  drhox1,ximax1,n1)
			call CompactDefNP(LcN1, UcN1,     px1, dxd11,    dpx1,ximax1,n1)

			do ix=0,n1
			   dvx2(ix,iy)=  dvx1(ix)
			   dux2(ix,iy)=  dux1(ix)
			   dTx2(ix,iy)=  dTx1(ix)
			 drhox2(ix,iy)=drhox1(ix)
			   dpx2(ix,iy)=  dpx1(ix)
			end do
		end do

	!XMR部
!$omp parallel do private(ix,ixshift,vx1,ux1,Tx1,rhox1,px1,dvx1,dux1,dTx1,drhox1,dpx1) ! 7/7 
		do iy=m1,m1+m2
			do ix=n1+n2,nt
				ixshift=ix-(n1+n2)
 			   vx1(ixshift)=   v(ix, iy)
 			   ux1(ixshift)=   u(ix, iy)
 			   Tx1(ixshift)=Temp(ix, iy)
 		     rhox1(ixshift)= rho(ix, iy) 
			   px1(ixshift)=   p(ix, iy)
			end do

			call CompactDefNP(LcN3, UcN3,     vx1, dxd13,    dvx1,ximax3,n3)
			call CompactDefNP(LcN3, UcN3,     ux1, dxd13,    dux1,ximax3,n3)
			call CompactDefNP(LcN3, UcN3,     Tx1, dxd13,    dTx1,ximax3,n3)
			call CompactDefNP(LcN3, UcN3,   rhox1, dxd13,  drhox1,ximax3,n3)
			call CompactDefNP(LcN3, UcN3,     px1, dxd13,    dpx1,ximax3,n3)

			do ix=n1+n2,nt
				ixshift=ix-(n1+n2)
			   dvx2(ix,iy)=  dvx1(ixshift)
			   dux2(ix,iy)=  dux1(ixshift)
			   dTx2(ix,iy)=  dTx1(ixshift)
			 drhox2(ix,iy)=drhox1(ixshift)
			   dpx2(ix,iy)=  dpx1(ixshift)
			end do
		end do

	!XUP部
!$omp parallel do private(ix,vx1,ux1,Tx1,rhox1,px1,dvx1,dux1,dTx1,drhox1,dpx1) ! 7/7 
		do iy=m1+m2+1,mt
			do ix=0,nt
 			   vx1(ix)=   v(ix, iy)
 			   ux1(ix)=   u(ix, iy)
 			   Tx1(ix)=Temp(ix, iy)
 		     rhox1(ix)= rho(ix, iy) 
			   px1(ix)=   p(ix, iy)
			end do

			call CompactDefNP(LcNT, UcNT,     vx1, dxd1,    dvx1,ximax,nt)
			call CompactDefNP(LcNT, UcNT,     ux1, dxd1,    dux1,ximax,nt)
			call CompactDefNP(LcNT, UcNT,     Tx1, dxd1,    dTx1,ximax,nt)
			call CompactDefNP(LcNT, UcNT,   rhox1, dxd1,  drhox1,ximax,nt)
			call CompactDefNP(LcNT, UcNT,     px1, dxd1,    dpx1,ximax,nt)

			do ix=0,nt
			   dvx2(ix,iy)=  dvx1(ix)
			   dux2(ix,iy)=  dux1(ix)
			   dTx2(ix,iy)=  dTx1(ix)
			 drhox2(ix,iy)=drhox1(ix)
			   dpx2(ix,iy)=  dpx1(ix)
			end do
		end do


 	     !dv/dy, du/dy, dT/dy, d(rho×v)/dy, dp/dy,d(ρ×u×v)/dy,d(ρ×v×v)/dy,d((e+p)×v)/dyの作成

	!YLH部
!$omp parallel do private(iy,vy1,uy1,Ty1,rhoy1,py1,dvy1,duy1,dTy1,drhoy1,dpy1)   ! 7/7
		do ix=0,n1-1																				!	0から？
			do iy=0,mt
			   vy1(iy)=   v(ix, iy)
			   uy1(iy)=   u(ix, iy)
			   Ty1(iy)=Temp(ix, iy)
 		     rhoy1(iy)= rho(ix, iy) 
			   py1(iy)=   p(ix, iy)
			end do

			call CompactDefNP(LcMT, UcMT,     vy1, dyd1,    dvy1,ytmax,mt)
			call CompactDefNP(LcMT, UcMT,     uy1, dyd1,    duy1,ytmax,mt)
			call CompactDefNP(LcMT, UcMT,     Ty1, dyd1,    dTy1,ytmax,mt)
			call CompactDefNP(LcMT, UcMT,   rhoy1, dyd1,  drhoy1,ytmax,mt)
			call CompactDefNP(LcMT, UcMT,     py1, dyd1,    dpy1,ytmax,mt)

			do iy=0,mt
			    dvy2(ix,iy)=  dvy1(iy)
			    duy2(ix,iy)=  duy1(iy)
			    dTy2(ix,iy)=  dTy1(iy)
			  drhoy2(ix,iy)=drhoy1(iy)
			    dpy2(ix,iy)=  dpy1(iy)
			end do
		end do

	!YML部
!$omp parallel do private(iy,vy1,uy1,Ty1,rhoy1,py1,dvy1,duy1,dTy1,drhoy1,dpy1)   ! 7/7
		do ix=n1,n1+n2
			do iy=0,m1
			   vy1(iy)=   v(ix, iy)
			   uy1(iy)=   u(ix, iy)
			   Ty1(iy)=Temp(ix, iy)
 		     rhoy1(iy)= rho(ix, iy) 
			   py1(iy)=   p(ix, iy)
			end do

			call CompactDefNP(LcM1, UcM1,     vy1, dyd11,    dvy1,ytmax1,m1)
			call CompactDefNP(LcM1, UcM1,     uy1, dyd11,    duy1,ytmax1,m1)
			call CompactDefNP(LcM1, UcM1,     Ty1, dyd11,    dTy1,ytmax1,m1)
			call CompactDefNP(LcM1, UcM1,   rhoy1, dyd11,  drhoy1,ytmax1,m1)
			call CompactDefNP(LcM1, UcM1,     py1, dyd11,    dpy1,ytmax1,m1)

			do iy=0,m1
			    dvy2(ix,iy)=  dvy1(iy)
			    duy2(ix,iy)=  duy1(iy)
			    dTy2(ix,iy)=  dTy1(iy)
			  drhoy2(ix,iy)=drhoy1(iy)
			    dpy2(ix,iy)=  dpy1(iy)
			end do
		end do


	!YMU部
!$omp parallel do private(iy,iyshift,vy1,uy1,Ty1,rhoy1,py1,dvy1,duy1,dTy1,drhoy1,dpy1)   ! 7/7
		do ix=n1,n1+n2
			do iy=m1+m2,mt
				iyshift=iy-(m1+m2)
			   vy1(iyshift)=   v(ix, iy)
			   uy1(iyshift)=   u(ix, iy)
			   Ty1(iyshift)=Temp(ix, iy)
 		     rhoy1(iyshift)= rho(ix, iy) 
			   py1(iyshift)=   p(ix, iy)
			end do

			call CompactDefNP(LcM3, UcM3,     vy1, dyd13,    dvy1,ytmax3,m3)
			call CompactDefNP(LcM3, UcM3,     uy1, dyd13,    duy1,ytmax3,m3)
			call CompactDefNP(LcM3, UcM3,     Ty1, dyd13,    dTy1,ytmax3,m3)
			call CompactDefNP(LcM3, UcM3,   rhoy1, dyd13,  drhoy1,ytmax3,m3)
			call CompactDefNP(LcM3, UcM3,     py1, dyd13,    dpy1,ytmax3,m3)

			do iy=m1+m2,mt
				iyshift=iy-(m1+m2)
			    dvy2(ix,iy)=  dvy1(iyshift)
			    duy2(ix,iy)=  duy1(iyshift)
			    dTy2(ix,iy)=  dTy1(iyshift)
			  drhoy2(ix,iy)=drhoy1(iyshift)
			    dpy2(ix,iy)=  dpy1(iyshift)
			end do
		end do


	!YRH部
!$omp parallel do private(iy,vy1,uy1,Ty1,rhoy1,py1,dvy1,duy1,dTy1,drhoy1,dpy1)   ! 7/7
		do ix=n1+n2+1,nt
			do iy=0,mt
			   vy1(iy)=   v(ix, iy)
			   uy1(iy)=   u(ix, iy)
			   Ty1(iy)=Temp(ix, iy)
 		     rhoy1(iy)= rho(ix, iy) 
			   py1(iy)=   p(ix, iy)
			end do

			call CompactDefNP(LcMT, UcMT,     vy1, dyd1,    dvy1,ytmax,mt)
			call CompactDefNP(LcMT, UcMT,     uy1, dyd1,    duy1,ytmax,mt)
			call CompactDefNP(LcMT, UcMT,     Ty1, dyd1,    dTy1,ytmax,mt)
			call CompactDefNP(LcMT, UcMT,   rhoy1, dyd1,  drhoy1,ytmax,mt)
			call CompactDefNP(LcMT, UcMT,     py1, dyd1,    dpy1,ytmax,mt)

			do iy=0,mt
			    dvy2(ix,iy)=  dvy1(iy)
			    duy2(ix,iy)=  duy1(iy)
			    dTy2(ix,iy)=  dTy1(iy)
			  drhoy2(ix,iy)=drhoy1(iy)
			    dpy2(ix,iy)=  dpy1(iy)
			end do
		end do




!全領域の物理量計算
!$omp parallel do private(i)  ! 7/7
	do j=0,mt
		do i=0,nt
			F(1,i,j)=QS(2,i,j)
			G(1,i,j)=QS(3,i,j)
				Q23(i,j)=-0.5d0/QS(1,i,j)*( QS(2,i,j)*QS(2,i,j)+QS(3,i,j)*QS(3,i,j) )

			F(2,i,j)=QS(2,i,j)*QS(2,i,j)/QS(1,i,j)+(gam-1.0d0)*(QS(4,i,j)+Q23(i,j))
			G(2,i,j)=QS(2,i,j)*QS(3,i,j)/QS(1,i,j)

			F(3,i,j)=G( 2,i,j)
			G(3,i,j)=QS(3,i,j)*QS(3,i,j)/QS(1,i,j)+(gam-1.0d0)*(QS(4,i,j)+Q23(i,j))

			F(4,i,j)=QS(2,i,j)/QS(1,i,j)*(gam*QS(4,i,j)+(gam-1.0d0)*Q23(i,j))
			G(4,i,j)=QS(3,i,j)/QS(1,i,j)*(gam*QS(4,i,j)+(gam-1.0d0)*Q23(i,j))

!			myu(i,j)=myu0*(Temp(i,j)**(2.0d0/3.0d0))/(Temp0**(2.0d0/3.0d0))
			myu(i,j)=(1.0d0+ST/T0)/(Temp(i,j)+ST/T0)*(Temp(i,j)**(3.0d0/2.0d0))		!7/7

!RHS(粘性項)の被微分項を作る			
			Fv(2,i,j)=2.0d0*myu(i,j)/(3.0d0*Re)*(2.0d0*dux2(i,j)-dvy2(i,j))
			Fv(3,i,j)=myu(i,j)/Re*(dvx2(i,j)+duy2(i,j))
			Fv(4,i,j)=u(i,j)*Fv(2,i,j)+v(i,j)*Fv(3,i,j)+myu(i,j)/((gam-1.0d0)*Re*Pr)*dTx2(i,j)
			Gv(2,i,j)=Fv(3,i,j)
			Gv(3,i,j)=2.0d0*myu(i,j)/(3.0d0*Re)*(2.0d0*dvy2(i,j)-dux2(i,j))
			Gv(4,i,j)=u(i,j)*Gv(2,i,j)+v(i,j)*Gv(3,i,j)+myu(i,j)/((gam-1.0d0)*Re*Pr)*dTy2(i,j)
		end do
	end do

!RHSのx方向微分、y方向微分を行う
!$omp parallel do private(iy,ix,ixshift,iyshift,Fv1,dFv1,Gv1,dGv1)   ! 5/18
	do k=1,4
!RHSのx方向微分

	!XLW部
		do iy=0,m1-1
			do ix=0,nt
			   Fv1(ix)=Fv(k,ix,iy)
			end do

			call CompactDefNP(LcNT, UcNT, Fv1,dxd1, dFv1,ximax,nt)
			do ix=0,nt
			   dFv3(k, ix, iy)=dFv1(ix)
			end do
		end do

	!XML部
		do iy=m1,m1+m2
			do ix=0,n1
			   Fv1(ix)=Fv(k,ix,iy)
			end do

			call CompactDefNP(LcN1, UcN1, Fv1,dxd11, dFv1,ximax1,n1)
			do ix=0,n1
			   dFv3(k, ix, iy)=dFv1(ix)
			end do
		end do

	!XMR部
		do iy=m1,m1+m2
			do ix=n1+n2,nt
				ixshift=ix-(n1+n2)
			   Fv1(ixshift)=Fv(k,ix,iy)
			end do

			call CompactDefNP(LcN3, UcN3, Fv1,dxd13, dFv1,ximax3,n3)
			do ix=n1+n2,nt
				ixshift=ix-(n1+n2)
			   dFv3(k, ix, iy)=dFv1(ixshift)
			end do
		end do

	!XUP部
		do iy=m1+m2+1,mt
			do ix=0,nt
			   Fv1(ix)=Fv(k,ix,iy)
			end do

			call CompactDefNP(LcNT, UcNT, Fv1,dxd1, dFv1,ximax,nt)
			do ix=0,nt
			   dFv3(k, ix, iy)=dFv1(ix)
			end do
		end do


!RHSのy方向微分

	!YLH部
		do ix=0,n1-1
			do iy=0,mt
			   Gv1(iy)=Gv(k,ix,iy)
			end do

			call CompactDefNP(LcMT, UcMT, Gv1,dyd1, dGv1,ytmax,mt)
			do iy=0,mt
			   dGv3(k, ix, iy)=dGv1(iy)
			end do
		end do


	!YML部
		do ix=n1,n1+n2
			do iy=0,m1
			   Gv1(iy)=Gv(k,ix,iy)
			end do

			call CompactDefNP(LcM1, UcM1, Gv1,dyd11, dGv1,ytmax1,m1)
			do iy=0,m1
			   dGv3(k, ix, iy)=dGv1(iy)
			end do
		end do



	!YMU部	
		do ix=n1,n1+n2
			do iy=m1+m2,mt
				iyshift=iy-(m1+m2)
			   Gv1(iyshift)=Gv(k,ix,iy)
			end do

			call CompactDefNP(LcM3, UcM3, Gv1,dyd13, dGv1,ytmax3,m3)
			do iy=m1+m2,mt
				iyshift=iy-(m1+m2)
			   dGv3(k, ix, iy)=dGv1(iyshift)
			end do
		end do



	!YRH部	
		do ix=n1+n2+1,nt
			do iy=0,mt
			   Gv1(iy)=Gv(k,ix,iy)
			end do

			call CompactDefNP(LcMT, UcMT, Gv1,dyd1, dGv1,ytmax,mt)
			do iy=0,mt
			   dGv3(k, ix, iy)=dGv1(iy)
			end do
		end do
	end do


!Lax-Friedrich
!$omp parallel do private(iy,ix,ixshift,iyshift, ff_F_p, ff_F_m, sig, Fp, Fm, ff_G_p, ff_G_m, Gp, Gm)   ! 5/18 
	do k=1,4
!x方向微分

	!XLW部
		do iy=0,m1-1
			do ix=0,nt
			   ff_F_p(ix)=-0.5d0*(F(k,ix,iy)+lam*QS(k,ix,iy))
			   ff_F_m(ix)=-0.5d0*(F(k,ix,iy)-lam*QS(k,ix,iy))
			end do

			sig = -0.25d0
			call CompactDefDCS_NP(LcNTdm, UcNTdm,ff_F_p,dxd1, Fp,ximax,sig,nt)
			sig = 0.25d0
			call CompactDefDCS_NP(LcNTdp, UcNTdp,ff_F_m,dxd1, Fm,ximax,sig,nt)

			do ix=0,nt
			   Fp3(k, ix, iy)=Fp(ix)
			   Fm3(k, ix, iy)=Fm(ix)
			  Fpm3(k, ix, iy)=Fp3(k, ix, iy)+Fm3(k, ix, iy)
			end do
		end do

	!XML部
		do iy=m1,m1+m2
			do ix=0,n1
			   ff_F_p(ix)=-0.5d0*(F(k,ix,iy)+lam*QS(k,ix,iy))
			   ff_F_m(ix)=-0.5d0*(F(k,ix,iy)-lam*QS(k,ix,iy))
			end do

			sig = -0.25d0
			call CompactDefDCS_NP(LcN1dm, UcN1dm,ff_F_p,dxd11, Fp,ximax1,sig,n1)
			sig = 0.25d0
			call CompactDefDCS_NP(LcN1dp, UcN1dp,ff_F_m,dxd11, Fm,ximax1,sig,n1)

			do ix=0,n1
			   Fp3(k, ix, iy)=Fp(ix)
			   Fm3(k, ix, iy)=Fm(ix)
			  Fpm3(k, ix, iy)=Fp3(k, ix, iy)+Fm3(k, ix, iy)
			end do
		end do

	!XMR部
		do iy=m1,m1+m2
			do ix=n1+n2,nt
				ixshift=ix-(n1+n2)
			   ff_F_p(ixshift)=-0.5d0*(F(k,ix,iy)+lam*QS(k,ix,iy))
			   ff_F_m(ixshift)=-0.5d0*(F(k,ix,iy)-lam*QS(k,ix,iy))
			end do

			sig = -0.25d0
			call CompactDefDCS_NP(LcN3dm, UcN3dm,ff_F_p,dxd13, Fp,ximax3,sig,n3)
			sig = 0.25d0
			call CompactDefDCS_NP(LcN3dp, UcN3dp,ff_F_m,dxd13, Fm,ximax3,sig,n3)

			do ix=n1+n2,nt
				ixshift=ix-(n1+n2)
			   Fp3(k, ix, iy)=Fp(ixshift)
			   Fm3(k, ix, iy)=Fm(ixshift)
			  Fpm3(k, ix, iy)=Fp3(k, ix, iy)+Fm3(k, ix, iy)
			end do
		end do

	!XUP部
		do iy=m1+m2+1,mt
			do ix=0,nt
			   ff_F_p(ix)=-0.5d0*(F(k,ix,iy)+lam*QS(k,ix,iy))
			   ff_F_m(ix)=-0.5d0*(F(k,ix,iy)-lam*QS(k,ix,iy))
			end do

			sig = -0.25d0
			call CompactDefDCS_NP(LcNTdm, UcNTdm,ff_F_p,dxd1, Fp,ximax,sig,nt)
			sig = 0.25d0
			call CompactDefDCS_NP(LcNTdp, UcNTdp,ff_F_m,dxd1, Fm,ximax,sig,nt)

			do ix=0,nt
			   Fp3(k, ix, iy)=Fp(ix)
			   Fm3(k, ix, iy)=Fm(ix)
			  Fpm3(k, ix, iy)=Fp3(k, ix, iy)+Fm3(k, ix, iy)
			end do
		end do

!y方向微分

	!YLH部
		do ix=0,n1-1
			do iy=0,mt
			   ff_G_p(iy)=-0.5d0*(G(k,ix,iy)+lam*QS(k,ix,iy))
			   ff_G_m(iy)=-0.5d0*(G(k,ix,iy)-lam*QS(k,ix,iy))
			end do

			sig = -0.25d0
			call CompactDefDCS_NP(LcMTdm, UcMTdm,ff_G_p,dyd1, Gp, ytmax, sig,mt)
			sig = 0.25d0
			call CompactDefDCS_NP(LcMTdp, UcMTdp,ff_G_m,dyd1, Gm, ytmax, sig,mt)

			do iy=0,mt
			   Gp3(k, ix, iy)=Gp(iy)
			   Gm3(k, ix, iy)=Gm(iy)
			  Gpm3(k, ix, iy)=Gp3(k, ix, iy)+Gm3(k, ix, iy)
			end do
		end do

	!YML部
		do ix=n1,n1+n2
			do iy=0,m1
			   ff_G_p(iy)=-0.5d0*(G(k,ix,iy)+lam*QS(k,ix,iy))
			   ff_G_m(iy)=-0.5d0*(G(k,ix,iy)-lam*QS(k,ix,iy))
			end do

			sig = -0.25d0
			call CompactDefDCS_NP(LcM1dm, UcM1dm,ff_G_p,dyd11, Gp, ytmax1, sig,m1)
			sig = 0.25d0
			call CompactDefDCS_NP(LcM1dp, UcM1dp,ff_G_m,dyd11, Gm, ytmax1, sig,m1)

			do iy=0,m1
			   Gp3(k, ix, iy)=Gp(iy)
			   Gm3(k, ix, iy)=Gm(iy)
			  Gpm3(k, ix, iy)=Gp3(k, ix, iy)+Gm3(k, ix, iy)
			end do
		end do


	!YMU部
		do ix=n1,n1+n2
			do iy=m1+m2,mt
				iyshift=iy-(m1+m2)
			   ff_G_p(iyshift)=-0.5d0*(G(k,ix,iy)+lam*QS(k,ix,iy))
			   ff_G_m(iyshift)=-0.5d0*(G(k,ix,iy)-lam*QS(k,ix,iy))
			end do

			sig = -0.25d0
			call CompactDefDCS_NP(LcM3dm, UcM3dm,ff_G_p,dyd13, Gp, ytmax3, sig,m3)
			sig = 0.25d0
			call CompactDefDCS_NP(LcM3dp, UcM3dp,ff_G_m,dyd13, Gm, ytmax3, sig,m3)

			do iy=m1+m2,mt
				iyshift=iy-(m1+m2)
			   Gp3(k, ix, iy)=Gp(iyshift)
			   Gm3(k, ix, iy)=Gm(iyshift)
			  Gpm3(k, ix, iy)=Gp3(k, ix, iy)+Gm3(k, ix, iy)
			end do
		end do


	!YRH部
		do ix=n1+n2+1,nt
			do iy=0,mt
			   ff_G_p(iy)=-0.5d0*(G(k,ix,iy)+lam*QS(k,ix,iy))
			   ff_G_m(iy)=-0.5d0*(G(k,ix,iy)-lam*QS(k,ix,iy))
			end do

			sig = -0.25d0
			call CompactDefDCS_NP(LcMTdm, UcMTdm,ff_G_p,dyd1, Gp, ytmax, sig,mt)
			sig = 0.25d0
			call CompactDefDCS_NP(LcMTdp, UcMTdp,ff_G_m,dyd1, Gm, ytmax, sig,mt)

			do iy=0,mt
			   Gp3(k, ix, iy)=Gp(iy)
			   Gm3(k, ix, iy)=Gm(iy)
			  Gpm3(k, ix, iy)=Gp3(k, ix, iy)+Gm3(k, ix, iy)
			end do
		end do
	end do

	!NSCBC用
	
	!粘性項の準備	
	!x方向
!		do j=0,mt
!			do i=0,nt
!!				Fv(4,i,j)=u(i,j)*Fv(2,i,j)+v(i,j)*Fv(3,i,j)-myu(i,j)/((gam-1.0d0)*Re*Pr)*dTx2(i,j)
!				Fv_0(4,i,j)=u(i,j)*Fv(2,i,j)   !+v(i,j)*Fv(3,i,j)
!			end do
!		end do

!X方向壁 Non-reflecting Outflow

!	ix=nt
!	XLW部右左壁
!!$omp parallel do private(ix, Fv1, dFv1, dFv3_1)  ! 5/18 
		do iy=0, m1-1
			do ix=0,nt
			   Fv1(ix)=u(ix,iy)*Fv(2,ix,iy)
			end do
			call CompactDefNP(LcNT, UcNT, Fv1,dxd1, dFv1,ximax,nt)
			   dFv3_1(4,  0, iy)=dFv1( 0)
			   dFv3_1(4, nt, iy)=dFv1(nt)

		end do

!	XMR部右壁
!!$omp parallel do private(ix, ixshift, Fv1, dFv1)  ! 5/18 
		do iy=m1,m1+m2
			do ix=n1+n2,nt
				ixshift=ix-(n1+n2)
			   Fv1(ixshift)=u(ix,iy)*Fv(2,ix,iy)
			end do
			call CompactDefNP(LcN3, UcN3, Fv1,dxd13, dFv1,ximax3,n3)
				ix=nt
				ixshift=ix-(n1+n2)
			   dFv3_1(4, ix, iy)=dFv1(ixshift)
		end do

!	XML部左壁
!!$omp parallel do private(ix, Fv1,dFv1, dFv3_1)   ! 5/18 
		do iy=m1, m1+m2
			do ix=0, n1
			   Fv1(ix)=u(ix,iy)*Fv(2,ix,iy)
			end do
			call CompactDefNP(LcN1, UcN1, Fv1,dxd11, dFv1,ximax1,n1)
			   dFv3_1(4, 0, iy)=dFv1(0)
		end do

!	XUP部右左壁
!!$omp parallel do private(ix, Fv1,dFv1, dFv3_1)   ! 5/18 
		do iy=m1+m2+1, mt
			do ix=0,nt
			   Fv1(ix)=u(ix,iy)*Fv(2,ix,iy)
			end do
			call CompactDefNP(LcNT, UcNT, Fv1,dxd1, dFv1,ximax,nt)
			   dFv3_1(4,  0, iy)=dFv1( 0)
			   dFv3_1(4, nt, iy)=dFv1(nt)
		end do

!	X方向 右左壁部
		do iy=0,mt
			   dFv3(4,  0, iy)=dFv3_1(4,  0, iy)   +dvx2( 0,iy)*Fv(3, 0,iy)	!05/29
!			   dFv3(4,  0, iy)=dFv3_1(4,  0, iy)   +dvx2( 0,iy)*Fv(3, 0,iy)
			   dFv3(3,  0, iy)=0.0d0		!x=0での無反射条件τ１２=０でFv(３,n,j)=０のため
			   dFv3(4, nt, iy)=dFv3_1(4, nt, iy)   +dvx2(nt,iy)*Fv(3,nt,iy)
			   dFv3(3, nt, iy)=0.0d0		!x=Nでの無反射条件τ１２=０でFv(３,n,j)=０のため
		end do

	ix=0
		do iy=0,mt
!			LLx1(ix,iy)=LK*(p(ix, iy)-p0)                                 !もとのﾌﾟﾛｸﾞﾗﾑでは!
			LLx1(ix,iy)=(u(ix, iy)-c0(ix, iy))*(dpx2(ix,iy)-rho(ix,iy)*c0(ix, iy)*dux2(ix,iy))			! 04/22 要チェック！

			LLx2(ix,iy)=u(ix, iy)*(c0(ix, iy)*c0(ix, iy)*drhox2(ix,iy)-dpx2(ix,iy))
			LLx3(ix,iy)=u(ix, iy)*dvx2(ix,iy)
!			LLx4(ix,iy)=(u(ix, iy)+c0(ix, iy))*(dpx2(ix,iy)+rho(ix,iy)*c0(ix, iy)*dux2(ix,iy))
			LLx4(ix,iy)=0.0d0

!			LLx1(ix,iy)=(QS(2,ix, iy)/QS(1,ix, iy)-c0(ix, iy))*(dpx2(ix,iy)-QS(1,ix,iy)*c0(ix, iy)*dux2(ix,iy))

			Dx1(ix,iy)=(LLx2(ix,iy)+0.5d0*(LLx1(ix,iy)+LLx4(ix,iy)))/(c0(ix, iy)*c0(ix, iy))
			Dx2(ix,iy)=(LLx4(ix,iy)+LLx1(ix,iy))*0.5d0
			Dx3(ix,iy)=(LLx4(ix,iy)-LLx1(ix,iy))/(2.0d0*rho(ix, iy)*c0(ix, iy))
			Dx4(ix,iy)=LLx3(ix,iy)	
			
		end do

		do iy=0,mt
		  	Fpm3(1,ix,iy)=-Dx1(ix,iy)
		  	Fpm3(2,ix,iy)=-u(ix, iy)*Dx1(ix,iy)-rho(ix,iy)*Dx3(ix,iy)
	 		Fpm3(3,ix,iy)=-v(ix, iy)*Dx1(ix,iy)-rho(ix,iy)*Dx4(ix,iy)
		  	Fpm3(4,ix,iy)=-(0.5d0*Dx1(ix,iy)*(u(ix,iy)*u(ix,iy)+v(ix,iy)*v(ix,iy))+Dx2(ix,iy)/(gam-1.0d0) &
		  					+rho(ix,iy)*u(ix,iy)*Dx3(ix,iy)+rho(ix,iy)*v(ix,iy)*Dx4(ix,iy))
		end do


	ix=nt
		do iy=0,mt
			LLx1(ix,iy)=0.0d0
!			LLx1(ix,iy)=LK*(p(ix, iy)-p0)                                 !もとのﾌﾟﾛｸﾞﾗﾑでは!

			LLx2(ix,iy)=u(ix, iy)*(c0(ix, iy)*c0(ix, iy)*drhox2(ix,iy)-dpx2(ix,iy))
			LLx3(ix,iy)=u(ix, iy)*dvx2(ix,iy)
			LLx4(ix,iy)=(u(ix, iy)+c0(ix, iy))*(dpx2(ix,iy)+rho(ix,iy)*c0(ix, iy)*dux2(ix,iy))

!			LLx1(ix,iy)=(QS(2,ix, iy)/QS(1,ix, iy)-c0(ix, iy))*(dpx2(ix,iy)-QS(1,ix,iy)*c0(ix, iy)*dux2(ix,iy))

			Dx1(ix,iy)=(LLx2(ix,iy)+0.5d0*(LLx1(ix,iy)+LLx4(ix,iy)))/(c0(ix, iy)*c0(ix, iy))
			Dx2(ix,iy)=(LLx4(ix,iy)+LLx1(ix,iy))*0.5d0
			Dx3(ix,iy)=(LLx4(ix,iy)-LLx1(ix,iy))/(2.0d0*rho(ix, iy)*c0(ix, iy))
			Dx4(ix,iy)=LLx3(ix,iy)	
			
		end do

		do iy=0,mt
		  	Fpm3(1,ix,iy)=-Dx1(ix,iy)
		  	Fpm3(2,ix,iy)=-u(ix, iy)*Dx1(ix,iy)-rho(ix,iy)*Dx3(ix,iy)
	 		Fpm3(3,ix,iy)=-v(ix, iy)*Dx1(ix,iy)-rho(ix,iy)*Dx4(ix,iy)
		  	Fpm3(4,ix,iy)=-(0.5d0*Dx1(ix,iy)*(u(ix,iy)*u(ix,iy)+v(ix,iy)*v(ix,iy))+Dx2(ix,iy)/(gam-1.0d0) &
		  					+rho(ix,iy)*u(ix,iy)*Dx3(ix,iy)+rho(ix,iy)*v(ix,iy)*Dx4(ix,iy))
		end do




!Y方向壁 Non-reflecting Outflow

!	iy=mt
!	YLH部下上壁
		do ix=0, n1-1
			do iy=0,mt
			   Gv1(iy)=v(ix,iy)*Gv(2,ix,iy)
			end do
			call CompactDefNP(LcMT, UcMT, Gv1,dyd1, dGv1,ytmax,mt)
			   dGv3_1(4, ix,  0)=dGv1( 0)
			   dGv3_1(4, ix, mt)=dGv1(mt)
		end do

!	YMU部上壁
		do ix=n1, n1+n2
			do iy=m1+m2, mt
				iyshift=iy-(m1+m2)
			   Gv1(iyshift)=v(ix,iy)*Gv(2,ix,iy)
			end do
			call CompactDefNP(LcM3, UcM3, Gv1,dyd13, dGv1,ytmax3,m3)
				iy=mt
				iyshift=iy-(m1+m2)
			   dGv3_1(4, ix, iy)=dGv1(iyshift)
		end do

!	YML部下壁
		do ix=n1, n1+n2
			do iy=0 ,m1
			   Gv1(iy)=v(ix,iy)*Gv(2,ix,iy)
			end do
			call CompactDefNP(LcM1, UcM1, Gv1,dyd11, dGv1,ytmax1,m1)
			   dGv3_1(4, ix,  0)=dGv1( 0)
		end do


!	YRH部下上壁
		do ix=n1+n2+1, nt
			do iy=0,mt
			   Gv1(iy)=v(ix,iy)*Gv(2,ix,iy)
			end do
			call CompactDefNP(LcMT, UcMT, Gv1,dyd1, dGv1,ytmax,mt)
			   dGv3_1(4, ix,  0)=dGv1( 0)
			   dGv3_1(4, ix, mt)=dGv1(mt)
		end do

!	Y方向 下上壁部
		do ix=0,nt
			   dGv3(4, ix,  0)=dGv3_1(4, ix,  0)   +duy2(ix,0)*Gv(2,ix,0)	!05/29
!			   dGv3(4, ix,  0)=dGv3_1(4, ix,  0)   +duy2(nt,iy)*Gv(2,nt,iy)
			   dGv3(2, ix,  0)=0.0d0		!y=0での無反射条件τ１２=０でGv(３,i,m)=０のため
			   dGv3(4, ix, mt)=dGv3_1(4, ix, mt)   +duy2(ix,mt)*Gv(2,ix,mt)
			   dGv3(2, ix, mt)=0.0d0		!y=Mでの無反射条件τ１２=０でGv(３,i,m)=０のため
		end do


		iy=0
			do ix=0,nt
!			LLy1(ix,iy)=0.0d0      !L1は超音速なので流入はなく0
!			LLy1(ix,iy)=LK*(p(ix, iy)-p0)                                 !もとのﾌﾟﾛｸﾞﾗﾑでは!
			LLy1(ix,iy)=(v(ix, iy)-c0(ix, iy))*(dpy2(ix,iy)-rho(ix,iy)*c0(ix, iy)*dvy2(ix,iy))

!ここから下要チェック
			LLy2(ix,iy)=v(ix, iy)*(c0(ix, iy)*c0(ix, iy)*drhoy2(ix,iy)-dpy2(ix,iy))
			LLy3(ix,iy)=v(ix, iy)*duy2(ix,iy)
!			LLy4(ix,iy)=(v(ix, iy)+c0(ix, iy))*(dpy2(ix,iy)+rho(ix,iy)*c0(ix, iy)*dvy2(ix,iy))
			LLy4(ix,iy)=0.0d0

!			LLy1(ix,iy)=(QS(2,ix, iy)/QS(1,ix, iy)-c0(ix, iy))*(dpx2(ix,iy)-QS(1,ix,iy)*c0(ix, iy)*dux2(ix,iy))


			Dy1(ix,iy)=(LLy2(ix,iy)+0.5d0*(LLy1(ix,iy)+LLy4(ix,iy)))/(c0(ix, iy)*c0(ix, iy))
			Dy2(ix,iy)=(LLy4(ix,iy)+LLy1(ix,iy))*0.5d0
			Dy3(ix,iy)=(LLy4(ix,iy)-LLy1(ix,iy))/(2.0d0*rho(ix, iy)*c0(ix, iy))
			Dy4(ix,iy)=LLy3(ix,iy)	
			
			end do

		do ix=0,nt
		  	Gpm3(1,ix,iy)=-Dy1(ix,iy)
		  	Gpm3(3,ix,iy)=-v(ix, iy)*Dy1(ix,iy)-rho(ix,iy)*Dy3(ix,iy)  !8/31
	 		Gpm3(2,ix,iy)=-u(ix, iy)*Dy1(ix,iy)-rho(ix,iy)*Dy4(ix,iy)   !-dpx2(ix,iy)  !8/31 !9/2
		  	Gpm3(4,ix,iy)=-(0.5d0*Dy1(ix,iy)*(u(ix,iy)*u(ix,iy)+v(ix,iy)*v(ix,iy))+Dy2(ix,iy)/(gam-1.0d0) &
		  					+rho(ix,iy)*v(ix,iy)*Dy3(ix,iy)+rho(ix,iy)*u(ix,iy)*Dy4(ix,iy))
		end do


	iy=mt
		do ix=0,nt
			LLy1(ix,iy)=0.0d0      !L1は超音速なので流入はなく0
!			LLy1(ix,iy)=LK*(p(ix, iy)-p0)                                 !もとのﾌﾟﾛｸﾞﾗﾑでは!

!ここから下要チェック
			LLy2(ix,iy)=v(ix, iy)*(c0(ix, iy)*c0(ix, iy)*drhoy2(ix,iy)-dpy2(ix,iy))
			LLy3(ix,iy)=v(ix, iy)*duy2(ix,iy)
			LLy4(ix,iy)=(v(ix, iy)+c0(ix, iy))*(dpy2(ix,iy)+rho(ix,iy)*c0(ix, iy)*dvy2(ix,iy))

!			LLy1(ix,iy)=(QS(2,ix, iy)/QS(1,ix, iy)-c0(ix, iy))*(dpx2(ix,iy)-QS(1,ix,iy)*c0(ix, iy)*dux2(ix,iy))

			Dy1(ix,iy)=(LLy2(ix,iy)+0.5d0*(LLy1(ix,iy)+LLy4(ix,iy)))/(c0(ix, iy)*c0(ix, iy))
			Dy2(ix,iy)=(LLy4(ix,iy)+LLy1(ix,iy))*0.5d0
			Dy3(ix,iy)=(LLy4(ix,iy)-LLy1(ix,iy))/(2.0d0*rho(ix, iy)*c0(ix, iy))
			Dy4(ix,iy)=LLy3(ix,iy)	
			
			end do

		do ix=0,nt
		  	Gpm3(1,ix,iy)=-Dy1(ix,iy)
		  	Gpm3(3,ix,iy)=-v(ix, iy)*Dy1(ix,iy)-rho(ix,iy)*Dy3(ix,iy)
	 		Gpm3(2,ix,iy)=-u(ix, iy)*Dy1(ix,iy)-rho(ix,iy)*Dy4(ix,iy)
		  	Gpm3(4,ix,iy)=-(0.5d0*Dy1(ix,iy)*(u(ix,iy)*u(ix,iy)+v(ix,iy)*v(ix,iy))+Dy2(ix,iy)/(gam-1.0d0) &
		  					+rho(ix,iy)*v(ix,iy)*Dy3(ix,iy)+rho(ix,iy)*u(ix,iy)*Dy4(ix,iy))
		end do


!四隅の扱い
		ix=0
		iy=0
		Gpm3(2, ix, iy)=0.0d0
		Gpm3(3, ix, iy)=0.0d0
		Fpm3(2, ix, iy)=-(u(ix,iy)*Dy1(ix,iy)+rho(ix,iy)*Dy4(ix,iy)+u(ix,iy)*Dx1(ix,iy)+rho(ix,iy)*Dx3(ix,iy))
		Fpm3(3, ix, iy)=-(v(ix,iy)*Dy1(ix,iy)+rho(ix,iy)*Dy3(ix,iy)+v(ix,iy)*Dx1(ix,iy)+rho(ix,iy)*Dx4(ix,iy))

		ix=nt
		iy=0
		Gpm3(2, ix, iy)=0.0d0
		Gpm3(3, ix, iy)=0.0d0
		Fpm3(2, ix, iy)=-(u(ix,iy)*Dy1(ix,iy)+rho(ix,iy)*Dy4(ix,iy)+u(ix,iy)*Dx1(ix,iy)+rho(ix,iy)*Dx3(ix,iy))
		Fpm3(3, ix, iy)=-(v(ix,iy)*Dy1(ix,iy)+rho(ix,iy)*Dy3(ix,iy)+v(ix,iy)*Dx1(ix,iy)+rho(ix,iy)*Dx4(ix,iy))
		
		ix=0
		iy=mt
		Gpm3(2, ix, iy)=0.0d0
		Gpm3(3, ix, iy)=0.0d0
		Fpm3(2, ix, iy)=-(u(ix,iy)*Dy1(ix,iy)+rho(ix,iy)*Dy4(ix,iy)+u(ix,iy)*Dx1(ix,iy)+rho(ix,iy)*Dx3(ix,iy))
		Fpm3(3, ix, iy)=-(v(ix,iy)*Dy1(ix,iy)+rho(ix,iy)*Dy3(ix,iy)+v(ix,iy)*Dx1(ix,iy)+rho(ix,iy)*Dx4(ix,iy))

		ix=nt
		iy=mt
		Gpm3(2, ix, iy)=0.0d0
		Gpm3(3, ix, iy)=0.0d0
		Fpm3(2, ix, iy)=-(u(ix,iy)*Dy1(ix,iy)+rho(ix,iy)*Dy4(ix,iy)+u(ix,iy)*Dx1(ix,iy)+rho(ix,iy)*Dx3(ix,iy))
		Fpm3(3, ix, iy)=-(v(ix,iy)*Dy1(ix,iy)+rho(ix,iy)*Dy3(ix,iy)+v(ix,iy)*Dx1(ix,iy)+rho(ix,iy)*Dx4(ix,iy))

!Isothermal No-Slip Wall
!左右方向壁
!LH壁
	ix=n1	
!!$omp parallel do private(k)
		do iy=m1,m1+m2
			u(ix,iy)=0.0d0
			v(ix,iy)=0.0d0
			LLx4(ix,iy)=(u(ix, iy)+c0(ix, iy))*(dpx2(ix,iy)+rho(ix,iy)*c0(ix, iy)*dux2(ix,iy))
			LLx3(ix,iy)=0.0d0
			LLx2(ix,iy)=0.0d0
			LLx1(ix,iy)=LLx4(ix,iy)
			Dx1(ix,iy)=(LLx2(ix,iy)+0.5d0*(LLx1(ix,iy)+LLx4(ix,iy)))/(c0(ix, iy)*c0(ix, iy))
		  	Fpm3(1,ix,iy)=-Dx1(ix,iy)
		  do k=2,4
		  	Fpm3(k, ix, iy)=0.0d0
		  	Gpm3(k, ix, iy)=0.0d0
			dFv3(k, ix, iy)=0.0d0
			dGv3(k, ix, iy)=0.0d0
		  end do
		  	Gpm3(1, ix, iy)=0.0d0	
			dFv3(1, ix, iy)=0.0d0
			dGv3(1, ix, iy)=0.0d0
		end do
!RH壁
	ix=n1+n2
!!$omp parallel do private(k)
		do iy=m1,m1+m2
			u(ix,iy)=0.0d0
			v(ix,iy)=0.0d0
			LLx1(ix,iy)=(u(ix, iy)-c0(ix, iy))*(dpx2(ix,iy)-rho(ix,iy)*c0(ix, iy)*dux2(ix,iy))
			LLx2(ix,iy)=0.0d0
			LLx3(ix,iy)=0.0d0
			LLx4(ix,iy)=LLx1(ix,iy)
			Dx1(ix,iy)=(LLx2(ix,iy)+0.5d0*(LLx1(ix,iy)+LLx4(ix,iy)))/(c0(ix, iy)*c0(ix, iy))
		  	Fpm3(1,ix,iy)=-Dx1(ix,iy)
		  do k=2,4
		  	Fpm3(k, ix, iy)=0.0d0
		  	Gpm3(k, ix, iy)=0.0d0
			dFv3(k, ix, iy)=0.0d0
			dGv3(k, ix, iy)=0.0d0
		  end do
		  	Gpm3(1, ix, iy)=0.0d0
			dFv3(1, ix, iy)=0.0d0
			dGv3(1, ix, iy)=0.0d0
		end do
!上下方向壁
!LWR壁
	iy=m1
!!$omp parallel do private(k)
		do ix=n1,n1+n2
			u(ix,iy)=0.0d0
			v(ix,iy)=0.0d0
			LLy4(ix,iy)=(v(ix, iy)+c0(ix, iy))*(dpy2(ix,iy)+rho(ix,iy)*c0(ix, iy)*dvy2(ix,iy))
			LLy3(ix,iy)=0.0d0
			LLy2(ix,iy)=0.0d0
			LLy1(ix,iy)=LLy4(ix,iy)
			Dy1(ix,iy)=(LLy2(ix,iy)+0.5d0*(LLy1(ix,iy)+LLy4(ix,iy)))/(c0(ix, iy)*c0(ix, iy))
		  	Gpm3(1,ix,iy)=-Dy1(ix,iy)
		  do k=2,4
		  	Fpm3(k, ix, iy)=0.0d0
		  	Gpm3(k, ix, iy)=0.0d0
			dFv3(k, ix, iy)=0.0d0
			dGv3(k, ix, iy)=0.0d0
		  end do
		  	Fpm3(1, ix, iy)=0.0d0
			dFv3(1, ix, iy)=0.0d0
			dGv3(1, ix, iy)=0.0d0
		end do

!UPR壁
	iy=m1+m2
!!$omp parallel do private(k)
		do ix=n1,n1+n2
			u(ix,iy)=0.0d0
			v(ix,iy)=0.0d0
			LLy1(ix,iy)=(v(ix, iy)-c0(ix, iy))*(dpy2(ix,iy)-rho(ix,iy)*c0(ix, iy)*dvy2(ix,iy))
			LLy2(ix,iy)=0.0d0
			LLy3(ix,iy)=0.0d0
			LLy4(ix,iy)=LLy1(ix,iy)
			Dy1(ix,iy)=(LLy2(ix,iy)+0.5d0*(LLy1(ix,iy)+LLy4(ix,iy)))/(c0(ix, iy)*c0(ix, iy))
		  	Gpm3(1,ix,iy)=-Dy1(ix,iy)
		  do k=2,4
		  	Fpm3(k, ix, iy)=0.0d0
		  	Gpm3(k, ix, iy)=0.0d0
			dFv3(k, ix, iy)=0.0d0
			dGv3(k, ix, iy)=0.0d0
		  end do
		  	Fpm3(1, ix, iy)=0.0d0
			dFv3(1, ix, iy)=0.0d0
			dGv3(1, ix, iy)=0.0d0
		end do

!四隅の扱い
!UPR壁とRH壁の交点
		ix=n1+n2
		iy=m1+m2
		Fpm3(1, ix, iy)=-Dx1(ix,iy)-Dy1(ix,iy)
		Gpm3(1, ix, iy)=0.0d0
!RH壁とLWR壁の交点
		ix=n1+n2
		iy=m1
		Fpm3(1, ix, iy)=-Dx1(ix,iy)-Dy1(ix,iy)
		Gpm3(1, ix, iy)=0.0d0
!LWR壁とLH壁の交点
		ix=n1
		iy=m1
		Fpm3(1, ix, iy)=-Dx1(ix,iy)-Dy1(ix,iy)
		Gpm3(1, ix, iy)=0.0d0
!LH壁とUPR壁の交点
		ix=n1
		iy=m1+m2
		Fpm3(1, ix, iy)=-Dx1(ix,iy)-Dy1(ix,iy)
		Gpm3(1, ix, iy)=0.0d0

!$omp parallel do private(ix,iy)   ! 5/18 
	do k=1,4
		do ix=0,nt
			do iy=0,mt
	   	   	sout(k, ix, iy) = Fpm3(k, ix, iy)+Gpm3(k, ix, iy)  +dFv3(k, ix, iy)+dGv3(k, ix, iy)
	   	 	end do
		end do
	end do


		
		
end subroutine RK_substep


end module sub_compactP



!********************module*******************************************
module para
  implicit  none
!  real(8),parameter ::	x0 = 0.d0,   &
!						xm = 2.d0,   &
!						y0 = 0.d0,   &
!						ym = 2.d0,   &
!  real(8),parameter ::	omg = 1.9d0, &
!						eps = 1.d-20
  real(8),parameter ::	omg = 1.4d0, &
						eps = 1.d-13

  integer,parameter :: itrmax = 1000000
  integer :: i,j,k,LL
  real(8) :: dx,dy,pi
end module para

module Poisson_equation
implicit none
contains

!*************solver SOR method***************************************
subroutine sor(f,u,np,mp,dxd1,ddxd1,dyd1,ddyd1, umax, vmax, Step)  !ポワソン方程式の右辺を2次元配列fで与えると、
							!解が２次元配列uとして返ってくる。配列サイズは0～mp。
							!何回イタレーションをやったかの数がStepで返ってくる。
use para
use globals
implicit none
integer :: n, np,mp, Step
real(8),dimension(0:np,0:mp) :: u,f,uk,d
real(8) :: max
real(8) :: dxd1(0:np),ddxd1(0:np), dyd1(0:mp), ddyd1(0:mp), Jac(0:np,0:mp), alpfa(0:np,0:mp), gamma(0:np,0:mp)
real(8) :: a(0:np,0:mp), b(0:np,0:mp), AA(0:np,0:mp), Amz(0:np,0:mp), Azm(0:np,0:mp), Apz(0:np,0:mp), Azp(0:np,0:mp)
real(8) :: umax, vmax


  dx = (umax-x0)/dble(np)
  dy = (vmax-y0)/dble(mp)

do j=0,mp
do i=0,np
	Jac(i,j)  =dxd1(i)*dyd1(j)
	alpfa(i,j)=dyd1(j)*dyd1(j)
	gamma(i,j)=dxd1(i)*dxd1(i)
	a(i,j)    =dyd1(j)*dyd1(j)*dyd1(j)*ddxd1(i)
	b(i,j)    =dxd1(i)*dxd1(i)*dxd1(i)*ddyd1(j)
	AA(i,j)   =2.0d0*(alpfa(i,j)+gamma(i,j))/(Jac(i,j)*Jac(i,j))
	Apz(i,j)  =alpfa(i,j)/(Jac(i,j)*Jac(i,j))-0.5d0*a(i,j)*dx/(Jac(i,j)*Jac(i,j)*Jac(i,j))
	Amz(i,j)  =alpfa(i,j)/(Jac(i,j)*Jac(i,j))+0.5d0*a(i,j)*dx/(Jac(i,j)*Jac(i,j)*Jac(i,j))
	Azp(i,j)  =gamma(i,j)/(Jac(i,j)*Jac(i,j))-0.5d0*b(i,j)*dx/(Jac(i,j)*Jac(i,j)*Jac(i,j))
	Azm(i,j)  =gamma(i,j)/(Jac(i,j)*Jac(i,j))+0.5d0*b(i,j)*dx/(Jac(i,j)*Jac(i,j)*Jac(i,j))
enddo
enddo

!bcon  境界条件
	do i = 0,np
		u(i,0) = 1.0d0/gam
		u(i,mp) = 1.0d0/gam
	end do

	do j = 0,mp
		u(0,j) = 1.0d0/gam
		u(np,j) = 1.0d0/gam
	end do

!icon  初期条件
	do j = 1,mp-1
		do i = 1,np-1
			u(i,j) = 1.0d0/gam
		end do
	end do


  do j = 0,mp
    do i = 0,np
      uk(i,j) = 1.0d0/gam
    end do
  end do
  
  n = 0
  
  do LL = 1,itrmax
     
     n = n+1
!     write(*,*) n
     do j = 1,mp-1
       do k = 1,np-1
         uk(k,j) = u(k,j)
       end do
     end do
     
     do j = 1, mp-1
       do i = 1,np-1
			u(i,j) = u(i,j)+omg/AA(i,j)*(Amz(i,j)*u(i-1,j)+Azm(i,j)*u(i,j-1)+Apz(i,j)*u(i+1,j)+Azp(i,j)*u(i,j+1) &
			-dx*dx*f(i,j)-AA(i,j)*u(i,j))

			d(i,j) = (u(i,j)-uk(i,j))*(u(i,j)-uk(i,j))
       end do
     end do
     
     max = 0.0
     
     do j = 1,mp-1
       do k = 1,np-1
         if(max < d(k,j)) then
           max = d(k,j)
         end if
       end do
     end do
        
     if(max < eps) exit


  end do
  Step=n
  return
  
end subroutine sor
!*********************************************************************

end module Poisson_equation

!ここまではサブルーチン

!以下、メインのプログラム

program Fluid_cos_2phi
use globals
use sub_compactP
!use DCS_Subroutine
use Poisson_equation
implicit none

!!!--パラメータの宣言---------------------------------------------------------------------------------------------------------
!	real(8)            :: dt=0.006d0		!=0.006                                !これにRをかけて時間幅
	real(8)            :: dt=0.004d0		!=0.006                                !これにRをかけて時間幅
	integer, parameter :: End_Time_Step = 48000					!最大時間ステップ
	integer, parameter :: End_Time_Step2 =18000					!最大時間ステップ
!	integer, parameter :: End_Time_Step =  100					!最大時間ステップ					←テスト用
!	integer, parameter :: End_Time_Step2 =  10					!最大時間ステップ					←テスト用

!Output用
	character(6) :: digit, digit2, digit3, digit4
	integer :: Timestep

!RK用
integer  j, k, i, xi, yi, jt
integer  ix, iy
integer  iset, iset2

!二重極先行計算用
integer  jt_pre, MaxjI_pre
real(8)  ,allocatable ::Pressure_source3_pre(:,:)
real(8) PsL

!Lighthill Tensor用
integer  jarea, LiTens_Num, LiTens_Zone


integer  LT_Zone_LH_i,LT_Zone_RH_i,LT_Zone_LWR_i,LT_Zone_UPR_i

integer, parameter  ::  n   =  635    ! xをN分割 n=n1+n2+n3
integer, parameter  ::  n1  =  145    ! 箱の左から左端
integer, parameter  ::  n2  =  160    ! 箱の両端（X方向）
integer, parameter  ::  n3L =  200   ! 箱の右から細分領域右端
integer, parameter  ::  n3R =  130   ! 細分領域右端から右端
integer, parameter  ::  n3  =  330    ! 箱の右から右端 n3=n3L+n3R

integer, parameter  ::  m   =  635    ! xをN分割 m=m1+m2+m3
integer, parameter  ::  m1  =  145    ! 箱の下から下端
integer, parameter  ::  m2  =  160    ! 箱の両端（Y方向）
integer, parameter  ::  m3D =  200    ! 箱の上から細分領域上限
integer, parameter  ::  m3U =  130    ! 細分領域上限から上端
integer, parameter  ::  m3  =  330    ! 箱の上から上端 m3=m3L+m3U


!integer, parameter  ::  n  =  405    ! xをN分割 n=n1+n2+n3
!integer, parameter  ::  n1 =  145    ! 箱の左から左端
!integer, parameter  ::  n2 =   80    ! 箱の両端（X方向）
!integer, parameter  ::  n3 =  180    ! 箱の右から右端
!integer, parameter  ::  m  =  405    ! xをN分割 m=m1+m2+m3
!integer, parameter  ::  m1 =  145    ! 箱の下から下端
!integer, parameter  ::  m2 =   80    ! 箱の両端（Y方向）
!integer, parameter  ::  m3 =  180    ! 箱の上から上端


!integer, parameter  ::  MaxShift  =  240     !箱の上辺と右辺の和 n2+m2

!integer, parameter  ::  MeasurP_iX =  325     ! 観測点のX座標(刻み整数で)
!integer, parameter  ::  MeasurP_iY =   51     ! 観測点のX座標(刻み整数で)
real(8), parameter   ::  MeasurDistanceR =5.51d0     ! 観測点までの距離

real(8), parameter  ::  LT_Zone_LH  =   13.9d0   ! Lighthillテンソル計算領域 箱の左端X座標
real(8), parameter  ::  LT_Zone_RH  =   14.8d0   ! Lighthillテンソル計算領域 箱の右端X座標
real(8), parameter  ::  LT_Zone_LWR =   13.9d0  ! Lighthillテンソル計算領域 箱の下端Y座標
real(8), parameter  ::  LT_Zone_UPR =   14.4d0 ! Lighthillテンソル計算領域 箱の上端Y座標


!箱の大きさ
real(8), parameter :: Lx = 0.1d0                     !箱の左右の大きさ
real(8), parameter :: Ly = 0.1d0                     !箱の上下の大きさ

!箱の外の大きさ
real(8), parameter :: DL  = 14.0d0                     !箱の左の余白
real(8), parameter :: DR1 =  0.1d0                     !箱の右の余白
real(8), parameter :: DR2 = 15.9d0                     !箱の右の余白
real(8), parameter :: DR  = 16.0d0                     !箱の右の余白
real(8), parameter :: DD  = 14.0d0                     !箱の下の余白
real(8), parameter :: DU1 =  0.1d0                     !箱の上の余白
real(8), parameter :: DU2 = 15.9d0                     !箱の上の余白
real(8), parameter :: DU  = 16.0d0                     !箱の上の余白


!箱の中央から渦中心までの距離
real(8), parameter ::  Box2Vortex_X = 0.60d0                     !水平方向
real(8), parameter ::  Box2Vortex_Y = 0.15d0                     !高さ方向

integer, parameter      ::Write_Start_No1 =      0     ! tの書き出しを開始するステップ
integer, parameter      ::Write_Stop_No1  =   7200     ! tの書き出しを停止するステップ
integer, parameter      ::No1_set         =    400     ! tの書き出しをset置き
!	integer, parameter      ::No1_set      =     5   ! tの書き出しをset置き								←テスト用

integer, parameter      ::Write_Start_No2 =   7200     ! tの書き出しを開始するステップ
integer, parameter      ::Write_Stop_No2  =  18000     ! tの書き出しを停止するステップ
integer, parameter      ::No2_set         =     50     ! tの書き出しをset置き

integer, parameter      ::Write_Start_No3 =  18000     ! tの書き出しを開始するステップ
integer, parameter      ::Write_Stop_No3  =  20000     ! tの書き出しを停止するステップ
integer, parameter      ::No3_set         =    400     ! tの書き出しをset置き

integer, parameter      ::Write_Start_No4 =1170000     ! tの書き出しを開始するステップ
integer, parameter      ::Write_Stop_No4  =2170000     ! tの書き出しを停止するステップ
integer, parameter      ::No4_set         =    200     ! tの書き出しをset置き

integer, parameter      ::   set2         = 50000    ! グラフ用の書き出しファイルをset2置きに出力
!integer, parameter      ::   set2         =   20000   ! グラフ用の書き出しファイルをset2置きに出力			←テスト用
!	integer, parameter      ::   set2         =       5   ! グラフ用の書き出しファイルをset2置きに出力			←テスト用

integer, parameter      ::  Source_No     =      4     ! グラフを書き出す音源の各辺の数
integer, parameter      ::  Source_Sum_No     =  2     ! グラフを書き出す音源のX方向分とY方向分の2つ
integer, parameter      ::  LiTens_Source_No  = 17     ! Lighthillテンソルを分割して計算する音源エリアの数
integer, parameter      :: Monitor_Angle_every =  15     ! グラフを書き出すモニター点を、この角度ずつ（単位は度）
!integer, parameter      :: Monitor_Angle_every =  180    ! グラフを書き出すモニター点を、この角度ずつ（単位は度）	←テスト用
!integer, parameter      ::  Graph_No	  =      4     ! グラフを書き出すモニター点の数

integer, parameter      :: Graph_data_No_per_File =  8     ! ファイル内に書き出すデータ群の数	←修正箇所!!

integer, parameter      ::SumInitial              =100     !時間軸に-∞から積分する際に始める数字		←修正箇所!!
integer, parameter      ::SumRange                =5000    !時間軸で-∞から積分する代わりにこの分だけ積分する数字	←修正箇所!!
integer, parameter      ::Represenative_monitor_angle = 11 !代表して観察する観測点の角度（Monitor_Angle_everyの倍数）

integer                 ::   set ,set3          ! tの書き出しをset置き
integer                 ::   SumStart          !時間積算の開始

integer  Graph_No, jI, MaxjI, Output_File_No
!integer  MeasurP_iX(0:Graph_No) , MeasurP_iY(0:Graph_No)
!real(8)  Xd_dummyM,Xd_dummyP,	MeasurP_ix_R(0:Graph_No),	 MeasurP_iy_R(0:Graph_No), &
!		 						MeasurP_rx(0:Graph_No),		 MeasurP_ry(0:Graph_No)

integer, allocatable :: MeasurP_iX(:) , MeasurP_iY(:)
real(8)                 Xd_dummyM,Xd_dummyP, Yd_dummyM,Yd_dummyP

!音源用
!real(8) ::    pd2(0:n, 0:m, 0:n )=0.0d0	!pd2(0:n, 0:m, 0:MaxShift )のMaxShiftはnとmの大きい方
!real(8)		MeasurP_rx, MeasurP_ry, Dist2MP, Dummy
real(8), allocatable :: MeasurP_ix_R(:), MeasurP_iy_R(:), MeasurP_rx(:), MeasurP_ry(:),&
		 				Pressure_source_X(:), Pressure_source_Y(:), Ax(:,:), Ay(:,:),t_dist(:,:) ,&
		 				 Axx(:,:),Axy(:,:),Ayy(:,:), t_dist_quadro(:,:,:)
Integer		STEP2MP, shift, MaxShift
Integer ::  MinSTEP2MP=10000000		!十分に大きな数字を入れておく
Integer :: Delay_Step_ZONE(0:n,0:m), Delay_Step_UPR(0:n), Delay_Step_RH(0:m), Delay_Step_LWR(0:n), Delay_Step_LH(0:m)
real(8)		Ox, Oy, dtR, DummySQRT
integer  ixshift, iyshift

!CS用
integer  Nx, ia,  iw

real(8) ::     LcNP(0:n, 0:n)=0.0d0       !配列宣言　LcNP(m)
real(8) ::     UcNP(0:n, 0:n)=0.0d0       !配列宣言　UcNP(m)


real(8) ::     LcN1(  0:  n1, 0:  1)=0.0d0       !配列宣言　LcNP(m)
real(8) ::     LcN2(  0:  n2, 0:  1)=0.0d0       !配列宣言　LcNP(m)
real(8) ::     LcN3(  0:  n3, 0:  1)=0.0d0       !配列宣言　LcNP(m)
real(8) ::     LcNT(  0:   n, 0:  1)=0.0d0       !配列宣言　LcNP(m)
real(8) ::     LcM1(  0:  m1, 0:  1)=0.0d0       !配列宣言　LcNP(m)
real(8) ::     LcM2(  0:  m2, 0:  1)=0.0d0       !配列宣言　LcNP(m)
real(8) ::     LcM3(  0:  m3, 0:  1)=0.0d0       !配列宣言　LcNP(m)
real(8) ::     LcMT(  0:   m, 0:  1)=0.0d0       !配列宣言　LcNP(m)
real(8) ::     UcN1(  0:  n1, 0:  1)=0.0d0       !配列宣言　UcNP(m)
real(8) ::     UcN2(  0:  n2, 0:  1)=0.0d0       !配列宣言　UcNP(m)
real(8) ::     UcN3(  0:  n3, 0:  1)=0.0d0       !配列宣言　UcNP(m)
real(8) ::     UcNT(  0:   n, 0:  1)=0.0d0       !配列宣言　UcNP(m)
real(8) ::     UcM1(  0:  m1, 0:  1)=0.0d0       !配列宣言　UcNP(m)
real(8) ::     UcM2(  0:  m2, 0:  1)=0.0d0       !配列宣言　UcNP(m)
real(8) ::     UcM3(  0:  m3, 0:  1)=0.0d0       !配列宣言　UcNP(m)
real(8) ::     UcMT(  0:   m, 0:  1)=0.0d0       !配列宣言　UcNP(m)

	real(8) ::     LcN1dp(  0:  n1, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcN2dp(  0:  n2, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcN3dp(  0:  n3, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcNTdp(  0:   n, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcM1dp(  0:  m1, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcM2dp(  0:  m2, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcM3dp(  0:  m3, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcMTdp(  0:   m, 0:  1)=0.0d0       !配列宣言　LcNP(m)

	real(8) ::     LcN1dm(  0:  n1, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcN2dm(  0:  n2, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcN3dm(  0:  n3, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcNTdm(  0:   n, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcM1dm(  0:  m1, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcM2dm(  0:  m2, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcM3dm(  0:  m3, 0:  1)=0.0d0       !配列宣言　LcNP(m)
	real(8) ::     LcMTdm(  0:   m, 0:  1)=0.0d0       !配列宣言　LcNP(m)

	real(8) ::     UcN1dp(  0:  n1, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcN2dp(  0:  n2, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcN3dp(  0:  n3, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcNTdp(  0:   n, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcM1dp(  0:  m1, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcM2dp(  0:  m2, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcM3dp(  0:  m3, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcMTdp(  0:   m, 0:  1)=0.0d0       !配列宣言　UcNP(m)

	real(8) ::     UcN1dm(  0:  n1, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcN2dm(  0:  n2, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcN3dm(  0:  n3, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcNTdm(  0:   n, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcM1dm(  0:  m1, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcM2dm(  0:  m2, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcM3dm(  0:  m3, 0:  1)=0.0d0       !配列宣言　UcNP(m)
	real(8) ::     UcMTdm(  0:   m, 0:  1)=0.0d0       !配列宣言　UcNP(m)

!ゼロ割チェック用
real(8) ::     c0(0:n, 0:m)=0.0d0, sqrt_zero=0.0d0
Integer RK_Step

!RK用
real(8) ::     rho(0:n, 0:m)=0.0d0
real(8) ::       u(0:n, 0:m)=0.0d0
real(8) ::       v(0:n, 0:m)=0.0d0
real(8) ::      r1(0:n, 0:m)=0.0d0
real(8) ::      r2(0:n, 0:m)=0.0d0
real(8) ::       e(0:n, 0:m)=0.0d0
real(8) ::       p(0:n, 0:m)=0.0d0
real(8) ::     pm1(0:n, 0:m)=0.0d0
real(8) ::     pd1(0:n, 0:m)=0.0d0
real(8) ::   pd1su(0:n, 0:m)=0.0d0
real(8) ::   pd1sv(0:n, 0:m)=0.0d0
real(8) ::  Stress(0:n, 0:m)=0.0d0
real(8) ::Stressm1(0:n, 0:m)=0.0d0
real(8) ::Stressd1(0:n, 0:m)=0.0d0


real(8) ::    rhom1(0:n, 0:m)=0.0d0
real(8) ::      um1(0:n, 0:m)=0.0d0
real(8) ::      vm1(0:n, 0:m)=0.0d0
real(8) ::    rhom2(0:n, 0:m)=0.0d0
real(8) ::      um2(0:n, 0:m)=0.0d0
real(8) ::      vm2(0:n, 0:m)=0.0d0

real(8) ::LiTens(   1:3,0:n,0:m)=0.0d0
real(8) ::LiTensm1( 1:3,0:n,0:m)=0.0d0
real(8) ::LiTensm2( 1:3,0:n,0:m)=0.0d0
real(8) ::LiTensm3( 1:3,0:n,0:m)=0.0d0
real(8) ::LiTensm4( 1:3,0:n,0:m)=0.0d0
real(8) ::LiTens_d2(1:3,0:n,0:m)=0.0d0

real(8) ::AxxLiTens_d2(0:n,0:m)=0.0d0
real(8) ::AxyLiTens_d2(0:n,0:m)=0.0d0
real(8) ::AyyLiTens_d2(0:n,0:m)=0.0d0

real(8) ,allocatable :: pd2(:,:,:)
real(8) ,allocatable :: Coeff_V(:,:,:),Coeff_UPR(:,:),Coeff_RH(:,:),Coeff_LWR(:,:),Coeff_LH(:,:)

!real(8) ::    LiTens1(1:3, 0:n, 0:m)=0.0d0
!real(8) ::    LiTens2(1:3, 0:n, 0:m)=0.0d0

real(8) :: Q(1:4, 0:n, 0:m)=0.0d0
real(8) :: F(1:4, 0:n, 0:m)=0.0d0
real(8) :: G(1:4, 0:n, 0:m)=0.0d0
real(8) ::    Q23(0:n, 0:m)=0.0d0
real(8) :: ff_F_p(0:n       )=0.0d0
real(8) :: ff_F_m(0:n       )=0.0d0
real(8) ::     Fp(0:n       )=0.0d0
real(8) ::     Fm(0:n       )=0.0d0
real(8) ::Fp3(1:4,0:n, 0:m)=0.0d0
real(8) ::Fm3(1:4,0:n, 0:m)=0.0d0
real(8) :: ff_G_p(       0:m)=0.0d0
real(8) :: ff_G_m(       0:m)=0.0d0
real(8) ::     Gp(       0:m)=0.0d0
real(8) ::     Gm(       0:m)=0.0d0
real(8) ::Gp3(1:4,0:n, 0:m)=0.0d0
real(8) ::Gm3(1:4,0:n, 0:m)=0.0d0
real(8) ::s_in(1:4,0:n, 0:m)=0.0d0
real(8) :: s1(1:4,0:n, 0:m)=0.0d0
real(8) :: s2(1:4,0:n, 0:m)=0.0d0
real(8) :: s3(1:4,0:n, 0:m)=0.0d0
real(8) :: s4(1:4,0:n, 0:m)=0.0d0


!real(8)  gam, lam
real(8)  lam
real(8)  sk

!Omega用
real(8) ::Omega(0:n, 0:m)=0.0d0
real(8) ::dvx2(0:n, 0:m)=0.0d0
real(8) ::duy2(0:n, 0:m)=0.0d0
real(8) ::vx1(0:n), uy1(0:m), dvx1(0:n), duy1(0:m)

!Div用
real(8) ::Div(0:n, 0:m)=0.0d0
real(8) ::dux2(0:n, 0:m)=0.0d0
real(8) ::dvy2(0:n, 0:m)=0.0d0
real(8) ::ux1(0:n), vy1(0:m), dux1(0:n), dvy1(0:m)



real(8) xRK,dx,t
real(8) yRK,dy

!RK用
!real(8) ::     x(0:n, 0:m)=0            !配列宣言　x(n) を0リセット
!real(8) ::     f(0:n, 0:m)=0            !配列宣言　f(n) を0リセット
real(8) xr
real(8) yr
real(8) Pi
real(8) ixReal, iyReal
real(8) :: QS(1:4, 0:n, 0:m)=0.0d0

!Grasho Vortex用
real(8) LR, LR2, SR1, SR2, u0, ur1,ur2, a0 !, rho01, rho02 
real(8) x1, x2, y1, y2

!Teylor Vortex用
real(8) Mv, Rv

!Oseen Vortex用
real(8) ::u1(0:n,0:m)=0.0d0,u1x1(0:n)=0.0d0, du1x1(0:n)=0.0d0, du1x2(0:n, 0:m)=0.0d0
!real(8) ::					u1y1(0:m)=0.0d0, du1y1(0:m)=0.0d0, du1y2(0:n, 0:m)=0.0d0
real(8) ::v1(0:n,0:m)=0.0d0,v1x1(0:n)=0.0d0, dv1x1(0:n)=0.0d0, dv1x2(0:n, 0:m)=0.0d0
!real(8) ::					v1y1(0:m)=0.0d0, dv1y1(0:m)=0.0d0, dv1y2(0:n, 0:m)=0.0d0
!real(8) ::uv1(0:n,0:m)=0.0d0
!
real(8) ::u2(0:n,0:m)=0.0d0,u2x1(0:n)=0.0d0, du2x1(0:n)=0.0d0, du2x2(0:n, 0:m)=0.0d0
!real(8) ::					u2y1(0:m)=0.0d0, du2y1(0:m)=0.0d0, du2y2(0:n, 0:m)=0.0d0
real(8) ::v2(0:n,0:m)=0.0d0,v2x1(0:n)=0.0d0, dv2x1(0:n)=0.0d0, dv2x2(0:n, 0:m)=0.0d0
!real(8) ::					v2y1(0:m)=0.0d0, dv2y1(0:m)=0.0d0, dv2y2(0:n, 0:m)=0.0d0
!real(8) ::uv2(0:n,0:m)=0.0d0


real(8) ::u3x1(0:n)=0.0d0, du3x1(0:n)=0.0d0, du3x2(0:n, 0:m)=0.0d0
real(8) ::u3y1(0:m)=0.0d0, du3y1(0:m)=0.0d0, du3y2(0:n, 0:m)=0.0d0
real(8) ::v3x1(0:n)=0.0d0, dv3x1(0:n)=0.0d0, dv3x2(0:n, 0:m)=0.0d0
real(8) ::v3y1(0:m)=0.0d0, dv3y1(0:m)=0.0d0, dv3y2(0:n, 0:m)=0.0d0


!real(8) ::rho1(0:n,0:m)=0.0d0, rho2(0:n,0:m)=0.0d0, pr1(0:n,0:m)=0.0d0, pr2(0:n,0:m)=0.0d0, xp(0:n,0:m)=0.0d0, yp(0:n,0:m)=0.0d0
real(8) :: AA1(0:n, 0:m)=0.0d0, RHS1(0:n,0:m)=0.0d0
real(8) :: AA2(0:n, 0:m)=0.0d0, RHS2(0:n,0:m)=0.0d0
real(8) :: AA3(0:n, 0:m)=0.0d0, RHS3(0:n,0:m)=0.0d0
real(8) ::uv(0:n,0:m)=0.0d0
real(8) Circ,rho0

!双子渦用
real(8) alpfa, beta
real(8) x, y

!粘性計算用
real(8) Pr
real(8) :: Temp(   0:n, 0:m)
real(8) :: myu(    0:n, 0:m)

!不等長格子用
real(8) Defdx, defdy, Defdxi, defdyt, Defdxi1, defdyt1, Defdxi3, defdyt3, xmax, xmax1, xmax3, ymax, ymax1, ymax3
real(8) xsi, yta,  Ximax, Ximax1,Ximax2,Ximax3, Ytmax, Ytmax1,Ytmax2, Ytmax3
real(8) dXmin, dYmin, dXmax3, dXmin3, dYmax3, dYmin3, kappa, RDefix, RDefiy
Integer(8) Defix, Defiy
real(8) :: Defx(0:n), Defy(0:m), xd(0:n), yd(0:m), dxd1(0:n), dyd1(0:m), ddxd1(0:n), ddyd1(0:m)
real(8) :: xd1(0:n1),xd3(0:n3),yd1(0:m1),yd3(0:m3),dxd11(0:n1), dxd13(0:n3),dyd11(0:m1), dyd13(0:m3)
!real(8) :: AA=(1.05d0)**103-1.05d0
real(8) :: aa0, Ratio  !aa0 ⊿x、Ratio 格子伸張の拡大比
real(8) :: Xia, Xib, Ytaa, Ytab, Xa, Xb, Ya, Yb

!スポンジ領域用
real(8) absorb_sigma
real(8) :: Q_standard(1:4, 0:n, 0:m)=0.0d0
integer :: mn=0, ixmn=0, iymn=0
integer :: Step1=0, Step2=0, Step3=0


real(8) :: sgm, wa
real(8) :: sig

!時間モニター用
!real(8) :: pressure_monitor1(1:Graph_No, 0:End_Time_Step)=0.0d0
real(8), allocatable :: pressure_monitor1(:,:)

real(8) :: Pressure_source1(1:Source_No, 0:End_Time_Step)=0.0d0
real(8) :: Pressure_source2(1:Source_Sum_No, 0:End_Time_Step)=0.0d0
Integer :: point_i
real(8) :: LiTens_source1(1:LiTens_Source_No,1:3,0:End_Time_Step)=0.0d0
real(8)  ,allocatable :: Pressure_source_LiTens_Every(:,:,:), Pressure_source_LiTens_Every_Pict(:,:)
!real(8) ,allocatable:: LiTens_source1(:,:,:)
real(8) :: LiTens_source2(1:3,0:End_Time_Step)=0.0d0
real(8) :: Stress_source1(1:Source_No, 0:End_Time_Step)=0.0d0
real(8) :: Stress_source2(1:Source_Sum_No, 0:End_Time_Step)=0.0d0
real(8)  ,allocatable :: CTR(:,:)
real(8), allocatable :: Pressure_source_LiTens(:,:), Pressure_source_LiTens_pre(:)
real(8)  ,allocatable :: Pressure_source_XY(:,:), Pressure_source3(:,:)
real(8)  ,allocatable :: Stress_source_XY(:,:), Stress_source3(:,:)

!四重極計算用
real(8) LT_Zone_rLH, LT_Zone_rRH, LT_Zone_rLWR, LT_Zone_rUPR
!エラー検知
real(8)		daP, daRho	!圧力の平方根内を正に

			!継続用スタート
!再スタート用
Integer :: ns=0
Integer :: Start_Time_Step=0
Integer  i1, i2
real(8)  xs,ys,zs,z1,z2,z3,z4,z5
!real(8)  z,z1,z2,z3,z4,z5
			!継続用エンド

!多点モニタ点用
real(8) angle_rad
Integer k_Graph_No, i_LiTens_Zone
Graph_No=360/Monitor_Angle_every
allocate(   MeasurP_iX(0:Graph_No), MeasurP_iY(0:Graph_No), MeasurP_ix_R(0:Graph_No), MeasurP_iy_R(0:Graph_No), &
			MeasurP_rx(0:Graph_No), MeasurP_ry(0:Graph_No), pressure_monitor1(0:Graph_No,0:End_Time_Step) , &
!
!			Pressure_source_LiTens_Every(LT_Zone_LH_i:LT_Zone_RH_i,LT_Zone_LWR_i:LT_Zone_UPR_i,0:End_Time_Step),&
!
			Pressure_source_LiTens(0:Graph_No-1,0:End_Time_Step),&
			Pressure_source_LiTens_pre(0:End_Time_Step),&
			Ax(0:LiTens_Source_No,0:Graph_No-1), Ay(0:LiTens_Source_No,0:Graph_No-1),&
!		 	Axx(LT_Zone_LH_i:LT_Zone_RH_i,LT_Zone_LWR_i:LT_Zone_UPR_i),Axy(LT_Zone_LH_i:LT_Zone_RH_i,LT_Zone_LWR_i:LT_Zone_UPR_i),&
!		 	Ayy(LT_Zone_LH_i:LT_Zone_RH_i,LT_Zone_LWR_i:LT_Zone_UPR_i),&
		 	t_dist(0:LiTens_Source_No,0:Graph_No-1),&
			CTR(0:LiTens_Source_No,1:2),&
			Pressure_source_XY(0:Graph_No,0:End_Time_Step),Pressure_source3(0:Graph_No,0:End_Time_Step),&
			Pressure_source3_pre(0:Graph_No,0:End_Time_Step),&
			Stress_source_XY(0:Graph_No,0:End_Time_Step),Stress_source3(0:Graph_No,0:End_Time_Step))


!四重極分先行計算
	do jt_pre=0,End_Time_Step
	Pressure_source_LiTens_pre(jt_pre)=0.0d0
	end do

!

Nx=n   !---------------------------------------------------------------

!格子伸張用

!Ratio=1.05d0
LR=0.1d0  !渦直径

dt=0.5d0*LR*dt
Ximax=dble(n)
Ximax1=dble(n1)
Ximax2=dble(n2)
Ximax3=dble(n3)

Ytmax=dble(m)
Ytmax1=dble(m1)
Ytmax2=dble(m2)
Ytmax3=dble(m3)



xmax=DL+Lx+DR1+DR2

ymax=DD+Ly+DU1+DU2



!刻み用
t=0.0d0

xRK = 0.0d0
yRK = 0.0d0


Pi=acos(-1.d0)
Defdx = xmax/dble(n)
Defdy = ymax/dble(m)
Defdxi=ximax/dble(n)
Defdyt=ytmax/dble(m)
Defdxi1=ximax/dble(n1)
Defdyt1=ytmax/dble(m1)
Defdxi3=ximax/dble(n3)
Defdyt3=ytmax/dble(m3)


!gam=1.4d0



Pr=0.7d0

!双子渦用
alpfa=0.0d0  !双子渦の渦中心距離の半分(X)
beta =0.5d0*LR     !双子渦の渦中心距離の半分(Y)

!x1=0.375d0  *xmax+0.002d0+alpfa
!y1=0.337d0  *ymax+0.002d0-beta
!x2=0.375d0  *xmax+0.002d0+alpfa
!y2=0.337d0  *ymax+0.002d0+beta

x1=DL+0.5d0*Lx+Box2Vortex_X+0.002d0+alpfa
y1=DD+0.5d0*Ly+Box2Vortex_Y+0.002d0-beta
x2=DL+0.5d0*Lx+Box2Vortex_X+0.002d0+alpfa
y2=DD+0.5d0*Ly+Box2Vortex_Y+0.002d0+beta
!LR=0.4d0
!LR2=LR/2.0d0
rho0=1.0d0
a0=dsqrt(gam)
u0=0.0d0

!p0=1.0d0

!Mv=0.6d0  !周速
!Mv=0.19d0  !7/21
!Mv=0.095d0  !渦進行速度Ma0.05
!Mv=0.1825d0  !渦進行速度Ma0.1
!Mv=0.365d0  !渦進行速度Ma0.2
Mv=0.547458d0  !渦進行速度Ma0.3
Rv=0.5d0*LR  !渦半径を渦間距離の半分と定義する
!Rv=alpfa !渦半径

Circ=2.0d0*Pi/0.7d0*Mv*Rv	!循環

absorb_sigma=0.0d0
!absorb_sigma=10.0d0/dble(n-mn)


!LU分解行列の作成
	call LUNPcomp(  n1,  LcN1,  UcN1)
	call LUNPcomp(  n2,  LcN2,  UcN2)
	call LUNPcomp(  n3,  LcN3,  UcN3)
	call LUNPcomp(   n,  LcNT,  UcNT)
	call LUNPcomp(  m1,  LcM1,  UcM1)
	call LUNPcomp(  m2,  LcM2,  UcM2)
	call LUNPcomp(  m3,  LcM3,  UcM3)
	call LUNPcomp(   m,  LcMT,  UcMT)

sig= 0.25
	call LU_DCS_NPcomp(  n1,  LcN1dp,  UcN1dp, sig)
	call LU_DCS_NPcomp(  n2,  LcN2dp,  UcN2dp, sig)
	call LU_DCS_NPcomp(  n3,  LcN3dp,  UcN3dp, sig)
	call LU_DCS_NPcomp(   n,  LcNTdp,  UcNTdp, sig)
	call LU_DCS_NPcomp(  m1,  LcM1dp,  UcM1dp, sig)
	call LU_DCS_NPcomp(  m2,  LcM2dp,  UcM2dp, sig)
	call LU_DCS_NPcomp(  m3,  LcM3dp,  UcM3dp, sig)
	call LU_DCS_NPcomp(   m,  LcMTdp,  UcMTdp, sig)

sig=-0.25
	call LU_DCS_NPcomp(  n1,  LcN1dm,  UcN1dm, sig)
	call LU_DCS_NPcomp(  n2,  LcN2dm,  UcN2dm, sig)
	call LU_DCS_NPcomp(  n3,  LcN3dm,  UcN3dm, sig)
	call LU_DCS_NPcomp(   n,  LcNTdm,  UcNTdm, sig)
	call LU_DCS_NPcomp(  m1,  LcM1dm,  UcM1dm, sig)
	call LU_DCS_NPcomp(  m2,  LcM2dm,  UcM2dm, sig)
	call LU_DCS_NPcomp(  m3,  LcM3dm,  UcM3dm, sig)
	call LU_DCS_NPcomp(   m,  LcMTdm,  UcMTdm, sig)




! dx/dξを計算しておく

!	ゾーン②
			do Defix=n1+1,n1+n2
!				RDefix=dble(Defix)
				xd(Defix)=DL+Lx/Ximax2*(dble(Defix)-Ximax1)
			enddo

!	ゾーン①
		dXmin=Lx/Ximax2
		call Newton_min(dXmin, DL, n1, kappa)
			do Defix=0,n1
!				RDefix=dble(Defix)
				xd(Defix)=DL-DL*(exp(kappa*(Ximax1-dble(Defix))/Ximax1)-1.0d0)/(exp(kappa)-1.0d0)
			enddo		

!	ゾーン③左（右へ減少）
		dXmax3=xd(n1+n2)-xd(n1+n2-1)
		call Newton_max(dXmax3, DR1, n3L, kappa)
			do Defix=n1+n2+1,n1+n2+n3L
!				RDefix=dble(Defix)
				xd(Defix)=(DL+Lx)+DR1-DR1*(exp(kappa*(dble(n3L)-(dble(Defix)-dble(n1+n2)))/(dble(n3L)))-1.0d0)/(exp(kappa)-1.0d0)
			enddo

!!	ゾーン③左（右へ増加）
!		dXmin3=xd(n1+n2)-xd(n1+n2-1)
!		call Newton_min(dXmin3, DR1, n3L, kappa)
!			do Defix=n1+n2+1,n1+n2+n3L
!!				RDefix=dble(Defix)
!				xd(Defix)=DL+Lx+DR1*(exp(kappa*(dble(Defix)-(dble(n1+n2)))/(dble(n3L)))-1.0d0)/(exp(kappa)-1.0d0)
!			enddo


!	ゾーン③右
!		dXmin=Lx/Ximax2
		dXmin=xd(n1+n2+n3L)-xd(n1+n2+n3L-1)
		call Newton_min(dXmin, DR2, n3R, kappa)
			do Defix=n1+n2+n3L+1,n
!				RDefix=dble(Defix)
				xd(Defix)=DL+Lx+DR1+DR2*(exp(kappa*(dble(Defix)-(dble(n1+n2+n3L)))/(dble(n3R)))-1.0d0)/(exp(kappa)-1.0d0)
			enddo

		call CompactNP(LcNT, UcNT,xd, dxd1,ximax,n)

		do Defix=0,n1
			 xd1(Defix)= xd(Defix)
			dxd11(Defix)=dxd1(Defix)
		enddo

		do Defix=0,n3
			i=Defix+n1+n2
			 xd3(Defix)= xd(i)
			dxd13(Defix)=dxd1(i)
		enddo



! d2x/dξ2を計算しておく
		call CompactNP(LcNT, UcNT,dxd1, ddxd1,ximax,n)


! dy/dηを計算しておく

!	ゾーン②
			do Defiy=m1+1,m1+m2
!				RDefiy=dble(Defiy)
				yd(Defiy)=DD+Ly/Ytmax2*(dble(Defiy)-Ytmax1)
			enddo

!	ゾーン①
		dYmin=Ly/Ytmax2
		call Newton_min(dYmin, DD, m1, kappa)
			do Defiy=0,m1
!				RDefiy=dble(Defiy)
				yd(Defiy)=DD-DD*(exp(kappa*(Ytmax1-dble(Defiy))/Ytmax1)-1.0d0)/(exp(kappa)-1.0d0)
			enddo		

!	ゾーン③下（上へ減少）
		dYmax3=yd(m1+m2)-yd(m1+m2-1)
		call Newton_max(dYmax3, DU1, m3D, kappa)
			do Defiy=m1+m2+1,m1+m2+m3D
!				RDefiy=dble(Defiy)
				yd(Defiy)=(DD+Ly)+DU1-DU1*(exp(kappa*(dble(m3D)-(dble(Defiy)-dble(m1+m2)))/(dble(m3D)))-1.0d0)/(exp(kappa)-1.0d0)
			enddo

!!	ゾーン③下（上へ増加）
!		dYmin3=yd(m1+m2)-yd(m1+m2-1)
!		call Newton_min(dYmin3, DU1, m3D, kappa)
!			do Defiy=m1+m2+1,m1+m2+m3D
!!				RDefiy=dble(Defiy)
!				yd(Defiy)=DD+Ly+DU1*(exp(kappa*(dble(Defiy)-(dble(m1+m2)))/(dble(m3D)))-1.0d0)/(exp(kappa)-1.0d0)
!			enddo


!	ゾーン③上
!		dYmin=Ly/Ytmax2
		dYmin=yd(m1+m2+m3D)-yd(m1+m2+m3D-1)
		call Newton_min(dYmin, DU2, m3U, kappa)
			do Defiy=m1+m2+m3D+1,m
!				RDefiy=dble(Defiy)
				yd(Defiy)=DD+Ly+DU1+DU2*(exp(kappa*(dble(Defiy)-(dble(m1+m2+m3D)))/(dble(m3U)))-1.0d0)/(exp(kappa)-1.0d0)
			enddo

		call CompactNP(LcMT, UcMT,yd, dyd1,ytmax,m)

		do Defiy=0,m1
			 yd1(Defiy)= yd(Defiy)
			dyd11(Defiy)=dyd1(Defiy)
		enddo

		do Defiy=0,m3
			i=Defiy+m1+m2
			 yd3(Defiy)= yd(i)
			dyd13(Defiy)=dyd1(i)
		enddo



! d2y/dη2を計算しておく
		call CompactNP(LcNT, UcNT,dyd1, ddyd1,ytmax,m)


			open(70,file='X_result/X_result.txt')
		   write(70,* )  '  '
		do Defiy=0, m
			do Defix=0, n

				x=xd(Defix)
				y=yd(Defiy)

		        write(70,7000 )  Defix, Defiy, x, y
		   7000 format(2(I10.5), 2(E20.11) )
		     enddo
		write(70,* )     
		enddo
			close(70)


!測定点の決定
angle_rad=dble(Monitor_Angle_every/360.0d0)

do k=0, Graph_No-1		
		MeasurP_ix_R(k) =DL+0.5d0*Lx+MeasurDistanceR*dcos(dble(k)*2.0d0*Pi*angle_rad)
		MeasurP_iy_R(k) =DD+0.5d0*Ly+MeasurDistanceR*dsin(dble(k)*2.0d0*Pi*angle_rad)
end do

!測定点との差
!MeasurP_rx=xd(MeasurP_ix_R(k))
!MeasurP_ry=yd(MeasurP_iy_R(k))


!測定値(MeasurP_iX_R)近くのメッシュ点(MeasurP_rx)を算出
	write(*, '("Graph_No ",I2.0)') Graph_No
			open(80,file='Monitoring points/Monitoring point.txt')
		   write(80,* )  'Point_k ','MeasurP_ix_R ','ix_R_Converted_to ','MeasurP_ix ',&
		   						    'MeasurP_iy_R ','iy_R_Converted_to ','MeasurP_iy'
	do k=0, Graph_No-1

			do i=0, n
				!測定値近くのメッシュ点を算出
					Xd_dummyM=xd(i  )-MeasurP_ix_R(k)
					Xd_dummyP=xd(i+1)-MeasurP_ix_R(k)
				if (Xd_dummyP >= 0) then
					if (-Xd_dummyM < Xd_dummyP) then 
							MeasurP_ix(k)=i
					exit
					else 
							MeasurP_ix(k)=i+1
					exit
					end if
				else
				end if	
			end do
			MeasurP_rx(k)=xd(MeasurP_ix(k))
!	write(*,'("  Converted to     ",E12.5)') MeasurP_rx(k)
!	write(*,'("  MeasurP_ix(  ",I2.0,")= ",I3.0)') k,MeasurP_ix(k) 

!	write(*,'("  MeasurP_iy_R(",I2.0,")=",E12.5)') k,MeasurP_iy_R(k) 
			do j=0, m
				!測定値近くのメッシュ点を算出
					Yd_dummyM=yd(j  )-MeasurP_iy_R(k)
					Yd_dummyP=yd(j+1)-MeasurP_iy_R(k)
				if (Yd_dummyP >= 0) then
					if (-Yd_dummyM < Yd_dummyP) then 
							MeasurP_iy(k)=j
					exit
					else 
							MeasurP_iy(k)=j+1
					exit
					end if
				else
				end if	
			end do
			MeasurP_ry(k)=yd(MeasurP_iy(k))
!	write(*,'("  Converted to     ",E12.5)') MeasurP_ry(k)
!	write(*,*) "k=" ,k
!	write(*,'("  MeasurP_iy(  ",I2.0,")= ",I3.0)') k,MeasurP_iy(k) 

		        write(80,8000 )  k,MeasurP_ix_R(k),MeasurP_rx(k),MeasurP_ix(k), MeasurP_iy_R(k), MeasurP_ry(k), MeasurP_iy(k)
		   8000 format(I6.1, 2(E14.5), I6.1, 2(E14.5), I6.1)
!	write(80,*)
	end do
			close(80)

!LT計算領域の左側X座標(LT_Zone_LH)近くのメッシュ点(LT_Zone_rLH)を算出
			open(81,file='Monitoring points/LH_Zone_LH.txt')


			!LT Zone左端のメッシュ点を算出
			do i=0, n
					Xd_dummyM=xd(i  )-LT_Zone_LH
					Xd_dummyP=xd(i+1)-LT_Zone_LH
				if (Xd_dummyP >= 0) then
					if (-Xd_dummyM < Xd_dummyP) then 
							LT_Zone_LH_i=i
					exit
					else 
							LT_Zone_LH_i=i+1
					exit
					end if
				else
				end if	
			end do
			LT_Zone_rLH=xd(LT_Zone_LH_i)

			!LT Zone右端のメッシュ点を算出
			do i=0, n
					Xd_dummyM=xd(i  )-LT_Zone_RH
					Xd_dummyP=xd(i+1)-LT_Zone_RH
				if (Xd_dummyP >= 0) then
					if (-Xd_dummyM < Xd_dummyP) then 
							LT_Zone_RH_i=i
					exit
					else 
							LT_Zone_RH_i=i+1
					exit
					end if
				else
				end if	
			end do
			LT_Zone_rRH=xd(LT_Zone_RH_i)

			!LT Zone下端のメッシュ点を算出
			do j=0, m
				!測定値近くのメッシュ点を算出
					Yd_dummyM=yd(j  )-LT_Zone_LWR
					Yd_dummyP=yd(j+1)-LT_Zone_LWR
				if (Yd_dummyP >= 0) then
					if (-Yd_dummyM < Yd_dummyP) then 
							LT_Zone_LWR_i=j
					exit
					else 
							LT_Zone_LWR_i=j+1
					exit
					end if
				else
				end if	
			end do
			LT_Zone_rLWR=yd(LT_Zone_LWR_i)

			!LT Zone上端のメッシュ点を算出
			do j=0, m
				!測定値近くのメッシュ点を算出
					Yd_dummyM=yd(j  )-LT_Zone_UPR
					Yd_dummyP=yd(j+1)-LT_Zone_UPR
				if (Yd_dummyP >= 0) then
					if (-Yd_dummyM < Yd_dummyP) then 
							LT_Zone_UPR_i=j
					exit
					else 
							LT_Zone_UPR_i=j+1
					exit
					end if
				else
				end if	
			end do
			LT_Zone_rUPR=yd(LT_Zone_UPR_i)



!	write(*,'("  Converted to     ",E12.5)') MeasurP_ry(k)
!	write(*,'("  MeasurP_iy(  ",I2.0,")= ",I3.0)') k,MeasurP_iy(k) 

	write(*,'("LT_Zone_LH  ",E12.5,"  LH_Converted_to ",E12.5," Mesh Point at LH ",I6.1)')LT_Zone_LH,LT_Zone_rLH,LT_Zone_LH_i
	write(*,'("LT_Zone_RH  ",E12.5,"  RH_Converted_to ",E12.5," Mesh Point at RH ",I6.1)')LT_Zone_RH,LT_Zone_rRH,LT_Zone_RH_i
	write(*,'("LT_Zone_LWR ",E12.5," LWR_Converted_to ",E12.5," Mesh Point at LWR",I6.1)')LT_Zone_LWR,LT_Zone_rLWR,LT_Zone_LWR_i
	write(*,'("LT_Zone_UPR ",E12.5," UPR_Converted_to ",E12.5," Mesh Point at UPR",I6.1)')LT_Zone_UPR,LT_Zone_rUPR,LT_Zone_UPR_i

!		        write(81,8100 )  k,MeasurP_ix_R(k),MeasurP_rx(k),MeasurP_ix(k), MeasurP_iy_R(k), MeasurP_ry(k), MeasurP_iy(k)
!		   8100 format(I6.1, 2(E14.5), I6.1, 2(E14.5), I6.1)
!	write(80,*)

			close(81)

!!allocate(  	Axx(LT_Zone_LH_i:LT_Zone_RH_i,LT_Zone_LWR_i:LT_Zone_UPR_i),Axy(LT_Zone_LH_i:LT_Zone_RH_i,LT_Zone_LWR_i:LT_Zone_UPR_i!!),&
!!		 	Ayy(LT_Zone_LH_i:LT_Zone_RH_i,LT_Zone_LWR_i:LT_Zone_UPR_i),&
!!		 	t_dist_quadro(LT_Zone_LH_i:LT_Zone_RH_i,LT_Zone_LWR_i:LT_Zone_UPR_i,0:Graph_No-1)  )

allocate(  	Axx(0:LT_Zone_RH_i,0:LT_Zone_UPR_i),Axy(0:LT_Zone_RH_i,0:LT_Zone_UPR_i),&
		 	Ayy(0:LT_Zone_RH_i,0:LT_Zone_UPR_i),&
		 	t_dist_quadro(0:LT_Zone_RH_i,0:LT_Zone_UPR_i,0:Graph_No-1),&
			Pressure_source_LiTens_Every(LT_Zone_LH_i:LT_Zone_RH_i,LT_Zone_LWR_i:LT_Zone_UPR_i,0:End_Time_Step),& 
			Pressure_source_LiTens_Every_Pict(LT_Zone_LH_i:LT_Zone_RH_i,LT_Zone_LWR_i:LT_Zone_UPR_i) )
!箱の中心座標
	Ox=DL+0.5d0*Lx
	Oy=DD+0.5d0*Ly
	write(*,'("  BOX CTR(X,Y)=( ",E11.3,","E11.3,")")') Ox, Oy
	CTR(0,1)=Ox
	CTR(0,2)=Oy
	CTR(1,1)=13.625d0
	CTR(1,2)=14.325d0
	CTR(2,1)=13.875d0
	CTR(2,2)=14.325d0
	CTR(3,1)=14.200d0
	CTR(3,2)=14.325d0
	CTR(4,1)=14.375d0
	CTR(4,2)=14.325d0
	CTR(5,1)=13.625d0
	CTR(5,2)=14.150d0
	CTR(6,1)=13.875d0
	CTR(6,2)=14.150d0
	CTR(7,1)=14.200d0
	CTR(7,2)=14.150d0
	CTR(8,1)=14.375d0
	CTR(8,2)=14.150d0
	CTR(9,1)=13.750d0
	CTR(9,2)=13.850d0
	CTR(10,1)=14.003d0
	CTR(10,2)=13.850d0
	CTR(11,1)=14.050d0
	CTR(11,2)=14.325d0
	CTR(12,1)=14.050d0
	CTR(12,2)=14.150d0
	CTR(13,1)=14.200d0
	CTR(13,2)=14.050d0
	CTR(14,1)=14.375d0
	CTR(14,2)=14.050d0
	CTR(15,1)=13.625d0
	CTR(15,2)=14.050d0
	CTR(16,1)=13.875d0
	CTR(16,2)=14.050d0
	CTR(17,1)=14.050d0
	CTR(17,2)=13.850d0

k=Represenative_monitor_angle
do ix=LT_Zone_LH_i, LT_Zone_RH_i
do iy=LT_Zone_LWR_i,LT_Zone_UPR_i

	!測定値近くのメッシュ点から、正確な測定点の座標を算出
		MeasurP_rx(k)=xd(MeasurP_ix(k))	!最終的な測定点のX座標
		MeasurP_ry(k)=yd(MeasurP_iy(k))	!最終的な測定点のY座標

	Axx(ix,iy)=(MeasurP_rx(k)-xd(ix))*(MeasurP_rx(k)-xd(ix)) &
		/  (((MeasurP_rx(k)-xd(ix))*(MeasurP_rx(k)-xd(ix))+(MeasurP_ry(k)-yd(iy))*(MeasurP_ry(k)-yd(iy))) &
		*dsqrt(dsqrt((MeasurP_rx(k)-xd(ix))*(MeasurP_rx(k)-xd(ix))+(MeasurP_ry(k)-yd(iy))*(MeasurP_ry(k)-yd(iy))))&
		*((dsqrt(2.0d0))*(dsqrt(2.0d0))*(dsqrt(2.0d0)))*Pi)

	Axy(ix,iy)=(MeasurP_rx(k)-xd(ix))*(MeasurP_ry(k)-yd(iy)) &
		/  (((MeasurP_rx(k)-xd(ix))*(MeasurP_rx(k)-xd(ix))+(MeasurP_ry(k)-yd(iy))*(MeasurP_ry(k)-yd(iy))) &
		*dsqrt(dsqrt((MeasurP_rx(k)-xd(ix))*(MeasurP_rx(k)-xd(ix))+(MeasurP_ry(k)-yd(iy))*(MeasurP_ry(k)-yd(iy))))&
		*((dsqrt(2.0d0))*(dsqrt(2.0d0))*(dsqrt(2.0d0)))*Pi)

	Ayy(ix,iy)=(MeasurP_ry(k)-yd(iy))*(MeasurP_ry(k)-yd(iy)) &
		/     (((MeasurP_rx(k)-xd(ix))*(MeasurP_rx(k)-xd(ix))+(MeasurP_ry(k)-yd(iy))*(MeasurP_ry(k)-yd(iy))) &
		*dsqrt(dsqrt((MeasurP_rx(k)-xd(ix))*(MeasurP_rx(k)-xd(ix))+(MeasurP_ry(k)-yd(iy))*(MeasurP_ry(k)-yd(iy))))&
		*((dsqrt(2.0d0))*(dsqrt(2.0d0))*(dsqrt(2.0d0)))*Pi)

	t_dist_quadro(ix,iy,k)= dsqrt((MeasurP_rx(k)-xd(ix))*(MeasurP_rx(k)-xd(ix))+(MeasurP_ry(k)-yd(iy))*(MeasurP_ry(k)-yd(iy)))
end do
end do

do j=0, LiTens_Source_No

	Ay(j,k)=(MeasurP_ry(k)-CTR(j,2)) &
		/     (dsqrt((MeasurP_rx(k)-CTR(j,1))*(MeasurP_rx(k)-CTR(j,1))+(MeasurP_ry(k)-CTR(j,2))*(MeasurP_ry(k)-CTR(j,2))) &
		*dsqrt(dsqrt((MeasurP_rx(k)-CTR(j,1))*(MeasurP_rx(k)-CTR(j,1))+(MeasurP_ry(k)-CTR(j,2))*(MeasurP_ry(k)-CTR(j,2))))&
		*((dsqrt(2.0d0))*(dsqrt(2.0d0))*(dsqrt(2.0d0)))*Pi)

	Ax(j,k)=(MeasurP_rx(k)-CTR(j,1)) &
		/     (dsqrt((MeasurP_rx(k)-CTR(j,1))*(MeasurP_rx(k)-CTR(j,1))+(MeasurP_ry(k)-CTR(j,2))*(MeasurP_ry(k)-CTR(j,2))) &
		*dsqrt(dsqrt((MeasurP_rx(k)-CTR(j,1))*(MeasurP_rx(k)-CTR(j,1))+(MeasurP_ry(k)-CTR(j,2))*(MeasurP_ry(k)-CTR(j,2))))&
		*((dsqrt(2.0d0))*(dsqrt(2.0d0))*(dsqrt(2.0d0)))*Pi)

	t_dist(j,k)= dsqrt((MeasurP_rx(k)-CTR(j,1))*(MeasurP_rx(k)-CTR(j,1))+(MeasurP_ry(k)-CTR(j,2))*(MeasurP_ry(k)-CTR(j,2)))
	write(*,'("  t_dist(  ", I2.0 , "," , I2.0 , ")= " , E12.5)') j,k,t_dist(j,k) 

end do


	do Defiy=0, m
		do Defix=0, n
				xsi=Defdxi*Defix
				yta=Defdyt*Defiy

				x=xd(Defix)
				y=yd(Defiy)


				SR1=dsqrt((x-x1)*(x-x1)+(y-y1)*(y-y1))	!中心からの距離
				SR2=dsqrt((x-x2)*(x-x2)+(y-y2)*(y-y2))

				ur1=+Circ/(2.0d0*Pi*SR1)*(1.0d0-dexp((-5.0d0/4.0d0)*SR1*SR1/(Rv*Rv)))	!周速
				u1(Defix, Defiy)=-ur1*(y-y1)/SR1
				v1(Defix, Defiy)= ur1*(x-x1)/SR1
!				uv1(Defix, Defiy)=dsqrt(u1(Defix, Defiy)*u1(Defix, Defiy)+v1(Defix, Defiy)*v1(Defix, Defiy))

				ur2=-Circ/(2.0d0*Pi*SR2)*(1.0d0-dexp((-5.0d0/4.0d0)*SR2*SR2/(Rv*Rv)))	!周速
				u2(Defix, Defiy)=-ur2*(y-y2)/SR2
				v2(Defix, Defiy)= ur2*(x-x2)/SR2
!				uv2(Defix, Defiy)=dsqrt(u2(Defix, Defiy)*u2(Defix, Defiy)+v2(Defix, Defiy)*v2(Defix, Defiy))
				u(Defix, Defiy) = u1(Defix, Defiy) + u2(Defix, Defiy)
				v(Defix, Defiy) = v1(Defix, Defiy) + v2(Defix, Defiy)				
		enddo
	enddo



!ポアソン方程式の右辺計算用
!!1の渦
!		!du1/dx
!			do Defiy=0,m
!				do Defix=0,n
!		 		   u1x1(Defix)=u1(Defix,Defiy)
!				end do
!!				call CompactNP(LcNP, UcNP, u1x1,du1x1,n)
!			call CompactDefNP(LcNT, UcNT, u1x1, dxd1, du1x1,Ximax, n)
!				do Defix=0,n
!					du1x2(Defix,Defiy)=du1x1(Defix)
!				end do
!			end do
!		!du1/dy
!			do Defix=0,n
!				do Defiy=0,m
!				   u1y1(Defiy)=u1(Defix,Defiy)
!				end do
!!				call CompactNP(LcNP, UcNP, u1y1,du1y1,m)
!			call CompactDefNP(LcMT, UcMT, u1y1, dyd1, du1y1,Ytmax, m)
!				do Defiy=0,m
!					du1y2(Defix,Defiy)=du1y1(Defiy)
!				end do
!			end do
!		!dv1/dx
!			do Defiy=0,m
!				do Defix=0,n
!		 		   v1x1(Defix)=v1(Defix,Defiy)
!				end do
!!				call CompactNP(LcNP, UcNP, v1x1,dv1x1,n)
!			call CompactDefNP(LcNT, UcNT, v1x1, dxd1, dv1x1,Ximax, n)
!				do Defix=0,n
!					dv1x2(Defix,Defiy)=dv1x1(Defix)
!				end do
!			end do
!		!dv1/dy
!			do Defix=0,n
!				do Defiy=0,m
!				   v1y1(Defiy)=v1(Defix,Defiy)
!				end do
!!				call CompactNP(LcNP, UcNP, v1y1,dv1y1,m)
!			call CompactDefNP(LcMT, UcMT, v1y1, dyd1, dv1y1,Ytmax, m)
!				do Defiy=0,m
!					dv1y2(Defix,Defiy)=dv1y1(Defiy)
!				end do
!			end do
!
!!2の渦
!		!du2/dx
!			do Defiy=0,m
!				do Defix=0,n
!		 		   u2x1(Defix)=u2(Defix,Defiy)
!				end do
!		!		call CompactNP(LcNP, UcNP, u2x1,du2x1,n)
!				call CompactDefNP(LcNT, UcNT, u2x1, dxd1, du2x1,Ximax, n)
!				do Defix=0,n
!				du2x2(Defix,Defiy)=du2x1(Defix)
!				end do
!			end do
!		!du2/dy
!			do Defix=0,n
!				do Defiy=0,m
!				   u2y1(Defiy)=u2(Defix,Defiy)
!				end do
!		!		call CompactNP(LcNP, UcNP, u2y1,du2y1,m)
!				call CompactDefNP(LcMT, UcMT, u2y1, dyd1, du2y1,Ytmax, m)
!				do Defiy=0,m
!				du2y2(Defix,Defiy)=du2y1(Defiy)
!				end do
!			end do
!		!dv2/dx
!			do Defiy=0,m
!				do Defix=0,n
!		 		   v2x1(Defix)=v2(Defix,Defiy)
!				end do
!		!		call CompactNP(LcNP, UcNP, v2x1,dv2x1,n)
!				call CompactDefNP(LcNT, UcNT, v2x1, dxd1, dv2x1,Ximax, n)
!				do Defix=0,n
!				dv2x2(Defix,Defiy)=dv2x1(Defix)
!				end do
!			end do
!		!dv2/dy
!			do Defix=0,n
!				do Defiy=0,m
!				   v2y1(Defiy)=v2(Defix,Defiy)
!				end do
!		!		call CompactNP(LcNP, UcNP, v2y1,dv2y1,m)
!				call CompactDefNP(LcMT, UcMT, v2y1, dyd1, dv2y1,Ytmax, m)
!				do Defiy=0,m
!				dv2y2(Defix,Defiy)=dv2y1(Defiy)
!				end do
!			end do

!3の渦
		!du3/dx
			do Defiy=0,m
				do Defix=0,n
		 		   u3x1(Defix)=u(Defix,Defiy)
				end do
			call CompactDefNP(LcNT, UcNT, u3x1, dxd1, du3x1,Ximax, n)
				do Defix=0,n
					du3x2(Defix,Defiy)=du3x1(Defix)
				end do
			end do
		!du3/dy
			do Defix=0,n
				do Defiy=0,m
				   u3y1(Defiy)=u(Defix,Defiy)
				end do
			call CompactDefNP(LcMT, UcMT, u3y1, dyd1, du3y1,Ytmax, m)
				do Defiy=0,m
					du3y2(Defix,Defiy)=du3y1(Defiy)
				end do
			end do
		!dv3/dx
			do Defiy=0,m
				do Defix=0,n
		 		   v3x1(Defix)=v(Defix,Defiy)
				end do
			call CompactDefNP(LcNT, UcNT, v3x1, dxd1, dv3x1,Ximax, n)
				do Defix=0,n
					dv3x2(Defix,Defiy)=dv3x1(Defix)
				end do
			end do
		!dv3/dy
			do Defix=0,n
				do Defiy=0,m
				   v3y1(Defiy)=v(Defix,Defiy)
				end do
			call CompactDefNP(LcMT, UcMT, v3y1, dyd1, dv3y1,Ytmax, m)
				do Defiy=0,m
					dv3y2(Defix,Defiy)=dv3y1(Defiy)
				end do
			end do


!ポアソン方程式の右辺計算
		do Defiy=0,m
			do Defix=0,n
!			RHS1(Defix,Defiy)=-(gam-1.0d0)/gam*(du1x2(Defix,Defiy)*du1x2(Defix,Defiy)	&
!								+2.0d0*du1y2(Defix,Defiy)*dv1x2(Defix,Defiy)+dv1y2(Defix,Defiy)*dv1y2(Defix,Defiy))
!			RHS2(Defix,Defiy)=-(gam-1.0d0)/gam*(du2x2(Defix,Defiy)*du2x2(Defix,Defiy)	&
!								+2.0d0*du2y2(Defix,Defiy)*dv2x2(Defix,Defiy)+dv2y2(Defix,Defiy)*dv2y2(Defix,Defiy))
			RHS3(Defix,Defiy)=-(gam-1.0d0)/gam*(du3x2(Defix,Defiy)*du3x2(Defix,Defiy)	&
								+2.0d0*du3y2(Defix,Defiy)*dv3x2(Defix,Defiy)+dv3y2(Defix,Defiy)*dv3y2(Defix,Defiy))
			end do
		end do




!ポアソン方程式の両辺計算
!	call sor(RHS1,AA1,n,dxd1,ddxd1,dyd1,ddyd1,Ximax,Step1)
!	call sor(RHS2,AA2,n,dxd1,ddxd1,dyd1,ddyd1,Ximax,Step2)
	call sor(RHS3,AA3,n,m,dxd1,ddxd1,dyd1,ddyd1,Ximax,Ytmax,Step3)

!		write(*,*)"Step1=",Step1
!		write(*,*)"Step2=",Step2
		write(*,*)"Step3=",Step3

	do yi=0, m
		do xi=0, n

!!		rho01=(gam*(rho0**(gam-1.0d0))*AA1(xi,yi)/(c0*c0))**(1.0d0/(gam-1.0d0))		!音速は？？？？
!		rho1(xi,yi)=(gam*(rho0**(gam-1.0d0))*AA1(xi,yi))**(1.0d0/(gam-1.0d0))
!		pr1(xi,yi)=AA1(xi,yi)*rho1(xi,yi)
!
!		rho2(xi,yi)=(gam*(rho0**(gam-1.0d0))*AA2(xi,yi))**(1.0d0/(gam-1.0d0))
!		pr2(xi,yi)=AA2(xi,yi)*rho2(xi,yi)		

		rho(xi,yi)=(gam*(rho0**(gam-1.0d0))*AA3(xi,yi))**(1.0d0/(gam-1.0d0))
		p(xi,yi)  =AA3(xi,yi)*rho(xi,yi)

		enddo
	enddo



do yi=0, m
do xi=0, n
		u(xi, yi)= u0	+u1(xi,yi)  +  u2(xi,yi)
		v(xi, yi)=   	+v1(xi,yi)  +  v2(xi,yi)
		e(xi, yi)=0.5d0*rho(xi, yi)*(u(xi, yi)*u(xi, yi)+v(xi, yi)*v(xi, yi))+p(xi, yi)/(gam-1.0d0) 

enddo
enddo



!rho,u,v,eをQ1～Q4へ
do j=0,m
	do i=0,n
	  Q(1,i,j)=rho(i,j)
	  Q(2,i,j)=Q(1,i,j)*u(i, j)
	  Q(3,i,j)=Q(1,i,j)*v(i, j)
	  Q(4,i,j)=e(i, j)

	  Q_standard(1,i,j)=rho(i,j)
	  Q_standard(2,i,j)=Q(1,i,j)*u(i, j)
	  Q_standard(3,i,j)=Q(1,i,j)*v(i, j)
	  Q_standard(4,i,j)=e(i, j)
	end do
end do




! t=0の初期値設定終了


xRK=0.0d0  !  ???
yRK=0.0d0  !  ???

iset=set  !初期値から書く
iset2=0


!iset=0  !初期値を書かない


do jt=0,End_Time_Step2            			!--------------------------------------------------------チェック

!do jt=19667,End_Time_Step            			!--------------------------------------------------------
!do jt=25000,End_Time_Step            			!--------------------------------------------------------

			!継続用エンド

		write(*,* )  jt, Q(1, 50,200), Q(1,200, 50)

!$omp parallel do private(iy)
	do ix=0, n
		do iy=0, m
			pm1(ix, iy) =p(  ix, iy)
			rho(ix, iy) =Q(1,ix, iy)
			  u(ix, iy) =Q(2,ix, iy)/Q(1,ix, iy)
			  v(ix, iy) =Q(3,ix, iy)/Q(1,ix, iy)
			  p(ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
		enddo
	enddo


!圧力モニター点0～Graph_No	!
do k=0, Graph_No-1
	ix=MeasurP_iX(k)
	iy=MeasurP_iY(k)
!			rho(ix, iy) =Q(1,ix, iy)
!			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
!			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
	Pressure_monitor1(k,jt)=p(ix,iy)
end do

	!Tempとmyuの作成
		!一般部
!$omp parallel do private(i, daP, daRho)
		do j=0,m
			do i=0,n
				Temp(i,j)=gam*(gam-1.0d0)/Q(1,i,j)*(Q(4,i,j)-0.5d0*Q(2,i,j)*Q(2,i,j)/Q(1,i,j) &
							-0.5d0*Q(3,i,j)*Q(3,i,j)/Q(1,i,j))
				myu(i,j)=(1.0d0+ST/T0)/(Temp(i,j)+ST/T0)*(Temp(i,j)**(3.0d0/2.0d0))

!		   	   c0(i,j)=dsqrt(gam*p(i, j)/rho(i, j))

				daP=dabs(p(i, j))
				daRho=dabs(rho(i, j))
!		   	   c0(i,j)=dsqrt(gam*daP/daRho)

		if( p(i,j)<0.0d0)then
		write(*,*) i,j,p(i,j),daP,"p for Temp Error hassei!"

		p(i,j)=daP		!要注意

		endif

		if( rho(i,j)<=0.0d0)then
		write(*,*) i,j, rho(i,j), daRho,"rho Error hassei!"

		rho(i,j)=daRho		!要注意
		QS(1,i,j)=daRho		!要注意

		endif

!		if( isnan(c0(i,j)))then
!		write(*,*)"c0 Error hassei!", RK_Step
!		write(*,*) i,j,daP, p(i, j), rho(i, j);stop
!		endif


			end do
		end do


	     !dv/dx, du/dxの作成
!各領域毎の微分設定
	!XLW部
!$omp parallel do private(ix,vx1,ux1,dvx1,dux1) ! 7/7 
		do iy=0,m1-1
			do ix=0,n
 			   vx1(ix)=   v(ix, iy)
 			   ux1(ix)=   u(ix, iy)
			end do

			call CompactDefNP(LcNT, UcNT,     vx1, dxd1,    dvx1,ximax,n)
			call CompactDefNP(LcNT, UcNT,     ux1, dxd1,    dux1,ximax,n)

			do ix=0,n
			   dvx2(ix,iy)=  dvx1(ix)
			   dux2(ix,iy)=  dux1(ix)
			end do
		end do

	!XML部
!$omp parallel do private(ix,vx1,ux1,dvx1,dux1) ! 7/7 
		do iy=m1,m1+m2
			do ix=0,n1
 			   vx1(ix)=   v(ix, iy)
 			   ux1(ix)=   u(ix, iy)
			end do

			call CompactDefNP(LcN1, UcN1,     vx1, dxd11,    dvx1,ximax1,n1)
			call CompactDefNP(LcN1, UcN1,     ux1, dxd11,    dux1,ximax1,n1)

			do ix=0,n1
			   dvx2(ix,iy)=  dvx1(ix)
			   dux2(ix,iy)=  dux1(ix)
			end do
		end do

	!XMR部
!$omp parallel do private(ix,ixshift,vx1,ux1,dvx1,dux1) ! 7/7 
		do iy=m1,m1+m2
			do ix=n1+n2,n
				ixshift=ix-(n1+n2)
 			   vx1(ixshift)=   v(ix, iy)
 			   ux1(ixshift)=   u(ix, iy)
			end do

			call CompactDefNP(LcN3, UcN3,     vx1, dxd13,    dvx1,ximax3,n3)
			call CompactDefNP(LcN3, UcN3,     ux1, dxd13,    dux1,ximax3,n3)

			do ix=n1+n2,n
				ixshift=ix-(n1+n2)
			   dvx2(ix,iy)=  dvx1(ixshift)
			   dux2(ix,iy)=  dux1(ixshift)
			end do
		end do

	!XUP部
!$omp parallel do private(ix,vx1,ux1,dvx1,dux1) ! 7/7 
		do iy=m1+m2+1,m
			do ix=0,n
 			   vx1(ix)=   v(ix, iy)
 			   ux1(ix)=   u(ix, iy)
			end do

			call CompactDefNP(LcNT, UcNT,     vx1, dxd1,    dvx1,ximax,n)
			call CompactDefNP(LcNT, UcNT,     ux1, dxd1,    dux1,ximax,n)

			do ix=0,n
			   dvx2(ix,iy)=  dvx1(ix)
			   dux2(ix,iy)=  dux1(ix)
			end do
		end do


 	     !dv/dy, du/dy, dT/dy, d(rho×v)/dy, dp/dy,d(ρ×u×v)/dy,d(ρ×v×v)/dy,d((e+p)×v)/dyの作成

	!YLH部
!$omp parallel do private(iy,vy1,uy1,dvy1,duy1)   ! 7/7
		do ix=0,n1-1																				!	0から？
			do iy=0,m
			   vy1(iy)=   v(ix, iy)
			   uy1(iy)=   u(ix, iy)
			end do

			call CompactDefNP(LcMT, UcMT,     vy1, dyd1,    dvy1,ytmax,m)
			call CompactDefNP(LcMT, UcMT,     uy1, dyd1,    duy1,ytmax,m)

			do iy=0,m
			    dvy2(ix,iy)=  dvy1(iy)
			    duy2(ix,iy)=  duy1(iy)
			end do
		end do

	!YML部
!$omp parallel do private(iy,vy1,uy1,dvy1,duy1)   ! 7/7
		do ix=n1,n1+n2
			do iy=0,m1
			   vy1(iy)=   v(ix, iy)
			   uy1(iy)=   u(ix, iy)
			end do

			call CompactDefNP(LcM1, UcM1,     vy1, dyd11,    dvy1,ytmax1,m1)
			call CompactDefNP(LcM1, UcM1,     uy1, dyd11,    duy1,ytmax1,m1)

			do iy=0,m1
			    dvy2(ix,iy)=  dvy1(iy)
			    duy2(ix,iy)=  duy1(iy)
			end do
		end do

	!YMU部
!$omp parallel do private(iy,iyshift,vy1,uy1,dvy1,duy1)   ! 7/7
		do ix=n1,n1+n2
			do iy=m1+m2,m
				iyshift=iy-(m1+m2)
			   vy1(iyshift)=   v(ix, iy)
			   uy1(iyshift)=   u(ix, iy)
			end do

			call CompactDefNP(LcM3, UcM3,     vy1, dyd13,    dvy1,ytmax3,m3)
			call CompactDefNP(LcM3, UcM3,     uy1, dyd13,    duy1,ytmax3,m3)

			do iy=m1+m2,m
				iyshift=iy-(m1+m2)
			    dvy2(ix,iy)=  dvy1(iyshift)
			    duy2(ix,iy)=  duy1(iyshift)
			end do
		end do

	!YRH部
!$omp parallel do private(iy,vy1,uy1,dvy1,duy1)   ! 7/7
		do ix=n1+n2+1,n
			do iy=0,m
			   vy1(iy)=   v(ix, iy)
			   uy1(iy)=   u(ix, iy)
			end do

			call CompactDefNP(LcMT, UcMT,     vy1, dyd1,    dvy1,ytmax,m)
			call CompactDefNP(LcMT, UcMT,     uy1, dyd1,    duy1,ytmax,m)

			do iy=0,m
			    dvy2(ix,iy)=  dvy1(iy)
			    duy2(ix,iy)=  duy1(iy)
			end do
		end do


!Source検討モニター点1
!境界上のpd1を積算
!	!UPR壁
		iy=m1+m2
!$omp parallel do private(shift) 
			do ix=n1+1, n1+n2-1
			Stressm1(ix, iy)=Stress(ix, iy)
!			Stress(ix, iy)=-2.0d0*myu(ix, iy)*dvy2(ix,iy)
			Stress(ix, iy)=-4.0d0/3.0d0*myu(ix, iy)*dvy2(ix,iy)
!			pm1(ix, iy) =p(  ix, iy)
			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
			uv( ix, iy)=dsqrt(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy))
			pd1(ix, iy)=(p(  	ix, iy)-pm1(	 ix, iy))/dt
	   Stressd1(ix, iy)=(Stress(ix, iy)-Stressm1(ix, iy))/dt
!				do shift=MaxShift,1,-1
!					pd2(ix,iy,shift)=pd2(ix,iy,shift-1)
!				end do
!					pd2(ix,iy,0)=pd1(ix, iy)
			end do
!	!RH-UPRコーナー
		iy=m1+m2
			ix=n1+n2
!			pm1(ix, iy) =p(  ix, iy)
			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
			uv(ix, iy)=dsqrt(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy))
			pd1(ix, iy)=(p(  ix, iy)-pm1(ix, iy))/dt

!	!RH壁
		ix=n1+n2
!$omp parallel do private(shift) 
			do iy=m1+1, m1+m2-1
			Stressm1(ix, iy)=Stress(ix, iy)
!			Stress(ix, iy)=-2.0d0*myu(ix, iy)*dux2(ix,iy)
			Stress(ix, iy)=-4.0d0/3.0d0*myu(ix, iy)*dux2(ix,iy)
!			pm1(ix, iy) =p(  ix, iy)
			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
			uv(ix, iy)=dsqrt(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy))
			pd1(ix, iy)=(p(  ix, iy)-pm1(ix, iy))/dt
	   Stressd1(ix, iy)=(Stress(ix, iy)-Stressm1(ix, iy))/dt
!				do shift=MaxShift,1,-1
!					pd2(ix,iy,shift)=pd2(ix,iy,shift-1)
!				end do
!					pd2(ix,iy,0)=pd1(ix, iy)
			end do

!	!RH-LWRコーナー
		iy=m1
			ix=n1+n2
!			pm1(ix, iy) =p(  ix, iy)
			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
			uv(ix, iy)=dsqrt(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy))
			pd1(ix, iy)=(p(  ix, iy)-pm1(ix, iy))/dt
			
	!LWR壁
		iy=m1
!$omp parallel do private(shift) 
			do ix=n1+1, n1+n2-1
			Stressm1(ix, iy)=Stress(ix, iy)
!			Stress(ix, iy)=-2.0d0*myu(ix, iy)*dvy2(ix,iy)
			Stress(ix, iy)=-4.0d0/3.0d0*myu(ix, iy)*dvy2(ix,iy)
!write(*,*) "Stressm1", Stressm1(ix,iy)
!			pm1(ix, iy) =p(  ix, iy)
			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
			uv(ix, iy)=dsqrt(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy))
			pd1(ix, iy)=(p(  ix, iy)-pm1(ix, iy))/dt
	   Stressd1(ix, iy)=(Stress(ix, iy)-Stressm1(ix, iy))/dt
!write(*,*) "Stressd1", Stressd1(ix,iy)
!				do shift=MaxShift,1,-1
!					pd2(ix,iy,shift)=pd2(ix,iy,shift-1)
!				end do
!					pd2(ix,iy,0)=pd1(ix, iy)
			end do

!	!LWR-LHコーナー
		iy=m1
			ix=n1
!			pm1(ix, iy) =p(  ix, iy)
			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
			uv(ix, iy)=dsqrt(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy))
			pd1(ix, iy)=(p(  ix, iy)-pm1(ix, iy))/dt

	!LH壁
		ix=n1
!$omp parallel do private(shift) 
			do iy=m1+1, m1+m2-1
			Stressm1(ix, iy)=Stress(ix, iy)
!			Stress(ix, iy)=-2.0d0*myu(ix, iy)*dux2(ix,iy)
			Stress(ix, iy)=-4.0d0/3.0d0*myu(ix, iy)*dux2(ix,iy)
!			pm1(ix, iy) =p(  ix, iy)
			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
			uv(ix, iy)=dsqrt(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy))
			pd1(ix, iy)=(p(  ix, iy)-pm1(ix, iy))/dt
	   Stressd1(ix, iy)=(Stress(ix, iy)-Stressm1(ix, iy))/dt


!write(*,*) 'p(n1, ',iy,')',  p(ix, iy)

!				do shift=MaxShift,1,-1
!					pd2(ix,iy,shift)=pd2(ix,iy,shift-1)
!				end do
!					pd2(ix,iy,0)=pd1(ix, iy)
			end do

!	!LH-UPRコーナー
		iy=m1+m2
			ix=n1
!			pm1(ix, iy) =p(  ix, iy)
			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
			uv(ix, iy)=dsqrt(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy))
			pd1(ix, iy)=(p(  ix, iy)-pm1(ix, iy))/dt


!壁面ソースでの時系列の発音を個々に記憶
	!UPR壁
	 Defiy=m1+m2
		Defix=n1
			Pressure_source1(1,jt)=Pressure_source1(1,jt)+0.7071d0*pd1(Defix,Defiy) * 0.5d0*(xd(Defix+1)-xd(Defix  ))
		do Defix=n1+1, n1+n2-1
			Pressure_source1(1,jt)=Pressure_source1(1,jt)+		   pd1(Defix,Defiy) * 0.5d0*(xd(Defix+1)-xd(Defix-1))
		 	  Stress_source1(1,jt)=  Stress_source1(1,jt)+	  Stressd1(Defix,Defiy) * 0.5d0*(xd(Defix+1)-xd(Defix-1))
		end do
		Defix=n1+n2
			Pressure_source1(1,jt)=Pressure_source1(1,jt)+0.7071d0*pd1(Defix,Defiy) * 0.5d0*(xd(Defix  )-xd(Defix-1))
	!RH壁
		Defix=n1+n2
		Defiy=m1
			Pressure_source1(2,jt)=Pressure_source1(2,jt)+0.7071d0*pd1(Defix,Defiy) * 0.5d0*(yd(Defiy+1)-yd(Defiy  ))
		do Defiy=m1+1, m1+m2-1
			Pressure_source1(2,jt)=Pressure_source1(2,jt)+		   pd1(Defix,Defiy) * 0.5d0*(yd(Defiy+1)-yd(Defiy-1))
		 	  Stress_source1(2,jt)=  Stress_source1(2,jt)+	  Stressd1(Defix,Defiy) * 0.5d0*(yd(Defix+1)-yd(Defix-1))
		end do
		Defiy=m1+m2
			Pressure_source1(2,jt)=Pressure_source1(2,jt)+0.7071d0*pd1(Defix,Defiy) * 0.5d0*(yd(Defiy  )-yd(Defiy-1))
	!LWR壁
		Defiy=m1
		Defix=n1
			Pressure_source1(3,jt)=Pressure_source1(3,jt)+0.7071d0*pd1(Defix,Defiy) * 0.5d0*(xd(Defix+1)-xd(Defix  ))
		do Defix=n1+1, n1+n2-1
			Pressure_source1(3,jt)=Pressure_source1(3,jt)+		   pd1(Defix,Defiy) * 0.5d0*(xd(Defix+1)-xd(Defix-1))
		 	  Stress_source1(3,jt)=  Stress_source1(3,jt)+	  Stressd1(Defix,Defiy) * 0.5d0*(xd(Defix+1)-xd(Defix-1))
		end do
		Defix=n1+n2
			Pressure_source1(3,jt)=Pressure_source1(3,jt)+0.7071d0*pd1(Defix,Defiy) * 0.5d0*(xd(Defix  )-xd(Defix-1))
	!LH壁
		Defix=n1
		Defiy=m1
			Pressure_source1(4,jt)=Pressure_source1(4,jt)+0.7071d0*pd1(Defix,Defiy) * 0.5d0*(yd(Defiy+1)-yd(Defiy  ))
		do Defiy=m1+1, m1+m2-1
			Pressure_source1(4,jt)=Pressure_source1(4,jt)+		   pd1(Defix,Defiy) * 0.5d0*(yd(Defiy+1)-yd(Defiy-1))
		 	  Stress_source1(4,jt)=  Stress_source1(4,jt)+	  Stressd1(Defix,Defiy) * 0.5d0*(yd(Defix+1)-yd(Defix-1))
		end do
		Defiy=m1+m2
			Pressure_source1(4,jt)=Pressure_source1(4,jt)+0.7071d0*pd1(Defix,Defiy) * 0.5d0*(yd(Defiy  )-yd(Defiy-1))

		Pressure_source2(1,jt)=Pressure_source1(1,jt)-Pressure_source1(3,jt)	!上下方向力の合算
		Pressure_source2(2,jt)=Pressure_source1(2,jt)-Pressure_source1(4,jt)	!左右方向力の合算
		  Stress_source2(1,jt)=  Stress_source1(1,jt)-  Stress_source1(3,jt)	!上下方向力の合算
		  Stress_source2(2,jt)=  Stress_source1(2,jt)-  Stress_source1(4,jt)	!左右方向力の合算
!write(*,*) 'Pressure_source1(1,', Pressure_source1(1,jt), Pressure_source1(3,jt)

!Lighthill項

do iy=0,m
	do ix=0,n
			rho(ix, iy) =Q(1,ix, iy)
			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
	end do
end do



	!Tempの作成
		!一般部
!$omp parallel do private(i, daP, daRho)
		do j=0,m
			do i=0,n
				Temp(i,j)=gam*(gam-1.0d0)/Q(1,i,j)*(Q(4,i,j)-0.5d0*Q(2,i,j)*Q(2,i,j)/Q(1,i,j) &
							-0.5d0*Q(3,i,j)*Q(3,i,j)/Q(1,i,j))
				myu(i,j)=(1.0d0+ST/T0)/(Temp(i,j)+ST/T0)*(Temp(i,j)**(3.0d0/2.0d0))

!		   	   c0(i,j)=dsqrt(gam*p(i, j)/rho(i, j))

				daP=dabs(p(i, j))
				daRho=dabs(rho(i, j))
!		   	   c0(i,j)=dsqrt(gam*daP/daRho)

		if( p(i,j)<0.0d0)then
		write(*,*) i,j,p(i,j),daP,"p for Temp Error hassei!"

		p(i,j)=daP		!要注意

		endif

		if( rho(i,j)<=0.0d0)then
		write(*,*) i,j, rho(i,j), daRho,"rho Error hassei!"

		rho(i,j)=daRho		!要注意
		QS(1,i,j)=daRho		!要注意

		endif

!		if( isnan(c0(i,j)))then
!		write(*,*)"c0 Error hassei!", RK_Step
!		write(*,*) i,j,daP, p(i, j), rho(i, j);stop
!		endif


			end do
		end do


	     !dv/dx, du/dxの作成
!各領域毎の微分設定
	!XLW部
!$omp parallel do private(ix,vx1,ux1,dvx1,dux1) ! 7/7 
		do iy=0,m1-1
			do ix=0,n
 			   vx1(ix)=   v(ix, iy)
 			   ux1(ix)=   u(ix, iy)
			end do

			call CompactDefNP(LcNT, UcNT,     vx1, dxd1,    dvx1,ximax,n)
			call CompactDefNP(LcNT, UcNT,     ux1, dxd1,    dux1,ximax,n)

			do ix=0,n
			   dvx2(ix,iy)=  dvx1(ix)
			   dux2(ix,iy)=  dux1(ix)
			end do
		end do

	!XML部
!$omp parallel do private(ix,vx1,ux1,dvx1,dux1) ! 7/7 
		do iy=m1,m1+m2
			do ix=0,n1
 			   vx1(ix)=   v(ix, iy)
 			   ux1(ix)=   u(ix, iy)
			end do

			call CompactDefNP(LcN1, UcN1,     vx1, dxd11,    dvx1,ximax1,n1)
			call CompactDefNP(LcN1, UcN1,     ux1, dxd11,    dux1,ximax1,n1)

			do ix=0,n1
			   dvx2(ix,iy)=  dvx1(ix)
			   dux2(ix,iy)=  dux1(ix)
			end do
		end do

	!XMR部
!$omp parallel do private(ix,ixshift,vx1,ux1,dvx1,dux1) ! 7/7 
		do iy=m1,m1+m2
			do ix=n1+n2,n
				ixshift=ix-(n1+n2)
 			   vx1(ixshift)=   v(ix, iy)
 			   ux1(ixshift)=   u(ix, iy)
			end do

			call CompactDefNP(LcN3, UcN3,     vx1, dxd13,    dvx1,ximax3,n3)
			call CompactDefNP(LcN3, UcN3,     ux1, dxd13,    dux1,ximax3,n3)

			do ix=n1+n2,n
				ixshift=ix-(n1+n2)
			   dvx2(ix,iy)=  dvx1(ixshift)
			   dux2(ix,iy)=  dux1(ixshift)
			end do
		end do

	!XUP部
!$omp parallel do private(ix,vx1,ux1,dvx1,dux1) ! 7/7 
		do iy=m1+m2+1,m
			do ix=0,n
 			   vx1(ix)=   v(ix, iy)
 			   ux1(ix)=   u(ix, iy)
			end do

			call CompactDefNP(LcNT, UcNT,     vx1, dxd1,    dvx1,ximax,n)
			call CompactDefNP(LcNT, UcNT,     ux1, dxd1,    dux1,ximax,n)

			do ix=0,n
			   dvx2(ix,iy)=  dvx1(ix)
			   dux2(ix,iy)=  dux1(ix)
			end do
		end do


 	     !dv/dy, du/dy, dT/dy, d(rho×v)/dy, dp/dy,d(ρ×u×v)/dy,d(ρ×v×v)/dy,d((e+p)×v)/dyの作成

	!YLH部
!$omp parallel do private(iy,vy1,uy1,dvy1,duy1)   ! 7/7
		do ix=0,n1-1																				!	0から？
			do iy=0,m
			   vy1(iy)=   v(ix, iy)
			   uy1(iy)=   u(ix, iy)
			end do

			call CompactDefNP(LcMT, UcMT,     vy1, dyd1,    dvy1,ytmax,m)
			call CompactDefNP(LcMT, UcMT,     uy1, dyd1,    duy1,ytmax,m)

			do iy=0,m
			    dvy2(ix,iy)=  dvy1(iy)
			    duy2(ix,iy)=  duy1(iy)
			end do
		end do

	!YML部
!$omp parallel do private(iy,vy1,uy1,dvy1,duy1)   ! 7/7
		do ix=n1,n1+n2
			do iy=0,m1
			   vy1(iy)=   v(ix, iy)
			   uy1(iy)=   u(ix, iy)
			end do

			call CompactDefNP(LcM1, UcM1,     vy1, dyd11,    dvy1,ytmax1,m1)
			call CompactDefNP(LcM1, UcM1,     uy1, dyd11,    duy1,ytmax1,m1)

			do iy=0,m1
			    dvy2(ix,iy)=  dvy1(iy)
			    duy2(ix,iy)=  duy1(iy)
			end do
		end do

	!YMU部
!$omp parallel do private(iy,iyshift,vy1,uy1,dvy1,duy1)   ! 7/7
		do ix=n1,n1+n2
			do iy=m1+m2,m
				iyshift=iy-(m1+m2)
			   vy1(iyshift)=   v(ix, iy)
			   uy1(iyshift)=   u(ix, iy)
			end do

			call CompactDefNP(LcM3, UcM3,     vy1, dyd13,    dvy1,ytmax3,m3)
			call CompactDefNP(LcM3, UcM3,     uy1, dyd13,    duy1,ytmax3,m3)

			do iy=m1+m2,m
				iyshift=iy-(m1+m2)
			    dvy2(ix,iy)=  dvy1(iyshift)
			    duy2(ix,iy)=  duy1(iyshift)
			end do
		end do

	!YRH部
!$omp parallel do private(iy,vy1,uy1,dvy1,duy1)   ! 7/7
		do ix=n1+n2+1,n
			do iy=0,m
			   vy1(iy)=   v(ix, iy)
			   uy1(iy)=   u(ix, iy)
			end do

			call CompactDefNP(LcMT, UcMT,     vy1, dyd1,    dvy1,ytmax,m)
			call CompactDefNP(LcMT, UcMT,     uy1, dyd1,    duy1,ytmax,m)

			do iy=0,m
			    dvy2(ix,iy)=  dvy1(iy)
			    duy2(ix,iy)=  duy1(iy)
			end do
		end do

!allocate ( LiTens(   1:3,0:n,0:m), LiTensm1( 1:3,0:n,0:m), LiTensm2( 1:3,0:n,0:m), LiTens_d2(1:3,0:n,0:m) )

!!領域内のTの２階時間微分
!$omp parallel do private(ix) 
	do iy=LT_Zone_LWR_i, LT_Zone_UPR_i
		do ix=LT_Zone_LH_i, LT_Zone_RH_i
			LiTensm4(1,ix,iy) =LiTensm3(1,ix,iy)
			LiTensm4(2,ix,iy) =LiTensm3(2,ix,iy)
			LiTensm4(3,ix,iy) =LiTensm3(3,ix,iy)
			LiTensm3(1,ix,iy) =LiTensm2(1,ix,iy)
			LiTensm3(2,ix,iy) =LiTensm2(2,ix,iy)
			LiTensm3(3,ix,iy) =LiTensm2(3,ix,iy)
			LiTensm2(1,ix,iy) =LiTensm1(1,ix,iy)
			LiTensm2(2,ix,iy) =LiTensm1(2,ix,iy)
			LiTensm2(3,ix,iy) =LiTensm1(3,ix,iy)
			LiTensm1(1,ix,iy) =LiTens(1,ix,iy)
			LiTensm1(2,ix,iy) =LiTens(2,ix,iy)
			LiTensm1(3,ix,iy) =LiTens(3,ix,iy)
			LiTens(1,ix,iy)=rho(ix,iy)*u(ix,iy)*u(ix,iy)+p(ix,iy)-rho(ix,iy)-2.0d0*myu(ix,iy)/3.0d0/Re*(2*dux2(ix,iy)-dvy2(ix,iy))
			LiTens(2,ix,iy)=rho(ix,iy)*u(ix,iy)*v(ix,iy)                          -myu(ix,iy)/      Re*(  dvx2(ix,iy)+duy2(ix,iy))
			LiTens(3,ix,iy)=rho(ix,iy)*v(ix,iy)*v(ix,iy)+p(ix,iy)-rho(ix,iy)-2.0d0*myu(ix,iy)/3.0d0/Re*(2*dvy2(ix,iy)-dux2(ix,iy))
!			LiTens_d2(1,ix,iy)=(LiTens(1,ix,iy)-2.0d0*LiTensm1(1,ix,iy)+LiTensm2(1,ix,iy))/dt/dt
!			LiTens_d2(2,ix,iy)=(LiTens(2,ix,iy)-2.0d0*LiTensm1(2,ix,iy)+LiTensm2(2,ix,iy))/dt/dt
!			LiTens_d2(3,ix,iy)=(LiTens(3,ix,iy)-2.0d0*LiTensm1(3,ix,iy)+LiTensm2(3,ix,iy))/dt/dt

			LiTens_d2(1,ix,iy)=(-LiTensm4(1,ix,iy)+16.0d0*LiTensm3(1,ix,iy)-30.0d0*LiTensm2(1,ix,iy)+16.0d0*LiTensm1(1,ix,iy) &
								-LiTens(1,ix,iy))/dt/dt/12.0d0
			LiTens_d2(2,ix,iy)=(-LiTensm4(2,ix,iy)+16.0d0*LiTensm3(2,ix,iy)-30.0d0*LiTensm2(2,ix,iy)+16.0d0*LiTensm1(2,ix,iy) &
								-LiTens(2,ix,iy))/dt/dt/12.0d0
			LiTens_d2(3,ix,iy)=(-LiTensm4(3,ix,iy)+16.0d0*LiTensm3(3,ix,iy)-30.0d0*LiTensm2(3,ix,iy)+16.0d0*LiTensm1(3,ix,iy) &
								-LiTens(3,ix,iy))/dt/dt/12.0d0
		end do
	end do
!deallocate ( LiTens, LiTensm1, LiTensm2, LiTens_d2 )

!!角柱内のTの２階時間微分をゼロ化
do ix=n1,n1+n2
do iy=m1,m1+m2
do k=1,3
			LiTens_d2(k,ix,iy)=0.0d0
end do
end do
end do




do k=0,Graph_No-1
	Pressure_source_XY(k,jt) = Ax(0,k)*Pressure_source2(2,jt)+Ay(0,k)*Pressure_source2(1,jt)
	Stress_source_XY(k,jt)   = Ax(0,k)*  Stress_source2(2,jt)+Ay(0,k)*  Stress_source2(1,jt)
!write(*,*) 'Ax(0,', k, ')', Ax(0,k), Ay(0,k), Stress_source_XY(k,jt), Stress_source2(2,jt), Stress_source2(1,jt)

end do


k=Represenative_monitor_angle

 do ix=LT_Zone_LH_i, LT_Zone_RH_i
 do iy=LT_Zone_LWR_i,LT_Zone_UPR_i
	Pressure_source_LiTens_Every(ix,iy,jt) = &
								(		 Axx(ix,iy) * LiTens_d2(1,ix,iy) &
								+2.0d0 * Axy(ix,iy) * LiTens_d2(2,ix,iy) &
							 			+Ayy(ix,iy) * LiTens_d2(3,ix,iy) )*0.25d0*(xd(ix+1)-xd(ix-1))*(yd(iy+1)-yd(iy-1))

!write(*,*) 'Pressure_source_LiTens_Every',i_LiTens_Zone, k,Pressure_source_LiTens_Every(i_LiTens_Zone,k,jt)
 end do
 end do

!	!使わないQS(k,i,j)領域をノンゼロにする
!		   do iy=m1+1, m1+m2-1
!			do ix=n1+1, n1+n2-1
!				Pressure_source_LiTens_Every(ix,iy,jt) = 1000.0d0
!				AxxLiTens_d2(ix,iy)=1000.0d0
!				AxyLiTens_d2(ix,iy)=1000.0d0
!				AyyLiTens_d2(ix,iy)=1000.0d0
!			end do
!		   end do



do k=0, Graph_No-1

			Pressure_source_LiTens(k,jt)=0.0d0
			Pressure_source3(k,jt)=0.0d0
			Stress_source3(k,jt)=0.0d0

!!!二重極分
!!		MaxjI=int((dt*jt-t_dist(0,k))/dt)
!!		do jI=SumInitial, MaxjI
!!!		do jI=1, MaxjI
!!			Pressure_source3(k,jt)=Pressure_source3(k,jt) &
!!						+Pressure_source_XY(k,jI)* 1.0d0/(SQRT(dt*jt-t_dist(0,k)-dt*jI)) * dt
!!			Stress_source3(k,jt)=Stress_source3(k,jt) &
!!						  +Stress_source_XY(k,jI)* 1.0d0/(SQRT(dt*jt-t_dist(0,k)-dt*jI)) * dt
!!!write(*,*) 'Stress_source3(,', k, ')', Stress_source_XY(k,jt), Stress_source3(k,jt),Pressure_source3(k,jt)
!!
!!		end do

!注；ここのPressure_source_LiTens_Everyのdtを中心差分の分だけ２つずらす必要あり。


!!四重極分
!	do i_LiTens_Zone=1,LiTens_Source_No
!		MaxjI=int((dt*jt-t_dist(i_LiTens_Zone,k))/dt)
!
!		!write(*,*) 'MaxjI',jt,i_LiTens_Zone, k,MaxjI
!		do jI=SumInitial, MaxjI
!!		do jI=1, MaxjI
!			Pressure_source_LiTens(k,jt)=Pressure_source_LiTens(k,jt) &
!						+Pressure_source_LiTens_Every(i_LiTens_Zone,k,jI)* 1.0d0/(SQRT(dt*jt-t_dist(i_LiTens_Zone,k)-dt*jI)) * dt
!		end do
!
!!		write(*,*) 'Pressure_source_LiTens',jt,i_LiTens_Zone, k,Pressure_source_LiTens(k,jt), Pressure_monitor1(k,jt)
!
!	end do

end do


!四重極分先行計算
k=Represenative_monitor_angle

		MaxjI_pre=jt

	if (MaxjI_pre < SumInitial+SumRange) then	! jt-SumRange < SumInitialのこと
		SumStart=SumInitial
	else
		SumStart=MaxjI_pre-SumRange
	endif


	do ix=LT_Zone_LH_i, LT_Zone_RH_i
	do iy=LT_Zone_LWR_i,LT_Zone_UPR_i
		jt_pre=jt+int(t_dist_quadro(ix,iy,k)/dt)
!		MaxjI=int((dt*jt-t_dist(i_LiTens_Zone,k))/dt)
		PsL=Pressure_source_LiTens_pre(jt_pre)
!$omp parallel do reduction(+: PsL)
		do jI=SumStart, MaxjI_pre-1
			PsL=PsL &
					+Pressure_source_LiTens_Every(ix,iy,jI)* 1.0d0/(SQRT(dt*(jt-jI))) * dt
		end do
!$omp end parallel do
		Pressure_source_LiTens_pre(jt_pre)=PsL
	end do
	end do

!deallocate ( LiTens_source1 )

	if (jt == Write_Start_No1) then
		iset=0
		 set=0
	endif

	if (jt>Write_Start_No1 .and. jt<=Write_Stop_No1) then
		set=No1_set
	endif


	if (jt == Write_Start_No2) then
		iset=0
		 set=0
	endif

	if (jt>Write_Start_No2 .and. jt<=Write_Stop_No2) then
		set=No2_set
	endif


	if (jt == Write_Start_No3) then
		iset=0
		 set=0
	endif

	if (jt>Write_Start_No3 .and. jt<=Write_Stop_No3) then
		set=No3_set
	endif


	if (jt == Write_Start_No4) then
		iset=0
		 set=0
	endif

	if (jt>Write_Start_No4 .and. jt<=Write_Stop_No4) then
		set=No4_set
	endif


	if ((jt>=Write_Start_No1 .and. jt<=Write_Stop_No1)  .or.   &
		(jt>=Write_Start_No2 .and. jt<=Write_Stop_No2)  .or.   &
		(jt>=Write_Start_No3 .and. jt<=Write_Stop_No3)  .or.   &
		(jt>=Write_Start_No4 .and. jt<=Write_Stop_No4)) then
	if (iset==set) then

!			write(*,*) 'Check1!'

	do ix=0, n
		do iy=0, m
			rho(ix, iy) =Q(1,ix, iy)
			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*rho(ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
			uv(ix, iy)=dsqrt(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy))
		enddo
	enddo

!		        write(*,* ) "u,v", u(50, 50), v(50, 50)


!境界をpd1sに
!	!UPR壁
		iy=m1+m2
			do ix=n1+1, n1+n2-1
				pd1su(ix, iy)= 0.0d0
				pd1sv(ix, iy)= pd1(ix, iy)
			end do

			ix=n1+n2
			iy=m1+m2
				pd1su(ix, iy)= 0.7071d0*pd1(ix, iy)
				pd1sv(ix, iy)= 0.7071d0*pd1(ix, iy)

!	!RH壁
		ix=n1+n2
			do iy=m1+1, m1+m2-1
				pd1su(ix, iy)= pd1(ix, iy)
				pd1sv(ix, iy)= 0.0d0
			end do

			ix=n1+n2
			iy=m1
				pd1su(ix, iy)= 0.7071d0*pd1(ix, iy)
				pd1sv(ix, iy)=-0.7071d0*pd1(ix, iy)

	!LWR壁
		iy=m1
			do ix=n1+1, n1+n2-1
				pd1su(ix, iy)= 0.0d0
				pd1sv(ix, iy)=-pd1(ix, iy)
			end do

			ix=n1
			iy=m1
				pd1su(ix, iy)=-0.7071d0*pd1(ix, iy)
				pd1sv(ix, iy)=-0.7071d0*pd1(ix, iy)

	!LH壁
		ix=n1
			do iy=m1+1, m1+m2-1
				pd1su(ix, iy)=-pd1(ix, iy)
				pd1sv(ix, iy)= 0.0d0
			end do

			ix=n1
			iy=m1+m2
				pd1su(ix, iy)=-0.7071d0*pd1(ix, iy)
				pd1sv(ix, iy)= 0.7071d0*pd1(ix, iy)

!!ωの算出
!
!      !ωの漸化式
!
!	do Defiy=0,m
!		do Defix=0,n
! 		   vx1(Defix)=v(Defix,Defiy)
!		end do
!		call CompactDefNP(LcNT, UcNT, vx1, dxd1, dvx1,Ximax, n)
!		do Defix=0,n
!		dvx2(Defix,Defiy)=dvx1(Defix)
!		end do
!	end do
!
!	do Defix=0,n
!		do Defiy=0,m
!		   uy1(Defiy)=u(Defix,Defiy)
!		end do
!		call CompactDefNP(LcMT, UcMT, uy1, dyd1, duy1,Ytmax, m)
!		do Defiy=0,m
!		duy2(Defix,Defiy)=duy1(Defiy)
!		end do
!	end do

	do Defiy=0,m
		do Defix=0,n
		Omega(Defix,Defiy)=dvx2(Defix,Defiy)-duy2(Defix,Defiy)
		end do
	end do

	!使わないOmega領域をノンゼロにする

		   do iy=m1+1, m1+m2-1
			do ix=n1+1, n1+n2-1
				Omega(ix,iy)=100.0d0	! 5にしておく？ 7/7
			end do
		   end do


!

!divの算出

!      !divの漸化式
!
!	do iy=0,m
!		do ix=0,n
! 		   ux1(ix)=u(ix,iy)
!		end do
!
!!			sgm = 0d0
!!			call bibun_DCS(ux1,dux1,L3,U3,sgm,n )
!
!!		call CompactNP(LcNP, UcNP, ux1,dux1,n)
!		call CompactDefNP(LcNP, UcNP, ux1, dxd1, dux1,Ximax, n)
!		do ix=0,n
!		dux2(ix,iy)=dux1(ix)
!		end do
!	end do
!
!	do ix=0,n
!		do iy=0,m
!		   vy1(iy)=v(ix,iy)
!		end do
!
!!			sgm = 0d0
!!			call bibun_DCS(vy1,dvy1,L3,U3,sgm,m )
!
!!		call CompactNP(LcNP, UcNP, vy1,dvy1,m)
!		call CompactDefNP(LcNP, UcNP, vy1, dyd1, dvy1,Ytmax, m)
!		do iy=0,m
!		dvy2(ix,iy)=dvy1(iy)
!		end do
!	end do
!
!	do iy=0,m
!		do ix=0,n
!		Div(ix,iy)=dux2(ix,iy)+dvy2(ix,iy)
!		end do
!	end do



			write(digit, '(I6.6)') jt


			open(10,file='Omega_result/' //trim(digit)//' .txt')
		   write(10,* )  '  '
		do Defiy=0, m
			do Defix=0, n
						xsi=Defdxi*Defix
						yta=Defdyt*Defiy

				x=xd(Defix)
				y=yd(Defiy)

		        write(10,1000 )  x, y, Omega(Defix, Defiy)
		   1000 format(2(E20.11), E20.11)
			enddo
		write(10,* )  '  '
		enddo
			close(10)


			open(20,file='P_result/' //trim(digit)//' .txt')
		   write(20,* )  '  '
		do Defiy=0, m
			do Defix=0, n
						xsi=Defdxi*Defix
						yta=Defdyt*Defiy

				x=xd(Defix)
				y=yd(Defiy)

		        write(20,2000 )  x, y, p(Defix, Defiy)
		   2000 format(2(E20.11), E20.11)
			enddo
		write(20,* )  '  '
		enddo
			close(20)


!			open(21,file='Pd1_result/' //trim(digit)//' .txt')
!		   write(21,* )  '  '
!		do Defiy=0, m
!			do Defix=0, n
!
!				x=xd(Defix)
!				y=yd(Defiy)
!
!		        write(21,2100 )  x, y, pd1(Defix, Defiy)
!		   2100 format(2(E20.11), E20.11)
!			enddo
!		write(21,* )  '  '
!		enddo
!			close(21)


!			open(22,file='Pd1s_result/' //trim(digit)//' .txt')
!		   write(22,* )  '  '
!		do Defiy=0, m
!			do Defix=0, n
!
!				x=xd(Defix)
!				y=yd(Defiy)
!
!		        write(22,2200 )  x, y, 0.0d0, pd1su(Defix, Defiy), pd1sv(Defix, Defiy), 0.0d0, 0.0d0
!		   2200 format(2(E20.11), E12.3, 2(E20.11), 2(E12.3))
!			enddo
!		write(22,* )  '  '
!		enddo
!			close(22)

!			open(40,file='Div_result/' //trim(digit)//' .txt')
!		   write(40,* )  '  '
!		do yi=0,m
!		     do xi=0,n
!			x=Defdx*xi
!			y=Defdy*yi  
!		        write(40,4000 )  x, y, Div(xi, yi)
!		   4000 format(2(E12.3), E22.9)
!		     enddo
!		   write(40,* )     
!		enddo
!			close(40)
!
			
!			open(50,file='rho_result/' //trim(digit)//' .txt')
!		   write(50,* )  '  '
!		do Defiy=0, m
!			do Defix=0, n
!						xsi=Defdxi*Defix
!						yta=Defdyt*Defiy
!
!				x=xd(Defix)
!				y=yd(Defiy)
!
!		        write(50,5000 )  x, y, rho(Defix, Defiy)
!		   5000 format(2(E20.11), E20.11)
!		     enddo
!		   write(50,* )     
!		enddo
!			close(50)
!
!
!			open(60,file='UV3D_result/' //trim(digit)//' .txt')
!		   write(60,* )  '  '
!		do Defiy=0, m
!			do Defix=0, n
!						xsi=Defdxi*Defix
!						yta=Defdyt*Defiy
!
!				x=xd(Defix)
!				y=yd(Defiy)
!
!		        write(60,6000 )  x, y, 0.0d0, u(Defix, Defiy), v(Defix, Defiy), 0.0d0, uv(Defix, Defiy)
!		   6000 format(2(E20.11), E12.3, 2(E20.11), E12.3, E20.11)
!		     enddo
!		write(60,* )     
!		enddo
!			close(60)


		open(320,file='LiTens_Every/' //trim(digit)//' .txt')
	!使わないQS(k,i,j)領域をノンゼロにする
			do Defiy=LT_Zone_LWR_i,LT_Zone_UPR_i
				do Defix=LT_Zone_LH_i, LT_Zone_RH_i
					Pressure_source_LiTens_Every_Pict(Defix,Defiy)=Pressure_source_LiTens_Every(Defix,Defiy,jt)
				enddo
			enddo

		   do iy=m1+1, m1+m2-1
			do ix=n1+1, n1+n2-1
				Pressure_source_LiTens_Every_Pict(ix,iy) = 1000.0d0
			end do
		   end do


		   write(320,* )  '  '
		do Defiy=LT_Zone_LWR_i,LT_Zone_UPR_i
			do Defix=LT_Zone_LH_i, LT_Zone_RH_i
						xsi=Defdxi*Defix
						yta=Defdyt*Defiy

				x=xd(Defix)
				y=yd(Defiy)

		        write(320,32000 )  x, y, Pressure_source_LiTens_Every_Pict(Defix,Defiy)
		   32000	format(2(E20.11),E20.11)
			enddo
		write(320,* )  '  '
		enddo
			close(320)


		!エラー検出
		if( isnan(Q( 1, 5, 5)))then
!			open( 99, file = "Error.txt"); write(99,*)"Error"; write(*,*)"Error" ;stop
!			close(99)
		write(*,*)"Error" ;stop
		endif
	
			
	  iset=0
	endif

	iset=iset+1
	endif


	if (iset2==set2) then

			write(digit2, '(I6.6)') jt

		do point_i=0,Graph_No-1
			write(digit3, '(I3.3)') point_i
!		open(100,file='Pressure_monitor/point'//Graph//' //trim(digit2)//' .txt')
		open(100,file='Pressure_monitor/' //trim('P')//trim(adjustl(digit3))//'F'//trim(digit2)//' .txt')
!		open(100,file='Pressure_monitor/' //trim(digit2)//' .txt')

			do j=jt-set2+1, jt
				write(100,10000 )  point_i,j, Pressure_monitor1(point_i,j)
		10000	format(I2,I8, E20.11)
			enddo
			write(100,* )  '  '

		close(100)
			enddo

		
!		open(210,file='Pressure_source1/' //trim(digit2)//' .txt')
!			do point_i=1,Source_No
!			do j=jt-set2+1, jt
!				write(210,21000 )  point_i,j, Pressure_source1(point_i,j)
!		21000	format(I2,I8, E20.11)
!			enddo
!			write(210,* )  '  '
!			enddo
!		close(210)

		open(220,file='Pressure_source2/' //trim(digit2)//' .txt')
			do point_i=1,Source_Sum_No
			do j=jt-set2+1, jt
				write(220,22000 )  point_i,j, Pressure_source2(point_i,j)
		22000	format(I2,I8, E20.11)
			enddo
			write(220,* )  '  '
			enddo
		close(220)


!!		open(320,file='LiTens_source2/' //trim(digit2)//' .txt')
!!			do i=1,LiTens_Source_No
!!			do k=0,Graph_No-1
!!			do j=jt-set2+1, jt
!!				write(320,32000 ) i, k, j, Pressure_source_LiTens_Every(i,k,j) ! LiTens_source2(LiTens_Num,j)
!!		32000	format(I5,I5,I8, E20.11)
!!			enddo
!!			write(320,* )  '  '
!!			enddo
!!			write(320,* )  '  '
!!			enddo
!!		close(320)


		open(330,file='LiTens_source3/' //trim(digit2)//' .txt')

			do k=0,Graph_No-1
			do j=jt-set2+1, jt
				write(330,33000 )  k, j, Pressure_source_LiTens(k,j)
		33000	format(I5,I8, E20.11)
			enddo
			write(330,* )  '  '
			enddo
		close(330)


		open(420,file='Stress_source2/' //trim(digit2)//' .txt')
			do point_i=1,Source_Sum_No
			do j=jt-set2+1, jt
				write(420,42000 )  point_i,j, Stress_source2(point_i,j)
		42000	format(I2,I8, E20.11)
			enddo
			write(420,* )  '  '
			enddo
		close(420)



	  iset2=0
	endif

	iset2=iset2+1

	!S1作成
		sk=0.0d0
		RK_Step=1
			Call RK_substep(Q,&
			LcN1,  UcN1,  LcN2,  UcN2,  LcN3,  UcN3,  LcNT,  UcNT,&
			LcM1,  UcM1,  LcM2,  UcM2,  LcM3,  UcM3,  LcMT,  UcMT,&
			LcN1dp,UcN1dp,LcN2dp,UcN2dp,LcN3dp,UcN3dp,LcNTdp,UcNTdp,&
			LcM1dp,UcM1dp,LcM2dp,UcM2dp,LcM3dp,UcM3dp,LcMTdp,UcMTdp,&
			LcN1dm,UcN1dm,LcN2dm,UcN2dm,LcN3dm,UcN3dm,LcNTdm,UcNTdm,&
			LcM1dm,UcM1dm,LcM2dm,UcM2dm,LcM3dm,UcM3dm,LcMTdm,UcMTdm,&
			n1,n2,n3,n,m1,m2,m3,m,dt,sk,s_in,s1,dxd1,dxd11,dxd13,dyd1,dyd11,dyd13,c0,RK_Step,&
			xmax,ximax,ximax1,ximax3,ymax,ytmax,ytmax1,ytmax3)	

	!S2作成
		sk=0.5d0
		RK_Step=2
			Call RK_substep(Q,&
			LcN1,  UcN1,  LcN2,  UcN2,  LcN3,  UcN3,  LcNT,  UcNT,&
			LcM1,  UcM1,  LcM2,  UcM2,  LcM3,  UcM3,  LcMT,  UcMT,&
			LcN1dp,UcN1dp,LcN2dp,UcN2dp,LcN3dp,UcN3dp,LcNTdp,UcNTdp,&
			LcM1dp,UcM1dp,LcM2dp,UcM2dp,LcM3dp,UcM3dp,LcMTdp,UcMTdp,&
			LcN1dm,UcN1dm,LcN2dm,UcN2dm,LcN3dm,UcN3dm,LcNTdm,UcNTdm,&
			LcM1dm,UcM1dm,LcM2dm,UcM2dm,LcM3dm,UcM3dm,LcMTdm,UcMTdm,&
			n1,n2,n3,n,m1,m2,m3,m,dt,sk,s1,s2,dxd1,dxd11,dxd13,dyd1,dyd11,dyd13,c0,RK_Step,&
			xmax,ximax,ximax1,ximax3,ymax,ytmax,ytmax1,ytmax3)	

	!S3作成
		sk=0.5d0
		RK_Step=3
			Call RK_substep(Q,&
			LcN1,  UcN1,  LcN2,  UcN2,  LcN3,  UcN3,  LcNT,  UcNT,&
			LcM1,  UcM1,  LcM2,  UcM2,  LcM3,  UcM3,  LcMT,  UcMT,&
			LcN1dp,UcN1dp,LcN2dp,UcN2dp,LcN3dp,UcN3dp,LcNTdp,UcNTdp,&
			LcM1dp,UcM1dp,LcM2dp,UcM2dp,LcM3dp,UcM3dp,LcMTdp,UcMTdp,&
			LcN1dm,UcN1dm,LcN2dm,UcN2dm,LcN3dm,UcN3dm,LcNTdm,UcNTdm,&
			LcM1dm,UcM1dm,LcM2dm,UcM2dm,LcM3dm,UcM3dm,LcMTdm,UcMTdm,&
			n1,n2,n3,n,m1,m2,m3,m,dt,sk,s2,s3,dxd1,dxd11,dxd13,dyd1,dyd11,dyd13,c0,RK_Step,&
			xmax,ximax,ximax1,ximax3,ymax,ytmax,ytmax1,ytmax3)	

	!S4作成
		sk=1.0d0
		RK_Step=4
			Call RK_substep(Q,&
			LcN1,  UcN1,  LcN2,  UcN2,  LcN3,  UcN3,  LcNT,  UcNT,&
			LcM1,  UcM1,  LcM2,  UcM2,  LcM3,  UcM3,  LcMT,  UcMT,&
			LcN1dp,UcN1dp,LcN2dp,UcN2dp,LcN3dp,UcN3dp,LcNTdp,UcNTdp,&
			LcM1dp,UcM1dp,LcM2dp,UcM2dp,LcM3dp,UcM3dp,LcMTdp,UcMTdp,&
			LcN1dm,UcN1dm,LcN2dm,UcN2dm,LcN3dm,UcN3dm,LcNTdm,UcNTdm,&
			LcM1dm,UcM1dm,LcM2dm,UcM2dm,LcM3dm,UcM3dm,LcMTdm,UcMTdm,&
			n1,n2,n3,n,m1,m2,m3,m,dt,sk,s3,s4,dxd1,dxd11,dxd13,dyd1,dyd11,dyd13,c0,RK_Step,&
			xmax,ximax,ximax1,ximax3,ymax,ytmax,ytmax1,ytmax3)	

!dp/dt用に、全ての変数を時間発展させる前に各原始変数を計算して確保しておく

!$omp parallel do private(ix)   ! 02/20
		do iy=LT_Zone_LWR_i, LT_Zone_UPR_i
			do ix=LT_Zone_LH_i, LT_Zone_RH_i
				rho(ix, iy) =Q(1,ix, iy)
				u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
				v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
				p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*Q(1,ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
			end do
		end do

!		iy=m1+m2
!			do ix=n1, n1+n2
!			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
!			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*Q(1,ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
!			enddo
!		ix=n1+n2
!			do iy=m1, m1+m2
!			u(  ix, iy) =Q(2,ix, iy)/ Q(1,ix, iy)
!			v(  ix, iy) =Q(3,ix, iy)/ Q(1,ix, iy)
!			p(  ix, iy) =(gam-1.0d0)*(Q(4,ix, iy)-0.5d0*Q(1,ix, iy)*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy)))
!			enddo

!全ての変数を時間発展させる
	do k=1,4
		do ix=0,n
			do iy=0,m

!	if (ix<mn .or. (n-mn)<ix) then
!		ixmn=abs((mn-ix))
!	endif
!
!	if (iy<mn .or. (m-mn)<iy) then
!		iymn=abs((mn-iy))
!	endif
	
	
	   	   	Q(k, ix, iy) = Q(k, ix, iy) + dt*( s1(k, ix, iy)+2.0d0*s2(k, ix, iy)+2.0d0*s3(k, ix, iy)+s4(k, ix, iy) )/6.0d0 !&
!	   	   					-absorb_sigma*(ixmn+iymn)*(Q(k, ix, iy)-Q_standard(k, ix, iy))
	   	 	end do
		end do
	end do




!Isothermal No-Slip Wall
!左右方向壁
!MLW壁
	ix=n1
		do iy=m1,m1+m2
			Q(2, ix, iy)=0.0d0
			Q(3, ix, iy)=0.0d0
			Q(4, ix, iy)=(Q(1, ix, iy)*T1/gam)/(gam-1.0d0)		!T1=壁面温度
		end do
!MRW壁
	ix=n1+n2
		do iy=m1,m1+m2
			Q(2, ix, iy)=0.0d0
			Q(3, ix, iy)=0.0d0
			Q(4, ix, iy)=(Q(1, ix, iy)*T1/gam)/(gam-1.0d0)		!T1=壁面温度
		end do

!上下方向壁
!MLW壁
	iy=m1
		do ix=n1, n1+n2
			Q(2, ix, iy)=0.0d0
			Q(3, ix, iy)=0.0d0
			Q(4, ix, iy)=(Q(1, ix, iy)*T1/gam)/(gam-1.0d0)		!T1=壁面温度
		end do
!MUW壁
	iy=m1+m2
		do ix=n1, n1+n2
			Q(2, ix, iy)=0.0d0
			Q(3, ix, iy)=0.0d0
			Q(4, ix, iy)=(Q(1, ix, iy)*T1/gam)/(gam-1.0d0)		!T1=壁面温度
		end do

	!使わないQS(k,i,j)領域をノンゼロにする
		do k=1,4,3
		   do iy=m1+1, m1+m2-1
			do ix=n1+1, n1+n2-1
				Q(k,ix,iy)=5.0d0	! 5にしておく？ 7/7
			end do
		   end do
		end do


		do k=2,3
		   do iy=m1+1, m1+m2-1
			do ix=n1+1, n1+n2-1
!				Q(k,ix,iy)=0.1d0	! 5にしておく？ 7/7
				Q(k,ix,iy)=0.0d0	! 05/01
			end do
		   end do
		end do

t  =dt*(jt+1)
enddo      

!!Output_File_No=ceiling(dble(LiTens_Source_No/LiTens_Source_No_per_File))
!!!Output_File_No=ceiling(dble(Graph_No/Graph_data_No_per_File))
!!!write(*,*) 'Output_File_No', Output_File_No
!!set3=0
!!do jt=1,Output_File_No
!!			write(digit4, '(I2.2)') jt
!!		open(300,file='LiTens_source2/' //trim('LiTens_Every')//trim(adjustl(digit4))//' .txt')
!!
!!			do i=set3,set3+LiTens_Source_No_per_File-1
!!			do k=0,Graph_No-1
!!			do j=0,End_Time_Step
!!				write(300,30000 )   i, k, j, Pressure_source_LiTens_Every(i,k,j)
!!		30000	format(I5,I5,I8, E20.11)
!!			enddo
!!			write(300,* )  '  '
!!			enddo
!!			write(300,* )  '  '
!!			enddo
!!set3=set3+LiTens_Source_No_per_File
!!		close(300)
!!
!!enddo


Output_File_No=ceiling(dble(Graph_No/Graph_data_No_per_File))
write(*,*) 'Output_File_No', Output_File_No
set3=0
do jt=1,Output_File_No
			write(digit4, '(I2.2)') jt
		open(300,file='LiTens_source3/' //trim('Pressure_source_LiTens_')//trim(adjustl(digit4))//' .txt')

			do point_i=set3,set3+Graph_data_No_per_File-1
			do j=0,End_Time_Step
			write(300,30000 )   point_i, j, Pressure_monitor1(point_i,j),&
								 Pressure_source3(point_i,j),Pressure_source3_pre(point_i,j), Pressure_source_LiTens_pre(j)
		30000	format(I5,I8, 3(E23.11),E23.11e4)
			enddo
			write(300,* )  '  '
			enddo
set3=set3+Graph_data_No_per_File
		close(300)

enddo


end program Fluid_cos_2phi