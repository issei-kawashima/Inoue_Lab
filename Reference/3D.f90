!吉野さんの乱流境界層の3次元コード
!第二不変量Qや、固有値・固有関数の使い方などが書いてあるから参考にできるはず
!あとは微分のsubroutineが細かくなっているなど(もしかしたら関数かも)計算高速化にも使えるかも？
!シューティング法など一部別手法で計算をしているので注意
module DCS5
!!!-----DCS and CCS parameter definition------------
!!!-----5th
	double precision :: alpc = 1.d0/3.d0, alpd= 1.d0
	double precision :: ac  = 14.d0/9.d0, ad  = 4.d0/9.d0, bc = 1.d0/9.d0, bd = 2.d0/9.d0
!!!-----3rd One-sided DCS (4th One-sided CCS)
	double precision :: alpb = 3.d0,  ab = -17.d0/6.d0, bb  = 1.5d0, cb  = 1.5d0, db  = -1.d0/6.d0
	double precision :: betab= 5.d0, ab2 = -3.5d0, bb2 = 0.5d0, cb2 = 3.5d0, db2 = -0.5d0
end module
!!!-------------------------------------------------------------------

program BoundaryLayer_3D
	!$ use omp_lib
	implicit none
!!!---------------------------------------------------------------------------------------------------------------------
!!!-------calclating area-----------------------------------------------------------------------------------------------
	integer,parameter :: idata=100                           				!number of output files
	integer,parameter :: Nx=200+1, Ny=100+1, Nz=80+1, Nt=10000         		!The number of grating points, TimeStep clock(Nt)
	double precision :: lx=200.d0, ly=50d0, lz=40.d0, dx, dy, dz            	!x,y,z range
	double precision :: lt=200.d0, dt, t										!time
	double precision,parameter :: pi=atan(1.0d0)*4.0d0
!!!
	double precision,parameter :: rhol=1.d0, ul=0.d0, vl=0.d0, wl=0.d0, pl=1.d0		!init. condition of ShockTube
	double precision,parameter :: zeta=1.d0, gamma=1.4d0, epsilon=0.2d0
	double precision,parameter :: Re=1000.d0, Ma=0.5d0, T1=117.d0, myu0=1.d0!822d-5
	double precision,parameter :: Pra=0.71d0, T0=20.d0+273.15d0
!!!
	integer :: i, j, k, l, hNx, num=1, timestep=0,num2,numt
	integer :: TIME_ARRAY(3)
!!!Main phisical quantities
	double precision :: p(1:Nx,1:Ny,1:Nz), Et(1:Nx,1:Ny,1:Nz), TT(1:Nx,1:Ny,1:Nz), rho(1:Nx,1:Ny,1:Nz)
	double precision :: x(1:Nx), y(1:Ny), z(1:Nz), myu(1:Nx,1:Ny,1:Nz),lambda(1:Nx, 1:Ny, 1:Nz)
	double precision :: Q1(1:Nx,1:Ny,1:Nz), Q2(1:Nx,1:Ny,1:Nz), Q3(1:Nx,1:Ny,1:Nz)
	double precision :: Q4(1:Nx,1:Ny,1:Nz), Q5(1:Nx,1:Ny,1:Nz), Q_ten(1:Nx,1:Ny,1:Nz)
!!!-------LU------!!!2016.07.30
	double precision :: L_xdp(1:Nx,1:Nx), U_xdp(1:Nx,1:Nx), L_xc(1:Nx,1:Nx), U_xc(1:Nx,1:Nx)
	double precision :: L_xdm(1:Nx,1:Nx), U_xdm(1:Nx,1:Nx)
	double precision :: L_ydp(1:Ny,1:Ny), U_ydp(1:Ny,1:Ny), L_yc(1:Ny,1:Ny), U_yc(1:Ny,1:Ny)
	double precision :: L_ydm(1:Ny,1:Ny), U_ydm(1:Ny,1:Ny)
	double precision :: L_zdp(1:Nz,1:Nz), U_zdp(1:Nz,1:Nz), L_zc(1:Nz,1:Nz), U_zc(1:Nz,1:Nz)
	double precision :: L_zdm(1:Nz,1:Nz), U_zdm(1:Nz,1:Nz)
!EE,RK4
	double precision :: Q1_next(1:Nx,1:Ny,1:Nz), Q2_next(1:Nx,1:Ny,1:Nz), Q3_next(1:Nx,1:Ny,1:Nz)
	double precision :: Q4_next(1:Nx,1:Ny,1:Nz), Q5_next(1:Nx,1:Ny,1:Nz)
	double precision :: Q1_0(1:Nx,1:Ny,1:Nz), Q2_0(1:Nx,1:Ny,1:Nz), Q3_0(1:Nx,1:Ny,1:Nz)
	double precision :: Q4_0(1:Nx,1:Ny,1:Nz), Q5_0(1:Nx,1:Ny,1:Nz)
	double precision :: Q1_1(1:Nx,1:Ny,1:Nz), Q2_1(1:Nx,1:Ny,1:Nz), Q3_1(1:Nx,1:Ny,1:Nz)
	double precision :: Q4_1(1:Nx,1:Ny,1:Nz), Q5_1(1:Nx,1:Ny,1:Nz)
	double precision :: Q1_2(1:Nx,1:Ny,1:Nz), Q2_2(1:Nx,1:Ny,1:Nz), Q3_2(1:Nx,1:Ny,1:Nz)
	double precision :: Q4_2(1:Nx,1:Ny,1:Nz), Q5_2(1:Nx,1:Ny,1:Nz)
	double precision :: Q1_3(1:Nx,1:Ny,1:Nz), Q2_3(1:Nx,1:Ny,1:Nz), Q3_3(1:Nx,1:Ny,1:Nz)
	double precision :: Q4_3(1:Nx,1:Ny,1:Nz), Q5_3(1:Nx,1:Ny,1:Nz)
	double precision :: Q1_4(1:Nx,1:Ny,1:Nz), Q2_4(1:Nx,1:Ny,1:Nz), Q3_4(1:Nx,1:Ny,1:Nz)
	double precision :: Q4_4(1:Nx,1:Ny,1:Nz), Q5_4(1:Nx,1:Ny,1:Nz)
!
	double precision :: k1_1(1:Nx,1:Ny,1:Nz), k2_1(1:Nx,1:Ny,1:Nz), k3_1(1:Nx,1:Ny,1:Nz)
	double precision :: k4_1(1:Nx,1:Ny,1:Nz), k5_1(1:Nx,1:Ny,1:Nz)
	double precision :: k1_2(1:Nx,1:Ny,1:Nz), k2_2(1:Nx,1:Ny,1:Nz), k3_2(1:Nx,1:Ny,1:Nz)
	double precision :: k4_2(1:Nx,1:Ny,1:Nz), k5_2(1:Nx,1:Ny,1:Nz)
	double precision :: k1_3(1:Nx,1:Ny,1:Nz), k2_3(1:Nx,1:Ny,1:Nz), k3_3(1:Nx,1:Ny,1:Nz)
	double precision :: k4_3(1:Nx,1:Ny,1:Nz), k5_3(1:Nx,1:Ny,1:Nz)
!Flux
	double precision :: Fx_1(1:Nx,1:Ny,1:Nz), Fx_2(1:Nx,1:Ny,1:Nz), Fx_3(1:Nx,1:Ny,1:Nz)
	double precision :: Fx_4(1:Nx,1:Ny,1:Nz), Fx_5(1:Nx,1:Ny,1:Nz)
	double precision :: Fy_1(1:Nx,1:Ny,1:Nz), Fy_2(1:Nx,1:Ny,1:Nz), Fy_3(1:Nx,1:Ny,1:Nz)
	double precision :: Fy_4(1:Nx,1:Ny,1:Nz), Fy_5(1:Nx,1:Ny,1:Nz)
	double precision :: Fz_1(1:Nx,1:Ny,1:Nz), Fz_2(1:Nx,1:Ny,1:Nz), Fz_3(1:Nx,1:Ny,1:Nz)
	double precision :: Fz_4(1:Nx,1:Ny,1:Nz), Fz_5(1:Nx,1:Ny,1:Nz)
!bibun
	double precision :: dFx_1(1:Nx,1:Ny,1:Nz), dFx_2(1:Nx,1:Ny,1:Nz), dFx_3(1:Nx,1:Ny,1:Nz)
	double precision :: dFx_4(1:Nx,1:Ny,1:Nz), dFx_5(1:Nx,1:Ny,1:Nz)
	double precision :: dFy_1(1:Nx,1:Ny,1:Nz), dFy_2(1:Nx,1:Ny,1:Nz), dFy_3(1:Nx,1:Ny,1:Nz)
	double precision :: dFy_4(1:Nx,1:Ny,1:Nz), dFy_5(1:Nx,1:Ny,1:Nz)
	double precision :: dFz_1(1:Nx,1:Ny,1:Nz), dFz_2(1:Nx,1:Ny,1:Nz), dFz_3(1:Nx,1:Ny,1:Nz)
	double precision :: dFz_4(1:Nx,1:Ny,1:Nz), dFz_5(1:Nx,1:Ny,1:Nz)
!visco
	double precision :: dFvx_2(1:Nx,1:Ny,1:Nz), dFvx_3(1:Nx,1:Ny,1:Nz), dFvx_4(1:Nx,1:Ny,1:Nz), dFvx_5(1:Nx,1:Ny,1:Nz)
	double precision :: dFvy_2(1:Nx,1:Ny,1:Nz), dFvy_3(1:Nx,1:Ny,1:Nz), dFvy_4(1:Nx,1:Ny,1:Nz), dFvy_5(1:Nx,1:Ny,1:Nz)
	double precision :: dFvz_2(1:Nx,1:Ny,1:Nz), dFvz_3(1:Nx,1:Ny,1:Nz), dFvz_4(1:Nx,1:Ny,1:Nz), dFvz_5(1:Nx,1:Ny,1:Nz)
	double precision :: rhorr(1:Nx,1:Ny),CFLx,CFLx_2,CFLy,CFLy_2,CFLz,CFLz_2
	double precision :: tau_xx(1:Nx,1:Ny,1:Nz), tau_xy(1:Nx,1:Ny,1:Nz), tau_xz(1:Nx,1:Ny,1:Nz)
	double precision :: tau_yy(1:Nx,1:Ny,1:Nz), tau_yz(1:Nx,1:Ny,1:Nz), tau_zz(1:Nx,1:Ny,1:Nz)
!!!
	double precision :: u(1:Nx, 1:Ny, 1:Nz), v(1:Nx, 1:Ny, 1:Nz), w(1:Nx, 1:Ny, 1:Nz)
	double precision :: dudx(1:Nx, 1:Ny, 1:Nz), dvdx(1:Nx, 1:Ny, 1:Nz), dwdx(1:Nx, 1:Ny, 1:Nz)
	double precision :: dudy(1:Nx, 1:Ny, 1:Nz), dvdy(1:Nx, 1:Ny, 1:Nz), dwdy(1:Nx, 1:Ny, 1:Nz)
	double precision :: dudz(1:Nx, 1:Ny, 1:Nz), dvdz(1:Nx, 1:Ny, 1:Nz), dwdz(1:Nx, 1:Ny, 1:Nz)
	double precision :: lambda_dx(1:Nx, 1:Ny, 1:Nz), lambda_dy(1:Nx, 1:Ny, 1:Nz), lambda_dz(1:Nx, 1:Ny, 1:Nz)
	double precision :: quex(1:Nx, 1:Ny, 1:Nz), quey(1:Nx, 1:Ny, 1:Nz), quez(1:Nx, 1:Ny, 1:Nz)
	double precision :: quex_dx(1:Nx, 1:Ny, 1:Nz), quey_dy(1:Nx, 1:Ny, 1:Nz), quez_dz(1:Nx, 1:Ny, 1:Nz)
	double precision :: tau_xx_dx(1:Nx, 1:Ny, 1:Nz), tau_yx_dx(1:Nx, 1:Ny, 1:Nz), tau_zx_dx(1:Nx, 1:Ny, 1:Nz)
	double precision :: tau_xy_dy(1:Nx, 1:Ny, 1:Nz), tau_yy_dy(1:Nx, 1:Ny, 1:Nz), tau_zy_dy(1:Nx, 1:Ny, 1:Nz)
	double precision :: tau_xz_dz(1:Nx, 1:Ny, 1:Nz), tau_yz_dz(1:Nx, 1:Ny, 1:Nz), tau_zz_dz(1:Nx, 1:Ny, 1:Nz)
	double precision :: TTx(1:Nx, 1:Ny, 1:Nz), TTy(1:Nx, 1:Ny, 1:Nz), TTz(1:Nx, 1:Ny, 1:Nz)
	double precision :: TTx2(1:Nx, 1:Ny, 1:Nz), TTy2(1:Nx, 1:Ny, 1:Nz), TTz2(1:Nx, 1:Ny, 1:Nz)
	double precision :: y_s(Ny),u_s(Ny),v_s(Ny),T_s(Ny),rho_s(Ny) !!shootings
	double precision :: x_s2(Nx,Ny),y_s2(Nx,Ny),u_s2(Nx,Ny),v_s2(Nx,Ny),T_s2(Nx,Ny),rho_s2(Nx,Ny)
!!!NSCBC
	double precision :: Ex_NSCBCin(1:5,1:Ny,1:Nz), Ex_NSCBCout(1:5,1:Ny,1:Nz), Ey_NSCBCin(1:5,1:Nx,1:Nz), Ey_NSCBCout(1:5,1:Nx,1:Nz)
!!expand
	double precision :: Cgs1,Cgs2,my(1:Ny)
!!!kakuran
	double precision :: Blank1(1:Ny), Blank2(1:Ny), y_u(1:Ny), u_func(1:Ny), u_func_r(1:Ny), u_func_i(1:Ny)
	double precision :: y_v(1:Ny), v_func(1:Ny), v_func_r(1:Ny), v_func_i(1:Ny)
	double precision :: u_kaku(1:Ny,1:Nz), v_kaku(1:Ny,1:Nz)
	double precision :: alpha=0.22d0, beta=0.d0, omegai=0.00188d0, omegar=0.07766d0, eigen_rate
!!!arg
	double precision :: arg_ufunc(1:Ny), amp_ufunc(1:Ny), arg_vfunc(1:Ny), amp_vfunc(1:Ny)
	integer :: TURB
!!!------------------------------------------------------------------------------------------------------------------
	character(len=13) :: Filename
	character(len=29) :: Filename_u, Filename_v
!!!!
	character(len=4) :: Ma_chara
	character(len=6) :: Re_chara
!	integer :: Re_chara
!	Re_chara1 = int(Re)     ! 実数から文字列への変換
!	write (Re_chara,'(i6.6)') Re_chara1     ! 実数から文字列への変換
!	write (Ma_chara,'(f4.2)') Ma     ! 実数から文字列への変換

	Filename = 'result000.prn'!!!!!!!!!!!!!!!!!!!!!フォルダ名
	Filename_u = 'M0.00Re000000a0.22b0.00ux.dat'
	Filename_v = 'M0.00Re000000a0.22b0.00vx.dat'

	write(Filename_u(2:5), '(f4.2)') Ma
	write(Filename_v(2:5), '(f4.2)') Ma
	write(Filename_u(8:13), '(i6.6)') int(Re)
	write(Filename_v(8:13), '(i6.6)') int(Re)

!!!
	dx = lx/dble(Nx-1)                                 				!x stride
	dy = ly/dble(Ny-1)                                 				!y stride
	dz = lz/dble(Nz-1)                                 				!y stride
	dt = lt/dble(Nt)                                 				!t stride
	CFLx=dt/dx
	CFLx_2=dt/(dx*dx)
	CFLy=dt/dy
	CFLy_2=dt/(dy*dy)
	CFLz=dt/dz
	CFLz_2=dt/(dz*dz)
	hNx=Nx/10
!!!-----CFL judgement--------
	write(*,*) 'CFLx=', CFLx
	write(*,*) 'CFL2x=', CFLx_2
	write(*,*) 'CFLy=', CFLy
	write(*,*) 'CFL2y=', CFLy_2
	write(*,*) 'CFLz=', CFLz
	write(*,*) 'CFL2z=', CFLz_2
	if(CFLx>=1 .or. CFLx_2>=1 .or. CFLy>=1 .or. CFLy_2>=1 .or. CFLz>=1 .or. CFLz_2>=1) stop 'CFL error'

!!!----------shooting read
	open(92, file='init.txt', status='old')
	do i=1,Nx
	do j=1,Ny
		read (92,*) x_s2(i,j), y_s2(i,j), u_s2(i,j), v_s2(i,j), T_s2(i,j), rho_s2(i,j)
	end do
	end do
	close(92)
!!!----------eigen func u
	open(93, file=Filename_u, status='old')
	do j=1,Ny
		read (93,*) Blank1(j), Blank2(j), y_u(j), u_func(j), u_func_r(j), u_func_i(j)
	end do
	close(93)
!!!----------eigen func u
	open(94, file=Filename_v, status='old')
	do j=1,Ny
		read (94,*) Blank1(j), Blank2(j), y_v(j), v_func(j), v_func_r(j), v_func_i(j)
	end do
	close(94)
!!!!!!!---------------------------------------------eigen func rate
	eigen_rate=0.01d0
	!$omp parallel do
	do j=1,Ny
	arg_ufunc(j) = dcos ( dasin( u_func_i(j)/u_func(j) ) )
	amp_ufunc(j) = eigen_rate*arg_ufunc(j) * u_func(j)
	arg_vfunc(j) = dcos ( dasin( v_func_i(j)/v_func(j) ) )
	amp_vfunc(j) = eigen_rate*arg_vfunc(j) * v_func(j)
	end do
	!$omp end parallel do
!!!--------------------------------------------------
	open(21, file='funccheck.dat', status='unknown', form='formatted')
	do j=1, Ny
		write(21, "(8(1x,g15.7))") Blank1(j), Blank2(j), y_u(j), u_func(j), arg_ufunc(j), amp_ufunc(j), dasin( u_func_i(j)/u_func(j) ), u_func_i(j)/u_func(j)
	end do
	close(21)
!!!-----init. condition setup-----------------
	!$omp parallel do
	do k=1, Nz
		do j=1, Ny
			do i=1, Nx
				p(i,j,k)	 = rho_s2(i,j)*T_s2(i,j)/gamma/Ma/Ma
			end do
		end do
	end do
	!$omp end parallel do

	!$omp parallel do
	do i=1,Nx
		x(i) = dx*dble(i-1)
	end do
	!$omp end parallel do

	!$omp parallel do
	do j=1,Ny
		y(j) = dy*dble(j-1)
!		my(j) = dy*dble(j-1)
	end do
	!$omp end parallel do

!!expand-------------------(tuning now)
!		Cgs1 = 0.05d0         ! : coefficient of grid stretching
!		do j=0,Ny
!			y(j)=dtanh(Cgs1*(-ly+j*dy))-dtanh(Cgs1*(-ly))
!		end do
!		Cgs2 = ly/y(Ny)   ! : coefficient of grid stretching
!		do j=0,Ny           ! : coordinate data of y-direction after grid stretcing
!			y(j)=Cgs2*(dtanh(Cgs1*(-ly+j*dy))-dtanh(Cgs1*(-ly)))
!		end do
!!!

	!$omp parallel do
	do k=1,Nz
		z(k) = dz*dble(k-1)
	end do
	!$omp end parallel do

	!$omp parallel do
	do k=1, Nz
		do j=1, Ny
			do i=1, Nx!hNx
				Et(i,j,k)    = 0.5d0*rho_s2(i,j)*(u_s2(i,j)*u_s2(i,j) + v_s2(i,j)*v_s2(i,j) + wl*wl) &
								&+ rho_s2(i,j)*T_s2(i,j)/gamma/Ma/Ma/(gamma - 1.d0)!!!
				Q1(i,j,k)   = rho_s2(i,j)!!!
				Q2(i,j,k)   = rho_s2(i,j)*u_s2(i,j)!!!
				Q3(i,j,k)   = rho_s2(i,j)*v_s2(i,j)!!!
				Q4(i,j,k)   = 0d0!rho_s2(i,j)*wl		!!!
				Q5(i,j,k)    = 0.5d0*rho_s2(i,j)*(u_s2(i,j)*u_s2(i,j) + v_s2(i,j)*v_s2(i,j) )&!+ wl*wl) &
								&+ rho_s2(i,j)*T_s2(i,j)/gamma/Ma/Ma/(gamma - 1.d0)!!!
!				TT(i,j,k) = T_s2(i,j)
				TT(i,j,k) = gamma*Ma*Ma/Q1(i,j,k)*(gamma - 1.d0)*( Q5(i,j,k) - 0.5d0/Q1(i,j,k)*( Q2(i,j,k)*Q2(i,j,k) &
						&+ Q3(i,j,k)*Q3(i,j,k) + Q4(i,j,k)*Q4(i,j,k)) )
				myu(i,j,k) =myu0*( T0 + T1 )/( T_s2(i,j) + T1)*T_s2(i,j)*dsqrt(T_s2(i,j))/( T0*dsqrt(T0) )
			end do
		end do
	end do
	!$omp end parallel do
	call visco_bibun()
   !$OMP parallel workshare
	dudy(1:Nx,1:Ny,1:Nz)=0.d0
	Q_ten(1:Nx,1:Ny,1:Nz)=0.d0
   !$OMP end parallel workshare
!!!-----init. condition writing---------------
	open(11, file=Filename, status='unknown', form='formatted')
	do k=1, Nz
		do j=1, Ny
			do i=1, Nx
!				write(11, "(11(1x,g15.7))") x(i), y(j), z(k), Q1(i,j,k), Q2(i,j,k)/Q1(i,j,k), &
!						&Q3(i,j,k)/Q1(i,j,k), Q4(i,j,k)/Q1(i,j,k), p(i,j,k), TT(i,j,k), dudy(i,j,k), Q_ten(i,j,k)
				write(11, "(10(1x,g15.7))") x(i), y(j), z(k), Q1(i,j,k), Q2(i,j,k)/Q1(i,j,k), &
						&Q3(i,j,k)/Q1(i,j,k), Q4(i,j,k)/Q1(i,j,k), p(i,j,k), TT(i,j,k), Q_ten(i,j,k)
!				write(11, "(9(1x,g15.7))") x(i), y(j), z(k), Q1(i,j,k), Q2(i,j,k)/Q1(i,j,k), Q3(i,j,k)/Q1(i,j,k), Q4(i,j,k)/Q1(i,j,k), p(i,j,k), TT(i,j,k)
			end do
!			write(11,*)
		end do
!		write(11,*)
!		write(11,*)
	end do
	close(11)
!!!-----------------------------------------------
	call BL_init()!!流入部
	open(11, file='result000_kaku.prn', status='unknown', form='formatted')
	do k=1, Nz
		do j=1, Ny
			do i=1, Nx
				write(11, "(10(1x,g15.7))") x(i), y(j), z(k), Q1(i,j,k), Q2(i,j,k)/Q1(i,j,k), &
						&Q3(i,j,k)/Q1(i,j,k), Q4(i,j,k)/Q1(i,j,k), p(i,j,k), TT(i,j,k), Q_ten(i,j,k)
			end do
		end do
	end do
	close(11)
!!!!-------LU Bunkai--------------------------------------
call itime(TIME_ARRAY)
write(*,*) 'LU calc',  TIME_ARRAY(1), TIME_ARRAY(2), TIME_ARRAY(3)
call LUbunkai(2, 2, Nx, 1, L_xdp, U_xdp, -0.25d0)
call LUbunkai(2, 2, Nx, 2, L_xdm, U_xdm,  0.25d0)
call LUbunkai(1, 2, Nx, 0, L_xc , U_xc ,  0.d0)
call LUbunkai(2, 2, Ny, 1, L_ydp, U_ydp, -0.25d0)
call LUbunkai(2, 2, Ny, 2, L_ydm, U_ydm,  0.25d0)
call LUbunkai(1, 2, Ny, 0, L_yc , U_yc ,  0.d0)
call LUbunkai(2, 1, Nz, 1, L_zdp, U_zdp, -0.25d0)
call LUbunkai(2, 1, Nz, 2, L_zdm, U_zdm,  0.25d0)
call LUbunkai(1, 1, Nz, 0, L_zc , U_zc ,  0.d0)
call itime(TIME_ARRAY)
write(*,*) 'LU OK',  TIME_ARRAY(1), TIME_ARRAY(2), TIME_ARRAY(3)
!!!--------------------------------------------------------------------------------
!!!----------------------------------------------
!┌─★┌─★┌─★┌─★┌─★┌─★┌─★┌─★
!│ｔ││ｉ││ｍ││ｅ││ｓ││ｔ││ｅ││ｐ│
!☆─┘☆─┘☆─┘☆─┘☆─┘☆─┘☆─┘☆─┘
!!!!---------------------------------------------
	do timestep=1, Nt
		t= dt * dble(timestep)
!		call BL_init()!!流入部
!		call euler()
		call RK4()
		call update()
		if(mod(timestep,Nt/idata)==0) then
			call output_files()
			num=num+1
		end if
		if(mod(timestep,10)==0) then
			call check()
		end if
		call NaN_check()
	end do
!!!===================================================================
!=====================================================================
!┌─★┌─★┌─★┌─★┌─★┌─★┌─★┌─★┌─★┌─★
!│ｓ││ｕ││ｂ││ｒ││ｏ││ｕ││ｔ││ｉ││ｎ││ｅ│
!☆─┘☆─┘☆─┘☆─┘☆─┘☆─┘☆─┘☆─┘☆─┘☆─┘
	contains
!=====================================================================
!!!-------------------------------------------------------
!!!-----1:X range(5th DCS)--------------
	subroutine LUbunkai(bib, b_cdtn, NN, chvel_pm, lm, um, sgm)
	use DCS5
	implicit none

	integer, intent(in)    :: bib!!CCSorDCSquality
	integer, intent(in)    :: NN!!グリッド数
	integer, intent(in)    :: b_cdtn!!周期ｏｒ片側
	integer, intent(in)    :: chvel_pm!!特性速度(0:CCS 1:正 2:負)
	integer :: ii, jj, kk

	double precision, intent(in) :: sgm
	double precision, intent(out) :: lm(1:NN,1:NN)!input
	double precision, intent(out) :: um(1:NN,1:NN)!input
	double precision,allocatable,dimension(:,:)  :: am
!!!
	double precision :: wa
 	allocate( am(NN,NN))


!!!-----(LU)-----------------------
!!!-----Matrix A definition------------
   !$OMP parallel workshare
	am(:,:)  = 0.d0
   !$OMP end parallel workshare
	am(1,1)  = 1.d0
!!--periodic-----------------------------------
	if(b_cdtn==1) then
		am(1,NN-1) = alpc*(1.d0 - alpd*sgm)
		am(NN,2)   = alpc*(1.d0 + alpd*sgm)
	!$omp parallel do
		do ii=2, NN
			am(ii,ii)   = 1.d0
			am(ii,ii-1) = alpc*(1.d0 - alpd*sgm)
			am(ii-1,ii) = alpc*(1.d0 + alpd*sgm)
		end do
	!$omp end parallel do
!!--One-Sided---------------------------------
	else if(b_cdtn==2) then
	!$omp parallel do
		do ii=1, NN-4
			am(ii+2,ii+1) = alpc*(1.d0 - alpd*sgm)
			am(ii+2,ii+3) = alpc*(1.d0 + alpd*sgm)
		end do
	!$omp end parallel do
	!$omp parallel do
		do ii=2, NN
			am(ii,ii)       = 1.d0
		end do
	!$omp end parallel do
!		if(chvel_pm==0 .or. chvel_pm==1) then
			am(1,2)       = alpb
!		else if(chvel_pm==2) then
!			am(1,2)       = betab
!		end if
!		if(chvel_pm==0 .or. chvel_pm==2) then
			am(NN,NN-1)   = alpb
!		else if(chvel_pm==1) then
!			am(NN,NN-1)   = betab
!		end if
		am(2,1)		  = 0.25d0!*(1.d0 - sgm)	!N=2,Nx-1 : 4thCCS(3rdDCS)
		am(2,3)		  = 0.25d0!*(1.d0 + sgm)
		am(NN-1,NN-2) = 0.25d0*(1.d0 - sgm)
		am(NN-1,NN)   = 0.25d0*(1.d0 + sgm)
	end if
!!!-----(LU分解)
!!!---LUを求める
	!$omp parallel do
	do ii=1,NN
		lm(ii,ii) = 1.d0	!!!-L行列の対角(i=j)成分
	end do
	!$omp end parallel do

	do jj=1, NN
		do ii=1, jj
			if(ii==1) then
				wa      = 0.d0
				um(ii,jj) = am(ii,jj)
			else
				do kk=1,ii-1
					wa = wa + lm(ii,kk)*um(kk,jj)
				end do
				um(ii,jj) = am(ii,jj) - wa
				wa = 0.d0
			end if
		end do

		do ii=jj+1, NN
			if(jj==1) then
				wa = 0.d0
				lm(ii,jj) = am(ii,jj)/um(jj,jj)
			else
				do kk=1,jj-1
					wa = wa + lm(ii,kk)*um(kk,jj)
				end do
				lm(ii,jj) = ( am(ii,jj) - wa )/um(jj,jj)
				wa = 0.d0
			end if
		end do
	end do
	end subroutine LUbunkai
!======================================================================================================
	subroutine bibun_DCS(bib, b_cdtn, NN, dn, sgm, lm, um, u_in, u_out,chvel_pm)!!!2016.07.06 chvel_pm追加
	use DCS5
	implicit none
	integer, intent(in)    :: bib!!CCSorDCSquality
	integer, intent(in)    :: NN!!グリッド数
	integer, intent(in)    :: b_cdtn!!周期ｏｒ片側
	integer, intent(in)    :: chvel_pm!!特性速度(0:CCS 1:正 2:負)
	integer :: ii, jj, kk

	double precision, intent(in) :: sgm
	double precision, intent(in) :: u_in(1:NN)!input
	double precision, intent(out) :: u_out(1:NN)!output
	double precision, intent(in) :: lm(1:NN,1:NN)!input
	double precision, intent(in) :: um(1:NN,1:NN)!input
	double precision,allocatable,dimension(:)  :: b, y1, y2
	double precision :: dn
!!!
	double precision :: wa

 	allocate( b(NN), y1(NN), y2(NN) )
!!!===========================================================================================================================
!!!-----vector b calclation (b_cdtn=1:periopdic，b_cdtn=2:one-sided scheme)
	!$omp parallel workshare
	b(1:NN)=0.d0
	!$omp end parallel workshare
!!--periodic-----------------------------------
	if(b_cdtn==1) then
		ii = 1
			b(1) = ( 0.5d0*ac*( u_in(1+1)-u_in(NN-1) ) +0.25d0*bc*( u_in(1+2)-u_in(NN-2) ) &
			      &+sgm*( ad*( u_in(1+1)-2.d0*u_in(1)+u_in(NN-1) ) + bd*( u_in(1+2)-2.d0*u_in(1)+u_in(NN-2) )*0.25d0 ) )/dn
!			b(1) = ( 0.5d0*ac*( u_in(1+1)-u_in(NN) ) +0.25d0*bc*( u_in(1+2)-u_in(NN-1) ) &
!			      &+sgm*( ad*( u_in(1+1)-2.d0*u_in(1)+u_in(NN) ) + bd*( u_in(1+2)-2.d0*u_in(1)+u_in(NN-1) )*0.25d0 ) )/dn
		ii = 2
			b(2) = ( 0.5d0*ac*( u_in(2+1)-u_in(2-1) ) +0.25d0*bc*( u_in(2+2)-u_in(NN-1) ) &
			      &+sgm*( ad*( u_in(2+1)-2.d0*u_in(2)+u_in(2-1) ) + bd*( u_in(2+2)-2.d0*u_in(2)+u_in(NN-1) )*0.25d0 ) )/dn
!			b(2) = ( 0.5d0*ac*( u_in(2+1)-u_in(2-1) ) +0.25d0*bc*( u_in(2+2)-u_in(NN) ) &
!			      &+sgm*( ad*( u_in(2+1)-2.d0*u_in(2)+u_in(2-1) ) + bd*( u_in(2+2)-2.d0*u_in(2)+u_in(NN) )*0.25d0 ) )/dn
!!--one-sided-----------------------------------
	else if(b_cdtn==2) then
		ii = 1
!		if(chvel_pm==0 .or. chvel_pm==1) then
			b(1) = ( ab*u_in(1) + bb*u_in(2) + cb*u_in(3) + db*u_in(4) )/dn
!		else if(chvel_pm==2) then
!			b(1) = ( ab2*u_in(1) + bb2*u_in(2) + cb2*u_in(3) + db2*u_in(4) )/dn
!		end if
		!!境界隣点
		ii = 2
			b(2) = ( 0.75d0*( u_in(2+1)-u_in(2-1) )  )/dn!+ sgm*( u_in(2+1)-2.d0*u_in(2)+u_in(2-1) )*0.5d0 )/dn
	end if !!!!(periodic)

	!$omp parallel do
	do ii=3,NN-2
		b(ii) = ( ac*( u_in(ii+1)-u_in(ii-1) )*0.5d0 +bc*( u_in(ii+2)-u_in(ii-2) )*0.25d0 &
		      &+sgm*( ad*( u_in(ii+1)-2.d0*u_in(ii)+u_in(ii-1) ) + bd*( u_in(ii+2)-2.d0*u_in(ii)+u_in(ii-2) )*0.25d0 ) )/dn
	end do
	!$omp end parallel do

	if(b_cdtn==1) then
		ii = NN-1
			b(NN-1) = ( 0.5d0*ac*( u_in(NN-1+1)-u_in(NN-1-1) ) +0.25d0*bc*( u_in(1+1)-u_in(NN-1-2) ) &
			      &+sgm*( ad*( u_in(NN-1+1)-2.d0*u_in(NN-1)+u_in(NN-1-1) ) + bd*( u_in(1+1)-2.d0*u_in(NN-1)+u_in(NN-1-2) )*0.25d0 ) )/dn
		ii = NN
			b(NN) = ( 0.5d0*ac*( u_in(1+1)-u_in(NN-1) ) +0.25d0*bc*( u_in(1+2)-u_in(NN-2) ) &
			      &+sgm*( ad*( u_in(1+1)-2.d0*u_in(NN)+u_in(NN-1) ) + bd*( u_in(1+2)-2.d0*u_in(NN)+u_in(NN-2) )*0.25d0 ) )/dn
!!--one-sided-----------------------------------
	else if(b_cdtn==2) then
		ii = NN-1    !!!!!!!修正
			b(NN-1) = ( 0.75d0*( u_in(NN-1+1)-u_in(NN-1-1) ) + sgm*( u_in(NN-1+1)-2.d0*u_in(NN-1)+u_in(NN-1-1) )*0.5d0 )/dn
		ii = NN
!		if(chvel_pm==0 .or. chvel_pm==2) then
				b(NN) =  -( ab*u_in(NN) + bb*u_in(NN-1) + cb*u_in(NN-2) + db*u_in(NN-3) )/dn
!		else if(chvel_pm==1) then
!				b(NN) =  -( ab2*u_in(NN) + bb2*u_in(NN-1) + cb2*u_in(NN-2) + db2*u_in(NN-3) )/dn
!		end if
	end if

!!!---前進代入
	ii = 1
		y1(ii) = b(ii)

	do ii=2, NN
		do jj=1, ii-1
			wa = wa + lm(ii,jj)*y1(jj)
		end do
		y1(ii) = b(ii) - wa
		wa = 0.d0
	end do
!!!---後退代入
	ii = NN
	y2(NN) = y1(NN)/um(NN,NN)

	do ii=NN-1, 1, -1
		do jj=ii+1, NN
			wa = wa + um(ii,jj)*y2(jj)
		end do
		y2(ii) = ( y1(ii) - wa )/um(ii,ii)
		wa = 0.d0
	end do
	!$OMP parallel workshare
	u_out(1:NN)=y2(1:NN)
	!$OMP end parallel workshare
 	deallocate( b, y1, y2)
	end subroutine bibun_DCS
!=====================================================================
!=====================================================================
	subroutine bunkatu(Qv_in1,Qv_in2,Qv_out,x_y_z)
	integer, intent(in)    :: x_y_z!!1:X, 2:Y, 3:Z
	integer :: iv,jv,kv
	double precision,intent(in) :: Qv_in1(1:Nx,1:Ny,1:Nz)
	double precision,intent(in) :: Qv_in2(1:Nx,1:Ny,1:Nz)
	double precision,intent(out) :: Qv_out(1:Nx,1:Ny,1:Nz)
	double precision :: Qv_plus(1:Nx,1:Ny,1:Nz),  Qv_minus(1:Nx,1:Ny,1:Nz)
	double precision :: dQv_plus(1:Nx,1:Ny,1:Nz),dQv_minus(1:Nx,1:Ny,1:Nz)
	!$OMP parallel workshare
	Qv_plus(1:Nx,1:Ny,1:Nz) = Qv_in1(1:Nx,1:Ny,1:Nz) + zeta *Qv_in2(1:Nx,1:Ny,1:Nz)
	Qv_minus(1:Nx,1:Ny,1:Nz) = Qv_in1(1:Nx,1:Ny,1:Nz) - zeta *Qv_in2(1:Nx,1:Ny,1:Nz)
	!$OMP end parallel workshare
	if(x_y_z==1) then
   !$OMP parallel do
		do kv=1,Nz
		do jv=1,Ny
			call bibun_DCS(2, 2, Nx, dx, -0.25d0, L_xdp, U_xdp, Qv_plus(:,jv,kv), dQv_plus(:,jv,kv),1)
			call bibun_DCS(2, 2, Nx, dx,  0.25d0, L_xdm, U_xdm, Qv_minus(:,jv,kv), dQv_minus(:,jv,kv),2)
		end do
		end do
   !$OMP end parallel do
	else if(x_y_z==2) then
   !$OMP parallel do
		do kv=1,Nz
		do iv=1,Nx
			call bibun_DCS(2, 2, Ny, dy, -0.25d0, L_ydp, U_ydp, Qv_plus(iv,:,kv), dQv_plus(iv,:,kv),1)
			call bibun_DCS(2, 2, Ny, dy, 0.25d0, L_ydm, U_ydm, Qv_minus(iv,:,kv), dQv_minus(iv,:,kv),2)
		end do
		end do
   !$OMP end parallel do
	else if(x_y_z==3) then
   !$OMP parallel do
		do jv=1,Ny
		do iv=1,Nx
			call bibun_DCS(2, 1, Nz, dz, -0.25d0, L_zdp, U_zdp, Qv_plus(iv,jv,:), dQv_plus(iv,jv,:),1)
			call bibun_DCS(2, 1, Nz, dz, 0.25d0, L_zdm, U_zdm, Qv_minus(iv,jv,:), dQv_minus(iv,jv,:),2)
		end do
		end do
   !$OMP end parallel do
	end if
   !$OMP parallel workshare
	Qv_out(1:Nx,1:Ny,1:Nz) = 0.5d0* ( dQv_plus(1:Nx,1:Ny,1:Nz) + dQv_minus(1:Nx,1:Ny,1:Nz) )
   !$OMP end parallel workshare
	end subroutine bunkatu
!=====================================================================
!=====================================================================
	subroutine CCS_3D(Qv_in1,Qv_out,x_y_z)
	integer, intent(in)    :: x_y_z!!1:X, 2:Y, 3:Z
	integer :: iv,jv,kv
	double precision,intent(in) :: Qv_in1(1:Nx,1:Ny,1:Nz)
	double precision,intent(out) :: Qv_out(1:Nx,1:Ny,1:Nz)
	if(x_y_z==1) then
   !$OMP parallel do
		do kv=1,Nz
		do jv=1,Ny
			call bibun_DCS(1, 2, Nx, dx, 0.d0, L_xc, U_xc, Qv_in1(:,jv,kv), Qv_out(:,jv,kv),0)
		end do
		end do
   !$OMP end parallel do
	else if(x_y_z==2) then
   !$OMP parallel do
		do kv=1,Nz
		do iv=1,Nx
			call bibun_DCS(1, 2, Ny, dy, 0.d0, L_yc, U_yc, Qv_in1(iv,:,kv), Qv_out(iv,:,kv),0)
		end do
		end do
   !$OMP end parallel do
	else if(x_y_z==3) then
   !$OMP parallel do
		do jv=1,Ny
		do iv=1,Nx
			call bibun_DCS(1, 1, Nz, dz, 0.d0, L_zc, U_zc, Qv_in1(iv,jv,:), Qv_out(iv,jv,:),0)
		end do
		end do
   !$OMP end parallel do
	end if
	end subroutine CCS_3D
!!!-----bounary layer initial condition-----------------------------------------------------------------------
	subroutine BL_init()
	!撹乱
	!$omp parallel do
	do k=1,Nz
	do j=1,Ny
		u_kaku(j,k) = u_s2(1,j)
		v_kaku(j,k) = v_s2(1,j)
	end do
	end do
	!$omp end parallel do
	!撹乱
	!$omp parallel do
	do k=1,Nz
	do j=1,Ny
		u_kaku(j,k) = u_s2(1,j) + amp_ufunc(j) * dcos( 2d0*pi/omegar/dble(lt)*(z(k)*beta - t) ) * dexp(2d0*pi/dble(lt)*omegar*t)
		v_kaku(j,k) = v_s2(1,j) + amp_vfunc(j) * dcos( 2d0*pi/omegar/dble(lt)*(z(k)*beta - t) ) * dexp(2d0*pi/dble(lt)*omegar*t)
	end do
	end do
	!$omp end parallel do
	!$omp parallel do
	do k=1,Nz
	do j=1,Ny
		Q1(1,j,k)   = rho_s2(1,j)!!!
		Q2(1,j,k)   = rho_s2(1,j)*u_kaku(j,k)!!!
		Q3(1,j,k)   = rho_s2(1,j)*v_kaku(j,k)!!!
		Q4(1,j,k)   = 0.d0!rho_s2(1,j)*wl		!!!
		Q5(1,j,k)   = 0.5d0*rho_s2(1,j)*(u_kaku(j,k)*u_kaku(j,k) + v_kaku(j,k)*v_kaku(j,k) ) &
							&+ rho_s2(1,j)*T_s2(1,j)/gamma/Ma/Ma/(gamma - 1.d0)!!!
	end do
	!!壁面
	do i=1,Nx
		Q2(i,1,k)   = 0.d0
		Q3(i,1,k)   = 0.d0
		Q4(i,1,k)   = 0.d0
		Q5(i,1,k)   = Q1(i,1,k)*TT(i,1,k)/gamma/Ma/Ma/(gamma - 1.d0)!!!
	end do
	end do
	!$omp end parallel do
	end subroutine BL_init
!!!-----5:calc tau------------------------------------------------------
!=====================================================================
!=====================================================================
!!!-----6:update------------------------------------------------------
	subroutine update()
	!$OMP parallel workshare
	u(1:Nx,1:Ny,1:Nz) = Q2(1:Nx,1:Ny,1:Nz)/Q1(1:Nx,1:Ny,1:Nz)
	v(1:Nx,1:Ny,1:Nz) = Q3(1:Nx,1:Ny,1:Nz)/Q1(1:Nx,1:Ny,1:Nz)
	w(1:Nx,1:Ny,1:Nz) = Q4(1:Nx,1:Ny,1:Nz)/Q1(1:Nx,1:Ny,1:Nz)
	p(1:Nx,1:Ny,1:Nz) = (gamma - 1.d0)*( Q5(1:Nx,1:Ny,1:Nz) - 0.5d0*Q1(1:Nx,1:Ny,1:Nz)*( u(1:Nx,1:Ny,1:Nz)*u(1:Nx,1:Ny,1:Nz) &
						&+ v(1:Nx,1:Ny,1:Nz)*v(1:Nx,1:Ny,1:Nz) + w(1:Nx,1:Ny,1:Nz)*w(1:Nx,1:Ny,1:Nz)) )
	Et(1:Nx,1:Ny,1:Nz)= Q5(1:Nx,1:Ny,1:Nz)
	TT(1:Nx,1:Ny,1:Nz) = gamma*Ma*Ma/Q1(1:Nx,1:Ny,1:Nz)*(gamma - 1.d0)*( Q5(1:Nx,1:Ny,1:Nz) &
						&- 0.5d0/Q1(1:Nx,1:Ny,1:Nz)*( Q2(1:Nx,1:Ny,1:Nz)*Q2(1:Nx,1:Ny,1:Nz) &
						&+ Q3(1:Nx,1:Ny,1:Nz)*Q3(1:Nx,1:Ny,1:Nz) + Q4(1:Nx,1:Ny,1:Nz)*Q4(1:Nx,1:Ny,1:Nz)) )
	myu(1:Nx,1:Ny,1:Nz) =myu0*( T0 + T1 )/( TT(1:Nx,1:Ny,1:Nz) + T1)*TT(1:Nx,1:Ny,1:Nz)*dsqrt(TT(1:Nx,1:Ny,1:Nz))/( T0*dsqrt(T0) )
!	myu(1:Nx,1:Ny,1:Nz) = ( TT(1:Nx,1:Ny,1:Nz) )**(2.d0/3.d0)!!
	!$OMP end parallel workshare
	end subroutine update
!=====================================================================
!=====================================================================
!!!-----7:visco_bibun and FxFyFz renew-----------------------------------
	subroutine visco_bibun()
	call CCS_3D(u,dudx,1)
	call CCS_3D(v,dvdx,1)
	call CCS_3D(w,dwdx,1)
	call CCS_3D(u,dudy,2)
	call CCS_3D(v,dvdy,2)
	call CCS_3D(w,dwdy,2)
	call CCS_3D(u,dudz,3)
	call CCS_3D(v,dvdz,3)
	call CCS_3D(w,dwdz,3)
	!$OMP parallel workshare
	lambda(1:Nx,1:Ny,1:Nz) = myu(1:Nx,1:Ny,1:Nz)/( (gamma-1.d0)*Pra*Ma*Ma )
	tau_xx(1:Nx,1:Ny,1:Nz) = myu(1:Nx,1:Ny,1:Nz)*(4.d0/3.d0*dudx(1:Nx,1:Ny,1:Nz) - 2.d0/3.d0*(dvdy(1:Nx,1:Ny,1:Nz) &
								&+ dwdz(1:Nx,1:Ny,1:Nz)))
	tau_xy(1:Nx,1:Ny,1:Nz) = myu(1:Nx,1:Ny,1:Nz)*(dudy(1:Nx,1:Ny,1:Nz) + dvdx(1:Nx,1:Ny,1:Nz))
	tau_xz(1:Nx,1:Ny,1:Nz) = myu(1:Nx,1:Ny,1:Nz)*(dudz(1:Nx,1:Ny,1:Nz) + dwdx(1:Nx,1:Ny,1:Nz))

	tau_yy(1:Nx,1:Ny,1:Nz) = myu(1:Nx,1:Ny,1:Nz)*(4.d0/3.d0*dvdy(1:Nx,1:Ny,1:Nz) - 2.d0/3.d0*(dudx(1:Nx,1:Ny,1:Nz) &
								&+ dwdz(1:Nx,1:Ny,1:Nz)))
	tau_yz(1:Nx,1:Ny,1:Nz) = myu(1:Nx,1:Ny,1:Nz)*(dvdz(1:Nx,1:Ny,1:Nz) + dwdy(1:Nx,1:Ny,1:Nz))
	tau_zz(1:Nx,1:Ny,1:Nz) = myu(1:Nx,1:Ny,1:Nz)*(4.d0/3.d0*dwdz(1:Nx,1:Ny,1:Nz) - 2.d0/3.d0*(dudx(1:Nx,1:Ny,1:Nz) &
								&+ dvdy(1:Nx,1:Ny,1:Nz)))
	!!第二不変量
	Q_ten(1:Nx,1:Ny,1:Nz) = dudx(1:Nx,1:Ny,1:Nz)*dvdy(1:Nx,1:Ny,1:Nz) + dvdy(1:Nx,1:Ny,1:Nz)*dwdz(1:Nx,1:Ny,1:Nz) + dwdz(1:Nx,1:Ny,1:Nz)*dudx(1:Nx,1:Ny,1:Nz) &
						&- dudz(1:Nx,1:Ny,1:Nz)*dwdx(1:Nx,1:Ny,1:Nz) - dwdy(1:Nx,1:Ny,1:Nz)*dvdz(1:Nx,1:Ny,1:Nz) + dvdx(1:Nx,1:Ny,1:Nz)*dudy(1:Nx,1:Ny,1:Nz)
	!$OMP end parallel workshare
!!!-----Calc tau and temp. (2nd derivative) :: in : u,v,myu,TT
!
	call CCS_3D(tau_xx,tau_xx_dx,1)
	call CCS_3D(tau_xy,tau_yx_dx,1)
	call CCS_3D(tau_xz,tau_zx_dx,1)
	call CCS_3D(tau_xy,tau_xy_dy,2)
	call CCS_3D(tau_yy,tau_yy_dy,2)
	call CCS_3D(tau_yz,tau_zy_dy,2)
	call CCS_3D(tau_xz,tau_xz_dz,3)
	call CCS_3D(tau_yz,tau_yz_dz,3)
	call CCS_3D(tau_zz,tau_zz_dz,3)
!
	call CCS_3D(TT,TTx,1)
	call CCS_3D(TT,TTy,2)
	call CCS_3D(TT,TTz,3)
!--------------------------------------------------
	!$OMP parallel workshare
	quex(1:Nx,1:Ny,1:Nz) = - ( lambda(1:Nx,1:Ny,1:Nz)*TTx(1:Nx,1:Ny,1:Nz) )
	quey(1:Nx,1:Ny,1:Nz) = - ( lambda(1:Nx,1:Ny,1:Nz)*TTy(1:Nx,1:Ny,1:Nz) )
	quez(1:Nx,1:Ny,1:Nz) = - ( lambda(1:Nx,1:Ny,1:Nz)*TTz(1:Nx,1:Ny,1:Nz) )
	!$OMP end parallel workshare
	call CCS_3D(quex,quex_dx,1)
	call CCS_3D(quey,quey_dy,2)
	call CCS_3D(quez,quez_dz,3)
!kyoukai-----------------------
!x
	!$OMP parallel workshare
	tau_yx_dx(Nx,1:Ny,1:Nz)=0.d0!!!!
	tau_zx_dx(Nx,1:Ny,1:Nz)=0.d0!!!!
	quex_dx(Nx,1:Ny,1:Nz)  =0.d0
!y
	tau_xy_dy(1:Nx,Ny,1:Nz)=0.d0!!!!
	tau_zy_dy(1:Nx,Ny,1:Nz)=0.d0!!!!
	quey_dy(1:Nx,Ny,1:Nz)  =0.d0
	!$OMP end parallel workshare
!--------------------------------------------------
	!$OMP parallel workshare
	!!visco term
	dFvx_2(1:Nx,1:Ny,1:Nz) = tau_xx_dx(1:Nx,1:Ny,1:Nz) /Re
	dFvx_3(1:Nx,1:Ny,1:Nz) = tau_yx_dx(1:Nx,1:Ny,1:Nz) /Re
	dFvx_4(1:Nx,1:Ny,1:Nz) = tau_zx_dx(1:Nx,1:Ny,1:Nz) /Re
	dFvx_5(1:Nx,1:Ny,1:Nz) = ( dudx(1:Nx,1:Ny,1:Nz) * tau_xx(1:Nx,1:Ny,1:Nz) + Q2(1:Nx,1:Ny,1:Nz) &
					&/ Q1(1:Nx,1:Ny,1:Nz) * tau_xx_dx(1:Nx,1:Ny,1:Nz) + dvdx(1:Nx,1:Ny,1:Nz) * tau_xy(1:Nx,1:Ny,1:Nz) &
					&+ Q3(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz) * tau_yx_dx(1:Nx,1:Ny,1:Nz) + dwdx(1:Nx,1:Ny,1:Nz) &
					&* tau_xz(1:Nx,1:Ny,1:Nz) + Q4(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz) * tau_zx_dx(1:Nx,1:Ny,1:Nz) &
					& - quex_dx(1:Nx,1:Ny,1:Nz) )/Re

	dFvy_2(1:Nx,1:Ny,1:Nz) = tau_xy_dy(1:Nx,1:Ny,1:Nz) /Re
	dFvy_3(1:Nx,1:Ny,1:Nz) = tau_yy_dy(1:Nx,1:Ny,1:Nz) /Re
	dFvy_4(1:Nx,1:Ny,1:Nz) = tau_zy_dy(1:Nx,1:Ny,1:Nz) /Re
	dFvy_5(1:Nx,1:Ny,1:Nz) = ( dudy(1:Nx,1:Ny,1:Nz) * tau_xy(1:Nx,1:Ny,1:Nz) + Q2(1:Nx,1:Ny,1:Nz) &
					&/ Q1(1:Nx,1:Ny,1:Nz) * tau_xy_dy(1:Nx,1:Ny,1:Nz) + dvdy(1:Nx,1:Ny,1:Nz) * tau_yy(1:Nx,1:Ny,1:Nz) &
					&+ Q3(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz) * tau_yy_dy(1:Nx,1:Ny,1:Nz) + dwdy(1:Nx,1:Ny,1:Nz) &
					&* tau_yz(1:Nx,1:Ny,1:Nz) + Q4(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz) * tau_zy_dy(1:Nx,1:Ny,1:Nz) &
					& - quey_dy(1:Nx,1:Ny,1:Nz) )/Re

	dFvz_2(1:Nx,1:Ny,1:Nz) = tau_xz_dz(1:Nx,1:Ny,1:Nz) /Re
	dFvz_3(1:Nx,1:Ny,1:Nz) = tau_yz_dz(1:Nx,1:Ny,1:Nz) /Re
	dFvz_4(1:Nx,1:Ny,1:Nz) = tau_zz_dz(1:Nx,1:Ny,1:Nz) /Re
	dFvz_5(1:Nx,1:Ny,1:Nz) = ( dudz(1:Nx,1:Ny,1:Nz) * tau_xz(1:Nx,1:Ny,1:Nz) + Q2(1:Nx,1:Ny,1:Nz) &
					&/ Q1(1:Nx,1:Ny,1:Nz) * tau_xz_dz(1:Nx,1:Ny,1:Nz) + dvdz(1:Nx,1:Ny,1:Nz) * tau_yz(1:Nx,1:Ny,1:Nz) &
					&+ Q3(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz) * tau_yz_dz(1:Nx,1:Ny,1:Nz) + dwdz(1:Nx,1:Ny,1:Nz) &
					&* tau_zz(1:Nx,1:Ny,1:Nz) + Q4(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz) * tau_zz_dz(1:Nx,1:Ny,1:Nz) &
					& - quez_dz(1:Nx,1:Ny,1:Nz) )/Re
	!$OMP end parallel workshare
	end subroutine visco_bibun
!=====================================================================
!=====================================================================
!!!------matome_iryu kou
	subroutine bibun_adv()
	double precision :: Fx1_2(1:Nx,1:Ny,1:Nz), Fx1_5(1:Nx,1:Ny,1:Nz), Fx2_2(1:Nx,1:Ny,1:Nz), Fx2_5(1:Nx,1:Ny,1:Nz)
	double precision :: dFx1_2(1:Nx,1:Ny,1:Nz), dFx1_5(1:Nx,1:Ny,1:Nz), dFx2_2(1:Nx,1:Ny,1:Nz), dFx2_5(1:Nx,1:Ny,1:Nz)
	double precision :: Fy1_3(1:Nx,1:Ny,1:Nz), Fy1_5(1:Nx,1:Ny,1:Nz), Fy2_3(1:Nx,1:Ny,1:Nz), Fy2_5(1:Nx,1:Ny,1:Nz)
	double precision :: dFy1_3(1:Nx,1:Ny,1:Nz), dFy1_5(1:Nx,1:Ny,1:Nz), dFy2_3(1:Nx,1:Ny,1:Nz), dFy2_5(1:Nx,1:Ny,1:Nz)
	double precision :: Fz1_4(1:Nx,1:Ny,1:Nz), Fz1_5(1:Nx,1:Ny,1:Nz), Fz2_4(1:Nx,1:Ny,1:Nz), Fz2_5(1:Nx,1:Ny,1:Nz)
	double precision :: dFz1_4(1:Nx,1:Ny,1:Nz), dFz1_5(1:Nx,1:Ny,1:Nz), dFz2_4(1:Nx,1:Ny,1:Nz), dFz2_5(1:Nx,1:Ny,1:Nz)
!

	!$OMP parallel workshare
	!!移流項の群をまとめて微分or項ごとに微分
	!!Fn_renew
	Fx_1(1:Nx,1:Ny,1:Nz) = Q2(1:Nx,1:Ny,1:Nz)
!	Fx_2(1:Nx,1:Ny,1:Nz) = Q2(1:Nx,1:Ny,1:Nz) * Q2(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz) + p(1:Nx,1:Ny,1:Nz)!!
	Fx_3(1:Nx,1:Ny,1:Nz) = Q2(1:Nx,1:Ny,1:Nz) * Q3(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)
	Fx_4(1:Nx,1:Ny,1:Nz) = Q2(1:Nx,1:Ny,1:Nz) * Q4(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)
!	Fx_5(1:Nx,1:Ny,1:Nz) = (Q5(1:Nx,1:Ny,1:Nz) + p(1:Nx,1:Ny,1:Nz)) * Q2(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)!!
	Fy_1(1:Nx,1:Ny,1:Nz) = Q3(1:Nx,1:Ny,1:Nz)
	Fy_2(1:Nx,1:Ny,1:Nz) = Q3(1:Nx,1:Ny,1:Nz) * Q2(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)
!	Fy_3(1:Nx,1:Ny,1:Nz) = Q3(1:Nx,1:Ny,1:Nz) * Q3(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz) + p(1:Nx,1:Ny,1:Nz)!!
	Fy_4(1:Nx,1:Ny,1:Nz) = Q3(1:Nx,1:Ny,1:Nz) * Q4(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)
!	Fy_5(1:Nx,1:Ny,1:Nz) = (Q5(1:Nx,1:Ny,1:Nz) + p(1:Nx,1:Ny,1:Nz)) * Q3(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)!!
	Fz_1(1:Nx,1:Ny,1:Nz) = Q4(1:Nx,1:Ny,1:Nz)
	Fz_2(1:Nx,1:Ny,1:Nz) = Q4(1:Nx,1:Ny,1:Nz) * Q2(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)
	Fz_3(1:Nx,1:Ny,1:Nz) = Q4(1:Nx,1:Ny,1:Nz) * Q3(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)
!	Fz_4(1:Nx,1:Ny,1:Nz) = Q4(1:Nx,1:Ny,1:Nz) * Q4(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz) + p(1:Nx,1:Ny,1:Nz)!!
!	Fz_5(1:Nx,1:Ny,1:Nz) = (Q5(1:Nx,1:Ny,1:Nz) + p(1:Nx,1:Ny,1:Nz)) * Q4(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)!!
	!$OMP end parallel workshare
	call bunkatu(Fx_1,Q1,dFx_1,1)
!	call bunkatu(Fx_2,Q2,dFx_2,1)
	call bunkatu(Fx_3,Q3,dFx_3,1)
	call bunkatu(Fx_4,Q4,dFx_4,1)
!	call bunkatu(Fx_5,Q5,dFx_5,1)
	call bunkatu(Fy_1,Q1,dFy_1,2)
	call bunkatu(Fy_2,Q2,dFy_2,2)
!	call bunkatu(Fy_3,Q3,dFy_3,2)
	call bunkatu(Fy_4,Q4,dFy_4,2)
!	call bunkatu(Fy_5,Q5,dFy_5,2)
	call bunkatu(Fz_1,Q1,dFz_1,3)
	call bunkatu(Fz_2,Q2,dFz_2,3)
	call bunkatu(Fz_3,Q3,dFz_3,3)
!	call bunkatu(Fz_4,Q4,dFz_4,3)
!	call bunkatu(Fz_5,Q5,dFz_5,3)
!!!----------------------------------
!!pressure - advection
	!$OMP parallel workshare
	Fx1_2(1:Nx,1:Ny,1:Nz) = Q2(1:Nx,1:Ny,1:Nz) * Q2(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)
	Fx2_2(1:Nx,1:Ny,1:Nz) = p(1:Nx,1:Ny,1:Nz)!!
	Fx1_5(1:Nx,1:Ny,1:Nz) = Q5(1:Nx,1:Ny,1:Nz) * Q2(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)!!
	Fx2_5(1:Nx,1:Ny,1:Nz) = p(1:Nx,1:Ny,1:Nz)  * Q2(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)!!
	Fy1_3(1:Nx,1:Ny,1:Nz) = Q3(1:Nx,1:Ny,1:Nz) * Q3(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)
	Fy2_3(1:Nx,1:Ny,1:Nz) = p(1:Nx,1:Ny,1:Nz)!!
	Fy1_5(1:Nx,1:Ny,1:Nz) = Q5(1:Nx,1:Ny,1:Nz) * Q3(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)!!
	Fy2_5(1:Nx,1:Ny,1:Nz) = p(1:Nx,1:Ny,1:Nz)  * Q3(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)!!
	Fz1_4(1:Nx,1:Ny,1:Nz) = Q4(1:Nx,1:Ny,1:Nz) * Q4(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)
	Fz2_4(1:Nx,1:Ny,1:Nz) = p(1:Nx,1:Ny,1:Nz)!!
	Fz1_5(1:Nx,1:Ny,1:Nz) = Q5(1:Nx,1:Ny,1:Nz) * Q4(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)!!
	Fz2_5(1:Nx,1:Ny,1:Nz) = p(1:Nx,1:Ny,1:Nz)  * Q4(1:Nx,1:Ny,1:Nz) / Q1(1:Nx,1:Ny,1:Nz)!!
	!$OMP end parallel workshare
	call bunkatu(Fx1_2,Q2,dFx1_2,1)
	call bunkatu(Fx2_2,Q2,dFx2_2,1)
	call bunkatu(Fx1_5,Q5,dFx1_5,1)
	call bunkatu(Fx2_5,Q5,dFx2_5,1)
	call bunkatu(Fy1_3,Q3,dFy1_3,2)
	call bunkatu(Fy2_3,Q3,dFy2_3,2)
	call bunkatu(Fy1_5,Q5,dFy1_5,2)
	call bunkatu(Fy2_5,Q5,dFy2_5,2)
	call bunkatu(Fz1_4,Q4,dFz1_4,3)
	call bunkatu(Fz2_4,Q4,dFz2_4,3)
	call bunkatu(Fz1_5,Q5,dFz1_5,3)
	call bunkatu(Fz2_5,Q5,dFz2_5,3)
!!!------------------------------
!	!$OMP parallel workshare
	dFx_2(1:Nx,1:Ny,1:Nz) = dFx1_2(1:Nx,1:Ny,1:Nz) + dFx2_2(1:Nx,1:Ny,1:Nz)
	dFx_5(1:Nx,1:Ny,1:Nz) = dFx1_5(1:Nx,1:Ny,1:Nz) + dFx2_5(1:Nx,1:Ny,1:Nz)
	dFy_3(1:Nx,1:Ny,1:Nz) = dFy1_3(1:Nx,1:Ny,1:Nz) + dFy2_3(1:Nx,1:Ny,1:Nz)
	dFy_5(1:Nx,1:Ny,1:Nz) = dFy1_5(1:Nx,1:Ny,1:Nz) + dFy2_5(1:Nx,1:Ny,1:Nz)
	dFz_4(1:Nx,1:Ny,1:Nz) = dFz1_4(1:Nx,1:Ny,1:Nz) + dFz2_4(1:Nx,1:Ny,1:Nz)
	dFz_5(1:Nx,1:Ny,1:Nz) = dFz1_5(1:Nx,1:Ny,1:Nz) + dFz2_5(1:Nx,1:Ny,1:Nz)
!	!$OMP end parallel workshare

!!!!==============================================================================		Q1_0(i,j,k) = Q1(i,j,k)
!!!Calclation NSCBC
	call NSCBC()

	end subroutine bibun_adv
!=====================================================================
!=====================================================================
!!!-----8:NSCBC------------------------------------------------------
	subroutine NSCBC()
	double precision :: chvel1x_in(1:Ny,1:Nz), chvel2x_in(1:Ny,1:Nz), chvel3x_in(1:Ny,1:Nz)
	double precision :: chvel4x_in(1:Ny,1:Nz), chvel5x_in(1:Ny,1:Nz)
	double precision :: chvel1x_out(1:Ny,1:Nz), chvel2x_out(1:Ny,1:Nz), chvel3x_out(1:Ny,1:Nz)
	double precision :: chvel4x_out(1:Ny,1:Nz), chvel5x_out(1:Ny,1:Nz)
	double precision :: chvel1y_in(1:Nx,1:Nz), chvel2y_in(1:Nx,1:Nz), chvel3y_in(1:Nx,1:Nz)
	double precision :: chvel4y_in(1:Nx,1:Nz), chvel5y_in(1:Nx,1:Nz)
	double precision :: chvel1y_out(1:Nx,1:Nz), chvel2y_out(1:Nx,1:Nz), chvel3y_out(1:Nx,1:Nz)
	double precision :: chvel4y_out(1:Nx,1:Nz), chvel5y_out(1:Nx,1:Nz)
	double precision :: L1x_in(1:Ny,1:Nz), L2x_in(1:Ny,1:Nz), L3x_in(1:Ny,1:Nz), L4x_in(1:Ny,1:Nz), L5x_in(1:Ny,1:Nz)
	double precision :: L1x_out(1:Ny,1:Nz), L2x_out(1:Ny,1:Nz), L3x_out(1:Ny,1:Nz), L4x_out(1:Ny,1:Nz), L5x_out(1:Ny,1:Nz)
	double precision :: L1y_in(1:Nx,1:Nz), L2y_in(1:Nx,1:Nz), L3y_in(1:Nx,1:Nz), L4y_in(1:Nx,1:Nz), L5y_in(1:Nx,1:Nz)
	double precision :: L1y_out(1:Nx,1:Nz), L2y_out(1:Nx,1:Nz), L3y_out(1:Nx,1:Nz), L4y_out(1:Nx,1:Nz), L5y_out(1:Nx,1:Nz)
	double precision :: d1x_in(1:Ny,1:Nz), d2x_in(1:Ny,1:Nz), d3x_in(1:Ny,1:Nz), d4x_in(1:Ny,1:Nz), d5x_in(1:Ny,1:Nz)
	double precision :: d1x_out(1:Ny,1:Nz), d2x_out(1:Ny,1:Nz), d3x_out(1:Ny,1:Nz), d4x_out(1:Ny,1:Nz), d5x_out(1:Ny,1:Nz)
	double precision :: d1y_in(1:Nx,1:Nz), d2y_in(1:Nx,1:Nz), d3y_in(1:Nx,1:Nz), d4y_in(1:Nx,1:Nz), d5y_in(1:Nx,1:Nz)
	double precision :: d1y_out(1:Nx,1:Nz), d2y_out(1:Nx,1:Nz), d3y_out(1:Nx,1:Nz), d4y_out(1:Nx,1:Nz), d5y_out(1:Nx,1:Nz)
!!!
	double precision :: dpdx_in(1:Ny,1:Nz), dudx_in(1:Ny,1:Nz), drhodx_in(1:Ny,1:Nz), dvdx_in(1:Ny,1:Nz), dwdx_in(1:Ny,1:Nz)
	double precision :: dpdx_out(1:Ny,1:Nz),dudx_out(1:Ny,1:Nz),drhodx_out(1:Ny,1:Nz),dvdx_out(1:Ny,1:Nz),dwdx_out(1:Ny,1:Nz)
	double precision :: dpdy_in(1:Nx,1:Nz), dudy_in(1:Nx,1:Nz), drhody_in(1:Nx,1:Nz), dvdy_in(1:Nx,1:Nz), dwdy_in(1:Nx,1:Nz)
	double precision :: dpdy_out(1:Nx,1:Nz),dudy_out(1:Nx,1:Nz),drhody_out(1:Nx,1:Nz),dvdy_out(1:Nx,1:Nz),dwdy_out(1:Nx,1:Nz)
	double precision :: xx_x(1:Nx),dxx_dx(1:Nx), yy_y(1:Ny), dyy_dy(1:Ny), vel_sound(1:Nx,1:Ny,1:Nz), vel_sound2(1:Nx,1:Ny,1:Nz)
	double precision :: vel_sq(1:Nx,1:Ny,1:Nz), vel_sq2(1:Nx,1:Ny,1:Nz)
	double precision :: xx_xp(1:Nx),dxx_dxp(1:Nx), yy_yp(1:Ny), dyy_dyp(1:Ny), xx_xm(1:Nx),dxx_dxm(1:Nx), yy_ym(1:Ny), dyy_dym(1:Ny)

!call bibun_DCS(1, 2, Nx, dx, sgm, u_in, u_out)!CCS,katagawa,Nx,dx,sgm=0.d0,in,out
	call update()
	!$OMP parallel workshare
	vel_sound2(1:Nx,1:Ny,1:Nz)=dabs( gamma*(gamma - 1.d0)*( Q5(1:Nx,1:Ny,1:Nz) - 0.5d0*( Q2(1:Nx,1:Ny,1:Nz)*Q2(1:Nx,1:Ny,1:Nz) &
			+ Q3(1:Nx,1:Ny,1:Nz)*Q3(1:Nx,1:Ny,1:Nz) + Q4(1:Nx,1:Ny,1:Nz)*Q4(1:Nx,1:Ny,1:Nz) )/Q1(1:Nx,1:Ny,1:Nz) )/Q1(1:Nx,1:Ny,1:Nz) )
	!$OMP end parallel workshare
	!$OMP parallel workshare
	vel_sound(1:Nx,1:Ny,1:Nz)=dsqrt( vel_sound2(1:Nx,1:Ny,1:Nz) )
	!$OMP end parallel workshare
!do k=1,Nz
!	do j=1,Ny
!		do i=1,Nx
!			if (isnan(vel_sound(i,j,k))) then
!			write(*,"(4(1x,g15.7))") 'in NSCBC -  vel_sound NaN ',i, j, k
!			call output_files()
!			call finish()
!			stop
!			end if
!		end do
!	end do
!end do

!!!!
!!必要な成分をここで微分
!!rho_x
	!$OMP parallel do
do k=1,Nz
	do j=1,Ny
!		do i=1,Nx
!			xx_x(i)=Q1(i,j,k)
!		end  do
!$omp critical
			call bibun_DCS(1, 2, Nx, dx, 0.d0, L_xdp, U_xdp, Q1(:,j,k), dxx_dx,1)!!!DCS (←CCS)2016.07.06修正
		drhodx_in(j,k)=dxx_dx(1)
		drhodx_out(j,k)=dxx_dx(Nx)
!$omp end critical
!			if (isnan(dxx_dx(1)) .or. isnan(dxx_dx(Nx))) then
!			write(*,"(3(1x,g15.7))") 'drhodx NaN error',j,k
!			call output_files()
!			call finish()
!			stop
!			end if
!	end do
!end do
!!p_x
!do k=1,Nz
!	do j=1,Ny
!$omp critical
		do i=1,Nx
			xx_x(i)= (gamma - 1.d0)*( Q5(i,j,k) - 0.5d0*( Q2(i,j,k)*Q2(i,j,k) + Q3(i,j,k)*Q3(i,j,k) + Q4(i,j,k)*Q4(i,j,k))/Q1(i,j,k) )
		end  do
		call bibun_DCS(1, 2, Nx, dx, 0.d0, L_xdp, U_xdp, xx_x, dxx_dx,1)
!		write(*,"(4(1x,g15.7))") dxx_dxp(1), dxx_dxm(1), dpdx_in(j,k), dxx_dxp(Nx), dxx_dxm(Nx), dpdx_out(j,k)

		dpdx_in(j,k)=dxx_dx(1)
		dpdx_out(j,k)=dxx_dx(Nx)
!$omp end critical
!			if (isnan(dpdx_out(j,k))) then
!			write(*,"(3(1x,g15.7))") 'dpdx_out NaN error',j, k
!			call output_files()
!			call finish()
!			stop
!			end if
!			if (isnan(dpdx_in(j,k))) then
!			write(*,"(3(1x,g15.7))") 'dpdx_in NaN error', j, k
!			call output_files()
!			call finish()
!			stop
!			end if
!	end do
!end do
!!u_x
!do k=1,Nz
!	do j=1,Ny
!		do i=1,Nx
!			xx_x(i)=Q2(i,j,k)/Q1(i,j,k)
!		end  do
!$omp critical
		call bibun_DCS(1, 2, Nx, dx, 0.d0, L_xdp, U_xdp, Q2(:,j,k)/Q1(:,j,k), dxx_dx,1)
		dudx_in(j,k)=dxx_dx(1)
		dudx_out(j,k)=dxx_dx(Nx)
!$omp end critical
!			if (isnan(dxx_dx(1)) .or. isnan(dxx_dx(Nx))) then
!			write(*,"(3(1x,g15.7))") 'dudx NaN error', j, k
!			call output_files()
!			call finish()
!			stop
!			end if

!	end do
!end do
!!v_x
!do k=1,Nz
!	do j=1,Ny
!		do i=1,Nx
!			xx_x(i)=Q3(i,j,k)/Q1(i,j,k)
!		end  do
!$omp critical
		call bibun_DCS(1, 2, Nx, dx, 0.d0, L_xdp, U_xdp, Q3(:,j,k)/Q1(:,j,k), dxx_dx,1)
		dvdx_in(j,k)=dxx_dx(1)
		dvdx_out(j,k)=dxx_dx(Nx)
!$omp end critical
!			if (isnan(dxx_dx(1)) .or. isnan(dxx_dx(Nx))) then
!			write(*,"(3(1x,g15.7))") 'dvdx NaN error', j, k
!			call output_files()
!			call finish()
!			stop
!			end if

!	end do
!end do
!!w_x
!do k=1,Nz
!	do j=1,Ny
!		do i=1,Nx
!			xx_x(i)=Q4(i,j,k)/Q1(i,j,k)
!		end  do
!$omp critical
		call bibun_DCS(1, 2, Nx, dx, 0.d0, L_xdp, U_xdp, Q4(:,j,k)/Q1(:,j,k), dxx_dx,1)
		dwdx_in(j,k) =dxx_dx(1)
		dwdx_out(j,k)=dxx_dx(Nx)
!$omp end critical
!			if (isnan(dxx_dx(1)) .or. isnan(dxx_dx(Nx))) then
!			write(*,"(3(1x,g15.7))") 'dwdx NaN error', j, k
!			call output_files()
!			call finish()
!			stop
!			end if

	end do
end do
!$OMP end parallel do
!	write(*,*) 'NSCBC'
!	call NaN_check2()
!----y bibun-----
!!rho_y
	!$OMP parallel do
do k=1,Nz
	do i=1,Nx
!		do j=1,Ny
!			yy_y(j)=Q1(i,j,k)
!		end  do
!$omp critical
		call bibun_DCS(1, 2, Ny, dy, 0.d0, L_ydp, U_ydp, Q1(i,:,k), dyy_dy,1)
		drhody_in(i,k)=dyy_dy(1)
		drhody_out(i,k)=dyy_dy(Ny)
!$omp end critical
!			if (isnan(dyy_dy(1)) .or. isnan(dyy_dy(Ny))) then
!			write(*,"(3(1x,g15.7))") 'drhody NaN error',i, k
!			call output_files()
!			call finish()
!			stop
!			end if
!	end do
!end do
!!p_y
!do k=1,Nz
!	do i=1,Nx
!$omp critical
		do j=1,Ny
			yy_y(j)= (gamma - 1.d0)*( Q5(i,j,k) - 0.5d0/Q1(i,j,k)*( Q2(i,j,k)*Q2(i,j,k) + Q3(i,j,k)*Q3(i,j,k) + Q4(i,j,k)*Q4(i,j,k)) )
		end  do
		call bibun_DCS(1, 2, Ny, dy, 0.d0, L_ydp, U_ydp, yy_y, dyy_dy,1)
		dpdy_in(i,k)=dyy_dy(1)
		dpdy_out(i,k)=dyy_dy(Ny)
!$omp end critical
!			if (isnan(dyy_dy(1)) .or. isnan(dyy_dy(Ny))) then
!			write(*,"(3(1x,g15.7))") 'dpdy NaN error',i, k
!			call output_files()
!			call finish()
!			stop
!			end if
!	end do
!end do
!!u_x
!do k=1,Nz
!	do i=1,Nx
!		do j=1,Ny
!			yy_y(j)=Q2(i,j,k)/Q1(i,j,k)
!		end  do
!$omp critical
		call bibun_DCS(1, 2, Ny, dy, 0.d0, L_ydp, U_ydp, Q2(i,:,k)/Q1(i,:,k), dyy_dy,1)
		dudy_in(i,k)=dyy_dy(1)
		dudy_out(i,k)=dyy_dy(Ny)
!$omp end critical
!			if (isnan(dyy_dy(1)) .or. isnan(dyy_dy(Ny))) then
!			write(*,"(3(1x,g15.7))") 'dudy NaN error',i, k
!			call output_files()
!			call finish()
!			stop
!			end if
!	end do
!end do
!!v_x
!do k=1,Nz
!	do i=1,Nx
!		do j=1,Ny
!			yy_y(j)=Q3(i,j,k)/Q1(i,j,k)
!		end  do
!$omp critical
		call bibun_DCS(1, 2, Ny, dy, 0.d0, L_ydp, U_ydp, Q3(i,:,k)/Q1(i,:,k), dyy_dy,1)
		dvdy_in(i,k)=dyy_dy(1)
		dvdy_out(i,k)=dyy_dy(Ny)
!$omp end critical
!			if (isnan(dyy_dy(1)) .or. isnan(dyy_dy(Ny))) then
!			write(*,"(3(1x,g15.7))") 'dvdy NaN error',i, k
!			call output_files()
!			call finish()
!			stop
!			end if
!	end do
!end do
!!w_x
!do k=1,Nz
!	do i=1,Nx
!		do j=1,Ny
!			yy_y(j)=Q4(i,j,k)/Q1(i,j,k)
!		end  do
!$omp critical
		call bibun_DCS(1, 2, Ny, dy, 0.d0, L_ydp, U_ydp, Q4(i,:,k)/Q1(i,:,k), dyy_dy,1)
		dwdy_in(i,k)=dyy_dy(1)
		dwdy_out(i,k)=dyy_dy(Ny)
!$omp end critical
!			if (isnan(dyy_dy(1)) .or. isnan(dyy_dy(Ny))) then
!			write(*,"(3(1x,g15.7))") 'dwdy NaN error',i, k
!			call output_files()
!			call finish()
!			stop
!			end if
	end do
end do
	!$OMP end parallel do
!!!!

!x----------------------------------------------------------------------------------------------------------------------------
	i=1!!!!x=0流入境界----------------------------
	!$omp parallel do
	do k=1,Nz
		do j=1,Ny
			chvel1x_in(j,k) = Q2(1,j,k)/Q1(1,j,k) - vel_sound(1,j,k)
		end do
	end do
	!$omp end parallel do

		!amplitudes of c:haracteristic waves

	!$omp parallel do
	do k=1,Nz
		do j=2,Ny
		!!流入
		L1x_in(j,k) = chvel1x_in(j,k) * ( dpdx_in(j,k) - Q1(1,j,k) * vel_sound(1,j,k) * dudx(1,j,k) )
		L3x_in(j,k) = 0.d0!chvel3x_in(j,k) * dvdx(1,j,k)
		L4x_in(j,k) = 0.d0!chvel4x_in(j,k) * dwdx(1,j,k)
!!流入あり----------------------------------------------------------------------------------------------------------------------------
!!!流入部u時間変化する場合
		L5x_in(j,k) = chvel1x_in(j,k) * ( dpdx_in(j,k) - Q1(1,j,k) * vel_sound(1,j,k) * dudx(1,j,k) )  -2.d0 * Q1(1,j,k) * vel_sound(1,j,k) * amp_ufunc(j) * 2d0*pi/omegar/dble(lt) * &
				&dexp(2d0*pi/dble(lt)*omegar*t) *( omegar * dcos( 2d0*pi/omegar/dble(lt)*(z(k)*beta - t) ) + omegai * dcos( 2d0*pi/omegar/dble(lt)*(z(k)*beta - t) ) )
		L2x_in(j,k) = (gamma-1.d0) * ( chvel1x_in(j,k) * ( dpdx_in(j,k) - Q1(1,j,k) * vel_sound(1,j,k) * dudx(1,j,k) )  &
				& - Q1(1,j,k) * vel_sound(1,j,k) * amp_ufunc(j) * 2d0*pi/omegar/dble(lt) * &
				&dexp(2d0*pi/dble(lt)*omegar*t) *( omegar * dcos( 2d0*pi/omegar/dble(lt)*(z(k)*beta - t) ) + omegai * dcos( 2d0*pi/omegar/dble(lt)*(z(k)*beta - t) ) ) )&
				&+ Q1(1,j,k)*Q2(1,j,k)*vel_sound(1,j,k)*vel_sound(1,j,k)*gamma*(gamma-1.d0)*Ma*Ma&
				&/(gamma*Ma*Ma)/Q1(1,j,k)*(gamma - 1.d0)*( Q5(1,j,k) - 0.5d0*Q1(1,j,k)*( u_kaku(j,k)*u_kaku(j,k) + v_kaku(j,k)*v_kaku(j,k) ))&
				&*amp_ufunc(j) * 2d0*pi/omegar/dble(lt) * dexp(2d0*pi/dble(lt)*omegar*t) *( omegar * dcos( 2d0*pi/omegar/dble(lt)*(z(k)*beta - t) ) &
				&+ omegai * dcos( 2d0*pi/omegar/dble(lt)*(z(k)*beta - t) ) )
!!流入なし----------------------------------------------------------------------------------------------------------------------------
!		L2x_in(j,k) = chvel1x_in(j,k) * ( dpdx_in(j,k) - Q1(1,j,k) * vel_sound(1,j,k) * dudx(1,j,k) )*(gamma -1.d0)!chvel2x_in(j,k) * ( vel_sound(1,j,k) *vel_sound(1,j,k) * drhodx_in(j,k) - dpdx_in(j,k) )
!		L5x_in(j,k) = chvel1x_in(j,k) * ( dpdx_in(j,k) - Q1(1,j,k) * vel_sound(1,j,k) * dudx(1,j,k) )!chvel5x_in(j,k) * ( dpdx_in(j,k) + Q1(1,j,k) * vel_sound(1,j,k) * dudx(1,j,k) )
		end do
	end do
	!$omp end parallel do

		!vector d

	!$omp parallel do

	do k=1,Nz
		do j=1,Ny
		d1x_in(j,k) = gamma * L1x_in(j,k) /( vel_sound2(1,j,k) )!( L2x_in(j,k) + L1x_in(j,k) ) /( vel_sound2(1,j,k) )
!		d2x_in(j,k) = 0.d0!L1x_in(j,k)	!!!LODIから修正(NS，nergyは不要，連続のみ計算)
!		d3x_in(j,k) = 0.5d0 * ( L5x_in(j,k) - L1x_in(j,k) ) /(Q1(1,j,k) *vel_sound(1,j,k) )
!		d4x_in(j,k) = L3x_in(j,k)
!		d5x_in(j,k) = L4x_in(j,k)
		end do
	end do
	!$omp end parallel do


	!!right side boundary(x=Nx流出境界)-----------------------------
	!$omp parallel do

	do k=1,Nz
		do j=1,Ny
		chvel1x_out(j,k) = Q2(Nx,j,k)/Q1(Nx,j,k) - vel_sound(Nx,j,k) !!new param
		chvel2x_out(j,k) = Q2(Nx,j,k)/Q1(Nx,j,k)
		chvel3x_out(j,k) = Q2(Nx,j,k)/Q1(Nx,j,k)
		chvel4x_out(j,k) = Q2(Nx,j,k)/Q1(Nx,j,k)
		chvel5x_out(j,k) = Q2(Nx,j,k)/Q1(Nx,j,k) + vel_sound(Nx,j,k)
		end do
	end do
	!$omp end parallel do


		!amplitudes of c:haracteristic waves
	!$omp parallel do

	do k=1,Nz
		do j=1,Ny
		L1x_out(j,k) = 0.d0!chvel1x_out(j,k) * ( dpdx_out(j,k) - Q1(Nx,j,k) * vel_sound(Nx,j,k) * dudx(Nx,j,k) )!!!!
!		if( chvel1x_out(j,k) < 0.d0 ) L1x_out(j,k) = 0.d0
		L2x_out(j,k) = chvel2x_out(j,k) * ( vel_sound2(Nx,j,k) * drhodx_out(j,k) - dpdx_out(j,k) )
		if( chvel2x_out(j,k) < 0.d0 ) L2x_out(j,k) = 0.d0
		L3x_out(j,k) = chvel3x_out(j,k) * dvdx(Nx,j,k)
		if( chvel3x_out(j,k) < 0.d0 ) L3x_out(j,k) = 0.d0
		L4x_out(j,k) = chvel4x_out(j,k) * dwdx(Nx,j,k)
		if( chvel4x_out(j,k) < 0.d0 ) L4x_out(j,k) = 0.d0
		L5x_out(j,k) = chvel5x_out(j,k) * ( dpdx_out(j,k) + Q1(Nx,j,k) * vel_sound(Nx,j,k) * dudx(Nx,j,k) )
		if( chvel5x_out(j,k) < 0.d0 ) L5x_out(j,k) = 0.d0
		end do
	end do
	!$omp end parallel do

		!vector d
	!$omp parallel do

	do k=1,Nz
		do j=1,Ny
!		d1x_out(j,k) = ( L2x_out(j,k) + 0.5d0 *( L5x_out(j,k) + L1x_out(j,k) ) ) /( vel_sound(Nx,j,k)*vel_sound(Nx,j,k) )
!		d2x_out(j,k) = 0.5d0 *( L5x_out(j,k) + L1x_out(j,k) )
!		d3x_out(j,k) = 0.5d0 *( L5x_out(j,k) - L1x_out(j,k) ) /(Q1(Nx,j,k) *vel_sound(Nx,j,k) )
		d1x_out(j,k) = ( L2x_out(j,k) + 0.5d0 * L5x_out(j,k) ) /( vel_sound2(Nx,j,k) )
		d2x_out(j,k) = 0.5d0 * L5x_out(j,k)
		d3x_out(j,k) = 0.5d0 * L5x_out(j,k) /(Q1(Nx,j,k) *vel_sound(Nx,j,k) )
		d4x_out(j,k) = L3x_out(j,k)
		d5x_out(j,k) = L4x_out(j,k)
		end do
	end do
	!$omp end parallel do


!x-------------------------------------------

!y-------------------------------------------
		!characteristic velocity
		j=1!!y=0壁側境界
	!$omp parallel do

	do k=1,Nz
		do i=1,Nx
		chvel1y_in(i,k) = Q3(i,1,k)/Q1(i,1,k) - vel_sound(i,1,k) !!new param
!		chvel2y_in(i,k) = Q3(i,1,k)/Q1(i,1,k)
!		chvel3y_in(i,k) = Q3(i,1,k)/Q1(i,1,k)
!		chvel4y_in(i,k) = Q3(i,1,k)/Q1(i,1,k)
!		chvel5y_in(i,k) = Q3(i,1,k)/Q1(i,1,k) + vel_sound(i,1,k)
		end do
	end do
	!$omp end parallel do


		!amplitudes of c:haracteristic waves
	!$omp parallel do

	do k=1,Nz
		do i=1,Nx
		L1y_in(i,k) = chvel1y_in(i,k) * ( dpdy_in(i,k) - Q1(i,1,k) * vel_sound(i,1,k) * dvdy(i,1,k) )
		L2y_in(i,k) = 0.d0!chvel2y_in(i,k) * ( vel_sound(i,1,k) * vel_sound(i,1,k) * drhody_in(i,k) - dpdy_in(i,k) )
		L3y_in(i,k) = 0.d0!chvel3y_in(i,k) * dudy(i,1,k)
		L4y_in(i,k) = 0.d0!chvel4y_in(i,k) * dwdy(i,1,k)
		L5y_in(i,k) = chvel1y_in(i,k) * ( dpdy_in(i,k) - Q1(i,1,k) * vel_sound(i,1,k) * dvdy(i,1,k) )!chvel5y_in(i,k) * ( dpdy_in(i,k) + Q1(i,1,k) * vel_sound(i,1,k) * dvdy(i,1,k) )
		end do
	end do
	!$omp end parallel do

		!vector d
	!$omp parallel do

	do k=1,Nz
		do i=1,Nx
		d1y_in(i,k) = L1y_in(i,k) /( vel_sound2(i,1,k) )
		d2y_in(i,k) = L1y_in(i,k)
		d3y_in(i,k) = 0.d0!L3y_in(i,k)
		d4y_in(i,k) = 0.d0!0.5d0 *( L5y_in(i,k) - L1y_in(i,k) ) /(Q1(i,1,k) *vel_sound(i,1,k) )
		d5y_in(i,k) = 0.d0!L4y_in(i,k)
		end do
	end do
	!$omp end parallel do


	!!higher boundary(y=Ny流出境界)
	j=Ny
	!$omp parallel do

	do k=1,Nz
		do i=1,Nx
		chvel1y_out(i,k) = Q3(i,Ny,k)/Q1(i,Ny,k) - vel_sound(i,Ny,k) !!new param
		chvel2y_out(i,k) = Q3(i,Ny,k)/Q1(i,Ny,k)
		chvel3y_out(i,k) = Q3(i,Ny,k)/Q1(i,Ny,k)
		chvel4y_out(i,k) = Q3(i,Ny,k)/Q1(i,Ny,k)
		chvel5y_out(i,k) = Q3(i,Ny,k)/Q1(i,Ny,k) + vel_sound(i,Ny,k)
		end do
	end do
	!$omp end parallel do

		!amplitudes of c:haracteristic waves
	!$omp parallel do

	do k=1,Nz
		do i=1,Nx
		L1y_out(i,k) = 0.d0!chvel1y_out(i,k) * ( dpdy_out(i,k) - Q1(i,Ny,k) * vel_sound(i,Ny,k) * dvdy(i,Ny,k) )
		L2y_out(i,k) = chvel2y_out(i,k) * ( vel_sound2(i,Ny,k) * drhody_out(i,k) - dpdy_out(i,k) )
		if( chvel2y_out(i,k) < 0.d0 ) L2y_out(i,k) = 0.d0
		L3y_out(i,k) = chvel3y_out(i,k) * dudy(i,Ny,k)
		if( chvel3y_out(i,k) < 0.d0 ) L3y_out(i,k) = 0.d0
		L4y_out(i,k) = chvel4y_out(i,k) * dwdy(i,Ny,k)
		if( chvel4y_out(i,k) < 0.d0 ) L4y_out(i,k) = 0.d0
		L5y_out(i,k) = chvel5y_out(i,k) * ( dpdy_out(i,k) + Q1(i,Ny,k) * vel_sound(i,Ny,k) * dvdy(i,Ny,k) )
		if( chvel5y_out(i,k) < 0.d0 ) L5y_out(i,k) = 0.d0
		end do
	end do
	!$omp end parallel do

		!vector d
	!$omp parallel do

	do k=1,Nz
		do i=1,Nx
		d1y_out(i,k) = ( L2y_out(i,k) + 0.5d0 *( L5y_out(i,k) ) ) /( vel_sound2(i,Ny,k) )
		d2y_out(i,k) = 0.5d0 * L5y_out(i,k)
		d3y_out(i,k) = L3y_out(i,k)
		d4y_out(i,k) = 0.5d0 * L5y_out(i,k) /(Q1(i,Ny,k) *vel_sound(i,Ny,k) )
		d5y_out(i,k) = L4y_out(i,k)
		end do
	end do
	!$omp end parallel do


	!!NSCBC vector!!!====================================================================================================================================

	!X
	!$omp parallel do

	do k=1,Nz!!!!!!!0812shuusei
		do j=1,Ny
		Ex_NSCBCin(1,j,k) = d1x_in(j,k)
		Ex_NSCBCin(2,j,k) = 0.d0!Q2(1,j,k)/Q1(1,j,k) * d1x_in(j,k)! + Q1(1,j,k) * d3x_in(j,k)
		Ex_NSCBCin(3,j,k) = 0.d0!Q3(1,j,k)/Q1(1,j,k) * d1x_in(j,k)! + Q1(1,j,k) * d4x_in(j,k)
		Ex_NSCBCin(4,j,k) = 0.d0!Q4(1,j,k)/Q1(1,j,k) * d1x_in(j,k)! + Q1(1,j,k) * d5x_in(j,k)
		Ex_NSCBCin(5,j,k) = 0.d0!d1x_in(j,k) * ( Q2(1,j,k)/Q1(1,j,k) * Q2(1,j,k)/Q1(1,j,k) + Q3(1,j,k)/Q1(1,j,k) &
							!&* Q3(1,j,k)/Q1(1,j,k) + Q4(1,j,k)/Q1(1,j,k) * Q4(1,j,k)/Q1(1,j,k) )/2.d0 &
							!&+ d2x_in(j,k)/( gamma - 1.d0 )! + Q2(1,j,k) * d3x_in(j,k) + Q3(1,j,k) * d4x_in(j,k) + Q4(1,j,k) * d5x_in(j,k)

		Ex_NSCBCout(1,j,k) = d1x_out(j,k)
		Ex_NSCBCout(2,j,k) = Q2(Nx,j,k)/Q1(Nx,j,k) * d1x_out(j,k) + Q1(Nx,j,k) * d3x_out(j,k)
		Ex_NSCBCout(3,j,k) = Q3(Nx,j,k)/Q1(Nx,j,k) * d1x_out(j,k) + Q1(Nx,j,k) * d4x_out(j,k)
		Ex_NSCBCout(4,j,k) = Q4(Nx,j,k)/Q1(Nx,j,k) * d1x_in(j,k) + Q1(Nx,j,k) * d5x_in(j,k)
		Ex_NSCBCout(5,j,k) = d1x_out(j,k) * ( Q2(Nx,j,k)/Q1(Nx,j,k) * Q2(Nx,j,k)/Q1(Nx,j,k) + Q3(Nx,j,k)/Q1(Nx,j,k) &
							&* Q3(Nx,j,k)/Q1(Nx,j,k) + Q4(Nx,j,k)/Q1(Nx,j,k) * Q4(Nx,j,k)/Q1(Nx,j,k) )/2.d0 &
							&+ d2x_out(j,k)/( gamma - 1.d0 ) + Q2(Nx,j,k) * d3x_out(j,k) + Q3(Nx,j,k) * d4x_out(j,k) + Q4(Nx,j,k) * d5x_out(j,k)
		end do
	end do
	!$omp end parallel do

	!Y
	!$omp parallel do

	do k=1,Nz
		do i=1,Nx
		Ey_NSCBCin(1,i,k) = d1y_in(i,k)
		Ey_NSCBCin(2,i,k) = 0.d0!Q2(i,1,k)/Q1(i,1,k) * d1y_in(i,k)! + Q1(i,1,k) * d3y_in(i,k)
		Ey_NSCBCin(3,i,k) = 0.d0!Q3(i,1,k)/Q1(i,1,k) * d1y_in(i,k)! + Q1(i,1,k) * d4y_in(i,k)
		Ey_NSCBCin(4,i,k) = 0.d0!Q4(i,1,k)/Q1(i,1,k) * d1y_in(i,k)! + Q1(i,1,k) * d5y_in(i,k)
		Ey_NSCBCin(5,i,k) = 0.d0!d2y_in(i,k)/( gamma - 1.d0 )!d1y_in(i,k) * ( Q2(i,1,k)/Q1(i,1,k) * Q2(i,1,k)/Q1(i,1,k) + Q3(i,1,k)/Q1(i,1,k) &
							!&* Q3(i,1,k)/Q1(i,1,k) + Q4(i,1,k)/Q1(i,1,k) * Q4(i,1,k)/Q1(i,1,k) )/2.d0 &
							!&+ d2y_in(i,k)/( gamma - 1.d0 )! + Q2(i,1,k) * d3y_in(i,k) + Q3(i,1,k) * d4y_in(i,k) + Q4(i,1,k) * d5y_in(i,k)

		Ey_NSCBCout(1,i,k) = d1y_out(i,k)
		Ey_NSCBCout(2,i,k) = Q2(i,Ny,k)/Q1(i,Ny,k) * d1y_out(i,k) + Q1(i,Ny,k) * d3y_out(i,k)
		Ey_NSCBCout(3,i,k) = Q3(i,Ny,k)/Q1(i,Ny,k) * d1y_out(i,k) + Q1(i,Ny,k) * d4y_out(i,k)
		Ey_NSCBCout(4,i,k) = Q4(i,Ny,k)/Q1(i,Ny,k) * d1y_out(i,k) + Q1(i,Ny,k) * d5y_out(i,k)
		Ey_NSCBCout(5,i,k) = d1y_out(i,k) * ( Q2(i,Ny,k)/Q1(i,Ny,k) * Q2(i,Ny,k)/Q1(i,Ny,k) + Q3(i,Ny,k)/Q1(i,Ny,k) &
							&* Q3(i,Ny,k)/Q1(i,Ny,k) + Q4(i,Ny,k)/Q1(i,Ny,k) * Q4(i,Ny,k)/Q1(i,Ny,k) )/2.d0 &
							&+ d2y_in(i,k)/( gamma - 1.d0 ) + Q2(i,Ny,k) * d3y_out(i,k) + Q3(i,Ny,k) * d4y_out(i,k) + Q4(i,Ny,k) * d5y_out(i,k)
		end do
	end do
	!$omp end parallel do

!	write(*,*) 'NAN2'
!	call NaN_check2()
	end subroutine NSCBC
!=====================================================================
!=====================================================================
!!!-----9:RK4,euler------------------------------------------------------
	subroutine RK4()
	!Qの値を保存しておく
	!$omp parallel workshare
		Q1_0(1:Nx,1:Ny,1:Nz) = Q1(1:Nx,1:Ny,1:Nz)!!!!2016.07.21 なかったので追加
		Q2_0(1:Nx,1:Ny,1:Nz) = Q2(1:Nx,1:Ny,1:Nz)
		Q3_0(1:Nx,1:Ny,1:Nz) = Q3(1:Nx,1:Ny,1:Nz)
		Q4_0(1:Nx,1:Ny,1:Nz) = Q4(1:Nx,1:Ny,1:Nz)
		Q5_0(1:Nx,1:Ny,1:Nz) = Q5(1:Nx,1:Ny,1:Nz)
	!$omp end parallel workshare
!	!$omp parallel do
!	do k=1,Nz
!	do j=1,Ny
!	do i=1,Nx
!1		Q1_0(i,j,k) = Q1(i,j,k)!!!!2016.07.21 なかったので追加
!		Q2_0(i,j,k) = Q2(i,j,k)
!1		Q3_0(i,j,k) = Q3(i,j,k)
!		Q4_0(i,j,k) = Q4(i,j,k)
!		Q5_0(i,j,k) = Q5(i,j,k)
!	end do
!	end do
!	end do
!	!$omp end parallel do
	call update()
	call visco_bibun()
	call bibun_adv()

	!!!!-----RK step1-----!!!!!
	!$omp parallel do
	do k=1,Nz
	do j=1,Ny
	do i=1,Nx
	if(i==1 .and. j==1) then!(1,1,k)流入＊壁
		Q1_1(1,1,k) = - Ex_NSCBCin(1,1,k) - Ey_NSCBCin(1,1,k) - dFz_1(1,1,k)
		Q2_1(1,1,k) = 0d0!- Ex_NSCBCin(2,1,k) - Ey_NSCBCin(2,1,k) - dFz_2(1,1,k) + dFvx_2(1,1,k) + dFvy_2(1,1,k) + dFvz_2(1,1,k)
		Q3_1(1,1,k) = 0d0!- Ex_NSCBCin(3,1,k) - Ey_NSCBCin(3,1,k) - dFz_3(1,1,k) + dFvx_3(1,1,k) + dFvy_3(1,1,k) + dFvz_3(1,1,k)
		Q4_1(1,1,k) = 0d0!- Ex_NSCBCin(4,1,k) - Ey_NSCBCin(4,1,k) - dFz_4(1,1,k) + dFvx_4(1,1,k) + dFvy_4(1,1,k) + dFvz_4(1,1,k)
		Q5_1(1,1,k) = 0d0!- Ex_NSCBCin(5,1,k) - Ey_NSCBCin(5,1,k) - dFz_5(1,1,k) + dFvx_5(1,1,k) + dFvy_5(1,1,k) + dFvz_5(1,1,k)
	else if(i==1 .and. j==Ny) then!(1,Ny,k)流入＊流出
		Q1_1(1,Ny,k) = - Ex_NSCBCin(1,Ny,k) - Ey_NSCBCout(1,1,k) - dFz_1(1,Ny,k)
		Q2_1(1,Ny,k) = 0d0!- Ex_NSCBCin(2,Ny,k) - Ey_NSCBCout(2,1,k) - dFz_2(1,Ny,k) + dFvx_2(1,Ny,k) + dFvy_2(1,Ny,k) + dFvz_2(1,Ny,k)
		Q3_1(1,Ny,k) = 0d0!- Ex_NSCBCin(3,Ny,k) - Ey_NSCBCout(3,1,k) - dFz_3(1,Ny,k) + dFvx_3(1,Ny,k) + dFvy_3(1,Ny,k) + dFvz_3(1,Ny,k)
		Q4_1(1,Ny,k) = 0d0!- Ex_NSCBCin(4,Ny,k) - Ey_NSCBCout(4,1,k) - dFz_4(1,Ny,k) + dFvx_4(1,Ny,k) + dFvy_4(1,Ny,k) + dFvz_4(1,Ny,k)
		Q5_1(1,Ny,k) = 0d0!- Ex_NSCBCin(5,Ny,k) - Ey_NSCBCout(5,1,k) - dFz_5(1,Ny,k) + dFvx_5(1,Ny,k) + dFvy_5(1,Ny,k) + dFvz_5(1,Ny,k)
	else if(i==1 .and. j/=1 .and. j/=Ny) then!(1,2:Ny-1,k)流入
		Q1_1(1,j,k) = - Ex_NSCBCin(1,j,k) - dFy_1(1,j,k) - dFz_1(1,j,k)
		Q2_1(1,j,k) = 0d0!- Ex_NSCBCin(2,j,k) - dFy_2(1,j,k) - dFz_2(1,j,k) + dFvx_2(1,j,k) + dFvy_2(1,j,k) + dFvz_2(1,j,k)
		Q3_1(1,j,k) = 0d0!- Ex_NSCBCin(3,j,k) - dFy_3(1,j,k) - dFz_3(1,j,k) + dFvx_3(1,j,k) + dFvy_3(1,j,k) + dFvz_3(1,j,k)
		Q4_1(1,j,k) = 0d0!- Ex_NSCBCin(4,j,k) - dFy_4(1,j,k) - dFz_4(1,j,k) + dFvx_4(1,j,k) + dFvy_4(1,j,k) + dFvz_4(1,j,k)
		Q5_1(1,j,k) = 0d0!- Ex_NSCBCin(5,j,k) - dFy_5(1,j,k) - dFz_5(1,j,k) + dFvx_5(1,j,k) + dFvy_5(1,j,k) + dFvz_5(1,j,k)
	else if(i==Nx .and. j==1) then!(Nx,1,k)流出＊壁
		Q1_1(Nx,1,k) = - Ex_NSCBCout(1,1,k) - Ey_NSCBCin(1,Nx,k) - dFz_1(Nx,1,k)
		Q2_1(Nx,1,k) = 0d0!- Ex_NSCBCout(2,1,k) - Ey_NSCBCin(2,Nx,k) - dFz_2(Nx,1,k) + dFvx_2(Nx,1,k) + dFvy_2(Nx,1,k) + dFvz_2(Nx,1,k)
		Q3_1(Nx,1,k) = 0d0!- Ex_NSCBCout(3,1,k) - Ey_NSCBCin(3,Nx,k) - dFz_3(Nx,1,k) + dFvx_3(Nx,1,k) + dFvy_3(Nx,1,k) + dFvz_3(Nx,1,k)
		Q4_1(Nx,1,k) = 0d0!- Ex_NSCBCout(4,1,k) - Ey_NSCBCin(4,Nx,k) - dFz_4(Nx,1,k) + dFvx_4(Nx,1,k) + dFvy_4(Nx,1,k) + dFvz_4(Nx,1,k)
		Q5_1(Nx,1,k) = 0d0!- Ex_NSCBCout(5,1,k) - Ey_NSCBCin(5,Nx,k) - dFz_5(Nx,1,k) + dFvx_5(Nx,1,k) + dFvy_5(Nx,1,k) + dFvz_5(Nx,1,k)
	else if(i==Nx .and. j==Ny) then!(1,1,k)流出＊流出
		Q1_1(Nx,Ny,k) = - Ex_NSCBCout(1,Ny,k) - Ey_NSCBCout(1,Nx,k) - dFz_1(Nx,Ny,k)
		Q2_1(Nx,Ny,k) = - Ex_NSCBCout(2,Ny,k) - Ey_NSCBCout(2,Nx,k) - dFz_2(Nx,Ny,k) + dFvx_2(Nx,Ny,k) + dFvy_2(Nx,Ny,k) + dFvz_2(Nx,Ny,k)
		Q3_1(Nx,Ny,k) = - Ex_NSCBCout(3,Ny,k) - Ey_NSCBCout(3,Nx,k) - dFz_3(Nx,Ny,k) + dFvx_3(Nx,Ny,k) + dFvy_3(Nx,Ny,k) + dFvz_3(Nx,Ny,k)
		Q4_1(Nx,Ny,k) = - Ex_NSCBCout(4,Ny,k) - Ey_NSCBCout(4,Nx,k) - dFz_4(Nx,Ny,k) + dFvx_4(Nx,Ny,k) + dFvy_4(Nx,Ny,k) + dFvz_4(Nx,Ny,k)
		Q5_1(Nx,Ny,k) = - Ex_NSCBCout(5,Ny,k) - Ey_NSCBCout(5,Nx,k) - dFz_5(Nx,Ny,k) + dFvx_5(Nx,Ny,k) + dFvy_5(Nx,Ny,k) + dFvz_5(Nx,Ny,k)
	else if(i==Nx .and. j/=1 .and. j/=Ny) then!(Nx,2:Ny-1,k)流出
		Q1_1(Nx,j,k) = - Ex_NSCBCout(1,j,k) - dFy_1(Nx,j,k) - dFz_1(Nx,j,k)
		Q2_1(Nx,j,k) = - Ex_NSCBCout(2,j,k) - dFy_2(Nx,j,k) - dFz_2(Nx,j,k) + dFvx_2(Nx,j,k) + dFvy_2(Nx,j,k) + dFvz_2(Nx,j,k)
		Q3_1(Nx,j,k) = - Ex_NSCBCout(3,j,k) - dFy_3(Nx,j,k) - dFz_3(Nx,j,k) + dFvx_3(Nx,j,k) + dFvy_3(Nx,j,k) + dFvz_3(Nx,j,k)
		Q4_1(Nx,j,k) = - Ex_NSCBCout(4,j,k) - dFy_4(Nx,j,k) - dFz_4(Nx,j,k) + dFvx_4(Nx,j,k) + dFvy_4(Nx,j,k) + dFvz_4(Nx,j,k)
		Q5_1(Nx,j,k) = - Ex_NSCBCout(5,j,k) - dFy_5(Nx,j,k) - dFz_5(Nx,j,k) + dFvx_5(Nx,j,k) + dFvy_5(Nx,j,k) + dFvz_5(Nx,j,k)
	else if(i/=1 .and. i/=Nx .and. j==1) then!壁
		Q1_1(i,1,k) = - dFx_1(i,1,k) - Ey_NSCBCin(1,i,k) - dFz_1(i,1,k)
		Q2_1(i,1,k) = 0d0!- dFx_2(i,1,k) - Ey_NSCBCin(2,i,k) - dFz_2(i,1,k) + dFvx_2(i,1,k) + dFvy_2(i,1,k) + dFvz_2(i,1,k)
		Q3_1(i,1,k) = 0d0!- dFx_3(i,1,k) - Ey_NSCBCin(3,i,k) - dFz_3(i,1,k) + dFvx_3(i,1,k) + dFvy_3(i,1,k) + dFvz_3(i,1,k)
		Q4_1(i,1,k) = 0d0!- dFx_4(i,1,k) - Ey_NSCBCin(4,i,k) - dFz_4(i,1,k) + dFvx_4(i,1,k) + dFvy_4(i,1,k) + dFvz_4(i,1,k)
		Q5_1(i,1,k) = 0d0!- dFx_5(i,1,k) - Ey_NSCBCin(5,i,k) - dFz_5(i,1,k) + dFvx_5(i,1,k) + dFvy_5(i,1,k) + dFvz_5(i,1,k)
	else if(i/=1 .and. i/=Nx .and. j==Ny) then!流出
		Q1_1(i,Ny,k) = - dFx_1(i,Ny,k) - Ey_NSCBCout(1,i,k) - dFz_1(i,Ny,k)
		Q2_1(i,Ny,k) = - dFx_2(i,Ny,k) - Ey_NSCBCout(2,i,k) - dFz_2(i,Ny,k) + dFvx_2(i,Ny,k) + dFvy_2(i,Ny,k) + dFvz_2(i,Ny,k)
		Q3_1(i,Ny,k) = - dFx_3(i,Ny,k) - Ey_NSCBCout(3,i,k) - dFz_3(i,Ny,k) + dFvx_3(i,Ny,k) + dFvy_3(i,Ny,k) + dFvz_3(i,Ny,k)
		Q4_1(i,Ny,k) = - dFx_4(i,Ny,k) - Ey_NSCBCout(4,i,k) - dFz_4(i,Ny,k) + dFvx_4(i,Ny,k) + dFvy_4(i,Ny,k) + dFvz_4(i,Ny,k)
		Q5_1(i,Ny,k) = - dFx_5(i,Ny,k) - Ey_NSCBCout(5,i,k) - dFz_5(i,Ny,k) + dFvx_5(i,Ny,k) + dFvy_5(i,Ny,k) + dFvz_5(i,Ny,k)
	else
		Q1_1(i,j,k) = - dFx_1(i,j,k) - dFy_1(i,j,k) - dFz_1(i,j,k)
		Q2_1(i,j,k) = - dFx_2(i,j,k) - dFy_2(i,j,k) - dFz_2(i,j,k) + dFvx_2(i,j,k) + dFvy_2(i,j,k) + dFvz_2(i,j,k)
		Q3_1(i,j,k) = - dFx_3(i,j,k) - dFy_3(i,j,k) - dFz_3(i,j,k) + dFvx_3(i,j,k) + dFvy_3(i,j,k) + dFvz_3(i,j,k)
		Q4_1(i,j,k) = - dFx_4(i,j,k) - dFy_4(i,j,k) - dFz_4(i,j,k) + dFvx_4(i,j,k) + dFvy_4(i,j,k) + dFvz_4(i,j,k)
		Q5_1(i,j,k) = - dFx_5(i,j,k) - dFy_5(i,j,k) - dFz_5(i,j,k) + dFvx_5(i,j,k) + dFvy_5(i,j,k) + dFvz_5(i,j,k)
	end if
		k1_1(i,j,k) = Q1(i,j,k) + 0.5d0*Q1_1(i,j,k)*dt
		k2_1(i,j,k) = Q2(i,j,k) + 0.5d0*Q2_1(i,j,k)*dt
		k3_1(i,j,k) = Q3(i,j,k) + 0.5d0*Q3_1(i,j,k)*dt
		k4_1(i,j,k) = Q4(i,j,k) + 0.5d0*Q4_1(i,j,k)*dt
		k5_1(i,j,k) = Q5(i,j,k) + 0.5d0*Q5_1(i,j,k)*dt

		Q1(i,j,k) = k1_1(i,j,k)
		Q2(i,j,k) = k2_1(i,j,k)
		Q3(i,j,k) = k3_1(i,j,k)
		Q4(i,j,k) = k4_1(i,j,k)
		Q5(i,j,k) = k5_1(i,j,k)
	end do
	end do
	end do
	!$omp end parallel do
!!!!=============================================================================
!!!!==============================================================================
!	num2=num
!	num=992
!	call output_files()
	call BL_init()!!流入部
	call update()
	call visco_bibun()
	call bibun_adv()
!	num=num2

	!2
	!$omp parallel do
	do k=1,Nz
	do j=1,Ny
	do i=1,Nx
	if(i==1 .and. j==1) then
		Q1_2(1,1,k) = - Ex_NSCBCin(1,1,k) - Ey_NSCBCin(1,1,k) - dFz_1(1,1,k)
		Q2_2(1,1,k) = 0d0!- Ex_NSCBCin(2,1,k) - Ey_NSCBCin(2,1,k) - dFz_2(1,1,k) + dFvx_2(1,1,k) + dFvy_2(1,1,k) + dFvz_2(1,1,k)
		Q3_2(1,1,k) = 0d0!- Ex_NSCBCin(3,1,k) - Ey_NSCBCin(3,1,k) - dFz_3(1,1,k) + dFvx_3(1,1,k) + dFvy_3(1,1,k) + dFvz_3(1,1,k)
		Q4_2(1,1,k) = 0d0!- Ex_NSCBCin(4,1,k) - Ey_NSCBCin(4,1,k) - dFz_4(1,1,k) + dFvx_4(1,1,k) + dFvy_4(1,1,k) + dFvz_4(1,1,k)
		Q5_2(1,1,k) = 0d0!- Ex_NSCBCin(5,1,k) - Ey_NSCBCin(5,1,k) - dFz_5(1,1,k) + dFvx_5(1,1,k) + dFvy_5(1,1,k) + dFvz_5(1,1,k)
	else if(i==1 .and. j==Ny) then
		Q1_2(1,Ny,k) = - Ex_NSCBCin(1,Ny,k) - Ey_NSCBCout(1,1,k) - dFz_1(1,Ny,k)
		Q2_2(1,Ny,k) = 0d0!- Ex_NSCBCin(2,Ny,k) - Ey_NSCBCout(2,1,k) - dFz_2(1,Ny,k) + dFvx_2(1,Ny,k) + dFvy_2(1,Ny,k) + dFvz_2(1,Ny,k)
		Q3_2(1,Ny,k) = 0d0!- Ex_NSCBCin(3,Ny,k) - Ey_NSCBCout(3,1,k) - dFz_3(1,Ny,k) + dFvx_3(1,Ny,k) + dFvy_3(1,Ny,k) + dFvz_3(1,Ny,k)
		Q4_2(1,Ny,k) = 0d0!- Ex_NSCBCin(4,Ny,k) - Ey_NSCBCout(4,1,k) - dFz_4(1,Ny,k) + dFvx_4(1,Ny,k) + dFvy_4(1,Ny,k) + dFvz_4(1,Ny,k)
		Q5_2(1,Ny,k) = 0d0!- Ex_NSCBCin(5,Ny,k) - Ey_NSCBCout(5,1,k) - dFz_5(1,Ny,k) + dFvx_5(1,Ny,k) + dFvy_5(1,Ny,k) + dFvz_5(1,Ny,k)
	else if(i==1 .and. j/=1 .and. j/=Ny) then
		Q1_2(1,j,k) = - Ex_NSCBCin(1,j,k) - dFy_1(1,j,k) - dFz_1(1,j,k)
		Q2_2(1,j,k) = 0d0!- Ex_NSCBCin(2,j,k) - dFy_2(1,j,k) - dFz_2(1,j,k) + dFvx_2(1,j,k) + dFvy_2(1,j,k) + dFvz_2(1,j,k)
		Q3_2(1,j,k) = 0d0!- Ex_NSCBCin(3,j,k) - dFy_3(1,j,k) - dFz_3(1,j,k) + dFvx_3(1,j,k) + dFvy_3(1,j,k) + dFvz_3(1,j,k)
		Q4_2(1,j,k) = 0d0!- Ex_NSCBCin(4,j,k) - dFy_4(1,j,k) - dFz_4(1,j,k) + dFvx_4(1,j,k) + dFvy_4(1,j,k) + dFvz_4(1,j,k)
		Q5_2(1,j,k) = 0d0!- Ex_NSCBCin(5,j,k) - dFy_5(1,j,k) - dFz_5(1,j,k) + dFvx_5(1,j,k) + dFvy_5(1,j,k) + dFvz_5(1,j,k)
	else if(i==Nx .and. j==1) then
		Q1_2(Nx,1,k) = - Ex_NSCBCout(1,1,k) - Ey_NSCBCin(1,Nx,k) - dFz_1(Nx,1,k)
		Q2_2(Nx,1,k) = 0d0!- Ex_NSCBCout(2,1,k) - Ey_NSCBCin(2,Nx,k) - dFz_2(Nx,1,k) + dFvx_2(Nx,1,k) + dFvy_2(Nx,1,k) + dFvz_2(Nx,1,k)
		Q3_2(Nx,1,k) = 0d0!- Ex_NSCBCout(3,1,k) - Ey_NSCBCin(3,Nx,k) - dFz_3(Nx,1,k) + dFvx_3(Nx,1,k) + dFvy_3(Nx,1,k) + dFvz_3(Nx,1,k)
		Q4_2(Nx,1,k) = 0d0!- Ex_NSCBCout(4,1,k) - Ey_NSCBCin(4,Nx,k) - dFz_4(Nx,1,k) + dFvx_4(Nx,1,k) + dFvy_4(Nx,1,k) + dFvz_4(Nx,1,k)
		Q5_2(Nx,1,k) = 0d0!- Ex_NSCBCout(5,1,k) - Ey_NSCBCin(5,Nx,k) - dFz_5(Nx,1,k) + dFvx_5(Nx,1,k) + dFvy_5(Nx,1,k) + dFvz_5(Nx,1,k)
	else if(i==Nx .and. j==Ny) then
		Q1_2(Nx,Ny,k) = - Ex_NSCBCout(1,Ny,k) - Ey_NSCBCout(1,Nx,k) - dFz_1(Nx,Ny,k)
		Q2_2(Nx,Ny,k) = - Ex_NSCBCout(2,Ny,k) - Ey_NSCBCout(2,Nx,k) - dFz_2(Nx,Ny,k) + dFvx_2(Nx,Ny,k) + dFvy_2(Nx,Ny,k) + dFvz_2(Nx,Ny,k)
		Q3_2(Nx,Ny,k) = - Ex_NSCBCout(3,Ny,k) - Ey_NSCBCout(3,Nx,k) - dFz_3(Nx,Ny,k) + dFvx_3(Nx,Ny,k) + dFvy_3(Nx,Ny,k) + dFvz_3(Nx,Ny,k)
		Q4_2(Nx,Ny,k) = - Ex_NSCBCout(4,Ny,k) - Ey_NSCBCout(4,Nx,k) - dFz_4(Nx,Ny,k) + dFvx_4(Nx,Ny,k) + dFvy_4(Nx,Ny,k) + dFvz_4(Nx,Ny,k)
		Q5_2(Nx,Ny,k) = - Ex_NSCBCout(5,Ny,k) - Ey_NSCBCout(5,Nx,k) - dFz_5(Nx,Ny,k) + dFvx_5(Nx,Ny,k) + dFvy_5(Nx,Ny,k) + dFvz_5(Nx,Ny,k)
	else if(i==Nx .and. j/=1 .and. j/=Ny) then
		Q1_2(Nx,j,k) = - Ex_NSCBCout(1,j,k) - dFy_1(Nx,j,k) - dFz_1(Nx,j,k)
		Q2_2(Nx,j,k) = - Ex_NSCBCout(2,j,k) - dFy_2(Nx,j,k) - dFz_2(Nx,j,k) + dFvx_2(Nx,j,k) + dFvy_2(Nx,j,k) + dFvz_2(Nx,j,k)
		Q3_2(Nx,j,k) = - Ex_NSCBCout(3,j,k) - dFy_3(Nx,j,k) - dFz_3(Nx,j,k) + dFvx_3(Nx,j,k) + dFvy_3(Nx,j,k) + dFvz_3(Nx,j,k)
		Q4_2(Nx,j,k) = - Ex_NSCBCout(4,j,k) - dFy_4(Nx,j,k) - dFz_4(Nx,j,k) + dFvx_4(Nx,j,k) + dFvy_4(Nx,j,k) + dFvz_4(Nx,j,k)
		Q5_2(Nx,j,k) = - Ex_NSCBCout(5,j,k) - dFy_5(Nx,j,k) - dFz_5(Nx,j,k) + dFvx_5(Nx,j,k) + dFvy_5(Nx,j,k) + dFvz_5(Nx,j,k)
	else if(i/=1 .and. i/=Nx .and. j==1) then
		Q1_2(i,1,k) = - dFx_1(i,1,k) - Ey_NSCBCin(1,i,k) - dFz_1(i,1,k)
		Q2_2(i,1,k) = 0d0!- dFx_2(i,1,k) - Ey_NSCBCin(2,i,k) - dFz_2(i,1,k) + dFvx_2(i,1,k) + dFvy_2(i,1,k) + dFvz_2(i,1,k)
		Q3_2(i,1,k) = 0d0!- dFx_3(i,1,k) - Ey_NSCBCin(3,i,k) - dFz_3(i,1,k) + dFvx_3(i,1,k) + dFvy_3(i,1,k) + dFvz_3(i,1,k)
		Q4_2(i,1,k) = 0d0!- dFx_4(i,1,k) - Ey_NSCBCin(4,i,k) - dFz_4(i,1,k) + dFvx_4(i,1,k) + dFvy_4(i,1,k) + dFvz_4(i,1,k)
		Q5_2(i,1,k) = 0d0!- dFx_5(i,1,k) - Ey_NSCBCin(5,i,k) - dFz_5(i,1,k) + dFvx_5(i,1,k) + dFvy_5(i,1,k) + dFvz_5(i,1,k)
	else if(i/=1 .and. i/=Nx .and. j==Ny) then
		Q1_2(i,Ny,k) = - dFx_1(i,Ny,k) - Ey_NSCBCout(1,i,k) - dFz_1(i,Ny,k)
		Q2_2(i,Ny,k) = - dFx_2(i,Ny,k) - Ey_NSCBCout(2,i,k) - dFz_2(i,Ny,k) + dFvx_2(i,Ny,k) + dFvy_2(i,Ny,k) + dFvz_2(i,Ny,k)
		Q3_2(i,Ny,k) = - dFx_3(i,Ny,k) - Ey_NSCBCout(3,i,k) - dFz_3(i,Ny,k) + dFvx_3(i,Ny,k) + dFvy_3(i,Ny,k) + dFvz_3(i,Ny,k)
		Q4_2(i,Ny,k) = - dFx_4(i,Ny,k) - Ey_NSCBCout(4,i,k) - dFz_4(i,Ny,k) + dFvx_4(i,Ny,k) + dFvy_4(i,Ny,k) + dFvz_4(i,Ny,k)
		Q5_2(i,Ny,k) = - dFx_5(i,Ny,k) - Ey_NSCBCout(5,i,k) - dFz_5(i,Ny,k) + dFvx_5(i,Ny,k) + dFvy_5(i,Ny,k) + dFvz_5(i,Ny,k)
	else
		Q1_2(i,j,k) = - dFx_1(i,j,k) - dFy_1(i,j,k) - dFz_1(i,j,k)
		Q2_2(i,j,k) = - dFx_2(i,j,k) - dFy_2(i,j,k) - dFz_2(i,j,k) + dFvx_2(i,j,k) + dFvy_2(i,j,k) + dFvz_2(i,j,k)
		Q3_2(i,j,k) = - dFx_3(i,j,k) - dFy_3(i,j,k) - dFz_3(i,j,k) + dFvx_3(i,j,k) + dFvy_3(i,j,k) + dFvz_3(i,j,k)
		Q4_2(i,j,k) = - dFx_4(i,j,k) - dFy_4(i,j,k) - dFz_4(i,j,k) + dFvx_4(i,j,k) + dFvy_4(i,j,k) + dFvz_4(i,j,k)
		Q5_2(i,j,k) = - dFx_5(i,j,k) - dFy_5(i,j,k) - dFz_5(i,j,k) + dFvx_5(i,j,k) + dFvy_5(i,j,k) + dFvz_5(i,j,k)
	end if
		k1_2(i,j,k) = Q1_0(i,j,k) + 0.5d0*Q1_2(i,j,k)*dt
		k2_2(i,j,k) = Q2_0(i,j,k) + 0.5d0*Q2_2(i,j,k)*dt
		k3_2(i,j,k) = Q3_0(i,j,k) + 0.5d0*Q3_2(i,j,k)*dt
		k4_2(i,j,k) = Q4_0(i,j,k) + 0.5d0*Q4_2(i,j,k)*dt
		k5_2(i,j,k) = Q5_0(i,j,k) + 0.5d0*Q5_2(i,j,k)*dt

		Q1(i,j,k) = k1_2(i,j,k)
		Q2(i,j,k) = k2_2(i,j,k)
		Q3(i,j,k) = k3_2(i,j,k)
		Q4(i,j,k) = k4_2(i,j,k)
		Q5(i,j,k) = k5_2(i,j,k)
	end do
	end do
	end do
	!$omp end parallel do
!!!!=============================================================================
!!!!==============================================================================
!	open(95, file='NSCBC.csv', status='unknown', form='formatted')
!	do k=1,Nz
!	do j=1,Ny
!		do i=1,Nx
!			write(95, "(10(1x,g15.7))") i,j,k,dFx_2(i,j,k),dFy_2(i,j,k),dFz_2(i,j,k), Ex_NSCBCin(2,j,k), Ex_NSCBCout(2,j,k), Ey_NSCBCin(2,i,k), Ey_NSCBCout(2,i,k)
!		end do
!	end do
!	end do
!	close(95)
!write(*,*) 'Check_RKstep2'

!	num2=num
!	num=993
!	call output_files()
	call BL_init()!!流入部
	call update()
	call visco_bibun()
	call bibun_adv()
!		write(*,*) 'CheckRK22'
!		call NaN_check()
!	num=num2

	!3
	!$omp parallel do
	do k=1,Nz
	do j=1,Ny
	do i=1,Nx
!		Q1_3(i,j,k) = - dFx_1(i,j,k) - dFy_1(i,j,k) - dFz_1(i,j,k)
!		Q2_3(i,j,k) = - dFx_2(i,j,k) - dFy_2(i,j,k) - dFz_2(i,j,k) + dFvx_2(i,j,k) + dFvy_2(i,j,k) + dFvz_2(i,j,k)
!		Q3_3(i,j,k) = - dFx_3(i,j,k) - dFy_3(i,j,k) - dFz_3(i,j,k) + dFvx_3(i,j,k) + dFvy_3(i,j,k) + dFvz_3(i,j,k)
!		Q4_3(i,j,k) = - dFx_4(i,j,k) - dFy_4(i,j,k) - dFz_4(i,j,k) + dFvx_4(i,j,k) + dFvy_4(i,j,k) + dFvz_4(i,j,k)
!		Q5_3(i,j,k) = - dFx_5(i,j,k) - dFy_5(i,j,k) - dFz_5(i,j,k) + dFvx_5(i,j,k) + dFvy_5(i,j,k) + dFvz_5(i,j,k)
	if(i==1 .and. j==1) then
		Q1_3(1,1,k) = - Ex_NSCBCin(1,1,k) - Ey_NSCBCin(1,1,k) - dFz_1(1,1,k)
		Q2_3(1,1,k) = 0d0!- Ex_NSCBCin(2,1,k) - Ey_NSCBCin(2,1,k) - dFz_2(1,1,k) + dFvx_2(1,1,k) + dFvy_2(1,1,k) + dFvz_2(1,1,k)
		Q3_3(1,1,k) = 0d0!- Ex_NSCBCin(3,1,k) - Ey_NSCBCin(3,1,k) - dFz_3(1,1,k) + dFvx_3(1,1,k) + dFvy_3(1,1,k) + dFvz_3(1,1,k)
		Q4_3(1,1,k) = 0d0!- Ex_NSCBCin(4,1,k) - Ey_NSCBCin(4,1,k) - dFz_4(1,1,k) + dFvx_4(1,1,k) + dFvy_4(1,1,k) + dFvz_4(1,1,k)
		Q5_3(1,1,k) = 0d0!- Ex_NSCBCin(5,1,k) - Ey_NSCBCin(5,1,k) - dFz_5(1,1,k) + dFvx_5(1,1,k) + dFvy_5(1,1,k) + dFvz_5(1,1,k)
	else if(i==1 .and. j==Ny) then
		Q1_3(1,Ny,k) = - Ex_NSCBCin(1,Ny,k) - Ey_NSCBCout(1,1,k) - dFz_1(1,Ny,k)
		Q2_3(1,Ny,k) = 0d0!- Ex_NSCBCin(2,Ny,k) - Ey_NSCBCout(2,1,k) - dFz_2(1,Ny,k) + dFvx_2(1,Ny,k) + dFvy_2(1,Ny,k) + dFvz_2(1,Ny,k)
		Q3_3(1,Ny,k) = 0d0!- Ex_NSCBCin(3,Ny,k) - Ey_NSCBCout(3,1,k) - dFz_3(1,Ny,k) + dFvx_3(1,Ny,k) + dFvy_3(1,Ny,k) + dFvz_3(1,Ny,k)
		Q4_3(1,Ny,k) = 0d0!- Ex_NSCBCin(4,Ny,k) - Ey_NSCBCout(4,1,k) - dFz_4(1,Ny,k) + dFvx_4(1,Ny,k) + dFvy_4(1,Ny,k) + dFvz_4(1,Ny,k)
		Q5_3(1,Ny,k) = 0d0!- Ex_NSCBCin(5,Ny,k) - Ey_NSCBCout(5,1,k) - dFz_5(1,Ny,k) + dFvx_5(1,Ny,k) + dFvy_5(1,Ny,k) + dFvz_5(1,Ny,k)
	else if(i==1 .and. j/=1 .and. j/=Ny) then
		Q1_3(1,j,k) = - Ex_NSCBCin(1,j,k) - dFy_1(1,j,k) - dFz_1(1,j,k)
		Q2_3(1,j,k) = 0d0!- Ex_NSCBCin(2,j,k) - dFy_2(1,j,k) - dFz_2(1,j,k) + dFvx_2(1,j,k) + dFvy_2(1,j,k) + dFvz_2(1,j,k)
		Q3_3(1,j,k) = 0d0!- Ex_NSCBCin(3,j,k) - dFy_3(1,j,k) - dFz_3(1,j,k) + dFvx_3(1,j,k) + dFvy_3(1,j,k) + dFvz_3(1,j,k)
		Q4_3(1,j,k) = 0d0!- Ex_NSCBCin(4,j,k) - dFy_4(1,j,k) - dFz_4(1,j,k) + dFvx_4(1,j,k) + dFvy_4(1,j,k) + dFvz_4(1,j,k)
		Q5_3(1,j,k) = 0d0!- Ex_NSCBCin(5,j,k) - dFy_5(1,j,k) - dFz_5(1,j,k) + dFvx_5(1,j,k) + dFvy_5(1,j,k) + dFvz_5(1,j,k)
	else if(i==Nx .and. j==1) then
		Q1_3(Nx,1,k) = - Ex_NSCBCout(1,1,k) - Ey_NSCBCin(1,Nx,k) - dFz_1(Nx,1,k)
		Q2_3(Nx,1,k) = 0d0!- Ex_NSCBCout(2,1,k) - Ey_NSCBCin(2,Nx,k) - dFz_2(Nx,1,k) + dFvx_2(Nx,1,k) + dFvy_2(Nx,1,k) + dFvz_2(Nx,1,k)
		Q3_3(Nx,1,k) = 0d0!- Ex_NSCBCout(3,1,k) - Ey_NSCBCin(3,Nx,k) - dFz_3(Nx,1,k) + dFvx_3(Nx,1,k) + dFvy_3(Nx,1,k) + dFvz_3(Nx,1,k)
		Q4_3(Nx,1,k) = 0d0!- Ex_NSCBCout(4,1,k) - Ey_NSCBCin(4,Nx,k) - dFz_4(Nx,1,k) + dFvx_4(Nx,1,k) + dFvy_4(Nx,1,k) + dFvz_4(Nx,1,k)
		Q5_3(Nx,1,k) = 0d0!- Ex_NSCBCout(5,1,k) - Ey_NSCBCin(5,Nx,k) - dFz_5(Nx,1,k) + dFvx_5(Nx,1,k) + dFvy_5(Nx,1,k) + dFvz_5(Nx,1,k)
	else if(i==Nx .and. j==Ny) then
		Q1_3(Nx,Ny,k) = - Ex_NSCBCout(1,Ny,k) - Ey_NSCBCout(1,Nx,k) - dFz_1(Nx,Ny,k)
		Q2_3(Nx,Ny,k) = - Ex_NSCBCout(2,Ny,k) - Ey_NSCBCout(2,Nx,k) - dFz_2(Nx,Ny,k) + dFvx_2(Nx,Ny,k) + dFvy_2(Nx,Ny,k) + dFvz_2(Nx,Ny,k)
		Q3_3(Nx,Ny,k) = - Ex_NSCBCout(3,Ny,k) - Ey_NSCBCout(3,Nx,k) - dFz_3(Nx,Ny,k) + dFvx_3(Nx,Ny,k) + dFvy_3(Nx,Ny,k) + dFvz_3(Nx,Ny,k)
		Q4_3(Nx,Ny,k) = - Ex_NSCBCout(4,Ny,k) - Ey_NSCBCout(4,Nx,k) - dFz_4(Nx,Ny,k) + dFvx_4(Nx,Ny,k) + dFvy_4(Nx,Ny,k) + dFvz_4(Nx,Ny,k)
		Q5_3(Nx,Ny,k) = - Ex_NSCBCout(5,Ny,k) - Ey_NSCBCout(5,Nx,k) - dFz_5(Nx,Ny,k) + dFvx_5(Nx,Ny,k) + dFvy_5(Nx,Ny,k) + dFvz_5(Nx,Ny,k)
	else if(i==Nx .and. j/=1 .and. j/=Ny) then
		Q1_3(Nx,j,k) = - Ex_NSCBCout(1,j,k) - dFy_1(Nx,j,k) - dFz_1(Nx,j,k)
		Q2_3(Nx,j,k) = - Ex_NSCBCout(2,j,k) - dFy_2(Nx,j,k) - dFz_2(Nx,j,k) + dFvx_2(Nx,j,k) + dFvy_2(Nx,j,k) + dFvz_2(Nx,j,k)
		Q3_3(Nx,j,k) = - Ex_NSCBCout(3,j,k) - dFy_3(Nx,j,k) - dFz_3(Nx,j,k) + dFvx_3(Nx,j,k) + dFvy_3(Nx,j,k) + dFvz_3(Nx,j,k)
		Q4_3(Nx,j,k) = - Ex_NSCBCout(4,j,k) - dFy_4(Nx,j,k) - dFz_4(Nx,j,k) + dFvx_4(Nx,j,k) + dFvy_4(Nx,j,k) + dFvz_4(Nx,j,k)
		Q5_3(Nx,j,k) = - Ex_NSCBCout(5,j,k) - dFy_5(Nx,j,k) - dFz_5(Nx,j,k) + dFvx_5(Nx,j,k) + dFvy_5(Nx,j,k) + dFvz_5(Nx,j,k)
	else if(i/=1 .and. i/=Nx .and. j==1) then
		Q1_3(i,1,k) = - dFx_1(i,1,k) - Ey_NSCBCin(1,i,k) - dFz_1(i,1,k)
		Q2_3(i,1,k) = 0d0!- dFx_2(i,1,k) - Ey_NSCBCin(2,i,k) - dFz_2(i,1,k) + dFvx_2(i,1,k) + dFvy_2(i,1,k) + dFvz_2(i,1,k)
		Q3_3(i,1,k) = 0d0!- dFx_3(i,1,k) - Ey_NSCBCin(3,i,k) - dFz_3(i,1,k) + dFvx_3(i,1,k) + dFvy_3(i,1,k) + dFvz_3(i,1,k)
		Q4_3(i,1,k) = 0d0!- dFx_4(i,1,k) - Ey_NSCBCin(4,i,k) - dFz_4(i,1,k) + dFvx_4(i,1,k) + dFvy_4(i,1,k) + dFvz_4(i,1,k)
		Q5_3(i,1,k) = 0d0!- dFx_5(i,1,k) - Ey_NSCBCin(5,i,k) - dFz_5(i,1,k) + dFvx_5(i,1,k) + dFvy_5(i,1,k) + dFvz_5(i,1,k)
	else if(i/=1 .and. i/=Nx .and. j==Ny) then
		Q1_3(i,Ny,k) = - dFx_1(i,Ny,k) - Ey_NSCBCout(1,i,k) - dFz_1(i,Ny,k)
		Q2_3(i,Ny,k) = - dFx_2(i,Ny,k) - Ey_NSCBCout(2,i,k) - dFz_2(i,Ny,k) + dFvx_2(i,Ny,k) + dFvy_2(i,Ny,k) + dFvz_2(i,Ny,k)
		Q3_3(i,Ny,k) = - dFx_3(i,Ny,k) - Ey_NSCBCout(3,i,k) - dFz_3(i,Ny,k) + dFvx_3(i,Ny,k) + dFvy_3(i,Ny,k) + dFvz_3(i,Ny,k)
		Q4_3(i,Ny,k) = - dFx_4(i,Ny,k) - Ey_NSCBCout(4,i,k) - dFz_4(i,Ny,k) + dFvx_4(i,Ny,k) + dFvy_4(i,Ny,k) + dFvz_4(i,Ny,k)
		Q5_3(i,Ny,k) = - dFx_5(i,Ny,k) - Ey_NSCBCout(5,i,k) - dFz_5(i,Ny,k) + dFvx_5(i,Ny,k) + dFvy_5(i,Ny,k) + dFvz_5(i,Ny,k)
	else
		Q1_3(i,j,k) = - dFx_1(i,j,k) - dFy_1(i,j,k) - dFz_1(i,j,k)
		Q2_3(i,j,k) = - dFx_2(i,j,k) - dFy_2(i,j,k) - dFz_2(i,j,k) + dFvx_2(i,j,k) + dFvy_2(i,j,k) + dFvz_2(i,j,k)
		Q3_3(i,j,k) = - dFx_3(i,j,k) - dFy_3(i,j,k) - dFz_3(i,j,k) + dFvx_3(i,j,k) + dFvy_3(i,j,k) + dFvz_3(i,j,k)
		Q4_3(i,j,k) = - dFx_4(i,j,k) - dFy_4(i,j,k) - dFz_4(i,j,k) + dFvx_4(i,j,k) + dFvy_4(i,j,k) + dFvz_4(i,j,k)
		Q5_3(i,j,k) = - dFx_5(i,j,k) - dFy_5(i,j,k) - dFz_5(i,j,k) + dFvx_5(i,j,k) + dFvy_5(i,j,k) + dFvz_5(i,j,k)
	end if

		k1_3(i,j,k) = Q1_0(i,j,k) + Q1_3(i,j,k)*dt
		k2_3(i,j,k) = Q2_0(i,j,k) + Q2_3(i,j,k)*dt
		k3_3(i,j,k) = Q3_0(i,j,k) + Q3_3(i,j,k)*dt
		k4_3(i,j,k) = Q4_0(i,j,k) + Q4_3(i,j,k)*dt
		k5_3(i,j,k) = Q5_0(i,j,k) + Q5_3(i,j,k)*dt

		Q1(i,j,k) = k1_3(i,j,k)
		Q2(i,j,k) = k2_3(i,j,k)
		Q3(i,j,k) = k3_3(i,j,k)
		Q4(i,j,k) = k4_3(i,j,k)
		Q5(i,j,k) = k5_3(i,j,k)
	end do
	end do
	end do
	!$omp end parallel do
!!!!=============================================================================
!!!!==============================================================================
!	num2=num
!	num=994
!	call output_files()

	call BL_init()!!流入部
	call update()
	call visco_bibun()
	call bibun_adv()
!		write(*,*) 'CheckRK23'
!		call NaN_check()
!	num=num2

	!4
	!$omp parallel do
	do k=1,Nz
	do j=1,Ny
	do i=1,Nx
	if(i==1 .and. j==1) then
		Q1_4(1,1,k) = - Ex_NSCBCin(1,1,k) - Ey_NSCBCin(1,1,k) - dFz_1(1,1,k)
		Q2_4(1,1,k) = 0d0!- Ex_NSCBCin(2,1,k) - Ey_NSCBCin(2,1,k) - dFz_2(1,1,k) + dFvx_2(1,1,k) + dFvy_2(1,1,k) + dFvz_2(1,1,k)
		Q3_4(1,1,k) = 0d0!- Ex_NSCBCin(3,1,k) - Ey_NSCBCin(3,1,k) - dFz_3(1,1,k) + dFvx_3(1,1,k) + dFvy_3(1,1,k) + dFvz_3(1,1,k)
		Q4_4(1,1,k) = 0d0!- Ex_NSCBCin(4,1,k) - Ey_NSCBCin(4,1,k) - dFz_4(1,1,k) + dFvx_4(1,1,k) + dFvy_4(1,1,k) + dFvz_4(1,1,k)
		Q5_4(1,1,k) = 0d0!- Ex_NSCBCin(5,1,k) - Ey_NSCBCin(5,1,k) - dFz_5(1,1,k) + dFvx_5(1,1,k) + dFvy_5(1,1,k) + dFvz_5(1,1,k)
	else if(i==1 .and. j==Ny) then
		Q1_4(1,Ny,k) = - Ex_NSCBCin(1,Ny,k) - Ey_NSCBCout(1,1,k) - dFz_1(1,Ny,k)
		Q2_4(1,Ny,k) = 0d0!- Ex_NSCBCin(2,Ny,k) - Ey_NSCBCout(2,1,k) - dFz_2(1,Ny,k) + dFvx_2(1,Ny,k) + dFvy_2(1,Ny,k) + dFvz_2(1,Ny,k)
		Q3_4(1,Ny,k) = 0d0!- Ex_NSCBCin(3,Ny,k) - Ey_NSCBCout(3,1,k) - dFz_3(1,Ny,k) + dFvx_3(1,Ny,k) + dFvy_3(1,Ny,k) + dFvz_3(1,Ny,k)
		Q4_4(1,Ny,k) = 0d0!- Ex_NSCBCin(4,Ny,k) - Ey_NSCBCout(4,1,k) - dFz_4(1,Ny,k) + dFvx_4(1,Ny,k) + dFvy_4(1,Ny,k) + dFvz_4(1,Ny,k)
		Q5_4(1,Ny,k) = 0d0!- Ex_NSCBCin(5,Ny,k) - Ey_NSCBCout(5,1,k) - dFz_5(1,Ny,k) + dFvx_5(1,Ny,k) + dFvy_5(1,Ny,k) + dFvz_5(1,Ny,k)
	else if(i==1 .and. j/=1 .and. j/=Ny) then
		Q1_4(1,j,k) = - Ex_NSCBCin(1,j,k) - dFy_1(1,j,k) - dFz_1(1,j,k)
		Q2_4(1,j,k) = 0d0!- Ex_NSCBCin(2,j,k) - dFy_2(1,j,k) - dFz_2(1,j,k) + dFvx_2(1,j,k) + dFvy_2(1,j,k) + dFvz_2(1,j,k)
		Q3_4(1,j,k) = 0d0!- Ex_NSCBCin(3,j,k) - dFy_3(1,j,k) - dFz_3(1,j,k) + dFvx_3(1,j,k) + dFvy_3(1,j,k) + dFvz_3(1,j,k)
		Q4_4(1,j,k) = 0d0!- Ex_NSCBCin(4,j,k) - dFy_4(1,j,k) - dFz_4(1,j,k) + dFvx_4(1,j,k) + dFvy_4(1,j,k) + dFvz_4(1,j,k)
		Q5_4(1,j,k) = 0d0!- Ex_NSCBCin(5,j,k) - dFy_5(1,j,k) - dFz_5(1,j,k) + dFvx_5(1,j,k) + dFvy_5(1,j,k) + dFvz_5(1,j,k)
	else if(i==Nx .and. j==1) then
		Q1_4(Nx,1,k) = - Ex_NSCBCout(1,1,k) - Ey_NSCBCin(1,Nx,k) - dFz_1(Nx,1,k)
		Q2_4(Nx,1,k) = 0d0!- Ex_NSCBCout(2,1,k) - Ey_NSCBCin(2,Nx,k) - dFz_2(Nx,1,k) + dFvx_2(Nx,1,k) + dFvy_2(Nx,1,k) + dFvz_2(Nx,1,k)
		Q3_4(Nx,1,k) = 0d0!- Ex_NSCBCout(3,1,k) - Ey_NSCBCin(3,Nx,k) - dFz_3(Nx,1,k) + dFvx_3(Nx,1,k) + dFvy_3(Nx,1,k) + dFvz_3(Nx,1,k)
		Q4_4(Nx,1,k) = 0d0!- Ex_NSCBCout(4,1,k) - Ey_NSCBCin(4,Nx,k) - dFz_4(Nx,1,k) + dFvx_4(Nx,1,k) + dFvy_4(Nx,1,k) + dFvz_4(Nx,1,k)
		Q5_4(Nx,1,k) = 0d0!- Ex_NSCBCout(5,1,k) - Ey_NSCBCin(5,Nx,k) - dFz_5(Nx,1,k) + dFvx_5(Nx,1,k) + dFvy_5(Nx,1,k) + dFvz_5(Nx,1,k)
	else if(i==Nx .and. j==Ny) then
		Q1_4(Nx,Ny,k) = - Ex_NSCBCout(1,Ny,k) - Ey_NSCBCout(1,Nx,k) - dFz_1(Nx,Ny,k)
		Q2_4(Nx,Ny,k) = - Ex_NSCBCout(2,Ny,k) - Ey_NSCBCout(2,Nx,k) - dFz_2(Nx,Ny,k) + dFvx_2(Nx,Ny,k) + dFvy_2(Nx,Ny,k) + dFvz_2(Nx,Ny,k)
		Q3_4(Nx,Ny,k) = - Ex_NSCBCout(3,Ny,k) - Ey_NSCBCout(3,Nx,k) - dFz_3(Nx,Ny,k) + dFvx_3(Nx,Ny,k) + dFvy_3(Nx,Ny,k) + dFvz_3(Nx,Ny,k)
		Q4_4(Nx,Ny,k) = - Ex_NSCBCout(4,Ny,k) - Ey_NSCBCout(4,Nx,k) - dFz_4(Nx,Ny,k) + dFvx_4(Nx,Ny,k) + dFvy_4(Nx,Ny,k) + dFvz_4(Nx,Ny,k)
		Q5_4(Nx,Ny,k) = - Ex_NSCBCout(5,Ny,k) - Ey_NSCBCout(5,Nx,k) - dFz_5(Nx,Ny,k) + dFvx_5(Nx,Ny,k) + dFvy_5(Nx,Ny,k) + dFvz_5(Nx,Ny,k)
	else if(i==Nx .and. j/=1 .and. j/=Ny) then
		Q1_4(Nx,j,k) = - Ex_NSCBCout(1,j,k) - dFy_1(Nx,j,k) - dFz_1(Nx,j,k)
		Q2_4(Nx,j,k) = - Ex_NSCBCout(2,j,k) - dFy_2(Nx,j,k) - dFz_2(Nx,j,k) + dFvx_2(Nx,j,k) + dFvy_2(Nx,j,k) + dFvz_2(Nx,j,k)
		Q3_4(Nx,j,k) = - Ex_NSCBCout(3,j,k) - dFy_3(Nx,j,k) - dFz_3(Nx,j,k) + dFvx_3(Nx,j,k) + dFvy_3(Nx,j,k) + dFvz_3(Nx,j,k)
		Q4_4(Nx,j,k) = - Ex_NSCBCout(4,j,k) - dFy_4(Nx,j,k) - dFz_4(Nx,j,k) + dFvx_4(Nx,j,k) + dFvy_4(Nx,j,k) + dFvz_4(Nx,j,k)
		Q5_4(Nx,j,k) = - Ex_NSCBCout(5,j,k) - dFy_5(Nx,j,k) - dFz_5(Nx,j,k) + dFvx_5(Nx,j,k) + dFvy_5(Nx,j,k) + dFvz_5(Nx,j,k)
	else if(i/=1 .and. i/=Nx .and. j==1) then
		Q1_4(i,1,k) = - dFx_1(i,1,k) - Ey_NSCBCin(1,i,k) - dFz_1(i,1,k)
		Q2_4(i,1,k) = 0d0!- dFx_2(i,1,k) - Ey_NSCBCin(2,i,k) - dFz_2(i,1,k) + dFvx_2(i,1,k) + dFvy_2(i,1,k) + dFvz_2(i,1,k)
		Q3_4(i,1,k) = 0d0!- dFx_3(i,1,k) - Ey_NSCBCin(3,i,k) - dFz_3(i,1,k) + dFvx_3(i,1,k) + dFvy_3(i,1,k) + dFvz_3(i,1,k)
		Q4_4(i,1,k) = 0d0!- dFx_4(i,1,k) - Ey_NSCBCin(4,i,k) - dFz_4(i,1,k) + dFvx_4(i,1,k) + dFvy_4(i,1,k) + dFvz_4(i,1,k)
		Q5_4(i,1,k) = 0d0!- dFx_5(i,1,k) - Ey_NSCBCin(5,i,k) - dFz_5(i,1,k) + dFvx_5(i,1,k) + dFvy_5(i,1,k) + dFvz_5(i,1,k)
	else if(i/=1 .and. i/=Nx .and. j==Ny) then
		Q1_4(i,Ny,k) = - dFx_1(i,Ny,k) - Ey_NSCBCout(1,i,k) - dFz_1(i,Ny,k)
		Q2_4(i,Ny,k) = - dFx_2(i,Ny,k) - Ey_NSCBCout(2,i,k) - dFz_2(i,Ny,k) + dFvx_2(i,Ny,k) + dFvy_2(i,Ny,k) + dFvz_2(i,Ny,k)
		Q3_4(i,Ny,k) = - dFx_3(i,Ny,k) - Ey_NSCBCout(3,i,k) - dFz_3(i,Ny,k) + dFvx_3(i,Ny,k) + dFvy_3(i,Ny,k) + dFvz_3(i,Ny,k)
		Q4_4(i,Ny,k) = - dFx_4(i,Ny,k) - Ey_NSCBCout(4,i,k) - dFz_4(i,Ny,k) + dFvx_4(i,Ny,k) + dFvy_4(i,Ny,k) + dFvz_4(i,Ny,k)
		Q5_4(i,Ny,k) = - dFx_5(i,Ny,k) - Ey_NSCBCout(5,i,k) - dFz_5(i,Ny,k) + dFvx_5(i,Ny,k) + dFvy_5(i,Ny,k) + dFvz_5(i,Ny,k)
	else
		Q1_4(i,j,k) = - dFx_1(i,j,k) - dFy_1(i,j,k) - dFz_1(i,j,k)
		Q2_4(i,j,k) = - dFx_2(i,j,k) - dFy_2(i,j,k) - dFz_2(i,j,k) + dFvx_2(i,j,k) + dFvy_2(i,j,k) + dFvz_2(i,j,k)
		Q3_4(i,j,k) = - dFx_3(i,j,k) - dFy_3(i,j,k) - dFz_3(i,j,k) + dFvx_3(i,j,k) + dFvy_3(i,j,k) + dFvz_3(i,j,k)
		Q4_4(i,j,k) = - dFx_4(i,j,k) - dFy_4(i,j,k) - dFz_4(i,j,k) + dFvx_4(i,j,k) + dFvy_4(i,j,k) + dFvz_4(i,j,k)
		Q5_4(i,j,k) = - dFx_5(i,j,k) - dFy_5(i,j,k) - dFz_5(i,j,k) + dFvx_5(i,j,k) + dFvy_5(i,j,k) + dFvz_5(i,j,k)
	end if
	end do
	end do
	end do
!!!!=============================================================================
	!$omp end parallel do
	!$omp parallel do
	do k=1,Nz
	do j=1,Ny
	do i=1,Nx
		Q1_next(i,j,k) = Q1_0(i,j,k) + ( Q1_1(i,j,k) + 2.d0*Q1_2(i,j,k) + 2.d0*Q1_3(i,j,k) + Q1_4(i,j,k))*dt /6.d0
		Q2_next(i,j,k) = Q2_0(i,j,k) + ( Q2_1(i,j,k) + 2.d0*Q2_2(i,j,k) + 2.d0*Q2_3(i,j,k) + Q2_4(i,j,k))*dt /6.d0
		Q3_next(i,j,k) = Q3_0(i,j,k) + ( Q3_1(i,j,k) + 2.d0*Q3_2(i,j,k) + 2.d0*Q3_3(i,j,k) + Q3_4(i,j,k))*dt /6.d0
		Q4_next(i,j,k) = Q4_0(i,j,k) + ( Q4_1(i,j,k) + 2.d0*Q4_2(i,j,k) + 2.d0*Q4_3(i,j,k) + Q4_4(i,j,k))*dt /6.d0
		Q5_next(i,j,k) = Q5_0(i,j,k) + ( Q5_1(i,j,k) + 2.d0*Q5_2(i,j,k) + 2.d0*Q5_3(i,j,k) + Q5_4(i,j,k))*dt /6.d0

		Q1(i,j,k) = Q1_next(i,j,k)
		Q2(i,j,k) = Q2_next(i,j,k)
		Q3(i,j,k) = Q3_next(i,j,k)
		Q4(i,j,k) = Q4_next(i,j,k)
		Q5(i,j,k) = Q5_next(i,j,k)

	end do
	end do
	end do
!!!!=============================================================================
	!$omp end parallel do
	call BL_init()!!流入部
	num2=num
	num=999
	call output_files()
	num=num2
	end subroutine RK4
!!!-------------------------
	subroutine euler()
	call update()
	call visco_bibun()
	!advective
	call bibun_adv()
	!
	!$omp parallel do
	do k=1,Nz
	do j=1,Ny
	do i=1,Nx
!xyzすべて反映
	if(i==1 .and. j==1) then
		Q1_4(1,1,k) =  Q1(1,1,k) + (- Ex_NSCBCin(1,1,k) - Ey_NSCBCin(1,1,k) - dFz_1(1,1,k) )*dt
		Q2_4(1,1,k) =  Q2(1,1,k) + 0d0!(- Ex_NSCBCin(2,1,k) - Ey_NSCBCin(2,1,k) - dFz_2(1,1,k) + dFvx_2(1,1,k) + dFvy_2(1,1,k) + dFvz_2(1,1,k) )*dt
		Q3_4(1,1,k) =  Q3(1,1,k) + 0d0!(- Ex_NSCBCin(3,1,k) - Ey_NSCBCin(3,1,k) - dFz_3(1,1,k) + dFvx_3(1,1,k) + dFvy_3(1,1,k) + dFvz_3(1,1,k) )*dt
		Q4_4(1,1,k) =  Q4(1,1,k) + 0d0!(- Ex_NSCBCin(4,1,k) - Ey_NSCBCin(4,1,k) - dFz_4(1,1,k) + dFvx_4(1,1,k) + dFvy_4(1,1,k) + dFvz_4(1,1,k) )*dt
		Q5_4(1,1,k) =  Q5(1,1,k) + 0d0!(- Ex_NSCBCin(5,1,k) - Ey_NSCBCin(5,1,k) - dFz_5(1,1,k) + dFvx_5(1,1,k) + dFvy_5(1,1,k) + dFvz_5(1,1,k) )*dt
	else if(i==1 .and. j==Ny) then
		Q1_4(1,Ny,k) =  Q1(1,Ny,k) + (- Ex_NSCBCin(1,Ny,k) - Ey_NSCBCout(1,1,k) - dFz_1(1,Ny,k) )*dt
		Q2_4(1,Ny,k) =  Q2(1,Ny,k) + 0d0!(- Ex_NSCBCin(2,Ny,k) - Ey_NSCBCout(2,1,k) - dFz_2(1,Ny,k) + dFvx_2(1,Ny,k) + dFvy_2(1,Ny,k) + dFvz_2(1,Ny,k) )*dt
		Q3_4(1,Ny,k) =  Q3(1,Ny,k) + 0d0!(- Ex_NSCBCin(3,Ny,k) - Ey_NSCBCout(3,1,k) - dFz_3(1,Ny,k) + dFvx_3(1,Ny,k) + dFvy_3(1,Ny,k) + dFvz_3(1,Ny,k) )*dt
		Q4_4(1,Ny,k) =  Q4(1,Ny,k) + 0d0!(- Ex_NSCBCin(4,Ny,k) - Ey_NSCBCout(4,1,k) - dFz_4(1,Ny,k) + dFvx_4(1,Ny,k) + dFvy_4(1,Ny,k) + dFvz_4(1,Ny,k) )*dt
		Q5_4(1,Ny,k) =  Q5(1,Ny,k) + 0d0!(- Ex_NSCBCin(5,Ny,k) - Ey_NSCBCout(5,1,k) - dFz_5(1,Ny,k) + dFvx_5(1,Ny,k) + dFvy_5(1,Ny,k) + dFvz_5(1,Ny,k) )*dt
	else if(i==1 .and. j/=1 .and. j/=Ny) then
		Q1_4(1,j,k) = Q1(1,j,k) + (- Ex_NSCBCin(1,j,k) - dFy_1(1,j,k) - dFz_1(1,j,k) )*dt
		Q2_4(1,j,k) = Q2(1,j,k) + 0d0!(- Ex_NSCBCin(2,j,k) - dFy_2(1,j,k) - dFz_2(1,j,k) + dFvx_2(1,j,k) + dFvy_2(1,j,k) + dFvz_2(1,j,k) )*dt
		Q3_4(1,j,k) = Q3(1,j,k) + 0d0!(- Ex_NSCBCin(3,j,k) - dFy_3(1,j,k) - dFz_3(1,j,k) + dFvx_3(1,j,k) + dFvy_3(1,j,k) + dFvz_3(1,j,k) )*dt
		Q4_4(1,j,k) = Q4(1,j,k) + 0d0!(- Ex_NSCBCin(4,j,k) - dFy_4(1,j,k) - dFz_4(1,j,k) + dFvx_4(1,j,k) + dFvy_4(1,j,k) + dFvz_4(1,j,k) )*dt
		Q5_4(1,j,k) = Q5(1,j,k) + 0d0!(- Ex_NSCBCin(5,j,k) - dFy_5(1,j,k) - dFz_5(1,j,k) + dFvx_5(1,j,k) + dFvy_5(1,j,k) + dFvz_5(1,j,k) )*dt
	else if(i==Nx .and. j==1) then
		Q1_4(Nx,1,k) = Q1(Nx,1,k) + (- Ex_NSCBCout(1,1,k) - Ey_NSCBCin(1,Nx,k) - dFz_1(Nx,1,k) )*dt
		Q2_4(Nx,1,k) = Q2(Nx,1,k) + 0d0!(- Ex_NSCBCout(2,1,k) - Ey_NSCBCin(2,Nx,k) - dFz_2(Nx,1,k) + dFvx_2(Nx,1,k) + dFvy_2(Nx,1,k) + dFvz_2(Nx,1,k) )*dt
		Q3_4(Nx,1,k) = Q3(Nx,1,k) + 0d0!(- Ex_NSCBCout(3,1,k) - Ey_NSCBCin(3,Nx,k) - dFz_3(Nx,1,k) + dFvx_3(Nx,1,k) + dFvy_3(Nx,1,k) + dFvz_3(Nx,1,k) )*dt
		Q4_4(Nx,1,k) = Q4(Nx,1,k) + 0d0!(- Ex_NSCBCout(4,1,k) - Ey_NSCBCin(4,Nx,k) - dFz_4(Nx,1,k) + dFvx_4(Nx,1,k) + dFvy_4(Nx,1,k) + dFvz_4(Nx,1,k) )*dt
		Q5_4(Nx,1,k) = Q5(Nx,1,k) + 0d0!(- Ex_NSCBCout(5,1,k) - Ey_NSCBCin(5,Nx,k) - dFz_5(Nx,1,k) + dFvx_5(Nx,1,k) + dFvy_5(Nx,1,k) + dFvz_5(Nx,1,k) )*dt
	else if(i==Nx .and. j==Ny) then
		Q1_4(Nx,Ny,k) = Q1(Nx,Ny,k) + (- Ex_NSCBCout(1,Ny,k) - Ey_NSCBCout(1,Nx,k) - dFz_1(Nx,Ny,k) )*dt
		Q2_4(Nx,Ny,k) = Q2(Nx,Ny,k) + (- Ex_NSCBCout(2,Ny,k) - Ey_NSCBCout(2,Nx,k) - dFz_2(Nx,Ny,k) + dFvx_2(Nx,Ny,k) + dFvy_2(Nx,Ny,k) + dFvz_2(Nx,Ny,k) )*dt
		Q3_4(Nx,Ny,k) = Q3(Nx,Ny,k) + (- Ex_NSCBCout(3,Ny,k) - Ey_NSCBCout(3,Nx,k) - dFz_3(Nx,Ny,k) + dFvx_3(Nx,Ny,k) + dFvy_3(Nx,Ny,k) + dFvz_3(Nx,Ny,k) )*dt
		Q4_4(Nx,Ny,k) = Q4(Nx,Ny,k) + (- Ex_NSCBCout(4,Ny,k) - Ey_NSCBCout(4,Nx,k) - dFz_4(Nx,Ny,k) + dFvx_4(Nx,Ny,k) + dFvy_4(Nx,Ny,k) + dFvz_4(Nx,Ny,k) )*dt
		Q5_4(Nx,Ny,k) = Q5(Nx,Ny,k) + (- Ex_NSCBCout(5,Ny,k) - Ey_NSCBCout(5,Nx,k) - dFz_5(Nx,Ny,k) + dFvx_5(Nx,Ny,k) + dFvy_5(Nx,Ny,k) + dFvz_5(Nx,Ny,k) )*dt
	else if(i==Nx .and. j/=1 .and. j/=Ny) then
		Q1_4(Nx,j,k) = Q1(Nx,j,k) + (- Ex_NSCBCout(1,j,k) - dFy_1(Nx,j,k) - dFz_1(Nx,j,k) )*dt
		Q2_4(Nx,j,k) = Q2(Nx,j,k) + (- Ex_NSCBCout(2,j,k) - dFy_2(Nx,j,k) - dFz_2(Nx,j,k) + dFvx_2(Nx,j,k) + dFvy_2(Nx,j,k) + dFvz_2(Nx,j,k) )*dt
		Q3_4(Nx,j,k) = Q3(Nx,j,k) + (- Ex_NSCBCout(3,j,k) - dFy_3(Nx,j,k) - dFz_3(Nx,j,k) + dFvx_3(Nx,j,k) + dFvy_3(Nx,j,k) + dFvz_3(Nx,j,k) )*dt
		Q4_4(Nx,j,k) = Q4(Nx,j,k) + (- Ex_NSCBCout(4,j,k) - dFy_4(Nx,j,k) - dFz_4(Nx,j,k) + dFvx_4(Nx,j,k) + dFvy_4(Nx,j,k) + dFvz_4(Nx,j,k) )*dt
		Q5_4(Nx,j,k) = Q5(Nx,j,k) + (- Ex_NSCBCout(5,j,k) - dFy_5(Nx,j,k) - dFz_5(Nx,j,k) + dFvx_5(Nx,j,k) + dFvy_5(Nx,j,k) + dFvz_5(Nx,j,k) )*dt
	else if(i/=1 .and. i/=Nx .and. j==1) then
		Q1_4(i,1,k) = Q1(i,1,k) + (- dFx_1(i,1,k) - Ey_NSCBCin(1,i,k) - dFz_1(i,1,k) )*dt
		Q2_4(i,1,k) = Q2(i,1,k) + 0d0!(- dFx_2(i,1,k) - Ey_NSCBCin(2,i,k) - dFz_2(i,1,k) + dFvx_2(i,1,k) + dFvy_2(i,1,k) + dFvz_2(i,1,k) )*dt
		Q3_4(i,1,k) = Q3(i,1,k) + 0d0!(- dFx_3(i,1,k) - Ey_NSCBCin(3,i,k) - dFz_3(i,1,k) + dFvx_3(i,1,k) + dFvy_3(i,1,k) + dFvz_3(i,1,k) )*dt
		Q4_4(i,1,k) = Q4(i,1,k) + 0d0!(- dFx_4(i,1,k) - Ey_NSCBCin(4,i,k) - dFz_4(i,1,k) + dFvx_4(i,1,k) + dFvy_4(i,1,k) + dFvz_4(i,1,k) )*dt
		Q5_4(i,1,k) = Q5(i,1,k) + 0d0!(- dFx_5(i,1,k) - Ey_NSCBCin(5,i,k) - dFz_5(i,1,k) + dFvx_5(i,1,k) + dFvy_5(i,1,k) + dFvz_5(i,1,k) )*dt
	else if(i/=1 .and. i/=Nx .and. j==Ny) then
		Q1_4(i,Ny,k) = Q1(i,Ny,k) + (- dFx_1(i,Ny,k) - Ey_NSCBCout(1,i,k) - dFz_1(i,Ny,k) )*dt
		Q2_4(i,Ny,k) = Q2(i,Ny,k) + (- dFx_2(i,Ny,k) - Ey_NSCBCout(2,i,k) - dFz_2(i,Ny,k) + dFvx_2(i,Ny,k) + dFvy_2(i,Ny,k) + dFvz_2(i,Ny,k) )*dt
		Q3_4(i,Ny,k) = Q3(i,Ny,k) + (- dFx_3(i,Ny,k) - Ey_NSCBCout(3,i,k) - dFz_3(i,Ny,k) + dFvx_3(i,Ny,k) + dFvy_3(i,Ny,k) + dFvz_3(i,Ny,k) )*dt
		Q4_4(i,Ny,k) = Q4(i,Ny,k) + (- dFx_4(i,Ny,k) - Ey_NSCBCout(4,i,k) - dFz_4(i,Ny,k) + dFvx_4(i,Ny,k) + dFvy_4(i,Ny,k) + dFvz_4(i,Ny,k) )*dt
		Q5_4(i,Ny,k) = Q5(i,Ny,k) + (- dFx_5(i,Ny,k) - Ey_NSCBCout(5,i,k) - dFz_5(i,Ny,k) + dFvx_5(i,Ny,k) + dFvy_5(i,Ny,k) + dFvz_5(i,Ny,k) )*dt
	else
		Q1_4(i,j,k) = Q1(i,j,k) + (- dFx_1(i,j,k) - dFy_1(i,j,k) - dFz_1(i,j,k) )*dt
		Q2_4(i,j,k) = Q2(i,j,k) + (- dFx_2(i,j,k) - dFy_2(i,j,k) - dFz_2(i,j,k) + dFvx_2(i,j,k) + dFvy_2(i,j,k) + dFvz_2(i,j,k) )*dt
		Q3_4(i,j,k) = Q3(i,j,k) + (- dFx_3(i,j,k) - dFy_3(i,j,k) - dFz_3(i,j,k) + dFvx_3(i,j,k) + dFvy_3(i,j,k) + dFvz_3(i,j,k) )*dt
		Q4_4(i,j,k) = Q4(i,j,k) + (- dFx_4(i,j,k) - dFy_4(i,j,k) - dFz_4(i,j,k) + dFvx_4(i,j,k) + dFvy_4(i,j,k) + dFvz_4(i,j,k) )*dt
		Q5_4(i,j,k) = Q5(i,j,k) + (- dFx_5(i,j,k) - dFy_5(i,j,k) - dFz_5(i,j,k) + dFvx_5(i,j,k) + dFvy_5(i,j,k) + dFvz_5(i,j,k) )*dt
	end if
		Q1(i,j,k) = Q1_4(i,j,k)
		Q2(i,j,k) = Q2_4(i,j,k)
		Q3(i,j,k) = Q3_4(i,j,k)
		Q4(i,j,k) = Q4_4(i,j,k)
		Q5(i,j,k) = Q5_4(i,j,k)

!	Q1(i,j,k) = Q1_next(i,j,k)
!	Q2(i,j,k) = Q2_next(i,j,k)
!	Q3(i,j,k) = Q3_next(i,j,k)
!	Q4(i,j,k) = Q4_next(i,j,k)
!	Q5(i,j,k) = Q5_next(i,j,k)
	end do
	end do
	end do
	!$omp end parallel do
!!!!=============================================================================
	call BL_init()!!流入部
	end subroutine euler
!=====================================================================
!=====================================================================
!!!-----CSVfile output----------------------
	subroutine output_files()
	write(Filename(7:9), '(i3.3)') num
	open(11, file=Filename, status='unknown', form='formatted')
	do k=1,Nz
		do j=1,Ny
			do i=1, Nx
!				write(11, "(11(1x,g15.7))") x(i), y(j), z(k), Q1(i,j,k), Q2(i,j,k)/Q1(i,j,k), &
!						&Q3(i,j,k)/Q1(i,j,k), Q4(i,j,k)/Q1(i,j,k), p(i,j,k), TT(i,j,k), dudy(i,j,k), Q_ten(i,j,k)
				write(11, "(11(1x,g15.7))") x(i), y(j), z(k), Q1(i,j,k), Q2(i,j,k)/Q1(i,j,k), &
						&Q3(i,j,k)/Q1(i,j,k), Q4(i,j,k)/Q1(i,j,k), p(i,j,k), TT(i,j,k), Q_ten(i,j,k)
			end do
!			write(11,*)
		end do
!	write(11,*)
!		write(11,*)
	end do
	close(11)
	end subroutine output_files
!!-------------------------------------------
	subroutine finish()
	open(12, file='error.csv', status='unknown', form='formatted')
		write(12, "(7(1x,g15.7))") 'error' , Nx, Ny, Nz, Nt, dt*timestep, timestep
	close(12)
	end subroutine finish
!!-----Timestep check-----------------------
	subroutine check()
		call itime(TIME_ARRAY)
		write(*,*) "t= ", dt*timestep, timestep, "(*'･ω･)v ＜ｱｳｱｳｱｰ OK",  TIME_ARRAY(1), TIME_ARRAY(2), TIME_ARRAY(3)
!   		write(6,*) TIME_ARRAY(1), TIME_ARRAY(2), TIME_ARRAY(3)
	end subroutine check
!!-----Not a Number judgement
	subroutine NaN_check()
	!$omp parallel do
	do k=1,Nz
	do j=1,Ny
		do i=1, Nx
			if (isnan(Q1(i,j,k))) then
			write(*,"(3(1x,g15.7))") 'rho NaN error',i, j, k
			call output_files()
			call finish()
			stop
			end if
			if (isnan(Q2(i,j,k))) then
			write(*,"(3(1x,g15.7))") 'u NaN error',i, j, k
			call output_files()
			call finish()
			stop
			end if
!			if (isnan(Q3(i,j,k))) then
!			write(*,"(3(1x,g15.7))") 'v NaN error',i, j, k
!			call output_files()
!			call finish()
!			stop
!			end if
			if (isnan(Q5(i,j,k))) then
			write(*,"(3(1x,g15.7))") 'Et NaN error',i, j, k
			call output_files()
			call finish()
			stop
			end if
		end do
	end do
	end do
	!$omp end parallel do
	end subroutine NaN_check
end program BoundaryLayer_3D
