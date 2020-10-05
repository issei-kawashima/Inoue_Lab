!角柱　z周期
!攪乱なしt10 → v=0.0001sin(2pit/9) t90 → w=0.0001sin(2pit/9 - piz/2) t100
! ifort -fast -o tsai.out 3d1.f90
! nohup time ./tsai.out > a.log 2>&1 < /dev/null &
module func
	implicit none
	double precision,parameter:: pi = acos(-1d0)
!!-------------------------------------------------------------------------------------------------------------------------------!!
!!--主要パラメータ---------------------------------------------------------------------------------------------------------------!!
!!-------------------------------------------------------------------------------------------------------------------------------!!
	double precision,parameter:: Re = 100d0				!! レイノルズ数
	double precision,parameter:: dt = 1d0 / 50d0		!! 時間刻み
	double precision,parameter:: xmax =  35d0			!! 最大x
	double precision,parameter:: xmin = -25d0			!! 最小x
	double precision,parameter:: ymax =  25d0			!! 最大y
	double precision,parameter:: ymin = -25d0			!! 最小y
	double precision,parameter:: zmax =   6d0			!! 最大z
	double precision,parameter:: zmin =   0d0			!! 最小z
	integer,parameter:: nx = 180						!! x方向の分割数 (格子点数はn+1)
	integer,parameter:: ny = 200						!! y方向の分割数
	integer,parameter:: nz = 48							!! z方向の分割数
	integer,parameter:: lp = 100 * int( 1d0 / dt)		!! 時間進行回数 (時間進行回数×dt＝無次元時間t)
	integer,parameter:: ln = 100						!! ファイル出力数
	double precision,parameter:: sig = 2.5d-1			!! DCSのσ
	double precision,parameter:: ma = 3d-1				!! マッハ数
	double precision,parameter:: ls = 0d0				!! NSCBC無反射流出条件のσ
!!-------------------------------------------------------------------------------------------------------------------------------!!
!!--定数-------------------------------------------------------------------------------------------------------------------------!!
!!-------------------------------------------------------------------------------------------------------------------------------!!
	double precision,parameter:: gam = 1.4d0			!! 比熱比
	double precision,parameter:: ze = 1d0				!! ラックス・フリードリッヒの流束分割の定数
	double precision,parameter:: Pr = 7.1d-1			!! プラントル数
	double precision,parameter:: s = 12d1				!! サザランド定数
	double precision,parameter:: r0 = 1d0				!! 密度ρ 初期値
	double precision,parameter:: u0 = 1d0				!! x速度u 初期値
	double precision,parameter:: t0 = 1d0				!! 温度T  初期値
!!-------------------------------------------------------------------------------------------------------------------------------!!
!!--プログラム動作制御-----------------------------------------------------------------------------------------------------------!!
!!-------------------------------------------------------------------------------------------------------------------------------!!
	integer,parameter:: sh = 1							!! 最初から計算を始める場合=0、途中から計算を始める場合=1
	integer,parameter:: lnp = 150						!! 途中から計算を始める場合は開始ファイル番号を指定
	integer,parameter:: ry = 0							!! 流入部が任意値流入条件の場合=0、無反射流入条件の場合=1
	integer,parameter:: ks = 1							!! 格子伸長あり=1、なし=0
	integer,parameter:: bf = 1							!! バッファー領域あり=1、なし=0
	double precision,parameter:: ay = 0d0				!! 任意流入時、y速度v = ay * sin( oy * t)
	double precision,parameter:: oy = 2d0/9d0 * pi		!!
	double precision,parameter:: az = 3d-1				!! 任意流入時、z速度w = az * sin( oz * t + kz * z)
	double precision,parameter:: oz = 2d0/9d0 * pi		!!
	double precision,parameter:: kz = - pi / 3d0		!!
!!-------------------------------------------------------------------------------------------------------------------------------!!
!!--従属変数など-----------------------------------------------------------------------------------------------------------------!!
!!-------------------------------------------------------------------------------------------------------------------------------!!
	integer,parameter:: nk = 21 !!微分用配列の数
	double precision,parameter:: dx = ( xmax - xmin) / dble( nx), dy = ( ymax - ymin) / dble( ny)
	double precision,parameter:: dz = ( zmax - zmin) / dble( nz)
	integer,parameter:: xk = nx/(xmax - xmin)*(-xmin) - int( 1d0 / dx) / ( 2 - ks)
	integer,parameter::xkk = nx/(xmax - xmin)*(-xmin) + int( 1d0 / dx) / ( 2 - ks) !! 角柱x方向始点,終点
	integer,parameter:: yk = ny/2 - int( 1d0 / dy) / ( 2 - ks), ykk = ny/2 + int( 1d0 / dy) / ( 2 - ks) !! 角柱y方向始点,終点
	double precision,parameter:: c0 = sqrt( t0) / ma, p0 = r0 / ( gam * ma * ma), tg0 = t0 / ( ma * ma * gam * ( gam - 1d0))
!!-------------------------------------------------------------------------------------------------------------------------------!!
!!--グローバル変数---------------------------------------------------------------------------------------------------------------!!
!!-------------------------------------------------------------------------------------------------------------------------------!!
	integer:: nn( nk+1), nA( nk)
	double precision:: tm = 0d0, sigA( nk), dlA( nk), c( nx+1, ny+1, nz), g( 5, nx+1, ny+1, nz), dg( 3, 5, nx+1, ny+1, nz)
	double precision:: f( 3, 5, nx+1, ny+1, nz), df( 3, 5, nx+1, ny+1, nz), fv( 3, 5, nx+1, ny+1, nz), dfv( 3, 5, nx+1, ny+1, nz)
	double precision,allocatable,dimension(:,:):: lu, lu2
	double precision:: dxx( nx+1), dyy( ny+1), x( nx+1), y( ny+1), z( nz+1), buix( nx+1), bukx( nx+1), buiy( ny+1), buky( ny+1)
!!-------------------------------------------------------------------------------------------------------------------------------!!
contains

function cal( q)
	implicit none
	integer:: i, j, k, m, o
	double precision,intent(in):: q(:,:,:,:)
	double precision,parameter:: tc = 2.9115d2, qtg( 5) = (/ r0, r0 * u0, 0d0, 0d0, p0/( gam-1d0) + r0*u0*u0*5d-1/)
	double precision,parameter:: h1( 2:4) = (/0d0,1d0,1d0/)
	integer:: an( 3), ai( 3), aj( 3), ar( 3), kr, pm, fba, fbb, dr, r1, r2
	double precision:: cal( 5, nx+1, ny+1, nz), uv, uw, d( 5), bc( 7, 4)
	double precision:: t( nx+1, ny+1, nz), dtt( 3, nx+1, ny+1, nz), myu( nx+1, ny+1, nz)
	double precision:: da( 3,10, nx+1, ny+1, nz), aa(10, nx+1, ny+1, nz), dq( 3, 5, nx+1, ny+1, nz)

	!g,f,t組み立て
	do k = 1, nz
		do j = 1, ny+1
			do i = 1, nx+1
				g( 1,i,j,k) = q( 1,i,j,k)
				g( 2,i,j,k) = q( 2,i,j,k) / q( 1,i,j,k)
				g( 3,i,j,k) = q( 3,i,j,k) / q( 1,i,j,k)
				g( 4,i,j,k) = q( 4,i,j,k) / q( 1,i,j,k)
				g( 5,i,j,k) = ( gam - 1d0) * ( q( 5,i,j,k)- 5d-1) / &
				q( 1,i,j,k) * ( q( 2,i,j,k) * q( 2,i,j,k) + q( 3,i,j,k) * q( 3,i,j,k) + q( 4,i,j,k) * q( 4,i,j,k))
			enddo
		enddo
	enddo
	do k = 1, nz
		do j = 1, ny+1
			do i = 1, nx+1
				f( 1, 1,i,j,k) = q( 2,i,j,k)
				f( 2, 1,i,j,k) = q( 3,i,j,k)
				f( 3, 1,i,j,k) = q( 4,i,j,k)
				f( 1, 2,i,j,k) = q( 2,i,j,k) * g( 2,i,j,k) + g( 5,i,j,k)
				f( 2, 2,i,j,k) = q( 2,i,j,k) * g( 3,i,j,k)
				f( 3, 2,i,j,k) = q( 2,i,j,k) * g( 4,i,j,k)
				f( 1, 3,i,j,k) = f( 2, 2,i,j,k)
				f( 2, 3,i,j,k) = q( 3,i,j,k) * g( 3,i,j,k) + g( 5,i,j,k)
				f( 3, 3,i,j,k) = q( 3,i,j,k) * g( 4,i,j,k)
				f( 1, 4,i,j,k) = f( 3, 2,i,j,k)
				f( 2, 4,i,j,k) = f( 3, 3,i,j,k)
				f( 3, 4,i,j,k) = q( 4,i,j,k) * g( 4,i,j,k) + g( 5,i,j,k)
				f( 1, 5,i,j,k) = ( q( 5,i,j,k) + g( 5,i,j,k)) * g( 2,i,j,k)
				f( 2, 5,i,j,k) = ( q( 5,i,j,k) + g( 5,i,j,k)) * g( 3,i,j,k)
				f( 3, 5,i,j,k) = ( q( 5,i,j,k) + g( 5,i,j,k)) * g( 4,i,j,k)
			enddo
		enddo
	enddo
	do k = 1, nz
		do j = 1, ny+1
			do i = 1, nx+1
				t(i,j,k) = gam * g( 5,i,j,k) / g( 1,i,j,k) * ma * ma
			enddo
		enddo
	enddo

	!等温条件設定
	t( 1, :, :) = t( 1,:,:) * dble( ry) + t0 * dble( 1-ry)
	t(  xk+1, yk+1:ykk+1, :) = t0
	t( xkk+1, yk+1:ykk+1, :) = t0
	t( xk+2:xkk,  yk+1, :) = t0
	t( xk+2:xkk, ykk+1, :) = t0

	!c,μの組み立て
	do k = 1, nz
		do j = 1, ny+1
			do i = 1, nx+1
				c(i,j,k) = sqrt( t(i,j,k)) / ma
				myu(i,j,k) = ( t(i,j,k) * c(i,j,k) * ma) * ( tc + s) / ( t(i,j,k) * tc + s)
			enddo
		enddo
	enddo

	!g,t,qの各方向微分
	aa(1:5 ,:,:,:) = g
	aa(  6 ,:,:,:) = t
	aa(7:10,:,:,:) = q(2:5,:,:,:)
	call bbn( da, aa)
	dg              = da(:,1:5 ,:,:,:)
	dtt             = da(:,  6 ,:,:,:)
	dq(:,  1,:,:,:) = da(:,  1 ,:,:,:)
	dq(:,2:5,:,:,:) = da(:,7:10,:,:,:)

	!fvの組み立て
	do k = 1, nz
		do j = 1, ny+1
			do i = 1, nx+1
				d( 1) = 2d0 / 3d0 * ( 2d0 * dg( 1, 2,i,j,k) - dg( 2, 3,i,j,k) - dg( 3, 4,i,j,k))
				d( 2) = 2d0 / 3d0 * ( 2d0 * dg( 2, 3,i,j,k) - dg( 3, 4,i,j,k) - dg( 1, 2,i,j,k))
				d( 3) = 2d0 / 3d0 * ( 2d0 * dg( 3, 4,i,j,k) - dg( 1, 2,i,j,k) - dg( 2, 3,i,j,k))
				fv( 1, 1,i,j,k) = d( 1) * g( 2,i,j,k)
				fv( 2, 1,i,j,k) = d( 2) * g( 3,i,j,k)
				fv( 3, 1,i,j,k) = d( 3) * g( 4,i,j,k)
				fv( 1, 2,i,j,k) = d( 1)
				fv( 2, 2,i,j,k) = dg( 1, 3,i,j,k) + dg( 2, 2,i,j,k)
				fv( 3, 2,i,j,k) = dg( 1, 4,i,j,k) + dg( 3, 2,i,j,k)
				fv( 1, 3,i,j,k) = fv( 2, 2,i,j,k)
				fv( 2, 3,i,j,k) = d( 2)
				fv( 3, 3,i,j,k) = dg( 2, 4,i,j,k) + dg( 3, 3,i,j,k)
				fv( 1, 4,i,j,k) = fv( 3, 2,i,j,k)
				fv( 2, 4,i,j,k) = fv( 3, 3,i,j,k)
				fv( 3, 4,i,j,k) = d( 3)
				fv( 1, 5,i,j,k) = fv( 1, 2,i,j,k) * g( 2,i,j,k) + fv( 1, 3,i,j,k) * g( 3,i,j,k) + fv( 1, 4,i,j,k) * g( 4,i,j,k)&
				+ dtt( 1,i,j,k) / ( ( gam - 1d0) * Pr * ma * ma)
				fv( 2, 5,i,j,k) = fv( 2, 2,i,j,k) * g( 2,i,j,k) + fv( 2, 3,i,j,k) * g( 3,i,j,k) + fv( 2, 4,i,j,k) * g( 4,i,j,k)&
				+ dtt( 2,i,j,k) / ( ( gam - 1d0) * Pr * ma * ma)
				fv( 3, 5,i,j,k) = fv( 3, 2,i,j,k) * g( 2,i,j,k) + fv( 3, 3,i,j,k) * g( 3,i,j,k) + fv( 3, 4,i,j,k) * g( 4,i,j,k)&
				+ dtt( 3,i,j,k) / ( ( gam - 1d0) * Pr * ma * ma)
			enddo
		enddo
	enddo
	do k = 1, nz
		do j = 1, ny+1
			do i = 1, nx+1
				do m = 1, 5
					do o = 1, 3
						fv( o,m,i,j,k) = fv( o,m,i,j,k) * myu(i,j,k) / Re
					enddo
				enddo
			enddo
		enddo
	enddo

	!移流項の微分
	do i = 1, 3
		f( i,:,:,:,:) = f( i,:,:,:,:) - ze * q
	enddo
	call bbn2( da(:,1:5,:,:,:), f, 1)
	do i = 1, 3
		f( i,:,:,:,:) = f( i,:,:,:,:) + 2d0 * ze * q
	enddo
	call bbn2( df, f, 2)
	df = 5d-1 * ( df + da(:,1:5,:,:,:))

	!粘性項の微分
	call bbn2( dfv, fv, 3)

	!x入口　流入一定
	df( 1, 1, 1,:,:) = dble( ry)   * df( 1, 1, 1,:,:)&
				   + dble( 1-ry) * ( u0 - c0) * ( dg( 1, 5, 1,:,:) / c0 - g( 1, 1,:,:) * dg( 1, 2, 1,:,:)) / c0 * gam

	!流出
	bc( 1, 1) = 1	 ; bc( 2, 1) = 1	; bc( 3, 1) = ny+1 ; bc( 4, 1) = 1	; bc( 5, 1) = nz   ; bc( 6, 1) = 1 ; bc( 7, 1) = 1
	bc( 1, 2) = nx+1 ; bc( 2, 2) = 1	; bc( 3, 2) = ny+1 ; bc( 4, 2) = 1	; bc( 5, 2) = nz   ; bc( 6, 2) = 2 ; bc( 7, 2) = 1
	bc( 1, 3) = 1	 ; bc( 2, 3) = 1	; bc( 3, 3) = nz   ; bc( 4, 3) = 1	; bc( 5, 3) = nx+1 ; bc( 6, 3) = 1 ; bc( 7, 3) = 2
	bc( 1, 4) = ny+1 ; bc( 2, 4) = 1	; bc( 3, 4) = nz   ; bc( 4, 4) = 1	; bc( 5, 4) = nx+1 ; bc( 6, 4) = 2 ; bc( 7, 4) = 2
	!場所			 ; 範囲1開始		; 範囲1終了		   ; 範囲2開始		; 範囲2終了		   ; 壁面方向1→2  ; 座標x1y2z3
	do k = 2-ry, 4
		dr = bc( 7, k) ; r1 = mod( dr, 3) + 1 ; r2 = mod( dr+1, 3) + 1
		pm = dble( 2 * bc( 6, k) - 3)
		fba =  4 * bc( 6, k) - 3
		fbb = -4 * bc( 6, k) + 9
		ar = 0; an = 0; ai = 0; aj = 0
		an( dr) = bc( 1, k)
		ai( r1) = 1
		aj( r2) = 1
		do j = bc( 4, k), bc( 5, k)
			do i = bc( 2, k), bc( 3, k)
				ar = an + ai * i + aj * j
				d(	fba) = ( g( dr+1, ar(1), ar(2), ar(3)) + c( ar(1), ar(2), ar(3)) * pm) * ( dg( dr, 5, ar(1), ar(2), ar(3))&
						+ g( 1, ar(1), ar(2), ar(3)) * c( ar(1), ar(2), ar(3)) * dg( dr, dr+1, ar(1), ar(2), ar(3)) * pm)
				d(	fbb) = 0d0
				d( dr+1) = g( dr+1, ar(1), ar(2), ar(3)) * ( c( ar(1), ar(2), ar(3)) * c( ar(1), ar(2), ar(3))&
				 * dg( dr, 1, ar(1), ar(2), ar(3)) - dg( dr, 5, ar(1), ar(2), ar(3)))
				d( r1+1) = g( dr+1, ar(1), ar(2), ar(3)) * dg( dr, r1+1, ar(1), ar(2), ar(3))
				d( r2+1) = g( dr+1, ar(1), ar(2), ar(3)) * dg( dr, r2+1, ar(1), ar(2), ar(3))
				d( 2) = d( 2) * h1( bc( 6, k) + dr)
				d( 3) = d( 3) * h1( bc( 6, k) + dr)
				d( 4) = d( 4) * h1( bc( 6, k) + dr)
				kr 	 = ( d(    5) - d( 1)) / ( c( ar(1), ar(2), ar(3)) * g( 1, ar(1), ar(2), ar(3))) * 5d-1
				d(    5) = ( d(    1) + d( 5)) * 5d-1
				d(    1) = ( d( dr+1) + d( 5)) / ( c( ar(1), ar(2), ar(3)) * c( ar(1), ar(2), ar(3)))
				d( dr+1) = kr
				df( dr, 1, ar(1), ar(2), ar(3)) = d( 1)
				df( dr, 2:4, ar(1), ar(2), ar(3)) = g( 2:4, ar(1), ar(2), ar(3)) * d( 1) + g( 1, ar(1), ar(2), ar(3)) * d( 2:4)
				df( dr, 5, ar(1), ar(2), ar(3))&
				= ( g( 2, ar(1), ar(2), ar(3)) **2 + g( 3, ar(1), ar(2), ar(3)) **2 + g( 4, ar(1), ar(2), ar(3)) **2) * d( 1)*5d-1&
				+ d( 5) / ( gam - 1d0)&
				+ g( 1, ar(1), ar(2), ar(3)) * ( g( 2, ar(1), ar(2), ar(3)) * d( 2)&
											   + g( 3, ar(1), ar(2), ar(3)) * d( 3)&
											   + g( 4, ar(1), ar(2), ar(3)) * d( 4))
				dfv( dr, r1+1, ar(1), ar(2), ar(3)) = 0d0
				dfv( dr, r2+1, ar(1), ar(2), ar(3)) = 0d0
				dfv( dr, 5, ar(1), ar(2), ar(3)) = dfv( dr,    1, ar(1), ar(2), ar(3))&
												  + fv( dr, r1+1, ar(1), ar(2), ar(3)) * dg( dr, r1+1, ar(1), ar(2), ar(3))&
												  + fv( dr, r2+1, ar(1), ar(2), ar(3)) * dg( dr, r2+1, ar(1), ar(2), ar(3))
			enddo
		enddo
	enddo

	!等温無滑り壁面
	bc( 1, 1) = xk+1 ; bc( 2, 1) = yk+1 ; bc( 3, 1) = ykk+1 ; bc( 4, 1) =    1 ; bc( 5, 1) =    nz ; bc( 6, 1) = 2 ; bc( 7, 1) = 1
	bc( 1, 2) = xkk+1; bc( 2, 2) = yk+1 ; bc( 3, 2) = ykk+1 ; bc( 4, 2) =    1 ; bc( 5, 2) =    nz ; bc( 6, 2) = 1 ; bc( 7, 2) = 1
	bc( 1, 3) = yk+1 ; bc( 2, 3) =    1 ; bc( 3, 3) =    nz ; bc( 4, 3) = xk+1 ; bc( 5, 3) = xkk+1 ; bc( 6, 3) = 2 ; bc( 7, 3) = 2
	bc( 1, 4) = ykk+1; bc( 2, 4) =    1 ; bc( 3, 4) =    nz ; bc( 4, 4) = xk+1 ; bc( 5, 4) = xkk+1 ; bc( 6, 4) = 1 ; bc( 7, 4) = 2
	!場所			 ; 範囲1開始 		; 範囲1終了			; 範囲2開始		   ; 範囲2終了		   ; 壁面方向1→2  ; 座標x1y2z3
	do k = 1, 4
		dr = bc( 7, k)
		pm = dble( 2 * bc( 6, k) - 3)
		ar = 0; an = 0; ai = 0; aj = 0
		an( dr) = bc( 1, k)
		ai( mod( dr, 3) + 1) = 1
		aj( mod( dr+1, 3) + 1) = 1
		do j = bc( 4, k), bc( 5, k)
			do i = bc( 2, k), bc( 3, k)
				ar = an + ai * i + aj * j
				df( dr, 1, ar(1), ar(2), ar(3)) = dg( dr,    5, ar(1), ar(2), ar(3)) / c0 * pm&
											    + dg( dr, dr+1, ar(1), ar(2), ar(3)) * g( 1, ar(1), ar(2), ar(3))
			enddo
		enddo
	enddo

	!流出のための値一時退避撤去
	dfv(:, 1,:,:,:) = 0d0

	!バッファー領域の処理
	do k = 1, nz
		do j = 1, ny+1
			do i = 1, nx+1
				 df( 1,:, i, j, k) =  df( 1,:, i, j, k) + buix( i) * dq( 1,:, i, j, k) + buiy( j) * dq( 2,:, i, j, k)
				dfv( 1,:, i, j, k) = dfv( 1,:, i, j, k) - (bukx( i) + buky( j)) * ( q(:, i, j, k) - qtg(:))
			enddo
		enddo
	enddo

	!時間進行
	do i = 1, nz
		do j = 1, ny+1
			do k = 1, nx+1
				cal(:, k, j, i) = q(:, k, j, i) + ( dfv( 1,:, k, j, i) + dfv( 2,:, k, j, i) + dfv( 3,:, k, j, i)&
												   - df( 1,:, k, j, i)  - df( 2,:, k, j, i)  - df( 3,:, k, j, i)) * dt
			enddo
		enddo
	enddo

	!入口流入
	uv = ay * sin( oy * ( tm + dt))
	do k = 1, nz
		do j = 1, ny+1
			uw = exp( -1d-1 * y( j) * y( j)) * az * sin( oz * ( tm + dt) + kz * z( k))
			cal( 2, 1, j, k) = dble( ry) * cal( 2, 1, j, k) + dble( 1-ry) * cal( 1, 1, j, k) * u0
			cal( 3, 1, j, k) = dble( ry) * cal( 3, 1, j, k) + dble( 1-ry) * cal( 1, 1, j, k) * uv
			cal( 4, 1, j, k) = dble( ry) * cal( 4, 1, j, k) + dble( 1-ry) * cal( 1, 1, j, k) * uw
			cal( 5, 1, j, k) = dble( ry) * cal( 5, 1, j, k) + dble( 1-ry) * cal( 1, 1, j, k) * ( tg0 + 5d-1*(u0**2+uv**2+uw**2))
		enddo
	enddo

	!等温無滑り壁面
	do k = 1, 4
		dr = bc( 7, k)
		ar = 0; an = 0; ai = 0; aj = 0
		an( dr) = bc( 1, k)
		ai( mod( dr, 3) + 1) = 1
		aj( mod( dr+1, 3) + 1) = 1
		do j = bc( 4, k), bc( 5, k)
			do i = bc( 2, k), bc( 3, k)
				ar = an + ai * i + aj * j
				cal( 2:4, ar(1), ar(2), ar(3)) = 0d0
				cal( 5, ar(1), ar(2), ar(3)) = cal( 1, ar(1), ar(2), ar(3)) * tg0
			enddo
		enddo
	enddo

end function

function nxcs( fn, num)
	implicit none
	integer:: i, j, k, m, n, num, fs1, fs2, fs3, nm
	double precision:: si, fn(:,:,:,:), dl
	double precision,allocatable,dimension(:,:,:,:):: nxcs

	si = sigA( num); n = nA( num); dl = dlA( num); nm = nn( num); fs1 = size( fn, 3); fs2 = size( fn, 4); fs3 = size( fn, 1)
	allocate( nxcs( fs3, n+1, fs1, fs2))
	do m = 1, fs3
	do k = 1, fs2
		do j = 1, fs1
			nxcs( m, 1, j, k) = ( ( -1.7d1 * fn( m, 1, j, k) - fn( m, 4, j, k)) / 6d0 + 1.5d0 * ( fn( m, 2, j, k) + fn( m, 3, j, k))) / dl
			nxcs( m, 2, j, k) = ( ( si + 1.5d0) * fn( m, 3, j, k) - si * 2d0 * fn( m, 2, j, k) + ( si - 1.5d0) * fn( m, 1, j, k)) / dl * 5d-1
			do i = 3, n-1
			    nxcs( m, i, j, k) = ( ( 1.25d-1 * si + 2.5d-1) * fn( m, i+2, j, k) + ( 4d0 * si + 7d0) * fn( m, i+1, j, k)&
			    - 8.25d0*si*fn( m, i, j, k) + ( 4d0 * si - 7d0)*fn( m, i-1, j, k) + ( 1.25d-1 * si - 2.5d-1)*fn( m, i-2, j, k))/9d0/dl
			enddo
			nxcs( m, n, j, k) = ( ( si + 1.5d0) * fn( m, n+1, j, k) - si * 2d0 * fn( m, n, j, k) + ( si - 1.5d0) * fn( m, n-1, j, k))/dl*5d-1
			nxcs( m, n+1, j, k) = ( ( 1.7d1 * fn( m, n+1, j, k)+ fn( m, n-2, j, k)) / 6d0 - 1.5d0 * ( fn( m, n, j, k) + fn( m, n-1, j, k)) ) / dl
			do i = 2, n+1
			    nxcs( m, i, j, k) = nxcs( m, i, j, k) - lu( 1, i+nm) * nxcs( m, i-1, j, k)
			enddo
			nxcs( m, n+1, j, k) = nxcs( m, n+1, j, k) / lu( 2, n+1+nm)
			do i = n, 1, -1
			    nxcs( m, i, j, k) = ( nxcs( m, i, j, k) - nxcs( m, i+1, j, k) * lu( 3, i+nm)) / lu( 2, i+nm)
			enddo
		enddo
	enddo
	enddo
end function

function nycs( fn, num)
	implicit none
	integer:: i, j, k, m, n, num, fs1, fs2, fs3, nm
	double precision:: si, fn(:,:,:,:), dl
	double precision,allocatable,dimension(:,:,:,:):: nycs

	si = sigA( num); n = nA( num); dl = dlA( num); nm = nn( num); fs1 = size( fn, 2); fs2 = size( fn, 4); fs3 = size( fn, 1)
	allocate( nycs( fs3, fs1, n+1, fs2))
	do m = 1, fs3
	do k = 1, fs2
		do j = 1, fs1
			nycs( m, j, 1, k) = ( ( -1.7d1 * fn( m, j, 1, k) - fn( m, j, 4, k)) / 6d0 + 1.5d0 * ( fn( m, j, 2, k) + fn( m, j, 3, k))) / dl
			nycs( m, j, 2, k) = ( ( si + 1.5d0) * fn( m, j, 3, k) - si * 2d0 * fn( m, j, 2, k) + ( si - 1.5d0) * fn( m, j, 1, k)) / dl * 5d-1
			do i = 3, n-1
			    nycs( m, j, i, k) = ( ( 1.25d-1 * si + 2.5d-1) * fn( m, j, i+2, k) + ( 4d0 * si + 7d0) * fn( m, j, i+1, k)&
			    - 8.25d0*si*fn( m, j, i, k) + ( 4d0 * si - 7d0)*fn( m, j, i-1, k) + ( 1.25d-1 * si - 2.5d-1)*fn( m, j, i-2, k))/9d0/dl
			enddo
			nycs( m, j, n, k) = ( ( si + 1.5d0) * fn( m, j, n+1, k) - si * 2d0 * fn( m, j, n, k) + ( si - 1.5d0) * fn( m, j, n-1, k))/dl*5d-1
			nycs( m, j, n+1, k) = ( ( 1.7d1 * fn( m, j, n+1, k)+ fn( m, j, n-2, k)) / 6d0 - 1.5d0 * ( fn( m, j, n, k) + fn( m, j, n-1, k)) ) / dl
			do i = 2, n+1
			    nycs( m, j, i, k) = nycs( m, j, i, k) - lu( 1, i+nm) * nycs( m, j, i-1, k)
			enddo
			nycs( m, j, n+1, k) = nycs( m, j, n+1, k) / lu( 2, n+1+nm)
			do i = n, 1, -1
			    nycs( m, j, i, k) = ( nycs( m, j, i, k) - nycs( m, j, i+1, k) * lu( 3, i+nm)) / lu( 2, i+nm)
			enddo
		enddo
	enddo
	enddo
end function

function nzcs( fn, num)
	implicit none
	integer:: i, j, k, m, n, num, fs1, fs2, fs3, nm
	double precision:: si, fn(:,:,:,:), dl
	double precision,allocatable,dimension(:,:,:,:):: nzcs

	si = sigA( num); n = nA( num); dl = dlA( num); nm = nn( num); fs1 = size( fn, 2); fs2 = size( fn, 3); fs3 = size( fn, 1)
	allocate( nzcs( fs3, fs1, fs2, n))
	do m = 1, fs3
	do k = 1, fs2
		do j = 1, fs1
			do i = 1, n
			    nzcs( m, j, k, i) = ( ( 1.25d-1 * si + 2.5d-1) * fn( m, j, k, mod( i+1  , n)+1)&
			    				 + ( 4d0     * si + 7d0   ) * fn( m, j, k, mod( i    , n)+1)&
								    - 8.25d0 * si           * fn( m, j, k,      i          )&
			    				 + ( 4d0     * si - 7d0   ) * fn( m, j, k, mod( i-2+n, n)+1)&
			    				 + ( 1.25d-1 * si - 2.5d-1) * fn( m, j, k, mod( i-3+n, n)+1)) / 9d0 / dl
			enddo
			do i = 2, n-1
			    nzcs( m, j, k, i) = nzcs( m, j, k, i) - lu2( 2, i+nm) * nzcs( m, j, k, i-1)
			enddo
			do i = 1, n-1
				nzcs( m, j, k, n) = nzcs( m, j, k, n) - lu2( 1, i+nm) * nzcs( m, j, k, i)
			enddo
			nzcs( m, j, k, n) = nzcs( m, j, k, n) / lu2( 3, n+nm)
			nzcs( m, j, k, n-1) = ( nzcs( m, j, k, n-1) - nzcs( m, j, k, n) * lu2( 4, n-1+nm)) / lu2( 3, n-1+nm)
			do i = n-2, 1, -1
			    nzcs( m, j, k, i) = ( nzcs( m, j, k, i) - nzcs( m, j, k, i+1) * lu2( 4, i+nm) - nzcs( m, j, k, n) * lu2( 5, i+nm))&
																											/ lu2( 3, i+nm)
			enddo
		enddo
	enddo
	enddo
end function

subroutine NLU( num)       ! LU分解
    implicit none
    integer:: h, k, i, j, n, num
    double precision:: si
    double precision,allocatable,dimension(:,:):: a

    si = sigA( num); n = nA( num)
    allocate( a( 3, n+1))
    a( 1, 1) = 0d0
    a( 2, 1) = 1d0
    a( 3, 1) = 3d0
    a( 1, 2) = 2.5d-1 * ( 1d0 - si)
	a( 2, 2) = 1d0
    a( 3, 2) = 2.5d-1 * ( 1d0 + si)
    do i = 3, n-1
        a( 1, i) = ( 1d0 - si) / 3d0
        a( 2, i) = 1d0
        a( 3, i) = ( 1d0 + si) / 3d0
    enddo
    a( 1, n) = 2.5d-1 * ( 1d0 - si)
	a( 2, n) = 1d0
    a( 3, n) = 2.5d-1 * ( 1d0 + si)
    a( 1, n+1) = 3d0
    a( 2, n+1) = 1d0
    a( 3, n+1) = 0d0
    do h = 2, n+1
    	a( 1, h) = a( 1, h) / a( 2, h-1)
		a( 2, h) = a( 2, h) - a( 1, h) * a( 3, h-1)
		a( 3, h) = a( 3, h)
	enddo
	lu( :, nn( num)+1 : nn( num+1)) = a(:,:)
	deallocate( a)
end subroutine

subroutine NLU2( num)       ! LU分解 周期条件
    implicit none
    integer:: h, k, i, j, n, num
    double precision:: si
    double precision,allocatable,dimension(:,:):: a, l, u, alu

    si = sigA( num); n = nA( num)
    allocate( a( n, n), l( n, n), u( n, n), alu( 5, n))
    a=0d0; l=0d0; u=0d0; alu=0d0;
    a( 1, 1) = 1d0
    a( 2, 1) = ( 1d0 - si) / 3d0
    a( n, 1) = ( 1d0 + si) / 3d0
    do i = 2, n-1
        a( i - 1, i) = ( 1d0 + si) / 3d0
        a( i    , i) = 1d0
        a( i + 1, i) = ( 1d0 - si) / 3d0
    enddo
    a( 1  , n) = ( 1d0 - si) / 3d0
    a( n-1, n) = ( 1d0 + si) / 3d0
    a( n  , n) = 1d0
    do h = 1, n
        do j = h, n
            u( h, j) = a( h, j)
            do k = 1, h-1
                u(  h, j) = u( h, j) - l( h, k) * u( k, j)
            enddo
        enddo
        do i = h, n
            l( i, h) = a( i, h)
            do k = 1, h-1
                l( i, h) = l( i, h) - l( i, k) * u( k, h)
            enddo
            l( i, h) = l( i, h) / u( h, h)
        enddo
    enddo
!    write(99,*)"L"
!    do i = 1, n
!    write(99,"(8f6.2)")l( i, 1:n)
!    enddo
!    write(99,*)"U"
!    do i = 1, n
!    write(99,"(8f6.2)")u( i, 1:n)
!    enddo
!    a=0d0
!    do i = 1, n
!    	do j = 1, n
!    		do k = 1, n
!    			a( i, j) = a( i, j) + l( i, k) * u( k , j)
!    		enddo
!    	enddo
!    enddo
!    write(99,*)"A"
!    do i = 1, n
!    write(99,"(8f6.2)")a( i, 1:n)
!    enddo
	do i = 1, n
		alu( 1, i) = l( n, i)
	enddo
	alu( 2, 1) = 0d0
	alu( 3, 1) = u( 1, 1)
	alu( 4, 1) = u( 1, 2)
	do i = 2, n-1
		alu( 2, i) = l( i, i-1)
		alu( 3, i) = u( i, i  )
		alu( 4, i) = u( i, i+1)
	enddo
	alu( 2, n) = l( n, n-1)
	alu( 3, n) = u( n, n  )
	alu( 4, n) = 0d0
	do i = 1, n
		alu( 5, i) = u( i, n)
	enddo
!    write(99,*)"ck"
!    do i = 1, n
!    write(99,"(8f6.2)")alu(1:5, i)
!    enddo
    lu2( :, nn( num)+1 : nn( num+1)) = alu(:,:)
end subroutine

!subroutine bbx( da, a)
!    implicit none
!    integer:: i
!    double precision:: a(:,:,:), da(:,:,:)
!	da(      :    ,     1:yk+1,:) = nxcs( a(      :    ,     1:yk+1,:), 3)
!	da(     1:xk+1,  yk+2:ykk ,:) = nxcs( a(     1:xk+1,  yk+2:ykk ,:), 6)
!	da( xkk+1:nx+1,  yk+2:ykk ,:) = nxcs( a( xkk+1:nx+1,  yk+2:ykk ,:), 9)
!	da(      :    , ykk+1:ny+1,:) = nxcs( a(      :    , ykk+1:ny+1,:), 3)
!	do i = 1, nx+1
!		da(  i,:,:) = da(  i,:,:) * dxx( i)
!	enddo
!end subroutine
!
!subroutine bby( da, a)
!    implicit none
!    integer:: i
!    double precision:: a(:,:,:), da(:,:,:)
!	da(    1:xk+1 ,      :    ,:) = nycs( a(    1:xk+1 ,      :    ,:),12)
!	da(  xk+2:xkk ,     1:yk+1,:) = nycs( a(  xk+2:xkk ,     1:yk+1,:),15)
!	da(  xk+2:xkk , ykk+1:ny+1,:) = nycs( a(  xk+2:xkk , ykk+1:ny+1,:),18)
!	da( xkk+1:nx+1,      :    ,:) = nycs( a( xkk+1:nx+1,      :    ,:),12)
!	do i = 1, ny+1
!		da(:, i,:) = da(:, i,:) * dyy( i)
!	enddo
!end subroutine
!
!subroutine bbz( da, a)
!    implicit none
!    double precision:: a(:,:,:), da(:,:,:)
!	da(    1:xk+1 ,      :    ,:) = nzcs( a(    1:xk+1 ,      :    ,:),21)
!	da(  xk+2:xkk ,     1:yk+1,:) = nzcs( a(  xk+2:xkk ,     1:yk+1,:),21)
!	da(  xk+2:xkk , ykk+1:ny+1,:) = nzcs( a(  xk+2:xkk , ykk+1:ny+1,:),21)
!	da( xkk+1:nx+1,      :    ,:) = nzcs( a( xkk+1:nx+1,      :    ,:),21)
!end subroutine

subroutine bbn( da, a)
	implicit none
	integer:: i
	double precision,intent(in):: a(:,:,:,:)
	double precision,intent(out):: da(:,:,:,:,:)
	da( 1,:,      :    ,     1:yk+1,:) = nxcs( a(:,      :    ,     1:yk+1,:), 3)
	da( 1,:,     1:xk+1,  yk+2:ykk ,:) = nxcs( a(:,     1:xk+1,  yk+2:ykk ,:), 6)
	da( 1,:, xkk+1:nx+1,  yk+2:ykk ,:) = nxcs( a(:, xkk+1:nx+1,  yk+2:ykk ,:), 9)
	da( 1,:,      :    , ykk+1:ny+1,:) = nxcs( a(:,      :    , ykk+1:ny+1,:), 3)
	da( 2,:,    1:xk+1 ,      :    ,:) = nycs( a(:,    1:xk+1 ,      :    ,:),12)
	da( 2,:,  xk+2:xkk ,     1:yk+1,:) = nycs( a(:,  xk+2:xkk ,     1:yk+1,:),15)
	da( 2,:,  xk+2:xkk , ykk+1:ny+1,:) = nycs( a(:,  xk+2:xkk , ykk+1:ny+1,:),18)
	da( 2,:, xkk+1:nx+1,      :    ,:) = nycs( a(:, xkk+1:nx+1,      :    ,:),12)
	da( 3,:,    1:xk+1 ,      :    ,:) = nzcs( a(:,    1:xk+1 ,      :    ,:),21)
	da( 3,:,  xk+2:xkk ,     1:yk+1,:) = nzcs( a(:,  xk+2:xkk ,     1:yk+1,:),21)
	da( 3,:,  xk+2:xkk , ykk+1:ny+1,:) = nzcs( a(:,  xk+2:xkk , ykk+1:ny+1,:),21)
	da( 3,:, xkk+1:nx+1,      :    ,:) = nzcs( a(:, xkk+1:nx+1,      :    ,:),21)
	do i = 1, nx+1
		da( 1, :, i,:,:) = da( 1, :, i,:,:) * dxx( i)
	enddo
	do i = 1, ny+1
		da( 2, :,:, i,:) = da( 2, :,:, i,:) * dyy( i)
	enddo
end subroutine

subroutine bbn2( da, a, nm)
	implicit none
	integer:: i
	integer,intent(in):: nm
	double precision,intent(in):: a(:,:,:,:,:)
	double precision,intent(out):: da(:,:,:,:,:)
	da( 1,:,      :    ,     1:yk+1,:) = nxcs( a( 1,:,      :    ,     1:yk+1,:), nm   )
	da( 1,:,     1:xk+1,  yk+2:ykk ,:) = nxcs( a( 1,:,     1:xk+1,  yk+2:ykk ,:), nm+3 )
	da( 1,:, xkk+1:nx+1,  yk+2:ykk ,:) = nxcs( a( 1,:, xkk+1:nx+1,  yk+2:ykk ,:), nm+6 )
	da( 1,:,      :    , ykk+1:ny+1,:) = nxcs( a( 1,:,      :    , ykk+1:ny+1,:), nm   )
	da( 2,:,    1:xk+1 ,      :    ,:) = nycs( a( 2,:,    1:xk+1 ,      :    ,:), nm+9 )
	da( 2,:,  xk+2:xkk ,     1:yk+1,:) = nycs( a( 2,:,  xk+2:xkk ,     1:yk+1,:), nm+12)
	da( 2,:,  xk+2:xkk , ykk+1:ny+1,:) = nycs( a( 2,:,  xk+2:xkk , ykk+1:ny+1,:), nm+15)
	da( 2,:, xkk+1:nx+1,      :    ,:) = nycs( a( 2,:, xkk+1:nx+1,      :    ,:), nm+9 )
	da( 3,:,    1:xk+1 ,      :    ,:) = nzcs( a( 3,:,    1:xk+1 ,      :    ,:), nm+18)
	da( 3,:,  xk+2:xkk ,     1:yk+1,:) = nzcs( a( 3,:,  xk+2:xkk ,     1:yk+1,:), nm+18)
	da( 3,:,  xk+2:xkk , ykk+1:ny+1,:) = nzcs( a( 3,:,  xk+2:xkk , ykk+1:ny+1,:), nm+18)
	da( 3,:, xkk+1:nx+1,      :    ,:) = nzcs( a( 3,:, xkk+1:nx+1,      :    ,:), nm+18)
	do i = 1, nx+1
		da( 1, :, i,:,:) = da( 1, :, i,:,:) * dxx( i)
	enddo
	do i = 1, ny+1
		da( 2, :,:, i,:) = da( 2, :,:, i,:) * dyy( i)
	enddo
end subroutine

!nscbc無滑り壁面条件　fb1→2　dr1x2y3z　k→kk
subroutine nsk( n, k1, kk1, k2, kk2, fb, dr)
	implicit none
	integer:: i, j, k, n, k1, kk1, k2, kk2, dr, ar( 3), an( 3), ai( 3), aj( 3), fb, pm

	pm = 2d0 * dble( fb) - 3d0
	ar = 0
	an = 0
	ai = 0
	aj = 0
	an( dr) = n
	ai( mod( dr, 3) + 1) = 1
	aj( mod( dr+1, 3) + 1) = 1
	do i = k1, kk1
		do j = k2, kk2
			ar = an + ai * i + aj * j
			df( dr, 1, ar(1), ar(2), ar(3)) = dg( dr,    5, ar(1), ar(2), ar(3)) / c0 * pm&
											+ dg( dr, dr+1, ar(1), ar(2), ar(3)) * g( 1, ar(1), ar(2), ar(3))
		enddo
	enddo
end subroutine

!nscbc無滑り壁面条件　dr1x2y3z　k→kk
subroutine nskk( n, k1, kk1, k2, kk2, cal, dr)
	implicit none
	integer:: i, j, k, n, k1, kk1, k2, kk2, dr, ar( 3), an( 3), ai( 3), aj( 3)
	double precision :: cal(:,:,:,:)
	ar = 0
	an = 0
	ai = 0
	aj = 0
	an( dr) = n
	ai( mod( dr, 3) + 1) = 1
	aj( mod( dr+1, 3) + 1) = 1
	do i = k1, kk1
		do j = k2, kk2
			ar = an + ai * i + aj * j
			do k = 2, 4
				cal( k, ar(1), ar(2), ar(3)) = 0d0
			enddo
			cal( 5, ar(1), ar(2), ar(3)) = cal( 1, ar(1), ar(2), ar(3)) * tg0
		enddo
	enddo
end subroutine

!nscbc流出条件　fb1→2　dr1x2y3z　k→kk
subroutine nsr( n, k1, kk1, k2, kk2, fb, dr)
	implicit none
	integer:: i, j, n, k1, kk1, k2, kk2, dr, ar( 3), an( 3), ai( 3), aj( 3), fb, pm
	double precision:: ll( 5), d( 5), kr
	double precision:: h11( 3) = (/0d0,1d0,1d0/), h12( 3) = (/1d0,0d0,0d0/)
	double precision:: h21( 3) = (/1d0,0d0,1d0/), h22( 3) = (/0d0,1d0,0d0/)
	double precision:: h31( 3) = (/1d0,1d0,0d0/), h32( 3) = (/0d0,0d0,1d0/)
	pm = 2d0 * dble( fb) - 3d0
	ar = 0
	an = 0
	ai = 0
	aj = 0
	an( dr) = n
	ai( mod( dr, 3) + 1) = 1
	aj( mod( dr+1, 3) + 1) = 1
	do i = k1, kk1
		do j = k2, kk2
			ar = an + ai * i + aj * j

			kr = ( g( dr+1, ar(1), ar(2), ar(3)) + c( ar(1), ar(2), ar(3)) * pm) * ( dg( dr, 5, ar(1), ar(2), ar(3))&
			+ g( 1, ar(1), ar(2), ar(3)) * c( ar(1), ar(2), ar(3)) * dg( dr, dr+1, ar(1), ar(2), ar(3)) * pm)
			ll( 1) = ( 2d0 - dble( fb)) * kr
			ll( 5) = ( dble( fb) - 1d0) * kr
			kr = c( ar(1), ar(2), ar(3))*c( ar(1), ar(2), ar(3))*dg( dr, 1, ar(1), ar(2), ar(3))-dg( dr, 5, ar(1), ar(2), ar(3))
			ll( 2) = g( dr+1, ar(1), ar(2), ar(3)) * ( dg( dr, 2, ar(1), ar(2), ar(3)) * h11( dr) + h12( dr) * kr)
			ll( 3) = g( dr+1, ar(1), ar(2), ar(3)) * ( dg( dr, 3, ar(1), ar(2), ar(3)) * h21( dr) + h22( dr) * kr)
			ll( 4) = g( dr+1, ar(1), ar(2), ar(3)) * ( dg( dr, 4, ar(1), ar(2), ar(3)) * h31( dr) + h32( dr) * kr)

			kr    = ( ll( 5) - ll( 1)) / ( c( ar(1), ar(2), ar(3)) * g( 1, ar(1), ar(2), ar(3))) * 5d-1
			d( 2) = ( ll( 1) + ll( 5)) * 5d-1
			d( 1) = ( ll( dr+1) +  d( 2)) / ( c( ar(1), ar(2), ar(3)) * c( ar(1), ar(2), ar(3)))
			d( 3) =   ll( 2) * h11( dr) + h12( dr) * kr
			d( 4) =   ll( 3) * h21( dr) + h22( dr) * kr
			d( 5) =   ll( 4) * h31( dr) + h32( dr) * kr

			df( dr, 1, ar(1), ar(2), ar(3)) = d( 1)
			df( dr, 2, ar(1), ar(2), ar(3)) = g( 2, ar(1), ar(2), ar(3)) * d( 1) + g( 1, ar(1), ar(2), ar(3)) * d( 3)
			df( dr, 3, ar(1), ar(2), ar(3)) = g( 3, ar(1), ar(2), ar(3)) * d( 1) + g( 1, ar(1), ar(2), ar(3)) * d( 4)
			df( dr, 4, ar(1), ar(2), ar(3)) = g( 4, ar(1), ar(2), ar(3)) * d( 1) + g( 1, ar(1), ar(2), ar(3)) * d( 5)
			df( dr, 5, ar(1), ar(2), ar(3)) =&
			( g( 2, ar(1), ar(2), ar(3)) **2 + g( 3, ar(1), ar(2), ar(3)) **2 + g( 4, ar(1), ar(2), ar(3)) **2) * d( 1) * 5d-1&
			+ d( 2) / ( gam - 1d0)&
			+ g( 1, ar(1), ar(2), ar(3)) * g( 2, ar(1), ar(2), ar(3)) * d( 3)&
			+ g( 1, ar(1), ar(2), ar(3)) * g( 3, ar(1), ar(2), ar(3)) * d( 4)&
			+ g( 1, ar(1), ar(2), ar(3)) * g( 4, ar(1), ar(2), ar(3)) * d( 5)

			dfv( dr, 2, ar(1), ar(2), ar(3)) = dfv( dr, 2, ar(1), ar(2), ar(3)) * h12( dr)
			dfv( dr, 3, ar(1), ar(2), ar(3)) = dfv( dr, 3, ar(1), ar(2), ar(3)) * h22( dr)
			dfv( dr, 4, ar(1), ar(2), ar(3)) = dfv( dr, 4, ar(1), ar(2), ar(3)) * h32( dr)
			dfv( dr, 5, ar(1), ar(2), ar(3)) = dfv( dr, 1, ar(1), ar(2), ar(3))&
											  + fv( dr, 2, ar(1), ar(2), ar(3)) * dg( dr, 2, ar(1), ar(2), ar(3)) * h11( dr)&
											  + fv( dr, 3, ar(1), ar(2), ar(3)) * dg( dr, 3, ar(1), ar(2), ar(3)) * h21( dr)&
											  + fv( dr, 4, ar(1), ar(2), ar(3)) * dg( dr, 4, ar(1), ar(2), ar(3)) * h31( dr)

		enddo
	enddo
end subroutine

subroutine wrixy( z, sta, num)
	implicit none
	integer:: z, i, j, k, sta, num

	do k = 0, num-1
		do i = 1, yk+1
			do j = 1, nx+1
				write( sta+k, "(2f10.4,e25.16)")x( j), y( i), g( k+1, j, i, z)
			enddo
			write( sta+k, *)
		enddo
		do i = yk+2, ykk
			do j = 1, xk+1
				write( sta+k, "(2f10.4,e25.16)")x( j), y( i), g( k+1, j, i, z)
			enddo
			do j = xk+2, xkk
				write( sta+k, "(2f10.4,i25)")x( j), y( i), 10
			enddo
			do j = xkk+1, nx+1
				write( sta+k, "(2f10.4,e25.16)")x( j), y( i), g( k+1, j, i, z)
			enddo
			write( sta+k, *)
		enddo
		do i = ykk+1, ny+1
			do j = 1, nx+1
				write( sta+k, "(2f10.4,e25.16)")x( j), y( i), g( k+1, j, i, z)
			enddo
			write( sta+k, *)
		enddo
	enddo
	do i = sta, sta+num-1
		close( i)
	enddo
end subroutine

subroutine wrixz( y, sta, num)
	implicit none
	integer:: y, i, j, k, sta, num

	do k = 0, num-1
		do i = 1, nz
			do j = 1, xk+1
				write( sta+k, "(2f10.4,e25.16)")x( j), z( i), g( k+1, j, y, i)
			enddo
			do j = xk+2, xkk
				write( sta+k, "(2f10.4,i25)")x( j), z( i), 10
			enddo
			do j = xkk+1, nx+1
				write( sta+k, "(2f10.4,e25.16)")x( j), z( i), g( k+1, j, y, i)
			enddo
			write( sta+k, *)
		enddo
		do j = 1, xk+1
			write( sta+k, "(2f10.4,e25.16)")x( j), z( nz+1), g( k+1, j, y, 1)
		enddo
		do j = xk+2, xkk
			write( sta+k, "(2f10.4,i25)")x( j), z( nz+1), 10
		enddo
		do j = xkk+1, nx+1
			write( sta+k, "(2f10.4,e25.16)")x( j), z( nz+1), g( k+1, j, y, 1)
		enddo
		write( sta+k, *)
	enddo
	do i = sta, sta+num-1
		close( i)
	enddo
end subroutine

subroutine wri3d( sta)
	implicit none
	integer:: i, j, k, sta

	do k = 1, nz
		do j = 1, yk+1
			do i = 1, nx+1
				write( sta, "(3f10.4,5e25.16)")x( i), y( j), z( k), g( 1,i,j,k), g( 2,i,j,k), g( 3,i,j,k),g( 4,i,j,k),g( 5,i,j,k)
			enddo
			write( sta, *)
		enddo
		do j = yk+2, ykk
			do i = 1, xk+1
				write( sta, "(3f10.4,5e25.16)")x( i), y( j), z( k), g( 1,i,j,k), g( 2,i,j,k), g( 3,i,j,k),g( 4,i,j,k),g( 5,i,j,k)
			enddo
			do i = xk+2, xkk
				write( sta, "(3f10.4,5i25)")x( i), y( j), z( k), 10, 10, 10, 10, 10
			enddo
			do i = xkk+1, nx+1
				write( sta, "(3f10.4,5e25.16)")x( i), y( j), z( k), g( 1,i,j,k), g( 2,i,j,k), g( 3,i,j,k),g( 4,i,j,k),g( 5,i,j,k)
			enddo
			write( sta, *)
		enddo
		do j = ykk+1, ny+1
			do i = 1, nx+1
				write( sta, "(3f10.4,5e25.16)")x( i), y( j), z( k), g( 1,i,j,k), g( 2,i,j,k), g( 3,i,j,k),g( 4,i,j,k),g( 5,i,j,k)
			enddo
			write( sta, *)
		enddo
	enddo
		do j = 1, yk+1
			do i = 1, nx+1
				write( sta, "(3f10.4,5e25.16)")x( i), y( j), z( nz+1),g( 1,i,j,1),g( 2,i,j,1),g( 3,i,j,1),g( 4,i,j,1),g( 5,i,j,1)
			enddo
			write( sta, *)
		enddo
		do j = yk+2, ykk
			do i = 1, xk+1
				write( sta, "(3f10.4,5e25.16)")x( i), y( j), z( nz+1),g( 1,i,j,1),g( 2,i,j,1),g( 3,i,j,1),g( 4,i,j,1),g( 5,i,j,1)
			enddo
			do i = xk+2, xkk
				write( sta, "(3f10.4,5i25)")x( i), y( j), z( nz+1), 10, 10, 10, 10, 10
			enddo
			do i = xkk+1, nx+1
				write( sta, "(3f10.4,5e25.16)")x( i), y( j), z( nz+1),g( 1,i,j,1),g( 2,i,j,1),g( 3,i,j,1),g( 4,i,j,1),g( 5,i,j,1)
			enddo
			write( sta, *)
		enddo
		do j = ykk+1, ny+1
			do i = 1, nx+1
				write( sta, "(3f10.4,5e25.16)")x( i), y( j), z( nz+1),g( 1,i,j,1),g( 2,i,j,1),g( 3,i,j,1),g( 4,i,j,1),g( 5,i,j,1)
			enddo
			write( sta, *)
		enddo
	close( sta)
end subroutine
end module

program compact
	use func
	implicit none
	integer:: i, j, k, mx1, mx2, my
	character:: fnr*16, fnu*16, fnv*16, fnw*16, fnp*16
	double precision::  q( 5, nx+1, ny+1, nz), al, be, ga, de, ep, et, th!, zz( nz), zzz( nz)

	!微分の準備
	sigA = (/ sig, -sig, 0d0, sig, -sig, 0d0, sig, -sig, 0d0,&
			  sig, -sig, 0d0, sig, -sig, 0d0, sig, -sig, 0d0,&
			  sig, -sig, 0d0/)
	dlA = (/ dx, dx, dx, dx, dx, dx, dx, dx, dx,&
			 dy, dy, dy, dy, dy, dy, dy, dy, dy,&
			 dz, dz, dz/)
	nA = (/ nx, nx, nx, xk, xk, xk, nx-xkk, nx-xkk, nx-xkk,&
			ny, ny, ny, yk, yk, yk, ny-ykk, ny-ykk, ny-ykk,&
			nz, nz, nz/)
	nn( 1) = 0
	do i = 2, nk-2
		nn( i) = nn( i-1) + nA( i-1)+1
	enddo
	do i = nk-1, nk+1
		nn( i) = nn( i-1) + nA( i-1)
	enddo
	allocate( lu( 3, nn( nk-2)), lu2( 5, nn( nk-2)+1 : nn( nk+1)))
	do i = 1, nk-3
		call NLU( i)
	enddo
	do i = nk-2, nk
		call NLU2( i)
	enddo

	!初期値設定
	if( sh == 1)then
		write (fnu , '("3d" , i3.3, ".dat")') lnp
		open( 25, file = fnu )
		do k = 1, nz
			do j = 1, ny+1
				do i = 1, nx+1
					read( 25, *) x( i), y( j), z( k), g( 1,i,j,k), g( 2,i,j,k), g( 3,i,j,k), g( 4,i,j,k), g( 5,i,j,k)
				enddo
				read( 25, *)
			enddo
		enddo
		close( 25)
	else
		do k = 1, nz
			do i = 1, ny+1
				do j = 1, nx+1
					g( 1, j, i, k) = r0
					g( 2, j, i, k) = u0
					g( 3, j, i, k) = 0d0
					g( 4, j, i, k) = 0d0
					g( 5, j, i, k) = p0
				enddo
			enddo
		enddo
	endif

	!x,y,z設定
	do i = 1, nx+1
		x( i) = dx * dble( i - 1) + xmin
	enddo
	do i = 1, ny+1
		y( i) = dy * dble( i - 1) + ymin
	enddo
	do i = 1, nz+1
		z( i) = dz * dble( i - 1) + zmin
	enddo

	!格子伸長
	if( ks == 0)then
		dxx=1d0;dyy=1d0
	else
		al = 5d0
		be = 4d0
		ep = 25d0
		et = 50d0
		th = 1d0/500d0
		ga = ((log(cosh(al*(ep-be))/cosh(al*(ep+be)))-2d0*et*log(cosh(al*(1d0-be))/cosh(al*(1d0+be))))*5d-1/al&
		 + th/3d0*(ep**3-2d0*et))/( 2d0 * et - ep)
		de = 2d0*log(cosh(al*(1d0-be))/cosh(al*(1d0+be)))*5d-1/al + 2d0*ga + 2d0/3d0*th
		do i = 1, nx+1
			dxx( i) = de / ( ga + th * x( i)**2 - 5d-1 * ( tanh( al*( x( i) + be)) - tanh( al*( x( i) - be))))!微分の逆数
			x( i) = ( ga * x( i) + th / 3d0 * x( i)**3 + 5d-1 / al * log( cosh( al*( x( i) - be)) / cosh( al*( x( i) + be))))/de
		enddo
		al = 5d0
		be = 4d0
		ep = 25d0
		et = 50d0
		th = 1d0/500d0
		ga = ((log(cosh(al*(ep-be))/cosh(al*(ep+be)))-2d0*et*log(cosh(al*(1d0-be))/cosh(al*(1d0+be))))*5d-1/al&
		 + th/3d0*(ep**3-2d0*et))/( 2d0 * et - ep)
		de = 2d0*log(cosh(al*(1d0-be))/cosh(al*(1d0+be)))*5d-1/al + 2d0*ga + 2d0/3d0*th
		do i = 1, ny+1
			dyy( i) = de / ( ga + th * y( i)**2 - 5d-1 * ( tanh( al*( y( i) + be)) - tanh( al*( y( i) - be))))!微分の逆数
			y( i) = ( ga * y( i) + th / 3d0 * y( i)**3 + 5d-1 / al * log( cosh( al*( y( i) - be)) / cosh( al*( y( i) + be))))/de
		enddo
	endif

	!バッファー領域の設定
	if( bf == 0)then
		buix=0d0;bukx=0d0;buiy=0d0;buky=0d0
	else
		mx1 = 10
		mx2 = 40
		my = 13
		al = 1.15d0
		be = 1d-2
		ga = 1.125d0
		do i = 1, nx+1
			buix( i) = al * c0 * ( dble( 1-ry)+ tanh( atanh( be / al - 1d0) / ( x( nx+1 - mx2) - x( nx+1)) * ( x( i) - x( nx+1)))&
								  - dble( ry) * tanh( atanh( be / al - 1d0) / ( x(    1 + mx1) - x(    1)) * ( x( i) - x(    1))))
		enddo
		do i = 1, mx1+1
			bukx( i) = dble( ry) * ga * c0 * ( ( x( mx1+1) - x( i)) / ( x( mx1+1) - x( 1)))**3
		enddo
		do i = mx1+2, nx - mx2
			bukx( i) = 0d0
		enddo
		do i = nx+1 - mx2, nx+1
			bukx( i) = ga * c0 * ( ( x( i) - x( nx+1 - mx2)) / ( x( nx+1) - x( nx+1 - mx2)))**3
		enddo
		do i = 1, ny+1
			buiy( i) = al * c0 * ( tanh( atanh( be / al - 1d0) / ( y( ny+1 - my) - y( ny+1)) * ( y( i) - y( ny+1)))&
								  -tanh( atanh( be / al - 1d0) / ( y(    1 + my) - y(    1)) * ( y( i) - y(    1))))
		enddo
		do i = 1, my+1
			buky( i) = ga * c0 * ( ( y( my+1) - y( i)) / ( y( my+1) - y( 1)))**3
		enddo
		do i = my+2, ny - my
			buky( i) = 0d0
		enddo
		do i = ny+1 - my, ny+1
			buky( i) = ga * c0 * ( ( y( i) - y( ny+1 - my)) / ( y( ny+1) - y( ny+1 - my)))**3
		enddo
	endif

	!微分確認用
!	open( 34, file = "bx.dat");open( 35, file = "by.dat")
!	do i = 1, nx+1
!		g( 1,i,:,:) = cos( pi * x( i) / 1d1)
!		g( 2,i,:,:) = sin( pi * x( i) / 1d1)
!	enddo
!	do i = 1, ny+1
!		g( 3,:,i,:) = cos( pi * y( i) / 1d1)
!		g( 4,:,i,:) = sin( pi * y( i) / 1d1)
!	enddo
!	call bbn( dg, g)
!	do i = 1, nx+1
!		write( 34, "(f15.4,4e25.16)") x( i), g( 1, i, 1, 1), g( 2, i, 1, 1), dg( 1, 1, i, 1, 1), dg( 1, 2, i, 1, 1)
!	enddo
!	do i = 1, ny+1
!		write( 35, "(f15.4,4e25.16)") y( i), g( 3, 1, i, 1), g( 4, 1, i, 1), dg( 2, 3, 1, i, 1), dg( 2, 4, 1, i, 1)
!	enddo
!	stop
!	do i = 1, nx+1
!		write( 34, "(4f10.4)") dble( i-1)/2d0-10d0, x( i), 1d0/dxx( i), x( i)-x( i-1)
!	enddo

	!g→q変換
	do k = 1, nz
		do j = 1, ny+1
			do i = 1, nx+1
				q( 1, i, j, k) = g( 1, i, j, k)
				q( 2, i, j, k) = g( 2, i, j, k) * g( 1, i, j, k)
				q( 3, i, j, k) = g( 3, i, j, k) * g( 1, i, j, k)
				q( 4, i, j, k) = g( 4, i, j, k) * g( 1, i, j, k)
				q( 5, i, j, k) = g( 5, i, j, k) / ( gam - 1d0) + 5d-1 * g( 1, i, j, k)&
				* ( g( 2, i, j, k) * g( 2, i, j, k) + g( 3, i, j, k) * g( 3, i, j, k) + g( 4, i, j, k) * g( 4, i, j, k))
			enddo
		enddo
	enddo

	!初期値書き出し
	if( sh == 0)then
!		open( 30, file = "xy2r000.dat")
!		open( 31, file = "xy2u000.dat")
!		open( 32, file = "xy2v000.dat")
!		open( 33, file = "xy2w000.dat")
!		open( 34, file = "xy2p000.dat")
!		call wrixy( nz/2+1, 30, 5)
!		open( 35, file = "xz10r000.dat")
!		open( 36, file = "xz10u000.dat")
!		open( 37, file = "xz10v000.dat")
!		open( 38, file = "xz10w000.dat")
!		open( 39, file = "xz10p000.dat")
!		call wrixz( ny/2+1, 35, 5)
		open( 50, file = "3d000.dat" )
		call wri3d( 50)
!		do i = 2, 4
!			call bbx( dg( 1, i,:,:,:), g( i,:,:,:))
!			call bby( dg( 2, i,:,:,:), g( i,:,:,:))
!			call bbz( dg( 3, i,:,:,:), g( i,:,:,:))
!		enddo
!		g( 1,:,:,:) = dg( 2, 4,:,:,:) - dg( 3, 3,:,:,:)
!		g( 2,:,:,:) = dg( 3, 2,:,:,:) - dg( 1, 4,:,:,:)
!		g( 3,:,:,:) = dg( 1, 3,:,:,:) - dg( 2, 2,:,:,:)
!		g( 4,:,:,:) = dg( 1, 2,:,:,:) * dg( 2, 3,:,:,:) + dg( 2, 3,:,:,:) * dg( 3, 4,:,:,:) + dg( 3, 4,:,:,:) * dg( 1, 2,:,:,:)&
!					- dg( 1, 3,:,:,:) * dg( 2, 2,:,:,:) - dg( 2, 4,:,:,:) * dg( 3, 3,:,:,:) - dg( 3, 2,:,:,:) * dg( 1, 4,:,:,:)
!		open( 30, file = "xy2uzx000.dat")
!		open( 31, file = "xy2uzy000.dat")
!		open( 32, file = "xy2uzz000.dat")
!		open( 33, file = "xy2q000.dat"  )
!		call wrixy( nz/2+1, 30, 4)
!		open( 35, file = "xz10uzx000.dat")
!		open( 36, file = "xz10uzy000.dat")
!		open( 37, file = "xz10uzz000.dat")
!		open( 38, file = "xz10q000.dat"  )
!		call wrixz( ny/2+1, 35, 4)
!		open( 50, file = "uz000.dat" )
!		call wri3d( 50)
	endif

	!時間ループ
	do m = 1 + lnp * sh, ln + lnp * sh
		do i = 1, lp/ln
			q = ( q + 2d0 * cal( 7.5d-1 * q + 2.5d-1 * cal( cal( q)))) / 3d0
			tm = tm + dt
		enddo

		!エラー検出
		if( isnan(q( 1, 5, 5, 5)))then
			open( 99, file = "err.dat"); write(99,*)"err"; write(*,*)"err" ;stop
		endif

		!q→g変換
		do k = 1, nz
			do j = 1, ny+1
				do i = 1, nx+1
					g( 1, i, j, k) = q( 1, i, j, k)
					g( 2, i, j, k) = q( 2, i, j, k) / q( 1, i, j, k)
					g( 3, i, j, k) = q( 3, i, j, k) / q( 1, i, j, k)
					g( 4, i, j, k) = q( 4, i, j, k) / q( 1, i, j, k)
					g( 5, i, j, k) = ( gam - 1d0) * ( q( 5, i, j, k) - 5d-1 / q( 1, i, j, k)&
					* ( q( 2, i, j, k) * q( 2, i, j, k) + q( 3, i, j, k) * q( 3, i, j, k) + q( 4, i, j, k) * q( 4, i, j, k)))
				enddo
			enddo
		enddo

		!値書き出し
!		write (fnr , '("xy2r" , i3.3, ".dat")') m
!		write (fnu , '("xy2u" , i3.3, ".dat")') m
!		write (fnv , '("xy2v" , i3.3, ".dat")') m
!		write (fnw , '("xy2w" , i3.3, ".dat")') m
!		write (fnp , '("xy2p" , i3.3, ".dat")') m
!		open( 30, file = fnr )
!		open( 31, file = fnu )
!		open( 32, file = fnv )
!		open( 33, file = fnw )
!		open( 34, file = fnp )
!		call wrixy( nz/2+1, 30, 5)
!		write (fnr , '("xz10r" , i3.3, ".dat")') m
!		write (fnu , '("xz10u" , i3.3, ".dat")') m
!		write (fnv , '("xz10v" , i3.3, ".dat")') m
!		write (fnw , '("xz10w" , i3.3, ".dat")') m
!		write (fnp , '("xz10p" , i3.3, ".dat")') m
!		open( 35, file = fnr )
!		open( 36, file = fnu )
!		open( 37, file = fnv )
!		open( 38, file = fnw )
!		open( 39, file = fnp )
!		call wrixz( ny/2+1, 35, 5)
		write (fnu , '("3d" , i3.3, ".dat")') m
		open( 50, file = fnu )
		call wri3d( 50)
!		!渦度
!		do i = 2, 4
!			call bbx( dg( 1, i,:,:,:), g( i,:,:,:))
!			call bby( dg( 2, i,:,:,:), g( i,:,:,:))
!			call bbz( dg( 3, i,:,:,:), g( i,:,:,:))
!		enddo
!		g( 1,:,:,:) = dg( 2, 4,:,:,:) - dg( 3, 3,:,:,:)
!		g( 2,:,:,:) = dg( 3, 2,:,:,:) - dg( 1, 4,:,:,:)
!		g( 3,:,:,:) = dg( 1, 3,:,:,:) - dg( 2, 2,:,:,:)
!		g( 4,:,:,:) = dg( 1, 2,:,:,:) * dg( 2, 3,:,:,:) + dg( 2, 3,:,:,:) * dg( 3, 4,:,:,:) + dg( 3, 4,:,:,:) * dg( 1, 2,:,:,:)&
!					- dg( 1, 3,:,:,:) * dg( 2, 2,:,:,:) - dg( 2, 4,:,:,:) * dg( 3, 3,:,:,:) - dg( 3, 2,:,:,:) * dg( 1, 4,:,:,:)
!		write (fnu , '("xy2uzx" , i3.3, ".dat")') m
!		write (fnv , '("xy2uzy" , i3.3, ".dat")') m
!		write (fnw , '("xy2uzz" , i3.3, ".dat")') m
!		write (fnp , '("xy2q"   , i3.3, ".dat")') m
!		open( 30, file = fnu )
!		open( 31, file = fnv )
!		open( 32, file = fnw )
!		open( 33, file = fnp )
!		call wrixy( nz/2+1, 30, 4)
!		write (fnu , '("xz10uzx" , i3.3, ".dat")') m
!		write (fnv , '("xz10uzy" , i3.3, ".dat")') m
!		write (fnw , '("xz10uzz" , i3.3, ".dat")') m
!		write (fnp , '("xz10q"   , i3.3, ".dat")') m
!		open( 35, file = fnu )
!		open( 36, file = fnv )
!		open( 37, file = fnw )
!		open( 38, file = fnp )
!		call wrixz( ny/2+1, 35, 4)
!		write (fnu , '("uz" , i3.3, ".dat")') m
!		open( 50, file = fnu )
!		call wri3d( 50)

	enddo
end program
