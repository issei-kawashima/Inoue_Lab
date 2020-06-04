!!!===============================================================================================================================
!!!超音速ジェット(攪乱なし)を解くプログラム(3次元)
!!!	x方向:NSCBC
!!!	y方向:NSCBC
!!!	z方向:PBC
!!!	サブルーチンでは、すべてNNの数字を読むが、x,y方向ではNNまで計算、z方向ではNN-1まで計算する。
!!!20181011現在 NSCBCまで改変
!!!===============================================================================================================================



module sub
	implicit none
	!!!計算パラメータの設定-------------------------------------------------------------------------------------------------------
!	integer, parameter :: NX=300,NY=200,NZ=30,NT=1000,Nout=NT,it=60
	integer, parameter :: NX=150,NY=100,NZ=15,NT=1000,Nout=NT,it=60
	double precision, parameter :: Re=3000d0,Ma=2d0,Pr=0.72d0,T0=291.15d0,Sc=120d0/T0
	double precision, parameter :: zeta=1d0,sigma=0.25d0,gamma=1.4d0
	double precision, parameter :: dx=60d0/dble(NX),dy=40d0/dble(NY),dz=6d0/dble(NZ),dt=1d0/dble(NT)
	double precision, parameter :: rho0=1d0,p0=1d0/(gamma*Ma*Ma),u0=0d0,u1=1d0,v0=0d0,w0=0d0
	double precision, parameter :: T_in=1.12d0,T_wall=1d0
	contains
	!!!サブルーチン---------------------------------------------------------------------------------------------------------------
	subroutine G_Q(G,Q)					!!!G→Q
		implicit none
		double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: G,Q
			Q(1,:,:,:)=G(1,:,:,:)
			Q(2,:,:,:)=G(1,:,:,:)*G(2,:,:,:)
			Q(3,:,:,:)=G(1,:,:,:)*G(3,:,:,:)
			Q(4,:,:,:)=G(1,:,:,:)*G(4,:,:,:)
			Q(5,:,:,:)=G(5,:,:,:)/(gamma-1d0)+G(1,:,:,:)*(G(2,:,:,:)*G(2,:,:,:)+G(3,:,:,:)*G(3,:,:,:)+	&
						G(4,:,:,:)*G(4,:,:,:))*0.5d0
	end subroutine G_Q

	subroutine Q_G(Q,G)					!!!Q→G
		implicit none
		double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: Q,G
			G(1,:,:,:)=Q(1,:,:,:)
			G(2,:,:,:)=Q(2,:,:,:)/Q(1,:,:,:)
			G(3,:,:,:)=Q(3,:,:,:)/Q(1,:,:,:)
			G(4,:,:,:)=Q(4,:,:,:)/Q(1,:,:,:)
			G(5,:,:,:)=(gamma-1d0)*(Q(5,:,:,:)-(Q(2,:,:,:)*Q(2,:,:,:)+Q(3,:,:,:)*Q(3,:,:,:)+	&
						Q(4,:,:,:)*Q(4,:,:,:))/(2d0*Q(1,:,:,:)))
	end subroutine Q_G

	subroutine Q_Fx_LF(Q,Fxp,Fxm)		!!!Q→Fx→Fxp,Fxm
		implicit none
		double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: Q,Fx,Fxp,Fxm
			Fx(1,:,:,:)=Q(2,:,:,:)
			Fx(2,:,:,:)=(Q(2,:,:,:)*Q(2,:,:,:))/Q(1,:,:,:)+(gamma-1d0)*(Q(5,:,:,:)-	&
						(Q(2,:,:,:)*Q(2,:,:,:)+Q(3,:,:,:)*Q(3,:,:,:)+Q(4,:,:,:)*Q(4,:,:,:))/(2d0*Q(1,:,:,:)))
			Fx(3,:,:,:)=(Q(2,:,:,:)*Q(3,:,:,:))/Q(1,:,:,:)
			Fx(4,:,:,:)=(Q(2,:,:,:)*Q(4,:,:,:))/Q(1,:,:,:)
			Fx(5,:,:,:)=(Q(5,:,:,:)+(gamma-1d0)*(Q(5,:,:,:)-(Q(2,:,:,:)*Q(2,:,:,:)+Q(3,:,:,:)*Q(3,:,:,:)+	&
						Q(4,:,:,:)*Q(4,:,:,:))/(2d0*Q(1,:,:,:))))*(Q(2,:,:,:)/Q(1,:,:,:))
			!!!Lax-Friedrich
			Fxp=(Fx+Q*zeta)*0.5d0
			Fxm=(Fx-Q*zeta)*0.5d0
	end subroutine Q_Fx_LF

	subroutine Q_Fy_LF(Q,Fyp,Fym)		!!!Q→Fy→Fyp,Fym
		implicit none
		double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: Q,Fy,Fyp,Fym
			Fy(1,:,:,:)=Q(3,:,:,:)
			Fy(2,:,:,:)=(Q(2,:,:,:)*Q(3,:,:,:))/Q(1,:,:,:)
			Fy(3,:,:,:)=(Q(3,:,:,:)*Q(3,:,:,:))/Q(1,:,:,:)+(gamma-1d0)*(Q(5,:,:,:)-	&
						(Q(2,:,:,:)*Q(2,:,:,:)+Q(3,:,:,:)*Q(3,:,:,:)+Q(4,:,:,:)*Q(4,:,:,:))/(2d0*Q(1,:,:,:)))
			Fy(4,:,:,:)=(Q(3,:,:,:)*Q(4,:,:,:))/Q(1,:,:,:)
			Fy(5,:,:,:)=(Q(5,:,:,:)+(gamma-1d0)*(Q(5,:,:,:)-(Q(2,:,:,:)*Q(2,:,:,:)+Q(3,:,:,:)*Q(3,:,:,:)+	&
						Q(4,:,:,:)*Q(4,:,:,:))/(2d0*Q(1,:,:,:))))*(Q(3,:,:,:)/Q(1,:,:,:))
			!!!Lax-Friedrich
			Fyp=(Fy+Q*zeta)*0.5d0
			Fym=(Fy-Q*zeta)*0.5d0
	end subroutine Q_Fy_LF

	subroutine Q_Fz_LF(Q,Fzp,Fzm)		!!!!!!Q→Fz→Fzp,Fzm
		implicit none
		double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: Q,Fz,Fzp,Fzm
			Fz(1,:,:,:)=Q(4,:,:,:)
			Fz(2,:,:,:)=(Q(2,:,:,:)*Q(4,:,:,:))/Q(1,:,:,:)
			Fz(3,:,:,:)=(Q(3,:,:,:)*Q(4,:,:,:))/Q(1,:,:,:)
			Fz(4,:,:,:)=(Q(4,:,:,:)*Q(4,:,:,:))/Q(1,:,:,:)+(gamma-1d0)*(Q(5,:,:,:)-	&
						(Q(2,:,:,:)*Q(2,:,:,:)+Q(3,:,:,:)*Q(3,:,:,:)+Q(4,:,:,:)*Q(4,:,:,:))/(2d0*Q(1,:,:,:)))
			Fz(5,:,:,:)=(Q(5,:,:,:)+(gamma-1d0)*(Q(5,:,:,:)-(Q(2,:,:,:)*Q(2,:,:,:)+Q(3,:,:,:)*Q(3,:,:,:)+	&
						Q(4,:,:,:)*Q(4,:,:,:))/(2d0*Q(1,:,:,:))))*(Q(4,:,:,:)/Q(1,:,:,:))
			!!!Lax-Friedrich
			Fzp=(Fz+Q*zeta)*0.5d0
			Fzm=(Fz-Q*zeta)*0.5d0
	end subroutine Q_Fz_LF

	subroutine Q_Fvx_Fvy_Fvz(G,dGx,dGy,dGz,Fvx,Fvy,Fvz,LU_C_NX,LU_C_NY,L_C_NZ,U_C_NZ,dydxi)	!!!G,dGx,dGy,dGz→mu,T→Fvx,Fvy,Fvz
		implicit none
		integer :: i,j,k
		double precision, dimension(0:NY) :: dydxi
		double precision, dimension(0:NX,0:NY,0:NZ-1) :: mu,T,dTx,dTy,dTz
		double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: G,dGx,dGy,dGz,Fvx,Fvy,Fvz
		double precision, dimension(-1:1,0:NX) :: LU_C_NX
		double precision, dimension(-1:1,0:NY) :: LU_C_NY
		double precision, dimension(0:NZ-1,0:NZ-1) :: L_C_NZ,U_C_NZ
			T(:,:,:)=(gamma*Ma*Ma*G(5,:,:,:))/G(1,:,:,:)
			mu(:,:,:)=(T(:,:,:)**(2d0/3d0))*((1d0+Sc)/(T(:,:,:)+Sc))
			do j=0,NY
				do k=0,NZ-1
					call CCS(0,NX,dx,T(:,j,k),dTx(:,j,k),LU_C_NX)
				end do
			end do
			do i=0,NX
				do k=0,NZ-1
					call CCS(0,NY,dy,T(i,:,k),dTy(i,:,k),LU_C_NY)
					dTy(i,:,k)=dTy(i,:,k)*dydxi(:)
				end do
			end do
			do i=0,NX
				do j=0,NY
					call CCS_PBC(0,NZ,dz,T(i,j,:),dTz(i,j,:),L_C_NZ,U_C_NZ)
				end do
			end do
			!!!Fvxの決定
			Fvx(1,:,:,:)=0d0
			Fvx(2,:,:,:)=((2d0*mu(:,:,:))/(3d0*Re))*(2d0*dGx(2,:,:,:)-dGy(3,:,:,:)-dGz(4,:,:,:))
			Fvx(3,:,:,:)=(mu(:,:,:)/Re)*(dGx(3,:,:,:)+dGy(2,:,:,:))
			Fvx(4,:,:,:)=(mu(:,:,:)/Re)*(dGx(4,:,:,:)+dGz(2,:,:,:))
			Fvx(5,:,:,:)=(((2d0*mu(:,:,:))/(3d0*Re))*(2d0*dGx(2,:,:,:)-dGy(3,:,:,:)-dGz(4,:,:,:)))*G(2,:,:,:)	&
						+((mu(:,:,:)/Re)*(dGx(3,:,:,:)+dGy(2,:,:,:)))*G(3,:,:,:)	&
						+((mu(:,:,:)/Re)*(dGx(4,:,:,:)+dGz(2,:,:,:)))*G(4,:,:,:)	&
						+mu(:,:,:)*dTx(:,:,:)/((gamma-1d0)*Pr*Ma*Ma*Re)
			!!!Fvyの決定
			Fvy(1,:,:,:)=0d0
			Fvy(2,:,:,:)=(mu(:,:,:)/Re)*(dGx(3,:,:,:)+dGy(2,:,:,:))
			Fvy(3,:,:,:)=((2d0*mu(:,:,:))/(3d0*Re))*(-dGx(2,:,:,:)+2d0*dGy(3,:,:,:)-dGz(4,:,:,:))
			Fvy(4,:,:,:)=(mu(:,:,:)/Re)*(dGz(3,:,:,:)+dGy(4,:,:,:))
			Fvy(5,:,:,:)=((mu(:,:,:)/Re)*(dGx(3,:,:,:)+dGy(2,:,:,:)))*G(2,:,:,:)	&
						+(((2d0*mu(:,:,:))/(3d0*Re))*(-dGx(2,:,:,:)+2d0*dGy(3,:,:,:)-dGz(4,:,:,:)))*G(3,:,:,:)	&
						+((mu(:,:,:)/Re)*(dGz(3,:,:,:)+dGy(4,:,:,:)))*G(4,:,:,:)	&
						+mu(:,:,:)*dTy(:,:,:)/((gamma-1d0)*Pr*Ma*Ma*Re)
			!!!Fvzの決定
			Fvz(1,:,:,:)=0d0
			Fvz(2,:,:,:)=(mu(:,:,:)/Re)*(dGx(4,:,:,:)+dGz(2,:,:,:))
			Fvz(3,:,:,:)=(mu(:,:,:)/Re)*(dGy(4,:,:,:)+dGz(3,:,:,:))
			Fvz(4,:,:,:)=((2d0*mu(:,:,:))/(3d0*Re))*(-dGx(2,:,:,:)-dGy(3,:,:,:)+2d0*dGz(4,:,:,:))
			Fvz(5,:,:,:)=((mu(:,:,:)/Re)*(dGx(4,:,:,:)+dGz(2,:,:,:)))*G(2,:,:,:)	&
						+((mu(:,:,:)/Re)*(dGy(4,:,:,:)+dGz(3,:,:,:)))*G(3,:,:,:)	&
						+(((2d0*mu(:,:,:))/(3d0*Re))*(-dGx(2,:,:,:)-dGy(3,:,:,:)+2d0*dGz(4,:,:,:)))*G(4,:,:,:)	&
						+mu(:,:,:)*dTz(:,:,:)/((gamma-1d0)*Pr*Ma*Ma*Re)
			!!!無反射流出条件による境界条件の追加
			!!!Fvx
			Fvx(5,NX,:,:)=(((2d0*mu(NX,:,:))/(3d0*Re))*(2d0*dGx(2,NX,:,:)-dGy(3,NX,:,:)-dGz(4,NX,:,:)))*G(2,NX,:,:)
			Fvx(5,:,0,:)=(((2d0*mu(:,0,:))/(3d0*Re))*(2d0*dGx(2,:,0,:)-dGy(3,:,0,:)-dGz(4,:,0,:)))*G(2,:,0,:)
			Fvx(5,:,NY,:)=(((2d0*mu(:,NY,:))/(3d0*Re))*(2d0*dGx(2,:,NY,:)-dGy(3,:,NY,:)-dGz(4,:,NY,:)))*G(2,:,NY,:)
			!!!Fvy
			Fvy(5,NX,:,:)=(((2d0*mu(NX,:,:))/(3d0*Re))*(-dGx(2,NX,:,:)+2d0*dGy(3,NX,:,:)-dGz(4,NX,:,:)))*G(3,NX,:,:)
			Fvy(5,:,0,:)=(((2d0*mu(:,0,:))/(3d0*Re))*(-dGx(2,:,0,:)+2d0*dGy(3,:,0,:)-dGz(4,:,0,:)))*G(3,:,0,:)
			Fvy(5,:,NY,:)=(((2d0*mu(:,NY,:))/(3d0*Re))*(-dGx(2,:,NY,:)+2d0*dGy(3,:,NY,:)-dGz(4,:,NY,:)))*G(3,:,NY,:)
			!!!Fvz
			Fvz(5,NX,:,:)=(((2d0*mu(NX,:,:))/(3d0*Re))*(-dGx(2,NX,:,:)-dGy(3,NX,:,:)+2d0*dGz(4,NX,:,:)))*G(4,NX,:,:)
			Fvz(5,:,0,:)=(((2d0*mu(:,0,:))/(3d0*Re))*(-dGx(2,:,0,:)-dGy(3,:,0,:)+2d0*dGz(4,:,0,:)))*G(4,:,0,:)
			Fvz(5,:,NY,:)=(((2d0*mu(:,NY,:))/(3d0*Re))*(-dGx(2,:,NY,:)-dGy(3,:,NY,:)+2d0*dGz(4,:,NY,:)))*G(4,:,NY,:)
	end subroutine Q_Fvx_Fvy_Fvz

	subroutine LU(fn,NN,LUcom,sigma)	!!!LU分解と圧縮
		implicit none
		integer :: i, j, k, n
		integer, intent(in) :: fn, NN
		double precision, intent(in) :: sigma
		double precision :: alpha_5, alpha_D5
		double precision :: alpha_3, alpha_D3
		double precision :: alpha_OS
		double precision :: sum
		double precision, intent(out) :: LUcom(-1:1,fn:NN)
		double precision, dimension(fn:NN,fn:NN) :: A, L, U
			!!!初期値および係数
			LUcom(:,:)=0d0
			A(:,:)=0d0
			L(:,:)=0d0
			U(:,:)=0d0
			sum=0d0
			alpha_5=1d0/3d0
			alpha_D5=1d0
			alpha_3=1d0/4d0
			alpha_D3=1d0
			alpha_OS=3d0
			!!!A行列の決定
			do i=fn,NN
				A(i,i)=1d0
			end do
			!!!片側CS
			A(fn,fn+1)=alpha_OS
			A(NN,NN-1)=alpha_OS
			!!!CCS4,DCS3
			A(fn+1,fn)=alpha_3*(1d0-sigma*alpha_D3)
			A(fn+1,fn+2)=alpha_3*(1d0+sigma*alpha_D3)
			A(NN-1,NN)=alpha_3*(1d0+sigma*alpha_D3)
			A(NN-1,NN-2)=alpha_3*(1d0-sigma*alpha_D3)
			!!!CCS6,DCS5
			do i=fn+2,NN-2
				A(i,i-1)=alpha_5*(1d0-sigma*alpha_D5)
				A(i,i+1)=alpha_5*(1d0+sigma*alpha_D5)
			end do
			!!!LU分解
			do j=fn,NN
				U(fn,j)=A(fn,j)
			end do
			do i=fn,NN
				L(i,fn)=A(i,fn)/U(fn,fn)
				L(i,i)=1d0
			end do
			do i=fn+1,NN
				do j=i,NN
					do k=fn,i-1
						sum=sum+L(i,k)*U(k,j)
					end do
					U(i,j)=A(i,j)-sum
					sum=0d0
				end do
				do j=i+1,NN
					do k=fn,i-1
						sum=sum+L(j,k)*U(k,i)
					end do
					L(j,i)=(A(j,i)-sum)/U(i,i)
					sum=0d0
				end do
			end do
			!!!LU圧縮
			!!!6次以上の精度を使用する場合, 5行必要になるため注意!
			LUcom(-1,fn)=1d0
			LUcom(0,NN)=U(NN,NN)
			LUcom(1,NN)=U(NN,NN)
			do i=fn+1,NN
				LUcom(-1,i)=L(i,i-1)
			end do
			do i=fn,NN-1
				do k=0,1
					LUcom(k,i)=U(i,i+k)
				end do
			end do
	end subroutine LU

	subroutine DCS(fn,NN,dh,F,dF,LUcom,sigma)				!!!fnからNNまで微分
		implicit none
		integer :: fn, NN, i, j, k
		double precision :: sigma, dh
		double precision, intent(in) :: LUcom(-1:1,fn:NN)
		double precision, dimension(fn:NN) :: B
		double precision, dimension(1:5,fn:NN) :: F, dF, RHS
			!!!DCS右辺の決定
			do i=1,5
				!!!片側CS
				RHS(i,fn)=(-17d0/6d0*F(i,fn)+3d0/2d0*F(i,fn+1)+3d0/2d0*F(i,fn+2)-1d0/6d0*F(i,fn+3))/dh
				RHS(i,NN)=(1d0/6d0*F(i,NN-3)-3d0/2d0*F(i,NN-2)-3d0/2d0*F(i,NN-1)+17d0/6d0*F(i,NN))/dh
				!!!3次精度DCS
				RHS(i,fn+1)=(3d0/4d0*(F(i,fn+2)-F(i,fn))+sigma*0.5d0*(F(i,fn)-2d0*F(i,fn+1)+F(i,fn+2)))/dh
				RHS(i,NN-1)=(3d0/4d0*(F(i,NN)-F(i,NN-2))+sigma*0.5d0*(F(i,NN-2)-2d0*F(i,NN-1)+F(i,NN)))/dh
				!!!5次精度DCS
				do j=fn+2,NN-2
					RHS(i,j)=((7d0/9d0*(F(i,j+1)-F(i,j-1))+1d0/36d0*(F(i,j+2)-F(i,j-2)))	&
							+sigma*(4d0/9d0*(F(i,j-1)-2d0*F(i,j)+F(i,j+1))+1d0/18d0*(F(i,j-2)-2d0*F(i,j)+F(i,j+2))))/dh
				end do
			end do
			!!!前進代入/後退代入
			do i=1,5
				!!!前進代入
				B(fn)=RHS(i,fn)
				do j=fn+1,NN
					B(j)=(RHS(i,j)-LUcom(-1,j)*B(j-1))
				end do
				!!!後退代入
				dF(i,NN)=B(NN)/LUcom(0,NN)
				do j=NN-1,fn,-1
					dF(i,j)=(B(j)-LUcom(1,j)*dF(i,j+1))/LUcom(0,j)
				end do
			end do
	end subroutine DCS

	subroutine CCS(fn,NN,dh,F,dF,LUcom)
		implicit none
		integer :: fn, NN, i, j, k
		double precision :: dh
		double precision, intent(in) :: LUcom(-1:1,fn:NN)
		double precision, dimension(fn:NN) :: B
		double precision, dimension(fn:NN) :: F, dF, RHS
			!!!CCS右辺の決定
			!!!片側CS
			RHS(fn)=(-17d0/6d0*F(fn)+3d0/2d0*F(fn+1)+3d0/2d0*F(fn+2)-1d0/6d0*F(fn+3))/dh
			RHS(NN)=(1d0/6d0*F(NN-3)-3d0/2d0*F(NN-2)-3d0/2d0*F(NN-1)+17d0/6d0*F(NN))/dh
			!!!4次精度CCS
			RHS(fn+1)=(3d0/4d0*(F(fn+2)-F(fn)))/dh!+sigma*0.5d0*(F(fn)-2d0*F(fn+1)+F(fn+2)))/dh
			RHS(NN-1)=(3d0/4d0*(F(NN)-F(NN-2)))/dh!+sigma*0.5d0*(F(NN-2)-2d0*F(NN-1)+F(NN)))/dh
			!!!6次精度CCS
			do j=fn+2,NN-2
				RHS(j)=((7d0/9d0*(F(j+1)-F(j-1))+1d0/36d0*(F(j+2)-F(j-2))))/dh	!&
						!+sigma*(4d0/9d0*(F(j-1)-2d0*F(j)+F(j+1))+1d0/18d0*(F(j-2)-2d0*F(j)+F(j+2))))/dh
			end do
			!!!前進代入/後退代入
			!!!前進代入
			B(fn)=RHS(fn)
			do j=fn+1,NN
				B(j)=(RHS(j)-LUcom(-1,j)*B(j-1))
			end do
			!!!後退代入
			dF(NN)=B(NN)/LUcom(0,NN)
			do j=NN-1,fn,-1
				dF(j)=(B(j)-LUcom(1,j)*dF(j+1))/LUcom(0,j)
			end do
	end subroutine CCS

	subroutine LU_PBC(fn,NN,L,U,sigma)	!!!LU分解(PBC)
		implicit none
		integer :: i, j, k, n
		integer :: im, ip
		integer, intent(in) :: fn, NN
		double precision, intent(in) :: sigma
		double precision :: alpha_5, alpha_D5
		double precision :: sum
		double precision, intent(out), dimension(fn:NN-1,fn:NN-1) :: L,U
		double precision, dimension(fn:NN-1,fn:NN-1) :: A
			!!!初期値および係数
			A(:,:)=0d0
			L(:,:)=0d0
			U(:,:)=0d0
			sum=0d0
			alpha_5=1d0/3d0
			alpha_D5=1d0
			!!!A行列の決定
			do i=0,NN-1
				im=mod(NN+i-1,NN)
				ip=mod(NN+i+1,NN)
				do j=0,NN-1
					if(j==i) then
						A(i,j)=1d0
					else if(j==im) then
						A(i,j)=alpha_5*(1d0-sigma*alpha_D5)
					else if(j==ip) then
						A(i,j)=alpha_5*(1d0+sigma*alpha_D5)
					endif
				end do
			end do
			!!!LU分解
			do j=fn,NN-1
				U(fn,j)=A(fn,j)
			end do
			do i=fn,NN-1
				L(i,fn)=A(i,fn)/U(fn,fn)
				L(i,i)=1d0
			end do
			do i=fn+1,NN-1
				do j=i,NN-1
					do k=fn,i-1
						sum=sum+L(i,k)*U(k,j)
					end do
					U(i,j)=A(i,j)-sum
					sum=0d0
				end do
				do j=i+1,NN-1
					do k=fn,i-1
						sum=sum+L(j,k)*U(k,i)
					end do
					L(j,i)=(A(j,i)-sum)/U(i,i)
					sum=0d0
				end do
			end do
	end subroutine LU_PBC

	subroutine DCS_PBC(fn,NN,dh,F,dF,L,U,sigma)				!!!fnからNN-1まで微分
		implicit none
		integer :: fn, NN, i, j, k
		integer :: jm2, jm1, jp1, jp2
		double precision :: sigma, dh, sum
		double precision, intent(in), dimension(fn:NN-1,fn:NN-1) :: L, U
		double precision, dimension(fn:NN-1) :: B
		double precision, dimension(1:5,fn:NN-1) :: F, dF, RHS
			!!!DCS右辺の決定
			do i=1,5
				do j=fn,NN-1
					jm1=mod(NN+j-1,NN)
					jp1=mod(NN+j+1,NN)
					jm2=mod(NN+j-2,NN)
					jp2=mod(NN+j+2,NN)
					RHS(i,j)=((7d0/9d0*(F(i,jp1)-F(i,jm1))+1d0/36d0*(F(i,jp2)-F(i,jm2)))	&
							+sigma*(4d0/9d0*(F(i,jm1)-2d0*F(i,j)+F(i,jp1))+1d0/18d0*(F(i,jm2)-2d0*F(i,j)+F(i,jp2))))/dh
				end do
			end do
			!!!前進代入/後退代入
			do i=1,5
				sum=0d0
				!!!前進代入
				B(fn)=RHS(i,fn)
				do j=fn+1,NN-2
					B(j)=(RHS(i,j)-L(j,j-1)*B(j-1))
				end do
				do j=fn,NN-2
					sum=sum+L(NN-1,j)*B(j)
				end do
				B(NN-1)=RHS(i,NN-1)-sum
				!!!後退代入
				dF(i,NN-1)=B(NN-1)/U(NN-1,NN-1)
				dF(i,NN-2)=(B(NN-2)-U(NN-2,NN-1)*dF(i,NN-1))/U(NN-2,NN-2)
				do j=NN-3,fn,-1
					dF(i,j)=(B(j)-U(j,j+1)*dF(i,j+1)-U(j,NN-1)*dF(i,NN-1))/U(j,j)
				end do
			end do
	end subroutine DCS_PBC

	subroutine CCS_PBC(fn,NN,dh,F,dF,L,U)
		implicit none
		integer :: fn, NN, i, j, k
		integer :: jm2, jm1, jp1, jp2
		double precision :: dh, sum
		double precision, intent(in), dimension(fn:NN-1,fn:NN-1) :: L, U
		double precision, dimension(fn:NN-1) :: B
		double precision, dimension(fn:NN-1) :: F, dF, RHS
			!!!CCS右辺の決定
			do j=fn,NN-1
				jm1=mod(NN+j-1,NN)
				jp1=mod(NN+j+1,NN)
				jm2=mod(NN+j-2,NN)
				jp2=mod(NN+j+2,NN)
				RHS(j)=((7d0/9d0*(F(jp1)-F(jm1))+1d0/36d0*(F(jp2)-F(jm2))))/dh	!&
						!+sigma*(4d0/9d0*(F(j-1)-2d0*F(j)+F(j+1))+1d0/18d0*(F(j-2)-2d0*F(j)+F(j+2))))/dh
			end do
			!!!前進代入/後退代入
			sum=0d0
			!!!前進代入
			B(fn)=RHS(fn)
			do j=fn+1,NN-2
				B(j)=(RHS(j)-L(j,j-1)*B(j-1))
			end do
			do j =fn,NN-1
				sum=sum+L(NN-1,j)*B(j)
			end do
			B(NN-1)=RHS(NN-1)-sum
			!!!後退代入
			dF(NN-1)=B(NN-1)/U(NN-1,NN-1)
			dF(NN-2)=(B(NN-2)-U(NN-2,NN-1)*dF(NN-1))/U(NN-2,NN-2)
			do j=NN-3,fn,-1
				dF(j)=(B(j)-U(j,j+1)*dF(j+1)-U(j,NN-1)*dF(NN-1))/U(j,j)
			end do
	end subroutine CCS_PBC

	subroutine NSCBC_x(G,dGx,dFx)
		implicit none
		double precision :: LX
		double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: G,dGx,dFx
		double precision, dimension(1:5,0:1,0:NY,0:NZ-1) :: L,d
		double precision, dimension(0:1,0:NY,0:NZ-1) :: c,Ma_NSCBC
			LX=dx*dble(NX)
			c(0,:,:)=sqrt(gamma*G(5,0,:,:)/G(1,0,:,:))
			c(1,:,:)=sqrt(gamma*G(5,NX,:,:)/G(1,NX,:,:))
!			Ma_NSCBC(0,:,:)=G(2,0,:,:)/c(0,:,:)
			Ma_NSCBC(1,:,:)=G(2,NX,:,:)/c(1,:,:)
			!!!L1～L5の設定/使う条件に合わせて変更すること
			!!!計算領域x負側
			!!!亜音速流入境界条件
			L(1,0,:,:)=(G(2,0,:,:)-c(0,:,:))*(-G(1,0,:,:)*c(0,:,:)*dGx(2,0,:,:)+dGx(5,0,:,:))
			L(5,0,:,:)=L(1,0,:,:)
			L(2,0,:,:)=0.5d0*(gamma-1d0)*(L(1,0,:,:)+L(5,0,:,:))
!			L(2,0,:,:)=G(2,0,:,:)*(c(0,:,:)*c(0,:,:)*dGx(1,0,:,:)-dGx(5,0,:,:))
!			L(3,0,:,:)=G(2,0,:,:)*dGx(3,0,:,:)
!			L(4,0,:,:)=G(2,0,:,:)*dGx(4,0,:,:)
!			L(5,0,:,:)=(G(2,0,:,:)+c(0,:,:))*(G(1,0,:,:)*c(0,:,:)*dGx(2,0,:,:)+dGx(5,0,:,:))
!			L(5,0,:,:)=(sigma*(1d0-Ma_NSCBC(0,:,:)*Ma_NSCBC(0,:,:))*c(0,:,:)*(G(5,0,:,:)-p_l))/LX
			!!!d行列の計算
			d(1,0,:,:)=(1d0/(c(0,:,:)*c(0,:,:)))*((L(1,0,:,:)+L(5,0,:,:))*0.5d0+L(2,0,:,:))
!			d(2,0,:,:)=(L(1,0,:,:)+L(5,0,:,:))*0.5d0
!			d(3,0,:,:)=(1d0/(c(0,:,:)*G(1,0,:,:)))*(-L(1,0,:,:)+L(5,0,:,:))*0.5d0
!			d(4,0,:,:)=L(3,0,:,:)
!			d(5,0,:,:)=L(4,0,:,:)
			!!!NSCBCの適用
			dFx(1,0,:,:)=d(1,0,:,:)
!			dFx(2,0,:,:)=G(2,0,:,:)*d(1,0,:,:)+G(1,0,:,:)*d(3,0,:,:)
!			dFx(3,0,:,:)=G(3,0,:,:)*d(1,0,:,:)+G(1,0,:,:)*d(4,0,:,:)
!			dFx(4,0,:,:)=G(4,0,:,:)*d(1,0,:,:)+G(1,0,:,:)*d(5,0,:,:)
!			dFx(5,0,:,:)=(G(2,0,:,:)*G(2,0,:,:)+G(3,0,:,:)*G(3,0,:,:)+G(4,0,:,:)*G(4,0,:,:))*d(1,0,:,:)*0.5d0	&
!						+d(2,0,:,:)/(gamma-1d0)	&
!						+G(1,0,:,:)*G(2,0,:,:)*d(3,0,:,:)	&
!						+G(1,0,:,:)*G(3,0,:,:)*d(4,0,:,:)	&
!						+G(1,0,:,:)*G(4,0,:,:)*d(5,0,:,:)
			!!!計算領域x正側
			!!!無反射流出条件
!			L(1,1,:,:)=(G(2,NX,:,:)-c(1,:,:))*(-G(1,NX,:,:)*c(1,:,:)*dGx(2,NX,:,:)+dGx(5,NX,:,:))
			L(1,1,:,:)=(sigma*(1d0-Ma_NSCBC(1,:,:)*Ma_NSCBC(1,:,:))*c(1,:,:)*(G(5,NX,:,:)-p0))/LX
			L(2,1,:,:)=G(2,NX,:,:)*(c(1,:,:)*c(1,:,:)*dGx(1,NX,:,:)-dGx(5,NX,:,:))
			L(3,1,:,:)=G(2,NX,:,:)*dGx(3,NX,:,:)
			L(4,1,:,:)=G(2,NX,:,:)*dGx(4,NX,:,:)
			L(5,1,:,:)=(G(2,NX,:,:)+c(1,:,:))*(G(1,NX,:,:)*c(1,:,:)*dGx(2,NX,:,:)+dGx(5,NX,:,:))
			!!!d行列の計算
			d(1,1,:,:)=(1d0/(c(1,:,:)*c(1,:,:)))*((L(1,1,:,:)+L(5,1,:,:))*0.5d0+L(2,1,:,:))
			d(2,1,:,:)=(L(1,1,:,:)+L(5,1,:,:))*0.5d0
			d(3,1,:,:)=(1d0/(c(1,:,:)*G(1,NX,:,:)))*(-L(1,1,:,:)+L(5,1,:,:))*0.5d0
			d(4,1,:,:)=L(3,1,:,:)
			d(5,1,:,:)=L(4,1,:,:)
			!!!NSCBCの適用
			dFx(1,NX,:,:)=d(1,1,:,:)
			dFx(2,NX,:,:)=G(2,NX,:,:)*d(1,1,:,:)+G(1,NX,:,:)*d(3,1,:,:)
			dFx(3,NX,:,:)=G(3,NX,:,:)*d(1,1,:,:)+G(1,NX,:,:)*d(4,1,:,:)
			dFx(4,NX,:,:)=G(4,NX,:,:)*d(1,1,:,:)+G(1,NX,:,:)*d(5,1,:,:)
			dFx(5,NX,:,:)=(G(2,NX,:,:)*G(2,NX,:,:)+G(3,NX,:,:)*G(3,NX,:,:)+G(4,NX,:,:)*G(4,NX,:,:))*d(1,1,:,:)*0.5d0	&
						+d(2,1,:,:)/(gamma-1d0)	&
						+G(1,NX,:,:)*G(2,NX,:,:)*d(3,1,:,:)	&
						+G(1,NX,:,:)*G(3,NX,:,:)*d(4,1,:,:)	&
						+G(1,NX,:,:)*G(4,NX,:,:)*d(5,1,:,:)
	end subroutine NSCBC_x

	subroutine NSCBC_y(G,dGy,dFy)
		implicit none
		double precision :: LY
		double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: G,dGy,dFy
		double precision, dimension(1:5,0:NX,0:1,0:NZ-1) :: L,d
		double precision, dimension(0:NX,0:1,0:NZ-1) :: c,Ma_NSCBC
			LY=dy*dble(NY)
			c(:,0,:)=sqrt(gamma*G(5,:,0,:)/G(1,:,0,:))
			c(:,1,:)=sqrt(gamma*G(5,:,NY,:)/G(1,:,NY,:))
			Ma_NSCBC(:,0,:)=G(2,:,0,:)/c(:,0,:)
			Ma_NSCBC(:,1,:)=G(2,:,NY,:)/c(:,1,:)
			!!!L1～L5の設定/使う条件に合わせて変更すること
			!!!計算領域y負側
			!!!等温滑りなし壁面条件
			L(1,:,0,:)=(G(3,:,0,:)-c(:,0,:))*(-G(1,:,0,:)*c(:,0,:)*dGy(3,:,0,:)+dGy(5,:,0,:))
			L(2,:,0,:)=G(3,:,0,:)*(c(:,0,:)*c(:,0,:)*dGy(1,:,0,:)-dGy(5,:,0,:))
!			L(2,:,0,:)=0d0
			L(3,:,0,:)=G(3,:,0,:)*dGy(2,:,0,:)
			L(4,:,0,:)=G(3,:,0,:)*dGy(4,:,0,:)
!			L(5,:,0,:)=(G(3,:,0,:)+c(:,0,:))*(G(1,:,0,:)*c(:,0,:)*dGy(3,:,0,:)+dGy(5,:,0,:))
			L(5,:,0,:)=(sigma*(1d0-Ma_NSCBC(:,0,:)*Ma_NSCBC(:,0,:))*c(:,0,:)*(G(5,:,0,:)-p0))/LY
!			L(5,:,0,:)=L(1,:,0,:)
			!!!d行列の計算
			d(1,:,0,:)=(1d0/(c(:,0,:)*c(:,0,:)))*((L(1,:,0,:)+L(5,:,0,:))*0.5d0+L(2,:,0,:))
			d(2,:,0,:)=(L(1,:,0,:)+L(5,:,0,:))*0.5d0
			d(3,:,0,:)=L(3,:,0,:)
			d(4,:,0,:)=(1d0/(c(:,0,:)*G(1,:,0,:)))*(-L(1,:,0,:)+L(5,:,0,:))*0.5d0
			d(5,:,0,:)=L(4,:,0,:)
			!!!NSCBCの適用
			dFy(1,:,0,:)=d(1,:,0,:)
			dFy(2,:,0,:)=G(2,:,0,:)*d(1,:,0,:)+G(1,:,0,:)*d(3,:,0,:)
			dFy(3,:,0,:)=G(3,:,0,:)*d(1,:,0,:)+G(1,:,0,:)*d(4,:,0,:)
			dFy(4,:,0,:)=G(4,:,0,:)*d(1,:,0,:)+G(1,:,0,:)*d(5,:,0,:)
			dFy(5,:,0,:)=(G(2,:,0,:)*G(2,:,0,:)+G(3,:,0,:)*G(3,:,0,:)+G(4,:,0,:)*G(4,:,0,:))*d(1,:,0,:)*0.5d0	&
						+d(2,:,0,:)/(gamma-1d0)	&
						+G(1,:,0,:)*G(2,:,0,:)*d(3,:,0,:)	&
						+G(1,:,0,:)*G(3,:,0,:)*d(4,:,0,:)	&
						+G(1,:,0,:)*G(4,:,0,:)*d(5,:,0,:)
			!!!計算領域y正側
			!!!等温滑りなし壁面条件
!			L(1,:,1,:)=(G(3,:,NY,:)-c(:,1,:))*(-G(1,:,NY,:)*c(:,1,:)*dGy(3,:,NY,:)+dGy(5,:,NY,:))
			L(1,:,1,:)=(sigma*(1d0-Ma_NSCBC(:,1,:)*Ma_NSCBC(:,1,:))*c(:,1,:)*(G(5,:,NY,:)-p0))/LY
			L(2,:,1,:)=G(3,:,NY,:)*(c(:,1,:)*c(:,1,:)*dGy(1,:,NY,:)-dGy(5,:,NY,:))
!			L(2,:,1,:)=0d0
			L(3,:,1,:)=G(3,:,NY,:)*dGy(2,:,NY,:)
			L(4,:,1,:)=G(3,:,NY,:)*dGy(4,:,NY,:)
			L(5,:,1,:)=(G(3,:,NY,:)+c(:,1,:))*(G(1,:,NY,:)*c(:,1,:)*dGy(3,:,NY,:)+dGy(5,:,NY,:))
!			L(1,:,1,:)=L(5,:,1,:)
			!!!d行列の計算
			d(1,:,1,:)=(1d0/(c(:,1,:)*c(:,1,:)))*((L(1,:,1,:)+L(5,:,1,:))*0.5d0+L(2,:,1,:))
			d(2,:,1,:)=(L(1,:,1,:)+L(5,:,1,:))*0.5d0
			d(3,:,1,:)=L(3,:,1,:)
			d(4,:,1,:)=(1d0/(c(:,1,:)*G(1,:,NY,:)))*(-L(1,:,1,:)+L(5,:,1,:))*0.5d0
			d(5,:,1,:)=L(4,:,1,:)
			!!!NSCBCの適用
			dFy(1,:,NY,:)=d(1,:,1,:)
			dFy(2,:,NY,:)=G(2,:,NY,:)*d(1,:,1,:)+G(1,:,NY,:)*d(3,:,1,:)
			dFy(3,:,NY,:)=G(3,:,NY,:)*d(1,:,1,:)+G(1,:,NY,:)*d(4,:,1,:)
			dFy(4,:,NY,:)=G(4,:,NY,:)*d(1,:,1,:)+G(1,:,NY,:)*d(5,:,1,:)
			dFy(5,:,NY,:)=(G(2,:,NY,:)*G(2,:,NY,:)+G(3,:,NY,:)*G(3,:,NY,:)+G(4,:,NY,:)*G(4,:,NY,:))*d(1,:,1,:)*0.5d0	&
						+d(2,:,1,:)/(gamma-1d0)	&
						+G(1,:,NY,:)*G(2,:,NY,:)*d(3,:,1,:)	&
						+G(1,:,NY,:)*G(3,:,NY,:)*d(4,:,1,:)	&
						+G(1,:,NY,:)*G(4,:,NY,:)*d(5,:,1,:)
	end subroutine NSCBC_y

	subroutine boundary_overwrite(Q,u_jet,xi)
		implicit none
		integer :: i,k
		double precision, dimension(0:NY) :: y,xi,u_jet
		double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: Q
			!!!亜音速流入境界条件
			Q(3,0,:,:)=0d0
			Q(4,0,:,:)=0d0
			do k=0,NZ-1
				do i=0,NY
					Q(2,0,i,k)=Q(1,0,i,k)*u_jet(i)
					Q(5,0,i,k)=(Q(1,0,i,k)*T_in/(gamma*Ma*Ma))/(gamma-1d0)+Q(1,0,i,k)*u_jet(i)*u_jet(i)*0.5d0
				end do
			end do
			!!!等温滑りなし壁面条件
!			Q(2,:,0,:)=0d0
!			Q(3,:,0,:)=0d0
!			Q(4,:,0,:)=0d0
!			Q(5,:,0,:)=(Q(1,:,0,:)*T_wall/(gamma*Ma*Ma))/(gamma-1d0)
!			Q(2,:,NY,:)=0d0
!			Q(3,:,NY,:)=0d0
!			Q(4,:,NY,:)=0d0
!			Q(5,:,NY,:)=(Q(1,:,NY,:)*T_wall/(gamma*Ma*Ma))/(gamma-1d0)
	end subroutine boundary_overwrite

	subroutine result(ts,G,T_k,xi)
		implicit none
		integer :: i,j,k
		integer, intent(in) :: ts
		double precision, intent(in), dimension(0:NY) :: xi
		double precision, intent(in), dimension(0:NX,0:NY,0:NZ-1) :: T_k
		double precision, intent(in), dimension(1:5,0:NX,0:NY,0:NZ-1) :: G
		character(len=24) charname
			charname='3dim0000000XXXYYYZZZ.txt'
			write(charname(5:11),"(i7.7)") ts
			write(charname(12:14),"(i3.3)") NX
			write(charname(15:17),"(i3.3)") NY
			write(charname(18:20),"(i3.3)") NZ
			open(10,file=charname,status='unknown',form='formatted')
			do i=0,NX
				do j =0,NY
!					write(10,"(2f10.4,6e25.16)")(i-NX/2)*dx,(j-NY/2)*dy,G(1,i,j,NZ/2),G(2,i,j,NZ/2),G(3,i,j,NZ/2),G(4,i,j,NZ/2),G(5,i,j,NZ/2),	&
!									T_k(i,j,NZ/2)
!					write(10,"(2f10.4,6e25.16)")(i-NX/2)*dx,(j-NY/2)*dy,G(1,i,j,0),G(2,i,j,0),G(3,i,j,0),G(4,i,j,0),G(5,i,j,0),	&
!									T_k(i,j,0)
					write(10,"(2f10.4,6e25.16)")(i-NX/2)*dx,xi(j),G(1,i,j,0),G(2,i,j,0),G(3,i,j,0),G(4,i,j,0),G(5,i,j,0),	&
									T_k(i,j,0)
				end do
				write(10,*)
			end do
			close(10)
	end subroutine result
end module sub

!!!メインプログラム
program main
	use sub
	implicit none
	integer :: h,i,j,k,ts,ii,jj,kk
	double precision :: cpu_t1,cpu_t2,t
	double precision, dimension(0:NY) :: y,dydxi,xi,u_jet
	double precision, dimension(0:NX,0:NY,0:NZ-1) :: T_k
	double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: G,dGx,dGy,dGz,Q,Q1,Q2,dQ,Fxp,Fxm,Fyp,Fym,Fzp,Fzm
	double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: dFxp,dFxm,dFyp,dFym,dFzp,dFzm,dFx,dFy,dFz
	double precision, dimension(1:5,0:NX,0:NY,0:NZ-1) :: Fvx,Fvy,Fvz,dFvx,dFvy,dFvz
	double precision, dimension(-1:1,0:NX) :: LU_D_p_NX,LU_D_m_NX,LU_C_NX
	double precision, dimension(-1:1,0:NY) :: LU_D_p_NY,LU_D_m_NY,LU_C_NY
	double precision, dimension(0:NZ-1,0:NZ-1) :: L_D_p_NZ,L_D_m_NZ,L_C_NZ,U_D_p_NZ,U_D_m_NZ,U_C_NZ
		!!!初期値設定
		G(1,:,:,:)=rho0
		G(2,:,:,:)=u0
		G(3,:,:,:)=v0
		G(4,:,:,:)=w0
		G(5,:,:,:)=p0
		!!!格子伸長の関数
		do i=0,NY
			y(i)=dy*dble(i)-20d0
			xi(i)=(1.7d0*y(i)-(1d0/10d0)*(-log(cosh(5d0*(y(i)-5d0)))+log(cosh(5d0*(y(i)+5d0)))))/1.4d0
			dydxi(i)=1.4d0/(1.7d0-(0.5d0*(-tanh(5d0*(y(i)-5d0))+tanh(5d0*(y(i)+5d0)))))
		end do
		!!!"top-hat型jet"
		do i=NY/2,NY
			u_jet(i)=(u1/2d0)*(1d0-tanh((12.5d0/4d0)*(xi(i)-(1d0/xi(i)))))
		end do
		do i=0,(NY/2)-1
			u_jet(i)=u_jet(NY-i)
		end do
		call G_Q(G,Q)
		!!!LU行列の作成と圧縮
		call LU(0,NX,LU_D_p_NX,-sigma)
		call LU(0,NY,LU_D_p_NY,-sigma)
		call LU_PBC(0,NZ,L_D_p_NZ,U_D_p_NZ,-sigma)
		call LU(0,NX,LU_D_m_NX,sigma)
		call LU(0,NY,LU_D_m_NY,sigma)
		call LU_PBC(0,NZ,L_D_m_NZ,U_D_m_NZ,sigma)
		call LU(0,NX,LU_C_NX,0d0)
		call LU(0,NY,LU_C_NY,0d0)
		call LU_PBC(0,NZ,L_C_NZ,U_C_NZ,0d0)
		call cpu_time(cpu_t1)
		T_k(:,:,:)=(gamma*Ma*Ma*G(5,:,:,:))/G(1,:,:,:)
		call result(0,G,T_k,xi)
		!!!以下計算ループ開始
		do ts=1,NT*it
			t=dt*dble(ts)
			!!!-------------------------------------------------------------------------------------------------------------------
			!!!1段目
			!!!-------------------------------------------------------------------------------------------------------------------
			!!!移流項
			call Q_Fx_LF(Q,Fxp,Fxm)
			call Q_Fy_LF(Q,Fyp,Fym)
			call Q_Fz_LF(Q,Fzp,Fzm)
			do j=0,NY
				do k=0,NZ-1
					call DCS(0,NX,dx,Fxp(:,:,j,k),dFxp(:,:,j,k),LU_D_p_NX,-sigma)
					call DCS(0,NX,dx,Fxm(:,:,j,k),dFxm(:,:,j,k),LU_D_m_NX,sigma)
				end do
			end do
			do i=0,NX
				do k=0,NZ-1
					call DCS(0,NY,dy,Fyp(:,i,:,k),dFyp(:,i,:,k),LU_D_p_NY,-sigma)
					call DCS(0,NY,dy,Fym(:,i,:,k),dFym(:,i,:,k),LU_D_m_NY,sigma)
					do h=1,5
						dFyp(h,i,:,k)=dFyp(h,i,:,k)*dydxi(:)
						dFym(h,i,:,k)=dFym(h,i,:,k)*dydxi(:)
					end do
				end do
			end do
			do i=0,NX
				do j=0,NY
					call DCS_PBC(0,NZ,dz,Fzp(:,i,j,:),dFzp(:,i,j,:),L_D_p_NZ,U_D_p_NZ,-sigma)
					call DCS_PBC(0,NZ,dz,Fzm(:,i,j,:),dFzm(:,i,j,:),L_D_m_NZ,U_D_m_NZ,sigma)
				end do
			end do
			dFx=dFxp+dFxm
			dFy=dFyp+dFym
			dFz=dFzp+dFzm
			!!!NSCBC
			do h=1,5
				do j=0,NY
					do k=0,NZ-1
						call CCS(0,NX,dx,G(h,:,j,k),dGx(h,:,j,k),LU_C_NX)
					end do
				end do
				do i =0,NX
					do k=0,NZ-1
						call CCS(0,NY,dy,G(h,i,:,k),dGy(h,i,:,k),LU_C_NY)
						dGy(h,i,:,k)=dGy(h,i,:,k)*dydxi(:)
					end do
				end do
				do i=0,NX
					do j=0,NY
						call CCS_PBC(0,NZ,dz,G(h,i,j,:),dGz(h,i,j,:),L_C_NZ,U_C_NZ)
					end do
				end do
			end do
			call NSCBC_x(G,dGx,dFx)
			call NSCBC_y(G,dGy,dFy)
			!!!拡散項
			call Q_Fvx_Fvy_Fvz(G,dGx,dGy,dGz,Fvx,Fvy,Fvz,LU_C_NX,LU_C_NY,L_C_NZ,U_C_NZ,dydxi)
			do h=1,5
				do j=0,NY
					do k=0,NZ-1
						call CCS(0,NX,dx,Fvx(h,:,j,k),dFvx(h,:,j,k),LU_C_NX)
					end do
				end do
				do i =0,NX
					do k=0,NZ-1
						call CCS(0,NY,dy,Fvy(h,i,:,k),dFvy(h,i,:,k),LU_C_NY)
						dFvy(h,i,:,k)=dFvy(h,i,:,k)*dydxi(:)
					end do
				end do
				do i=0,NX
					do j=0,NY
						call CCS_PBC(0,NZ,dz,Fvz(h,i,j,:),dFvz(h,i,j,:),L_C_NZ,U_C_NZ)
					end do
				end do
			end do
			!!!時間進展
			Q1=Q+dt*(dFvx+dFvy+dFvz-dFx-dFy-dFz)
!			Q1=Q+dt*(-dFx-dFy-dFz)
			!!!境界条件
			call boundary_overwrite(Q1,u_jet,xi)
			!!!-------------------------------------------------------------------------------------------------------------------
			!!!2段目
			!!!-------------------------------------------------------------------------------------------------------------------
			call Q_G(Q1,G)
			!!!移流項
			call Q_Fx_LF(Q1,Fxp,Fxm)
			call Q_Fy_LF(Q1,Fyp,Fym)
			call Q_Fz_LF(Q1,Fzp,Fzm)
			do j=0,NY
				do k=0,NZ-1
					call DCS(0,NX,dx,Fxp(:,:,j,k),dFxp(:,:,j,k),LU_D_p_NX,-sigma)
					call DCS(0,NX,dx,Fxm(:,:,j,k),dFxm(:,:,j,k),LU_D_m_NX,sigma)
				end do
			end do
			do i=0,NX
				do k=0,NZ-1
					call DCS(0,NY,dy,Fyp(:,i,:,k),dFyp(:,i,:,k),LU_D_p_NY,-sigma)
					call DCS(0,NY,dy,Fym(:,i,:,k),dFym(:,i,:,k),LU_D_m_NY,sigma)
					do h=1,5
						dFyp(h,i,:,k)=dFyp(h,i,:,k)*dydxi(:)
					    dFym(h,i,:,k)=dFym(h,i,:,k)*dydxi(:)
					end do
				end do
			end do
			do i=0,NX
				do j=0,NY
					call DCS_PBC(0,NZ,dz,Fzp(:,i,j,:),dFzp(:,i,j,:),L_D_p_NZ,U_D_p_NZ,-sigma)
					call DCS_PBC(0,NZ,dz,Fzm(:,i,j,:),dFzm(:,i,j,:),L_D_m_NZ,U_D_m_NZ,sigma)
				end do
			end do
			dFx=dFxp+dFxm
			dFy=dFyp+dFym
			dFz=dFzp+dFzm
			!!!NSCBC
			do h=1,5
				do j=0,NY
					do k=0,NZ-1
						call CCS(0,NX,dx,G(h,:,j,k),dGx(h,:,j,k),LU_C_NX)
					end do
				end do
				do i =0,NX
					do k=0,NZ-1
						call CCS(0,NY,dy,G(h,i,:,k),dGy(h,i,:,k),LU_C_NY)
						dGy(h,i,:,k)=dGy(h,i,:,k)*dydxi(:)
					end do
				end do
				do i=0,NX
					do j=0,NY
						call CCS_PBC(0,NZ,dz,G(h,i,j,:),dGz(h,i,j,:),L_C_NZ,U_C_NZ)
					end do
				end do
			end do
			call NSCBC_x(G,dGx,dFx)
			call NSCBC_y(G,dGy,dFy)
			!!!拡散項
			call Q_Fvx_Fvy_Fvz(G,dGx,dGy,dGz,Fvx,Fvy,Fvz,LU_C_NX,LU_C_NY,L_C_NZ,U_C_NZ,dydxi)
			do h=1,5
				do j=0,NY
					do k=0,NZ-1
						call CCS(0,NX,dx,Fvx(h,:,j,k),dFvx(h,:,j,k),LU_C_NX)
					end do
				end do
				do i =0,NX
					do k=0,NZ-1
						call CCS(0,NY,dy,Fvy(h,i,:,k),dFvy(h,i,:,k),LU_C_NY)
						dFvy(h,i,:,k)=dFvy(h,i,:,k)*dydxi(:)
					end do
				end do
				do i=0,NX
					do j=0,NY
						call CCS_PBC(0,NZ,dz,Fvz(h,i,j,:),dFvz(h,i,j,:),L_C_NZ,U_C_NZ)
					end do
				end do
			end do
			!!!時間進展
			Q2=(3d0/4d0)*Q+(1d0/4d0)*Q1+(dt/4d0)*(dFvx+dFvy+dFvz-dFx-dFy-dFz)
!			Q2=(3d0/4d0)*Q+(1d0/4d0)*Q1+(dt/4d0)*(-dFx-dFy-dFz)
			!!!境界条件
			call boundary_overwrite(Q2,u_jet,xi)
			!!!-------------------------------------------------------------------------------------------------------------------
			!!!3段目
			!!!-------------------------------------------------------------------------------------------------------------------
			call Q_G(Q2,G)
			!!!移流項
			call Q_Fx_LF(Q2,Fxp,Fxm)
			call Q_Fy_LF(Q2,Fyp,Fym)
			call Q_Fz_LF(Q2,Fzp,Fzm)
			do j=0,NY
				do k=0,NZ-1
					call DCS(0,NX,dx,Fxp(:,:,j,k),dFxp(:,:,j,k),LU_D_p_NX,-sigma)
					call DCS(0,NX,dx,Fxm(:,:,j,k),dFxm(:,:,j,k),LU_D_m_NX,sigma)
				end do
			end do
			do i=0,NX
				do k=0,NZ-1
					call DCS(0,NY,dy,Fyp(:,i,:,k),dFyp(:,i,:,k),LU_D_p_NY,-sigma)
					call DCS(0,NY,dy,Fym(:,i,:,k),dFym(:,i,:,k),LU_D_m_NY,sigma)
					do h=1,5
						dFyp(h,i,:,k)=dFyp(h,i,:,k)*dydxi(:)
						dFym(h,i,:,k)=dFym(h,i,:,k)*dydxi(:)
					end do
				end do
			end do
			do i=0,NX
				do j=0,NY
					call DCS_PBC(0,NZ,dz,Fzp(:,i,j,:),dFzp(:,i,j,:),L_D_p_NZ,U_D_p_NZ,-sigma)
					call DCS_PBC(0,NZ,dz,Fzm(:,i,j,:),dFzm(:,i,j,:),L_D_m_NZ,U_D_m_NZ,sigma)
				end do
			end do
			dFx=dFxp+dFxm
			dFy=dFyp+dFym
			dFz=dFzp+dFzm
			!!!NSCBC
			do h=1,5
				do j=0,NY
					do k=0,NZ-1
						call CCS(0,NX,dx,G(h,:,j,k),dGx(h,:,j,k),LU_C_NX)
					end do
				end do
				do i =0,NX
					do k=0,NZ-1
						call CCS(0,NY,dy,G(h,i,:,k),dGy(h,i,:,k),LU_C_NY)
						dGy(h,i,:,k)=dGy(h,i,:,k)*dydxi(:)
					end do
				end do
				do i=0,NX
					do j=0,NY
						call CCS_PBC(0,NZ,dz,G(h,i,j,:),dGz(h,i,j,:),L_C_NZ,U_C_NZ)
					end do
				end do
			end do
			call NSCBC_x(G,dGx,dFx)
			call NSCBC_y(G,dGy,dFy)
			!!!拡散項
			call Q_Fvx_Fvy_Fvz(G,dGx,dGy,dGz,Fvx,Fvy,Fvz,LU_C_NX,LU_C_NY,L_C_NZ,U_C_NZ,dydxi)
			do h=1,5
				do j=0,NY
					do k=0,NZ-1
						call CCS(0,NX,dx,Fvx(h,:,j,k),dFvx(h,:,j,k),LU_C_NX)
					end do
				end do
				do i =0,NX
					do k=0,NZ-1
						call CCS(0,NY,dy,Fvy(h,i,:,k),dFvy(h,i,:,k),LU_C_NY)
						dFvy(h,i,:,k)=dFvy(h,i,:,k)*dydxi(:)
					end do
				end do
				do i=0,NX
					do j=0,NY
						call CCS_PBC(0,NZ,dz,Fvz(h,i,j,:),dFvz(h,i,j,:),L_C_NZ,U_C_NZ)
					end do
				end do
			end do
			!!!時間進展
			dQ=(1d0/3d0)*Q+(2d0/3d0)*Q2+((2d0*dt)/3d0)*(dFvx+dFvy+dFvz-dFx-dFy-dFz)
!			dQ=(1d0/3d0)*Q+(2d0/3d0)*Q2+((2d0*dt)/3d0)*(-dFx-dFy-dFz)
			!!!境界条件
			call boundary_overwrite(dQ,u_jet,xi)
			call Q_G(dQ,G)
			write(*,*) ts
			!!!計算破産時の緊急出力
			do i=0,NX
				do j=0,NY
					do k=0,NZ-1
						if(isnan(dQ(1,i,j,k))) then
							write(*,*) "NaN"
							call Q_G(Q,G)
							T_k(:,:,:)=(gamma*Ma*Ma*G(5,:,:,:))/G(1,:,:,:)
							call result(ts,G,T_k,xi)
							stop
						end if
					end do
				end do
			end do
			!!!出力の条件分岐
			if(mod(ts,Nout)==0) then
				T_k(:,:,:)=(gamma*Ma*Ma*G(5,:,:,:))/G(1,:,:,:)
				call result(ts,G,T_k,xi)
				call cpu_time(cpu_t2)
				write(*,*) "cpu_time : ",cpu_t2-cpu_t1, "seconds"
			end if
			Q=dQ
		end do
end program main
