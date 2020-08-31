program average3D
integer ::  MX,MY,MZ,num,count,count2,count3
double precision :: ad,bd,cd,a,b,c,omega0,a3,ad3,Re
double precision :: alpha,beta,alphaD,betaD,alpha3,alphaD3
double precision :: ymax,ymin
double precision, allocatable:: frhoup(:,:,:)
double precision, allocatable:: uzu_W(:,:,:),uzuW_ave(:,:)
double precision, allocatable:: tauw(:),cf(:),AAuave_y(:,:),LUuave_y(:,:)
double precision, allocatable:: QRuave_y(:,:,:),guave_y(:,:,:)
double precision, allocatable:: x(:),y(:),z(:),ys(:),yys(:),ysy(:)
double precision, allocatable:: rhoup(:,:,:,:),Tuv(:,:,:,:)
double precision, allocatable :: grhoupx(:,:,:,:),grhoupy(:,:,:,:)
double precision, allocatable :: grhoupz(:,:,:,:)
double precision, allocatable :: rho_read(:,:),u_read(:,:),Rex(:)


num=5
MX=500 ;MY=100 ;MZ=20 ;MXt=100
!'分割時間dt'
dt=1.d0/200d0

allocate(frhoup(1:6,0:MX,0:MY))
allocate(uzu_W(0:MX,0:MY,0:MZ),uzuW_ave(0:MX,0:MY))
allocate(tauw(0:MX),cf(0:MX),AAuave_y(0:MY,0:MY),LUuave_y(0:MY,0:MY))
allocate(QRuave_y(1:6,0:MX,0:MY),guave_y(1:6,0:MX,0:MY))
allocate(x(0:MX),y(0:MY),z(0:MZ),ys(0:MY),yys(0:MY),ysy(0:MY))
allocate(rhoup(1:num,0:MX,0:MY,0:MZ),Tuv(1:num,0:MX,0:MY,0:MZ))
allocate(grhoupx(1:num,0:MX,0:MY,0:MZ),grhoupy(1:num,0:MX,0:MY,0:MZ))
allocate(grhoupz(1:num,0:MX,0:MY,0:MZ))
allocate(rho_read(0:MX,0:MY),u_read(0:MX,0:MY),Rex(0:MX))


tmax=250d0/dt
pi=dacos(-1.0d0)
c=1d0
a=14.d0/9.d0
b=1.d0/9.d0
c=0.d0
ad=4.d0/9.d0
bd=2.d0/9.d0
cd=0.d0
!cc=1.d0
alpha=1.d0/3.d0
beta=0.d0
alphaD=1.d0
betaD=0.d0
omegaNS=0.25d0
omega0=0.d0
omega1=0.25d0
omega2=-0.25d0
a3=1.5d0
ad3=0.5d0
alpha=1.d0/3.d0
alpha3=0.25d0
alphaD3=1.d0

xmax=100d0
xmin=0d0
ymax=15d0
ymin=0d0
zmax=2d0
zmin=0d0
hx=(xmax-xmin)/dble(MX)
hy=(ymax-ymin)/dble(MY)
hz=(zmax-zmin)/dble(MZ)
dQ1=0.d0
dQ2=0.d0
Q=0.d0

gam=1.4d0
zeta=1.d0

Re=1000.d0
Pr=1.d0
Ma=0.3d0

Pmugen=1.d0/gam/(Ma**2.d0)
Lx=(xmax-xmin)
Ly=(ymax-ymin)
Lz=(zmax-zmin)

count=0 ;count2=0 ;count3=0

a11=0.13d0


!!x,y,zの定義-----------------------------------------------------------------------------------------------------------------------
do i=0,MX
   x(i)=(dble(MX-i)*xmin+dble(i)*xmax)/dble(MX)
enddo
do j=0,MY
	ys(j)=(dble(MY-j)*ymin+dble(j)*ymax)/dble(MY)
enddo

y=ymin+(ymax - ymin)*dexp(-a11*(ymax - ymin - ys)) - (ymax - ymin - ys)*dexp(-a11*(ymax - ymin))  !格子伸長(ここでのyは資料式(3.49)のygsを指す)
y(0)=0.d0

yys=a11*(ymax - ymin)*dexp(-a11*(ymax - ymin - ys))+dexp(-a11*(ymax - ymin)) !上の式をyで微分した値dy/dygs)
ysy=1.d0/yys																 !逆数を取ったdygs/dy

do k=0,MZ
	z(k)=(dble(MZ-k)*zmin+dble(k)*zmax)/dble(MZ)
enddo
!!-----------------------------------------------------------------------------------------------------------------------------------
!!shooting法による初期値の導入-------------------------------------------------------------------------------------------------------
!open(21,file='r_2d.txt',status='old')
!open(22,file='u_2d.txt',status='old')
!open(25,file='rhooo.txt',status='unknown')   
!
!do i=0,MX
!	do j=0,MY
!	read(21,*) rho_read(i,j)
!	read(22,*) u_read(i,j)
!	write(25,*) rho_read(i,j)
!	enddo
!enddo
!
!close(21)
!close(22)
!close(25)
!!-----------------------------------------------------------------------------------------------------------------------------------
!!平均化するρ,u,v,w,P,Tの導入-------------------------------------------------------------------------------------------------------
open(501,file='3D_ρ.txt',status='old')
open(502,file='3D_u.txt',status='old')
open(503,file='3D_v.txt',status='old')
open(504,file='3D_w.txt',status='old')
open(505,file='3D_P.txt',status='old')
open(506,file='3D_T.txt',status='old')

do k=0,MZ
	do j=0,MY
		do i=0,MX
		read(501,*) rhoup(1,i,j,k)
		read(502,*) rhoup(2,i,j,k)
		read(503,*) rhoup(3,i,j,k)
		read(504,*) rhoup(4,i,j,k)
		read(505,*) rhoup(5,i,j,k)
		read(506,*) Tuv(5,i,j,k)
		enddo
	enddo
enddo
close(501)
close(502)
close(503)
close(504)
close(505)
close(506)		

!!------------------------------------------------------------------------------------------------------------------------------------
!!平均化開始--------------------------------------------------------------------------------------------------------------------------
call favre_ave(MX,MY,MZ,num,x,y,rhoup,Tuv,frhoup,count)

call uzu(MX,MY,MZ,num,x,y,grhoupx,grhoupy,grhoupz,uzu_W,uzuW_ave,count2)

call cf_tauw(MX,MY,MZ,x,y,Re,a,b,c,ad,bd,cd,omega0,alpha,beta,alphaD,betaD,alpha3,alphaD3,ymax,ymin,ysy&
			&,rho_read,u_read,frhoup,tauw,cf,count3)

!!------------------------------------------------------------------------------------------------------------------------------------
end program average3D

!!アンサンブル平均-----------------------------------------------------------------------
subroutine an_ave(MX,MY,MZ,num,rhoup,Tuv)
implicit none
integer ::  MX,MY,MZ,num
double precision :: sum_rho(0:MX,0:MY),sum_u(0:MX,0:MY),sum_v(0:MX,0:MY),sum_T(0:MX,0:MY)
double precision :: rhoup(1:num,0:MX,0:MY,0:MZ),Tuv(1:num,0:MX,0:MY,0:MZ)
double precision :: rho_ave(0:MX,0:MY),u_ave(0:MX,0:MY),v_ave(0:MX,0:MY),T_ave(0:MX,0:MY)
integer i,j,k

sum_rho=0.d0; sum_u=0d0; sum_v=0d0; sum_T=0d0

do k=0,MZ
sum_rho=sum_rho+rhoup(1,:,:,k)
sum_u=sum_u+rhoup(2,:,:,k)
sum_v=sum_v+rhoup(3,:,:,k)
sum_T=sum_T+Tuv(5,:,:,k)
enddo

rho_ave(:,:)=sum_rho(:,:)/MZ
u_ave(:,:)=sum_u(:,:)/MZ
v_ave(:,:)=sum_v(:,:)/MZ
T_ave(:,:)=sum_T(:,:)/MZ

open(531,file='3D_ρave.txt',status='unknown')
open(532,file='3D_uave.txt',status='unknown')
open(533,file='3D_vave.txt',status='unknown')
open(534,file='3D_Tave.txt',status='unknown')
do j=0,MY
	do i=0,MX
	write(531,*) rho_ave(i,j)
	write(532,*) u_ave(i,j)
	write(533,*) v_ave(i,j)
	write(534,*) T_ave(i,j)
	enddo
enddo

close(531)
close(532)
close(533)
close(534)

end subroutine an_ave





!!!ファーブル平均-----------------------------------------------------------------------------------------------------
subroutine favre_ave(MX,MY,MZ,num,x,y,rhoup,Tuv,frhoup,count)
implicit none
integer ::  MX,MY,MZ,num
double precision :: rho(0:MX,0:MY),sum_rho(0:MX,0:MY),sum_rhou(0:MX,0:MY),sum_rhov(0:MX,0:MY)
double precision :: sum_rhow(0:MX,0:MY),sum_rhoP(0:MX,0:MY),sum_rhoT(0:MX,0:MY)
double precision :: rhoup(1:num,0:MX,0:MY,0:MZ),Tuv(1:num,0:MX,0:MY,0:MZ)
double precision :: frhoup(1:6,0:MX,0:MY)
double precision :: x(0:MX),y(0:MY)
character(len=16) :: Chara
integer :: count,i,j,k

rho=0.d0
sum_rho=0.d0; sum_rhou=0d0; sum_rhov=0d0; sum_rhow=0d0
sum_rhop=0d0 ;sum_rhoT=0.d0

do k=1,MZ
rho=rho+rhoup(1,:,:,k)
sum_rho=sum_rho+rhoup(1,:,:,k)*rhoup(1,:,:,k)
sum_rhou=sum_rhou+rhoup(2,:,:,k)*rhoup(1,:,:,k)
sum_rhov=sum_rhov+rhoup(3,:,:,k)*rhoup(1,:,:,k)
sum_rhow=sum_rhow+rhoup(4,:,:,k)*rhoup(1,:,:,k)
sum_rhoP=sum_rhop+rhoup(5,:,:,k)*rhoup(1,:,:,k)
sum_rhoT=sum_rhoT+Tuv(5,:,:,k)*rhoup(1,:,:,k)
enddo

rho=rho/MZ
frhoup(1,:,:)=sum_rho/MZ
frhoup(2,:,:)=sum_rhou/MZ
frhoup(3,:,:)=sum_rhov/MZ
frhoup(4,:,:)=sum_rhow/MZ
frhoup(5,:,:)=sum_rhoP/MZ
frhoup(6,:,:)=sum_rhoT/MZ

frhoup(1,:,:)=frhoup(1,:,:)/rho(:,:)
frhoup(2,:,:)=frhoup(2,:,:)/rho(:,:)
frhoup(3,:,:)=frhoup(3,:,:)/rho(:,:)
frhoup(4,:,:)=frhoup(4,:,:)/rho(:,:)
frhoup(5,:,:)=frhoup(5,:,:)/rho(:,:)
frhoup(6,:,:)=frhoup(6,:,:)/rho(:,:)


count=count+1
write(Chara, '(i4.4)') count
open(41,file='ρ_heikin'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(42,file='u_heikin'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(43,file='v_heikin'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(44,file='w_heikin'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(45,file='P_heikin'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(46,file='T_heikin'//trim(Chara)//'.txt',status='unknown',form='formatted')


do j=0,MY
	do i=0,MX
	write(41,'(3f24.16)') x(i),y(j),frhoup(1,i,j)
	write(42,'(3f24.16)') x(i),y(j),frhoup(2,i,j)
	write(43,'(3f24.16)') x(i),y(j),frhoup(3,i,j)
	write(44,'(3f24.16)') x(i),y(j),frhoup(4,i,j)
	write(45,'(3f24.16)') x(i),y(j),frhoup(5,i,j)
	write(46,'(3f24.16)') x(i),y(j),frhoup(6,i,j)
	enddo
	write(41,*)
	write(42,*)
	write(43,*)
	write(44,*)
	write(45,*)
	write(46,*)	
enddo
	close(41)
	close(42)
	close(43)
	close(44)
	close(45)
	close(46)


end subroutine favre_ave

!!渦構造の可視化用--------------------------------------------------------------------------------
subroutine uzu(MX,MY,MZ,num,x,y,grhoupx,grhoupy,grhoupz,uzu_W,uzuW_ave,count2)
implicit none
integer ::  MX,MY,MZ,num
double precision :: grhoupx(1:num,0:MX,0:MY,0:MZ),grhoupy(1:num,0:MX,0:MY,0:MZ)
double precision :: grhoupz(1:num,0:MX,0:MY,0:MZ)
double precision :: x(0:MX),y(0:MY)
double precision :: uzu_W(0:MX,0:MY,0:MZ),uzuW_ave(0:MX,0:MY)
double precision :: sum_uzu(0:MX,0:MY)
character(len=16) :: Chara
integer :: count2,i,j,k


uzu_W(:,:,:)=grhoupx(2,:,:,:)*grhoupy(3,:,:,:)+grhoupy(3,:,:,:)*grhoupz(4,:,:,:)&
     &+grhoupz(4,:,:,:)*grhoupx(2,:,:,:)-(grhoupz(2,:,:,:)*grhoupx(4,:,:,:)&
     &+grhoupy(4,:,:,:)*grhoupz(3,:,:,:)+grhoupx(3,:,:,:)*grhoupy(2,:,:,:))

sum_uzu=0.d0
!!z方向の平均を取る---------------------------
do k=0,MZ-1
sum_uzu(:,:)=sum_uzu(:,:)+uzu_W(:,:,k)
enddo

uzuW_ave=sum_uzu/MZ
!!--------------------------------------------

count2=count2+1
write(Chara, '(i4.4)') count2
open(1000,file='uzu_W'//trim(Chara)//'.txt',status='unknown',form='formatted')
do j=0,MY
	do i=0,MX
	write(1000,'(3f24.16)') x(i),y(j),uzuW_ave(i,j)
	enddo
	write(1000,*)
enddo

close(1000)
end subroutine uzu

!!!サブルーチン左辺定義(Aの定義)---------------------------------------------------------------
subroutine Ax(MX,AA,alpha,beta,alphaD,betaD,omega,alpha3,alphaD3)
integer :: MX
double precision :: alpha,beta,alphaD,betaD,omega,alpha3,alphaD3
double precision :: AA(0:MX,0:MX)
integer i

AA=0.0d0

do i=0,MX
	AA(i,i)=1.0d0
enddo

	AA(0,1)=3.0d0                !片側差分
	AA(MX,MX-1)=3.0d0
	
	AA(1,0)=alpha3*(1.d0-alphaD3*omega)  !3次DCS
	AA(1,2)=alpha3*(1.d0+alphaD3*omega)
	AA(MX-1,MX)=alpha3*(1.d0+alphaD3*omega)
	AA(MX-1,MX-2)=alpha3*(1.d0-alphaD3*omega)

do i=2,MX-2								!5次DCS
	AA(i,i+1)=alpha*(1.d0+omega*alphaD)
	AA(i,i+2)=beta*(1.d0+omega*betaD)
	AA(i,i-1)=alpha*(1.d0-omega*alphaD)
	AA(i,i-2)=beta*(1.d0-omega*betaD)
enddo

end subroutine Ax

!!!サブルーチンLU分解 x＆y＆z方向共有 そして、LU圧縮---------------------------------------------------------------------------
subroutine LU(MX,AA,LUcom)
implicit none
integer :: MX
double precision :: AA(0:MX,0:MX)
double precision  :: L(0:MX,0:MX),U(0:MX,0:MX)
double precision :: LUcom(-2:2,0:MX)
double precision s1,s2
integer i,j,k
L=0.0d0               
U=0.0d0
do i=0,MX              !式開始
L(i,i)=1.0d0
enddo
	do i=0,MX		
		do j=0,MX
		s1=0.0d0
		s2=0.0d0
	   	if(i<=j) then        !上三角行列を作る(U)
	   		do k=0,i-1 		
	   		s1=s1+L(i,k)*U(k,j)
	   		enddo
	   	U(i,j)=AA(i,j)-s1
	   	else if(i>=j) then   !下三角行列を作る(L)
	   		do k=0,j-1
	   		s2=s2+L(i,k)*U(k,j)
	   		enddo
	   	L(i,j)=(AA(i,j)-s2)/U(j,j)
	   	endif
	   	enddo
	enddo
	
	do i = 0,MX                   !0を抜いて圧縮
		LUcom(0,i) = U(i,i)	
	enddo	
	do i = 1,MX-1
		LUcom(-1,i) = L(i,i-1)	
		LUcom(1,i) = U(i,i+1)
	enddo
	
	LUcom(-1,0)= 0
	LUcom(1,0)=U(0,1)
	LUcom(-1,MX)=L(MX,MX-1)
	LUcom(1,MX)=0

end subroutine LU

!!!サブルーチン右辺定義 y方向(行列x,y成分のみ) LUの前進後退------------------------------------------------------------------------
subroutine RHS2y(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega,ymax,ymin,QR,F,a3,ad3,LUcom,dF)
integer :: MX,MY,MZ,num
double precision :: ad,bd,cd,a,b,c,omega,a3,ad3
double precision :: f(1:num,0:MX,0:MY)
double precision :: QR(1:num,0:MX,0:MY)
double precision :: LUcom(-2:2,0:MY)
double precision :: Y(1:num,0:MX,0:MY)
double precision :: dF(1:num,0:MX,0:MY)
double precision :: ymax,ymin
double precision h
integer j,i

h=(ymax-ymin)/dble(MY)

	QR(:,:,0)=((-17.d0/6.d0)*f(:,:,0)+1.5d0*f(:,:,1)+1.5d0*f(:,:,2)-(1.0d0/6.d0)*f(:,:,3))/h
	QR(:,:,MY)=((1.d0/6.d0)*f(:,:,MY-3)-1.5d0*f(:,:,MY-2)-1.5d0*f(:,:,MY-1)+(17.d0/6.d0)*f(:,:,MY))/h


	QR(:,:,1)=a3*(-f(:,:,0)+f(:,:,2))/(2.d0*h)+omega*(ad3*(f(:,:,0)-2.d0*f(:,:,1)+f(:,:,2))/h)
	QR(:,:,MY-1)=a3*(-f(:,:,MY-2)+f(:,:,MY))/(2.d0*h)+omega*(ad3*(f(:,:,MY-2)-2.d0*f(:,:,MY-1)+f(:,:,MY))/h)

		do i=2,MY-2	
		QR(:,:,i)=a*(-f(:,:,i-1)+f(:,:,i+1))/(2.d0*h)+b*(-f(:,:,i-2)+f(:,:,i+2))/(4.d0*h)&
	     &+omega*(ad*(f(:,:,i-1)-2.d0*f(:,:,i)+f(:,:,i+1))/h+bd*(f(:,:,i-2)-2.d0*f(:,:,i)+f(:,:,i+2))/(4.d0*h))
		enddo 

!!  前進代入
	Y(:,:,0)=QR(:,:,0)
	do j=1,MY
		Y(:,:,j)=QR(:,:,j)-LUcom(-1,j) * Y(:,:,j-1)
	enddo
!!  後退代入
	dF(:,:,MY)=Y(:,:,MY)/LUcom(0,MY)
	dF(:,:,MY-1)=(Y(:,:,MY-1)-dF(:,:,MY)*LUcom(1,MY-1))/LUcom(0,MY-1)
	do j=MY-2,0,-1
		dF(:,:,j)=(Y(:,:,j)-LUcom(1,j)*dF(:,:,j+1)-LUcom(2,j)*dF(:,:,j+2))/LUcom(0,j)
	enddo

end subroutine RHS2y


!!せん断応力及び摩擦係数Cf----------------------------------------------------------------------------------
subroutine cf_tauw(MX,MY,MZ,x,y,Re,a,b,c,ad,bd,cd,omega,alpha,beta,alphaD,betaD,alpha3,alphaD3,ymax,ymin,ysy&
                  &,rho_read,u_read,frhoup,tauw,cf,count3)
implicit none
integer :: MX,MY,MZ
double precision :: ad,bd,cd,a,b,c,omega,a3,ad3
double precision :: alpha,beta,alphaD,betaD,alpha3,alphaD3
double precision :: ymax,ymin
double precision :: ysy(0:MY),rho_read(0:MX,0:MY),u_read(0:MX,0:MY),frhoup(1:6,0:MX,0:MY),x(0:MX),y(0:MY)
double precision :: tauw(0:MX),cf(0:MX)
double precision :: AAuave_y(0:MY,0:MY),LUuave_y(0:MY,0:MY)
double precision :: QRuave_y(1:6,0:MX,0:MY),guave_y(1:6,0:MX,0:MY)
double precision :: Re
double precision :: mu1
double precision :: Rex(0:MX),Cftheo(0:MX)
double precision :: uplus(0:MX,0:MY),yplus(0:MX,0:MY),Ut(0:MX)
character(len=16) :: Chara
integer :: i,j,o,count3

mu1=1.d0/Re

call Ax(MY,AAuave_y,alpha,beta,alphaD,betaD,0d0,alpha3,alphaD3)
call LU(MY,AAuave_y,LUuave_y)
call RHS2y(MX,MY,MZ,6,a,b,c,ad,bd,cd,0d0,ymax,ymin,QRuave_y,frhoup,a3,ad3,LUuave_y,guave_y)
do o=1,6
	do i=0,MX
	guave_y(o,i,:)=guave_y(o,i,:)*ysy(:)    !格子伸長を適用
	enddo
enddo

tauw=mu1*guave_y(2,:,0)       !!μ*du/dy

Ut=dsqrt(tauw/frhoup(1,:,0))

do j=0,MY
uplus(:,j)=frhoup(2,:,j)/Ut(:)
enddo

do i=0,MX
yplus(i,:)=(Ut(i)*y(:)*frhoup(1,i,:)/mu1)
enddo

open(231,file='bunsekiy.csv')
do i=0,MX
	write(231,*)'x=', x(i)
	do j=0,MY
	write(231,*) yplus(i,j),',',uplus(i,j)
	enddo
	write(231,*)
enddo
close(231)

cf(:)=tauw(:)/(0.5d0*rho_read(:,0)*(u_read(:,0)**2.d0))

do i=0,MX
Rex(i)=(Re**2.d0)/(1.721**2.d0)+Re*x(i)
enddo

Cftheo(:)=0.455d0/((dlog10(Rex(:)))**2.58d0)    !!理論解

count3=count3+1
write(Chara, '(i4.4)') count3
open(2000,file='tauw'//trim(Chara)//'.csv',status='unknown',form='formatted')
open(2001,file='cf'//trim(Chara)//'.csv',status='unknown',form='formatted')
	do i=0,MX
	write(2000,'(3f24.16)') x(i),tauw(i)
	write(2001,'(3f24.16)') Rex(i),cftheo(i),cf(i)
	enddo
close(2000)
close(2001)


end subroutine cf_tauw
