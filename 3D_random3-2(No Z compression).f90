


!!���i��A(1,1)�̂悤��1����X�^�[�g�����Ă������A�����ł�A(0,0)�X�^�[�g�ƂȂ��Ă���
!!�������AF��Q�̎���Q(1,0)�X�^�[�g�ł���B�܂�Q(1,0),Q(2,0),Q(3,0)�����ꂼ��̎��̍ŏ��̓_��\��



!!���L
!!�܂��T�u���[�`������(�{����j
!!num�̒ǉ�!!!!!!!!!
!!�T�u���[�`���E�Ӓ�` x,y,z�����I��
!!�T�u���[�`��A�s���LU������A��L,U�̓����̐���ς����ȊO�A���ɉ��������I��
!!�T�u���[�`���O�i��� x��y��z�����I��
!!�T�u���[�`��x��y��z���� F(FP��FM)�̒�`�I��
!!�T�u���[�`��!!!��,u,P���`����I��
!!�T�u���[�`��NSCBC��x��y��z�����I��(Lz�̒ǉ��y�т����͉�����)
!!T,u,v�̔����T�u���[�`���E�Ӓ�`�͏�L�̃T�u���[�`���E�Ӓ�`�ŕ₤
!!�T�u���[�`��du��dT����S�����̓��o����єS�x�ʂ̓��o x��y��z�����I��
!!�T�u���[�`��LU���� x��y��z�������L �����āALU���k�̏I��
!!�T�u���[�`���E�Ӓ�` ����y��z����  LU�O�i��ނ̏I��
!!�ȉ��{���J�n

!!hx��hy���T�u���[�`�����ɓK�p���邽�߂ɁA�T�u���[�`���̘g����xmax,xmin,ymax,ymin�̓���
!!Tuv��T��Ma**2�����Y�ꂽ
!!Pr=1.d0,Ma=3d0�ɕύX
!!Pmugen�̕ύX(1.d0/gam/(Ma**2.d0))
!!x��y������NSCBC�̏d�Ȃ镔���̒���
!!�f�B���N���ǉ�(������)(�����shooting�p����Ȃ�)
!!Re��1000

!!���̌�A�����ǖʂɊւ��邱�Ƃ��s��
!!P��Ma�o��
!!�f�B���N���ǉ�(�����ǖ�)
!!mu�̕ύX(Pr=1.d0������T�u���[�`���̂悤�Ȏ�)
!!T���t�@�C���ŏo��

!!�i�q�L����K�p����
!!�{���̃J�E���g��o����times�ɕύX
!!rho_read,u_read,v_read,P_read,T_read,w_read�̒ǉ�

!!�h���̃��[�h �ǉ�
!!�����ύX�E�i�q�_�ύX�E�������ԕύX
!!tmax�̕ύX
!!Q���`�����w_read��0�Ȃ̂ŁA�����Ȃ������B(�ł��邾��0�Ȃ͕̂ϐ��Ƃ��ď����K�v������)
!!y-z�ʂł̝�����MXt�̓o��irandom_3D��MX�̕ύX)
!!�����̃T�u���[�`���y�уR�[�����̓���
!!MX,MY�Ȃǂ����[�h�ɂ����A�l�����̂܂ܖ{���ɑł�����
!!u_aver��u_margin�𓊓��y�т��̃e�L�X�g�t�@�C���̍쐬��{���ɒǉ�(�s���_)
!!�{����sum_uz��0�N���A��Y��Ă���
!!�{����s�̕���(if���̒�)��s=s+1�Ƃ���
!!�h������������read�ł͂Ȃ�write�ł�������C�ԍ��ԈႦ����Cold�ł͂Ȃ�unknown�ɂ��Ă���
!!Q�̏����l�C�T�u���[�`�������������D(x=0�Œ�ōs�����)
!!TVD�̂��ꂼ��̒i��call kakuran_Q�𓱓������D�܂�if���̒���call kakuran_Q���O�ɏo����
!!���֐��̑�����~�߁C���ۂɑ��֐����h���Ɏ��t����
!!�f�B���N���s��������
!!�T�u���[�`��kakuran_Q��Q(5)��ω�������(�S��Q�ŕ\�킷�悤��)

!!3�����h���̊ԈႢ����!--------->z�����O�i��ނ̍s�񐔂�4�ɂȂ��Ă����Dnum�ɂ�������ƕς����D

!!���́C���Ϗ�̍쐬!!
!!�t�@�[�u�����ς̃T�u���[�`���쐬(�s��)
!!�t�@�[�u�����ς̃T�u���[�`����i,j�͈̔͂�1����MX,MY�������̂�0����ɕς���
!!�t�@�[�u�����ς̃T�u���[�`����close�������Ă����̂Œ���
!!�Q�\���̉����p�̃T�u���[�`���̍쐬�i�K���͂܂�)
!!�T�u���[�`���E�Ӓ�` y����(�s��x,y�����̂�) LU�̑O�i��ނ̒ǉ��̍쐬(���C�W���₹��f���͗p)
!!�T�u���[�`������f���͋y�і��C�W��Cf�̍쐬

!!�h���̕ύX(120�X�e�b�v���ƂɊh����1�i�߂Ă������C1000�X�e�b�v(5s����)�ɕύX)
!!80s��150s�ɕύX
!!�T�u���[�`���t�@�[���u�����ς̒���write(41,'(3f24.16)')�Ȃǂ�write(41,*)�ɂ���(���������1/21�̒i�K�ł͌v�Z���񂵂Ă��Ȃ�)
!!�T�u���[�`������f���͋y�і��C�W��Cf��cf�̕ύX(rho_read�̓���)
!!�T�u���[�`������f���͋y�і��C�W��Cf��tauw�̃t�@�C���o�͗p��open�Ȃǂ���ꂽ
!!�T�u���[�`������f���͋y�і��C�W��Cf��cf�̕ύX(u_read�̓���)
!!count,count2,count3��0���(�J�E���g�̂���) (1/22)
!!Rex�̓����C�����cf�̗��_�l�̓���(1/23)
!!�{����5s���Ƃ�3�����p�̃�,u,v,w,P,T�̃e�L�X�g�������(�����ǂݍ���ŕ��ω��ɂȂ��邽��)
!!count,count2,count3��0���(�J�E���g�̂���)���T�u���[�`���ł͂Ȃ��{���Ɉړ�(1/24)

!!call�y��subroutine�̍s��ς���ۂ�(,)�J���}�̈ʒu�ɒ���(1/29)(������)
!!tauw��1.4�{�s����(������)(1/29)


program DCS
implicit none
double precision, allocatable:: AAMx(:,:),LMx(:,:),UMx(:,:)
double precision, allocatable:: AAMy(:,:),LMy(:,:),UMy(:,:)
double precision hx,hy,hz,pi,s1,s2,s3,s4,omegaNS,omega,omega2,omega1,a,b,c,ad,bd,cd,alpha,beta,alphaD,betaD,cc
double precision, allocatable:: an(:),L6(:),x(:),y(:),z(:),Bound1(:),BoundMX(:)
double precision, allocatable:: Q(:,:,:,:),dQ1(:,:,:,:),dQ2(:,:,:,:)
double precision, allocatable:: Fx(:,:,:,:),FPx(:,:,:,:),FMx(:,:,:,:),Fy(:,:,:,:),FPy(:,:,:,:),FMy(:,:,:,:)
double precision, allocatable:: QRMx(:,:,:,:),QRPx(:,:,:,:),gMx(:,:,:,:),gPx(:,:,:,:)
double precision, allocatable:: QRMy(:,:,:,:),QRPy(:,:,:,:),gMy(:,:,:,:),gPy(:,:,:,:)
double precision, allocatable:: AAPx(:,:),UPx(:,:),LPx(:,:)
double precision, allocatable:: AAPy(:,:),UPy(:,:),LPy(:,:)
double precision, allocatable:: AATux(:,:),LTux(:,:),UTux(:,:)
double precision, allocatable:: AATuy(:,:),LTuy(:,:),UTuy(:,:)
double precision, allocatable:: Tuv(:,:,:,:),TuRHSx(:,:,:,:),TuRHSy(:,:,:,:),gTux(:,:,:,:),gTuy(:,:,:,:)
double precision, allocatable:: Tuvisx(:,:,:,:),Tuvisy(:,:,:,:)
double precision, allocatable:: TuRHS2x(:,:,:,:),TuRHS2y(:,:,:,:),ggTux(:,:,:,:),ggTuy(:,:,:,:)
double precision, allocatable:: gx(:,:,:,:),gy(:,:,:,:),rhoup(:,:,:,:),rhoupRHSx(:,:,:,:),grhoupx(:,:,:,:)
double precision, allocatable:: LNSx(:,:),UNSx(:,:),AANSx(:,:)
double precision, allocatable:: c0x(:,:),cMXx(:,:),Ma0x(:,:),MaMXx(:,:),NSLx(:,:,:,:),NSdx(:,:,:,:),mu(:,:,:)
double precision, allocatable:: AANSy(:,:),LNSy(:,:),UNSy(:,:)
double precision, allocatable:: rhoupRHSy(:,:,:,:),grhoupy(:,:,:,:)
double precision, allocatable:: c0y(:,:),cMXy(:,:),Ma0y(:,:),MaMXy(:,:),NSLy(:,:,:,:),NSdy(:,:,:,:)
double precision, allocatable:: LUMx(:,:),LUMy(:,:),LUPx(:,:),LUPy(:,:),LUNSx(:,:),LUNSy(:,:),LUTux(:,:),LUTuy(:,:)
double precision, allocatable:: LUMz(:,:),LUPz(:,:),LUNSz(:,:),LUTuz(:,:)
double precision, allocatable:: AAMz(:,:),AAPz(:,:),AANSz(:,:),AATuz(:,:)
double precision, allocatable:: Fz(:,:,:,:),FPz(:,:,:,:),FMz(:,:,:,:),QRMz(:,:,:,:),gMz(:,:,:,:)
double precision, allocatable:: QRPz(:,:,:,:),gPz(:,:,:,:),gz(:,:,:,:)
double precision, allocatable:: rhoupRHSz(:,:,:,:),grhoupz(:,:,:,:)
double precision, allocatable:: c0z(:,:),cMXz(:,:),Ma0z(:,:),MaMXz(:,:),NSLz(:,:,:,:),NSdz(:,:,:,:)
double precision, allocatable:: TuRHSz(:,:,:,:),gTuz(:,:,:,:),Tuvisz(:,:,:,:),TuRHS2z(:,:,:,:),ggTuz(:,:,:,:)
double precision, allocatable:: LMz(:,:),UMz(:,:)
double precision, allocatable:: LPz(:,:),UPz(:,:)
double precision, allocatable:: LNSz(:,:),UNSz(:,:)
double precision, allocatable:: LTuz(:,:),UTuz(:,:)
double precision, allocatable:: ys(:),rho(:,:,:),P(:,:,:),v(:,:,:),u(:,:,:),w(:,:,:),T(:,:,:),Sc(:,:,:)
double precision, allocatable:: yys(:),ysy(:)
double precision, allocatable:: rho_read(:,:),P_read(:,:),v_read(:,:),u_read(:,:),T_read(:,:),w_read(:,:)
double precision, allocatable:: kakuran_u(:,:,:),kakuran_v(:,:,:),kakuran_w(:,:,:)
double precision, allocatable:: u_ave(:,:),u_margin(:,:,:),sum_uz(:,:)
double precision, allocatable:: window(:,:)
double precision, allocatable:: frhoup(:,:,:)
double precision, allocatable:: uzu_W(:,:,:),uzuW_ave(:,:)
double precision, allocatable:: tauw(:),cf(:),AAuave_y(:,:),LUuave_y(:,:)
double precision, allocatable:: QRuave_y(:,:,:),guave_y(:,:,:)
double precision xmax,xmin,ymax,ymin,zmax,zmin
double precision a3,ad3,alpha3,alphaD3
integer MX,k,o,times
double precision delx,delt,d,dt,Lth,err6,LL6
double precision gam,zeta
double precision Re,Pr,Ma
integer N,i,j,im1,im2,im3,im4,ip1,ip2,ip3,im0,tmax
double precision omega0,Pmugen,Lx,Ly,Lz,a11
integer MY,filecount
integer MXt,s
character(len=16) Chara
double precision CT1,CT2
integer num,MZ,count,count2,count3

num=5

!write(*,*)'������x'
!read(*,*) MX
!write(*,*) '������y'
!read(*,*) MY
!write(*,*) '������z'
!read(*,*) MZ
!write(*,*) 'y-z�ʂł̝�����'
!read(*,*) MXt

MX=200 ;MY=70 ;MZ=50 ;MXt=100
!'��������dt'
dt=1.d0/125d0

allocate(AAMx(0:MX,0:MX),LMx(0:MX,0:MX),UMx(0:MX,0:MX))
allocate(AAMy(0:MY,0:MY),LMy(0:MY,0:MY),UMy(0:MY,0:MY))
allocate(x(0:MX),y(0:MY),z(0:MZ))
allocate(Q(1:num,0:MX,0:MY,0:MZ),dQ1(1:num,0:MX,0:MY,0:MZ),dQ2(1:num,0:MX,0:MY,0:MZ))
allocate(Fx(1:num,0:MX,0:MY,0:MZ),FMx(1:num,0:MX,0:MY,0:MZ),FPx(1:num,0:MX,0:MY,0:MZ))
allocate(Fy(1:num,0:MX,0:MY,0:MZ),FMy(1:num,0:MX,0:MY,0:MZ),FPy(1:num,0:MX,0:MY,0:MZ))
allocate(QRMx(1:num,0:MX,0:MY,0:MZ),QRPx(1:num,0:MX,0:MY,0:MZ),gMx(1:num,0:MX,0:MY,0:MZ))
allocate(gPx(1:num,0:MX,0:MY,0:MZ),AAPx(0:MX,0:MX))
allocate(QRMy(1:num,0:MX,0:MY,0:MZ),QRPy(1:num,0:MX,0:MY,0:MZ),gMy(1:num,0:MX,0:MY,0:MZ))
allocate(gPy(1:num,0:MX,0:MY,0:MZ))
allocate(LPx(0:MX,0:MX),UPx(0:MX,0:MX))
allocate(AAPy(0:MY,0:MY),LPy(0:MY,0:MY),UPy(0:MY,0:MY))
allocate(Bound1(-2:MX+2),BoundMX(-2:MX+2))
allocate(AATux(0:MX,0:MX),LTux(0:MX,0:MX),UTux(0:MX,0:MX))
allocate(AATuy(0:MY,0:MY),LTuy(0:MY,0:MY),UTuy(0:MY,0:MY))
allocate(Tuv(1:num,0:MX,0:MY,0:MZ),TuRHSx(1:num,0:MX,0:MY,0:MZ),TuRHSy(1:num,0:MX,0:MY,0:MZ))
allocate(gTux(1:num,0:MX,0:MY,0:MZ),gTuy(1:num,0:MX,0:MY,0:MZ))
allocate(TuRHS2x(1:num,0:MX,0:MY,0:MZ),TuRHS2y(1:num,0:MX,0:MY,0:MZ))
allocate(Tuvisx(1:num,0:MX,0:MY,0:MZ),Tuvisy(1:num,0:MX,0:MY,0:MZ))
allocate(ggTux(1:num,0:MX,0:MY,0:MZ),ggTuy(1:num,0:MX,0:MY,0:MZ))
allocate(gx(1:num,0:MX,0:MY,0:MZ),gy(1:num,0:MX,0:MY,0:MZ),rhoup(1:num,0:MX,0:MY,0:MZ))
allocate(rhoupRHSx(1:num,0:MX,0:MY,0:MZ),grhoupx(1:num,0:MX,0:MY,0:MZ))
allocate(AANSx(0:MX,0:MX),LNSx(0:MX,0:MX),UNSx(0:MX,0:MX))
allocate(NSLx(1:num,0:MX,0:MY,0:MZ),NSdx(1:num,0:MX,0:MY,0:MZ))
allocate(c0x(0:MY,0:MZ),cMXx(0:MY,0:MZ),Ma0x(0:MY,0:MZ),MaMXx(0:MY,0:MZ))
allocate(AANSy(0:MY,0:MY),LNSy(0:MY,0:MY),UNSy(0:MY,0:MY))
allocate(rhoupRHSy(1:num,0:MX,0:MY,0:MZ),grhoupy(1:num,0:MX,0:MY,0:MZ))
allocate(c0y(0:MX,0:MZ),cMXy(0:MX,0:MZ),Ma0y(0:MX,0:MZ),MaMXy(0:MX,0:MZ))
allocate(NSLy(1:num,0:MX,0:MY,0:MZ),NSdy(1:num,0:MX,0:MY,0:MZ),mu(0:MX,0:MY,0:MZ))
allocate(LUMx(-2:2,0:MX),LUMy(-2:2,0:MY),LUPx(-2:2,0:MX),LUPy(-2:2,0:MY),LUNSx(-2:2,0:MX))
allocate(LUNSy(-2:2,0:MY),LUTux(-2:2,0:MX),LUTuy(-2:2,0:MY))
allocate(LUMz(-2:2,0:MZ),LUPz(-2:2,0:MZ),LUNSz(-2:2,0:MZ),LUTuz(-2:2,0:MZ))
allocate(AAMz(0:MZ,0:MZ),AAPz(0:MZ,0:MZ),AANSz(0:MZ,0:MZ),AATuz(0:MZ,0:MZ))
allocate(Fz(1:num,0:MX,0:MY,0:MZ),FPz(1:num,0:MX,0:MY,0:MZ),FMz(1:num,0:MX,0:MY,0:MZ))
allocate(QRMz(1:num,0:MX,0:MY,0:MZ),gMz(1:num,0:MX,0:MY,0:MZ))
allocate(QRPz(1:num,0:MX,0:MY,0:MZ),gPz(1:num,0:MX,0:MY,0:MZ),gz(1:num,0:MX,0:MY,0:MZ))
allocate(rhoupRHSz(1:num,0:MX,0:MY,0:MZ),grhoupz(1:num,0:MX,0:MY,0:MZ))
allocate(c0z(0:MX,0:MY),cMXz(0:MX,0:MY),Ma0z(0:MX,0:MY),MaMXz(0:MX,0:MY))
allocate(NSLz(1:num,0:MX,0:MY,0:MZ),NSdz(1:num,0:MX,0:MY,0:MZ))
allocate(TuRHSz(1:num,0:MX,0:MY,0:MZ),gTuz(1:num,0:MX,0:MY,0:MZ),Tuvisz(1:num,0:MX,0:MY,0:MZ))
allocate(TuRHS2z(1:num,0:MX,0:MY,0:MZ),ggTuz(1:num,0:MX,0:MY,0:MZ))
allocate(LMz(0:MZ,0:MZ),UMz(0:MZ,0:MZ))
allocate(LPz(0:MZ,0:MZ),UPz(0:MZ,0:MZ))
allocate(LNSz(0:MZ,0:MZ),UNSz(0:MZ,0:MZ))
allocate(LTuz(0:MZ,0:MZ),UTuz(0:MZ,0:MZ))
allocate(ys(0:MY),rho(0:MX,0:MY,0:MZ),P(0:MX,0:MY,0:MZ),v(0:MX,0:MY,0:MZ),u(0:MX,0:MY,0:MZ))
allocate(w(0:MX,0:MY,0:MZ),T(0:MX,0:MY,0:MZ),Sc(0:MX,0:MY,0:MZ))
allocate(yys(0:MY),ysy(0:MY))
allocate(rho_read(0:MX,0:MY),P_read(0:MX,0:MY),v_read(0:MX,0:MY),u_read(0:MX,0:MY),T_read(0:MX,0:MY))
allocate(w_read(0:MX,0:MY))
allocate(kakuran_u(0:MXt,0:MY,0:MZ),kakuran_v(0:MXt,0:MY,0:MZ),kakuran_w(0:MXt,0:MY,0:MZ))
allocate(u_margin(0:MX,0:MY,0:MZ),u_ave(0:MX,0:MY),sum_uz(0:MX,0:MY))
allocate(window(0:MY,0:MZ))
allocate(frhoup(1:6,0:MX,0:MY))
allocate(uzu_W(0:MX,0:MY,0:MZ),uzuW_ave(0:MX,0:MY))
allocate(tauw(0:MX),cf(0:MX),AAuave_y(0:MY,0:MY),LUuave_y(0:MY,0:MY))
allocate(QRuave_y(1:6,0:MX,0:MY),guave_y(1:6,0:MX,0:MY))



tmax=1000d0/dt
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

xmax=150d0
xmin=0d0
ymax=8d0
ymin=0d0
zmax=12d0
zmin=0d0
hx=(xmax-xmin)/dble(MX)
hy=(ymax-ymin)/dble(MY)
hz=(zmax-zmin)/dble(MZ)
dQ1=0.d0
dQ2=0.d0
Q=0.d0

gam=1.4d0
zeta=1.d0

Sc(:,:,:)=120.d0/(273.15d0+T(:,:,:))  !!�̂���
Re=1150.d0
Pr=1.d0
Ma=0.3d0

Pmugen=1.d0/gam/(Ma**2.d0)
Lx=(xmax-xmin)
Ly=(ymax-ymin)
Lz=(zmax-zmin)

a11=0.13d0

QRMx=0.d0 ;QRMy=0.d0 ;QRPx=0.d0 ;QRPy=0.d0 ;gMx=0.d0 ;gMy=0.d0 ;gPx=0.d0 ;gPy=0.d0
AAMx=0.d0 ;AAMy=0.d0 ;AAPx=0.d0 ;AAPy=0.d0 ;UMx=0.d0 ;UMy=0.d0 ;UPx=0.d0 ;UPy=0.d0
LMx=0.d0 ;LMy=0.d0 ;LPx=0.d0 ;LPy=0d0 ;Q=0.d0 ;Fx=0.d0 ;Fy=0.d0 ;FPx=0.d0 ;FPy=0.d0 ;FMx=0.d0 ;FMy=0.d0

AATux=0.d0 ;AATuy=0.d0 ;LTux=0.d0 ;LTuy=0.d0 ;UTux=0.d0 ;UTuy=0.d0 ;Tuv=0.d0 ;gTux=0.d0 ;gTuy=0.d0
TuRHSx=0.d0 ;TuRHSy=0.d0 ;TuRHS2x=0.d0 ;TuRHS2y=0.d0 ;Tuvisx=0.d0 ;Tuvisy=0.d0 ;ggTux=0.d0 ;ggTuy=0.d0

gx=0.d0 ;gy=0.d0 ;rhoup=0.d0 ;rhoupRHSx=0.d0 ;grhoupx=0.d0 ;AANSx=0.d0 ;LNSx=0.d0 ;UNSx=0.d0 
c0x=0.d0 ;cMXx=0.d0 ;Ma0x=0.d0 ;MaMXx=0.d0 ;NSLx=0.d0 ;NSdx=0.d0 ;AANSy=0.d0 ;LNSy=0.d0 ;UNSy=0.d0
rhoupRHSy=0.d0 ;grhoupy=0.d0 ;c0y=0.d0 ;cMXy=0.d0 ;Ma0y=0.d0 ;MaMXy=0.d0 ;NSLy=0.d0 ;NSdy=0.d0 ;mu=0d0

LUMx=0.d0 ;LUMy=0.d0 ;LUPx=0.d0 ;LUPy=0.d0 ;LUNSx=0.d0 ;LUNSy=0.d0 ;LUTux=0.d0 ;LUTuy=0.d0
LUMz=0.d0 ;LUPz=0.d0 ;LUNSz=0.d0 ;LUTuz=0.d0 ;AAMz=0.d0 ;AAPz=0.d0 ;AANSz=0.d0 ;AATuz=0.d0 
Fz=0.d0 ;FPz=0.d0 ;FMz=0.d0 ;QRMz=0.d0 ;gMz=0.d0
QRPz=0.d0 ;gPz=0.d0 ;gz=0.d0 ;rhoupRHSz=0.d0 ;grhoupz=0.d0
c0z=0.d0 ;cMXz=0.d0 ;Ma0z=0.d0 ;MaMXz=0.d0 ;NSLz=0.d0 ;NSdz=0.d0
TuRHSz=0.d0 ;gTuz=0.d0 ;Tuvisz=0.d0 ;TuRHS2z=0.d0 ;ggTuz=0.d0

ys=0.d0 ;rho=0.d0 ;P=0.d0 ;v=0.d0 ;u=0.d0 ;w=0.d0 

LMz=0.d0 ;UMz=0.d0
LPz=0.d0 ;UPz=0.d0
LNSz=0.d0 ;UNSz=0.d0 
LTuz=0.d0 ;UTuz=0.d0

yys=0.d0 ;ysy=0.d0

rho_read=0.d0 ;P_read=0.d0 ;v_read=0.d0 ;u_read=0.d0 ;T_read=0.d0 ;w_read=0.d0

kakuran_u=0.d0 ;kakuran_v=0.d0 ;kakuran_w=0.d0

u_ave=0.d0 ;u_margin=0.d0 ;s=2 ;sum_uz=0.d0

window=0.d0 ;times=0.d0

count=0 ;count2=0 ;count3=0

filecount=0
call cpu_time(CT1)

!!x,y,z�̒�`-----------------------------------------------------------------------------------------------------------------------
do i=0,MX
   x(i)=(dble(MX-i)*xmin+dble(i)*xmax)/dble(MX)
enddo
do j=0,MY
	ys(j)=(dble(MY-j)*ymin+dble(j)*ymax)/dble(MY)
enddo

y=ymin+(ymax - ymin)*dexp(-a11*(ymax - ymin - ys)) - (ymax - ymin - ys)*dexp(-a11*(ymax - ymin))  !�i�q�L��(�����ł�y�͎�����(3.49)��ygs���w��)
y(0)=0.d0

yys=a11*(ymax - ymin)*dexp(-a11*(ymax - ymin - ys))+dexp(-a11*(ymax - ymin)) !��̎���y�Ŕ��������ldy/dygs)
ysy=1.d0/yys																 !�t���������dygs/dy

do k=0,MZ
	z(k)=(dble(MZ-k)*zmin+dble(k)*zmax)/dble(MZ)
enddo


!!�h���̓���-------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------
open(31,file='kakuran3D_u.txt',status='old')
open(32,file='kakuran3D_v.txt',status='old')
open(33,file='kakuran3D_w.txt',status='old')
open(34,file='kakuran_uuuu.txt',status='unknown') 

do i=0,MXt
	do j=0,MY
		do k=0,MZ-1
 		read(31,*) kakuran_u(i,j,k)
 		read(32,*) kakuran_v(i,j,k)
 		read(33,*) kakuran_w(i,j,k)
 		write(34,*) kakuran_u(i,j,k) !�m�F�p
 		enddo
 	enddo
enddo

close(31)
close(32)
close(33)
close(34)

do i=0,MXt
	do j=0,MY
	kakuran_u(i,j,MZ)=kakuran_u(i,j,0)
	kakuran_v(i,j,MZ)=kakuran_v(i,j,0)
	kakuran_w(i,j,MZ)=kakuran_w(i,j,0)
	enddo
enddo

do j=0,MY
window(j,:)=dtanh(3.d0*y(j))*dexp(-2d-5*(y(j)**6.d0))
enddo

do i=0,MXt
kakuran_u(i,:,:)=kakuran_u(i,:,:)*window(:,:)
kakuran_v(i,:,:)=kakuran_v(i,:,:)*window(:,:)
kakuran_w(i,:,:)=kakuran_w(i,:,:)*window(:,:)
enddo

open(100,file='kakkuran_kakunin.txt',status='unknown')
!open(101,file='window.csv',status='unknown')
do k=0,MZ
	do j=0,MY
	write(100,'(3f24.16)') y(j),z(k),kakuran_u(2,j,k)
!	write(101,'(3f24.16)') y(j),z(k),window(j,k)
	enddo
	write(100,*)
!	write(101,*)
enddo
close(100)
!close(101)

!kakuran_u=0.d0
!kakuran_v=0.d0
!kakuran_w=0.d0

write(*,*) kakuran_u(2,10,5)



!!-----------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------


write(Chara, '(i4.4)') filecount      !!t=0�̎��̃�,u,P���o�͂��邽��
open(11,file='��'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(12,file='u'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(13,file='v'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(14,file='w'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(15,file='P'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(16,file='T'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(17,file='udif'//trim(Chara)//'.txt',status='unknown',form='formatted')

!!Q�̒�`-------------------------------------------------------------------------------------
open(21,file='r_2d.txt',status='old')
open(22,file='u_2d.txt',status='old')
open(23,file='v_2d.txt',status='old')
open(24,file='T_2d.txt',status='old')
open(25,file='rhooo.txt',status='unknown')    !!!w�͌��

do i=0,MX
	do j=0,MY
	read(21,*) rho_read(i,j)
	read(22,*) u_read(i,j)
	read(23,*) v_read(i,j)
	read(24,*) T_read(i,j)
	write(25,*) rho_read(i,j)
	enddo
enddo

close(21)
close(22)
close(23)
close(24)
close(25)

!do i=0,MZ
!if(mod(i,2)==0) then 
!dFy(1,0,0,i)=1.4d0*dFy(1,0,0,i)
!endif
!end do

P_read(:,:)=rho_read(:,:)*T_read(:,:)/gam/(Ma**2.d0)  !P���ρAT���狁�߂�

do k=0,MZ                                             !�ǂݍ��񂾃�,u,v,P,T��3�����ɂ���
rho(:,:,k)=rho_read(:,:)
u(:,:,k)=u_read(:,:)
v(:,:,k)=v_read(:,:)
P(:,:,k)=P_read(:,:)
T(:,:,k)=T_read(:,:)
enddo

do i=0,MZ
if(mod(i,2)==0) then 
v(30,0,i)=0.3d0
endif
end do


i=0
do k=0,MZ                                              !�����ő��֐��̖������ʂ���(0<Ly<5��10<Ly<20�́A�������Ȃ��B������)
		do j=0,MY
!		if(j <= 79 .or. j > 112) then 
!			Q(1,i,j,k)=rho_read(i,j)
!			Q(2,i,j,k)=rho_read(i,j)*(u_read(i,j))
!			Q(3,i,j,k)=rho_read(i,j)*(v_read(i,j))
!			Q(4,i,j,k)=0.d0
!			Q(5,i,j,k)=P_read(i,j)/(gam-1.d0)+0.5d0*rho_read(i,j)&
!			           &*(u_read(i,j)*u_read(i,j)+v_read(i,j)*v_read(i,j)+w_read(i,j)*w_read(i,j))		
!		else if(j > 79 .and. j <= 112) then
			Q(1,i,j,k)=rho_read(i,j)
			Q(2,i,j,k)=rho_read(i,j)*(u_read(i,j)+kakuran_u(s,j,k))
			Q(3,i,j,k)=rho_read(i,j)*(v(i,j,k)+kakuran_v(s,j,k))
			Q(4,i,j,k)=rho_read(i,j)*(kakuran_w(s,j,k))
			Q(5,i,j,k)=P_read(i,j)/(gam-1.d0)+0.5d0*rho_read(i,j)&
			           &*((u_read(i,j)+kakuran_u(s,j,k))**2.d0+(v(i,j,k)+kakuran_v(s,j,k))**2.d0+(w_read(i,j)+kakuran_w(s,j,k))**2.d0)	
!		end if
		enddo
enddo

do k=0,MZ
	do i=1,MX
		do j=0,MY
			Q(1,i,j,k)=rho_read(i,j)
			Q(2,i,j,k)=rho_read(i,j)*u_read(i,j)
			Q(3,i,j,k)=rho_read(i,j)*v(i,j,k)
			Q(4,i,j,k)=rho_read(i,j)*w_read(i,j)
			Q(5,i,j,k)=P_read(i,j)/(gam-1.d0)+0.5d0*rho_read(i,j)&
			           &*(u_read(i,j)*u_read(i,j)+v(i,j,k)*v_read(i,j)+w_read(i,j)*w_read(i,j))
		enddo
	enddo
enddo


call Q_rhoup(MX,MY,MZ,num,Q,rhoup,gam)

!!u��z�����̕��ς����A���ۂ̑��x�Ƃ̍������---------------------------
do k=0,MZ-1
sum_uz(:,:)=sum_uz(:,:)+rhoup(2,:,:,k)
enddo

u_ave=sum_uz/MZ

do k=0,MZ
u_margin(:,:,k)=rhoup(2,:,:,k)-u_ave(:,:)
enddo
!!------------------------------------------------------------------------
	do j=0,MY				!!t=0�̎��̃�,u,P���o�͂��邽��
		do i=0,MX		
		write(11,'(3f24.16)') x(i),y(j),rhoup(1,i,j,5)
		write(12,'(3f24.16)') x(i),y(j),rhoup(2,i,j,5)
		write(13,'(3f24.16)') x(i),y(j),rhoup(3,i,j,6)
		write(14,'(3f24.16)') x(i),y(j),rhoup(4,i,j,5)
		write(15,'(3f24.16)') x(i),y(j),rhoup(5,i,j,5)
		write(16,'(3f24.16)') x(i),y(j),T_read(i,j)
		write(17,'(3f24.16)') x(i),y(j),u_margin(i,j,5)
		enddo
	write(11,*)
	write(12,*)
	write(13,*)
	write(14,*)
	write(15,*)
	write(16,*)
	write(17,*)
	enddo
	
	close(11)
	close(12)
	close(13)
	close(14)
	close(15)
	close(16)
	close(17)

!!A��LU�����̃T�u���[�`����call---------------------------------------------------------------------------------------
call Ax(MX,AAMx,alpha,beta,alphaD,betaD,omega1,alpha3,alphaD3)
call LU(MX,AAMx,LUMx)

call Ax(MY,AAMy,alpha,beta,alphaD,betaD,omega1,alpha3,alphaD3)
call LU(MY,AAMy,LUMy)

!call Ax(MZ,AAMz,alpha,beta,alphaD,betaD,omega1,alpha3,alphaD3)
call Ax_shuuki(MZ,AAMz,alpha,beta,alphaD,betaD,omega1) 
call LU_shuuki(MZ,AAMz,LMz,UMz)

call Ax(MX,AAPx,alpha,beta,alphaD,betaD,omega2,alpha3,alphaD3)
call LU(MX,AAPx,LUPx)

call Ax(MY,AAPy,alpha,beta,alphaD,betaD,omega2,alpha3,alphaD3)
call LU(MY,AAPy,LUPy)

!call Ax(MZ,AAPz,alpha,beta,alphaD,betaD,omega2,alpha3,alphaD3)
call Ax_shuuki(MZ,AAPz,alpha,beta,alphaD,betaD,omega2)
call LU_shuuki(MZ,AAPz,LPz,UPz)

call Ax(MX,AANSx,alpha,beta,alphaD,betaD,omega0,alpha3,alphaD3) !NSCBCx���������p��LU�y�ёO�i���
call LU(MX,AANSx,LUNSx)

call Ax(MY,AANSy,alpha,beta,alphaD,betaD,omega0,alpha3,alphaD3) !NSCBCy���������p��LU�y�ёO�i���
call LU(MY,AANSy,LUNSy)

!call Ax(MZ,AANSz,alpha,beta,alphaD,betaD,omega0,alpha3,alphaD3) !NSCBCz���������p��LU�y�ёO�i���
call Ax_shuuki(MZ,AANSz,alpha,beta,alphaD,betaD,omega0)
call LU_shuuki(MZ,AANSz,LNSz,UNSz)

call Ax(MX,AATux,alpha,beta,alphaD,betaD,omega0,alpha3,alphaD3)    !x�S�����̂��߂�LU����(CCS)
call LU(MX,AATux,LUTux)

call Ax(MY,AATuy,alpha,beta,alphaD,betaD,omega0,alpha3,alphaD3)    !y�S�����̂��߂�LU����(CCS)
call LU(MY,AATuy,LUTuy)

!call Ax(MZ,AATuz,alpha,beta,alphaD,betaD,omega0,alpha3,alphaD3)    !z�S�����̂��߂�LU����(CCS)
call Ax_shuuki(MZ,AATuz,alpha,beta,alphaD,betaD,omega0)
call LU_shuuki(MZ,AATuz,LTuz,UTuz)

write(*,*) tmax

!!�������߂Ă���----------------------------------------------------------------------------
do times=1,tmax+1

if(mod(times,1000)==0) then   !����������0.8���ƂɃJ�E���g(s��MXt�܂ł����悤�ȃJ�E���g)�@�@�@�����̓���
 s=s+1
 if(s >= 99) s=1
endif


call kakuran_Q(MX,MY,MZ,MXt,num,s,gam,Ma,Q,kakuran_u,kakuran_v,kakuran_w,times,u,v,w,T)


call FPMx(zeta,MX,MY,MZ,num,gam,Q,Fx,FPx,FMx)           !x���� F,FP,FM�̒�`���̃R�[��
call FPMy(zeta,MX,MY,MZ,num,gam,Q,Fy,FPy,FMy)           !y���� F,FP,FM�̒�`���̃R�[��
call FPMz(zeta,MX,MY,MZ,num,gam,Q,Fz,FPz,FMz)           !z���� F,FP,FM�̒�`���̃R�[��

call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega1,xmax,xmin,QRMx,FMx,a3,ad3,LUMx,gMx) !FM x�����̑O�i���

call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega1,ymax,ymin,QRMy,FMy,a3,ad3,LUMy,gMy) !FM y�����̑O�i���

call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega1,zmax,zmin,QRMz,FMz,a3,ad3) !FM z�����̑O�i���
call forebackz(MX,MY,MZ,num,LMz,UMz,gMz,QRMz) 

call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega2,xmax,xmin,QRPx,FPx,a3,ad3,LUPx,gPx) !FP x�����̑O�i���

call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega2,ymax,ymin,QRPy,FPy,a3,ad3,LUPy,gPy) !FP y�����̑O�i���

call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega2,zmax,zmin,QRPz,FPz,a3,ad3) !FM z�����̑O�i���
call forebackz(MX,MY,MZ,num,LPz,UPz,gPz,QRPz)

gx(:,:,:,:)=gPx(:,:,:,:)+gMx(:,:,:,:)
gy(:,:,:,:)=gPy(:,:,:,:)+gMy(:,:,:,:)
gz(:,:,:,:)=gPz(:,:,:,:)+gMz(:,:,:,:)

do o=1,num
	do i=0,MX
		do k=0,MZ
		gy(o,i,:,k)=gy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo

							                        !!!NSCBC�֘A(x����)��call�J�n
call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,xmax,xmin,rhoupRHSx,rhoup,a3,ad3,LUNSx,grhoupx) !!drhoup(��,u,P��x��������)�����߂�
call NSCBC_Lx(MX,MY,MZ,num,NSLx,grhoupx,rhoup,omegaNS,gam,c0x,cMXx,Ma0x,MaMXx,Pmugen,Lx)
call NSCBC_dx(MX,MY,MZ,num,c0x,cMXx,NSLx,NSdx,rhoup)
call NSCBC_dFx(MX,MY,MZ,num,NSdx,rhoup,gx,gam)                  !!!x��������NSCBC�I��

							                             !!!NSCBC�֘A(y����)��call�J�n
call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,ymax,ymin,rhoupRHSy,rhoup,a3,ad3,LUNSy,grhoupy) !!drhoup(��,u,P��y��������)�����߂�
do o=1,num
	do i=0,MX
		do k=0,MZ
		grhoupy(o,i,:,k)=grhoupy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo
call NSCBC_Ly(MX,MY,MZ,num,NSLy,grhoupy,rhoup,omegaNS,gam,c0y,cMXy,Ma0y,MaMXy,Pmugen,Ly)
call NSCBC_dy(MX,MY,MZ,num,c0y,cMXy,NSLy,NSdy,rhoup)
call NSCBC_dFy(MX,MY,MZ,num,NSdy,rhoup,gy,gam)                  !!!y��������NSCBC�I��

							                             !!!NSCBC�֘A(z����)��call�J�n
!call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,rhoupRHSz,rhoup,a3,ad3)
!call forebackz(MX,MY,MZ,num,LNSz,UNSz,grhoupz,rhoupRHSz)
!call NSCBC_Lz(MX,MY,MZ,num,NSLz,grhoupz,rhoup,omegaNS,gam,c0z,cMXz,Ma0z,MaMXz,Pmugen,Lz)
!call NSCBC_dz(MX,MY,MZ,num,c0z,cMXz,NSLz,NSdz,rhoup)
!call NSCBC_dFz(MX,MY,MZ,num,NSdz,rhoup,gz,gam)                  !!!z��������NSCBC�I��

call Q_uT(MX,MY,MZ,num,gam,Ma,Q,Tuv)     !��������S���������߂� x����
call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,xmax,xmin,TuRHSx,Tuv,a3,ad3,LUTux,gTux)

call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,ymax,ymin,TuRHSy,Tuv,a3,ad3,LUTuy,gTuy)  !y�����S����
do o=1,num
	do i=0,MX
		do k=0,MZ
		gTuy(o,i,:,k)=gTuy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo

call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,zmax,zmin,TuRHSz,Tuv,a3,ad3)  !z�����S����
call forebackz(MX,MY,MZ,num,LTuz,UTuz,gTuz,TuRHSz)

call viscosityx(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvisx,gTux,gTuy,gTuz,mu)  !x���� �S������2�K����
call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,xmax,xmin,TuRHS2x,Tuvisx,a3,ad3,LUTux,ggTux)

call viscosityy(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvisy,gTux,gTuy,gTuz,mu)  !y���� �S������2�K����
call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,ymax,ymin,TuRHS2y,Tuvisy,a3,ad3,LUTuy,ggTuy)
do o=1,num
	do i=0,MX
		do k=0,MZ
		ggTuy(o,i,:,k)=ggTuy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo

call viscosityz(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvisz,gTux,gTuy,gTuz,mu)  !z���� �S������2�K����
call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,zmax,zmin,TuRHS2z,Tuvisz,a3,ad3)
call forebackz(MX,MY,MZ,num,LTuz,UTuz,ggTuz,TuRHS2z)
call visout(MX,MY,MZ,num,ggTux,ggTuy,grhoupx,grhoupy,Tuvisx,Tuvisy,rhoup)


			dQ1(:,:,:,:)=Q(:,:,:,:)-1.d0*dt*(gx(:,:,:,:)+gy(:,:,:,:)+gz(:,:,:,:)&
			             &-(ggTux(:,:,:,:)+ggTuy(:,:,:,:)+ggTuz(:,:,:,:)))	

			!dQ1(:,:,:,:)=Q(:,:,:,:)-1.d0*dt*(gx(:,:,:,:)+gy(:,:,:,:)&
			!             &-(ggTux(:,:,:,:)+ggTuy(:,:,:,:)))	



call kakuran_Q(MX,MY,MZ,MXt,num,s,gam,Ma,dQ1,kakuran_u,kakuran_v,kakuran_w,times,u,v,w,T)


dQ1(2,:,0,:)=dQ1(1,:,0,:)*u(:,0,:)
dQ1(3,:,0,:)=dQ1(1,:,0,:)*v(:,0,:)    !!�f�B��������(�����ǖ�)
dQ1(4,:,0,:)=dQ1(1,:,0,:)*w(:,0,:) 
dQ1(5,:,0,:)=(dQ1(1,:,0,:)*T(:,0,:)/gam/(Ma**2.d0))/(gam-1.d0)+0.5d0*dQ1(1,:,0,:)&
             &*(u(:,0,:)*u(:,0,:)+v(:,0,:)*v(:,0,:)+w(:,0,:)*w(:,0,:))

!call NBound(MX,MY,Q)
call Q_rhoup(MX,MY,MZ,num,dQ1,rhoup,gam)


call FPMx(zeta,MX,MY,MZ,num,gam,dQ1,Fx,FPx,FMx)           !x���� F,FP,FM�̒�`���̃R�[��
call FPMy(zeta,MX,MY,MZ,num,gam,dQ1,Fy,FPy,FMy)           !y���� F,FP,FM�̒�`���̃R�[��
call FPMz(zeta,MX,MY,MZ,num,gam,dQ1,Fz,FPz,FMz)           !z���� F,FP,FM�̒�`���̃R�[��

call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega1,xmax,xmin,QRMx,FMx,a3,ad3,LUMx,gMx) !FM x�����̑O�i���

call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega1,ymax,ymin,QRMy,FMy,a3,ad3,LUMy,gMy) !FM y�����̑O�i���

call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega1,zmax,zmin,QRMz,FMz,a3,ad3) !FM z�����̑O�i���
call forebackz(MX,MY,MZ,num,LMz,UMz,gMz,QRMz) 

call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega2,xmax,xmin,QRPx,FPx,a3,ad3,LUPx,gPx) !FP x�����̑O�i���

call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega2,ymax,ymin,QRPy,FPy,a3,ad3,LUPy,gPy) !FP y�����̑O�i���

call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega2,zmax,zmin,QRPz,FPz,a3,ad3) !FM z�����̑O�i���
call forebackz(MX,MY,MZ,num,LPz,UPz,gPz,QRPz)

gx(:,:,:,:)=gPx(:,:,:,:)+gMx(:,:,:,:)
gy(:,:,:,:)=gPy(:,:,:,:)+gMy(:,:,:,:)
gz(:,:,:,:)=gPz(:,:,:,:)+gMz(:,:,:,:)

do o=1,num
	do i=0,MX
		do k=0,MZ
		gy(o,i,:,k)=gy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo

							                        !!!NSCBC�֘A(x����)��call�J�n
call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,xmax,xmin,rhoupRHSx,rhoup,a3,ad3,LUNSx,grhoupx) !!drhoup(��,u,P��x��������)�����߂�
call NSCBC_Lx(MX,MY,MZ,num,NSLx,grhoupx,rhoup,omegaNS,gam,c0x,cMXx,Ma0x,MaMXx,Pmugen,Lx)
call NSCBC_dx(MX,MY,MZ,num,c0x,cMXx,NSLx,NSdx,rhoup)
call NSCBC_dFx(MX,MY,MZ,num,NSdx,rhoup,gx,gam)                  !!!x��������NSCBC�I��

							                             !!!NSCBC�֘A(y����)��call�J�n
call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,ymax,ymin,rhoupRHSy,rhoup,a3,ad3,LUNSy,grhoupy) !!drhoup(��,u,P��y��������)�����߂�
do o=1,num
	do i=0,MX
		do k=0,MZ
		grhoupy(o,i,:,k)=grhoupy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo
call NSCBC_Ly(MX,MY,MZ,num,NSLy,grhoupy,rhoup,omegaNS,gam,c0y,cMXy,Ma0y,MaMXy,Pmugen,Ly)
call NSCBC_dy(MX,MY,MZ,num,c0y,cMXy,NSLy,NSdy,rhoup)
call NSCBC_dFy(MX,MY,MZ,num,NSdy,rhoup,gy,gam)                  !!!y��������NSCBC�I��

							                             !!!NSCBC�֘A(z����)��call�J�n
!call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,rhoupRHSz,rhoup,a3,ad3)
!call forebackz(MX,MY,MZ,num,LNSz,UNSz,grhoupz,rhoupRHSz)
!call NSCBC_Lz(MX,MY,MZ,num,NSLz,grhoupz,rhoup,omegaNS,gam,c0z,cMXz,Ma0z,MaMXz,Pmugen,Lz)
!call NSCBC_dz(MX,MY,MZ,num,c0z,cMXz,NSLz,NSdz,rhoup)
!call NSCBC_dFz(MX,MY,MZ,num,NSdz,rhoup,gz,gam)                  !!!z��������NSCBC�I��

call Q_uT(MX,MY,MZ,num,gam,Ma,dQ1,Tuv)     !��������S���������߂� x����
call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,xmax,xmin,TuRHSx,Tuv,a3,ad3,LUTux,gTux)

call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,ymax,ymin,TuRHSy,Tuv,a3,ad3,LUTuy,gTuy)  !y�����S����
do o=1,num
	do i=0,MX
		do k=0,MZ
		gTuy(o,i,:,k)=gTuy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo

call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,zmax,zmin,TuRHSz,Tuv,a3,ad3)  !z�����S����
call forebackz(MX,MY,MZ,num,LTuz,UTuz,gTuz,TuRHSz)

call viscosityx(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvisx,gTux,gTuy,gTuz,mu)  !x���� �S������2�K����
call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,xmax,xmin,TuRHS2x,Tuvisx,a3,ad3,LUTux,ggTux)

call viscosityy(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvisy,gTux,gTuy,gTuz,mu)  !y���� �S������2�K����
call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,ymax,ymin,TuRHS2y,Tuvisy,a3,ad3,LUTuy,ggTuy)
do o=1,num
	do i=0,MX
		do k=0,MZ
		ggTuy(o,i,:,k)=ggTuy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo

call viscosityz(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvisz,gTux,gTuy,gTuz,mu)  !z���� �S������2�K����
call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,zmax,zmin,TuRHS2z,Tuvisz,a3,ad3)
call forebackz(MX,MY,MZ,num,LTuz,UTuz,ggTuz,TuRHS2z)
call visout(MX,MY,MZ,num,ggTux,ggTuy,grhoupx,grhoupy,Tuvisx,Tuvisy,rhoup)

!			Q=Q-1.d0*dt*(gx+gy+gz)
			dQ2(:,:,:,:)=0.75d0*Q(:,:,:,:)+0.25d0*dQ1(:,:,:,:)-0.25d0*dt*(gx(:,:,:,:)+gy(:,:,:,:)+gz(:,:,:,:)&
			             &-(ggTux(:,:,:,:)+ggTuy(:,:,:,:)+ggTuz(:,:,:,:)))	

			!dQ2(:,:,:,:)=0.75d0*Q(:,:,:,:)+0.25d0*dQ1(:,:,:,:)-0.25d0*dt*(gx(:,:,:,:)+gy(:,:,:,:)&
			!             &-(ggTux(:,:,:,:)+ggTuy(:,:,:,:)))	


call kakuran_Q(MX,MY,MZ,MXt,num,s,gam,Ma,dQ2,kakuran_u,kakuran_v,kakuran_w,times,u,v,w,T)


dQ2(2,:,0,:)=dQ2(1,:,0,:)*u(:,0,:)
dQ2(3,:,0,:)=dQ2(1,:,0,:)*v(:,0,:)    !!�f�B��������(�����ǖ�)
dQ2(4,:,0,:)=dQ2(1,:,0,:)*w(:,0,:) 
dQ2(5,:,0,:)=(dQ2(1,:,0,:)*T(:,0,:)/gam/(Ma**2.d0))/(gam-1.d0)+0.5d0*dQ2(1,:,0,:)&
             &*(u(:,0,:)*u(:,0,:)+v(:,0,:)*v(:,0,:)+w(:,0,:)*w(:,0,:))

!call NBound(MX,MY,Q)
call Q_rhoup(MX,MY,MZ,num,dQ2,rhoup,gam)

			
call FPMx(zeta,MX,MY,MZ,num,gam,dQ2,Fx,FPx,FMx)           !x���� F,FP,FM�̒�`���̃R�[��
call FPMy(zeta,MX,MY,MZ,num,gam,dQ2,Fy,FPy,FMy)           !y���� F,FP,FM�̒�`���̃R�[��
call FPMz(zeta,MX,MY,MZ,num,gam,dQ2,Fz,FPz,FMz)           !z���� F,FP,FM�̒�`���̃R�[��

call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega1,xmax,xmin,QRMx,FMx,a3,ad3,LUMx,gMx) !FM x�����̑O�i���

call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega1,ymax,ymin,QRMy,FMy,a3,ad3,LUMy,gMy) !FM y�����̑O�i���

call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega1,zmax,zmin,QRMz,FMz,a3,ad3) !FM z�����̑O�i���
call forebackz(MX,MY,MZ,num,LMz,UMz,gMz,QRMz) 

call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega2,xmax,xmin,QRPx,FPx,a3,ad3,LUPx,gPx) !FP x�����̑O�i���

call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega2,ymax,ymin,QRPy,FPy,a3,ad3,LUPy,gPy) !FP y�����̑O�i���

call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega2,zmax,zmin,QRPz,FPz,a3,ad3) !FM z�����̑O�i���
call forebackz(MX,MY,MZ,num,LPz,UPz,gPz,QRPz)

gx(:,:,:,:)=gPx(:,:,:,:)+gMx(:,:,:,:)
gy(:,:,:,:)=gPy(:,:,:,:)+gMy(:,:,:,:)
gz(:,:,:,:)=gPz(:,:,:,:)+gMz(:,:,:,:)

do o=1,num
	do i=0,MX
		do k=0,MZ
		gy(o,i,:,k)=gy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo

							                        !!!NSCBC�֘A(x����)��call�J�n
call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,xmax,xmin,rhoupRHSx,rhoup,a3,ad3,LUNSx,grhoupx) !!drhoup(��,u,P��x��������)�����߂�
call NSCBC_Lx(MX,MY,MZ,num,NSLx,grhoupx,rhoup,omegaNS,gam,c0x,cMXx,Ma0x,MaMXx,Pmugen,Lx)
call NSCBC_dx(MX,MY,MZ,num,c0x,cMXx,NSLx,NSdx,rhoup)
call NSCBC_dFx(MX,MY,MZ,num,NSdx,rhoup,gx,gam)                  !!!x��������NSCBC�I��

							                             !!!NSCBC�֘A(y����)��call�J�n
call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,ymax,ymin,rhoupRHSy,rhoup,a3,ad3,LUNSy,grhoupy) !!drhoup(��,u,P��y��������)�����߂�
do o=1,num
	do i=0,MX
		do k=0,MZ
		grhoupy(o,i,:,k)=grhoupy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo
call NSCBC_Ly(MX,MY,MZ,num,NSLy,grhoupy,rhoup,omegaNS,gam,c0y,cMXy,Ma0y,MaMXy,Pmugen,Ly)
call NSCBC_dy(MX,MY,MZ,num,c0y,cMXy,NSLy,NSdy,rhoup)
call NSCBC_dFy(MX,MY,MZ,num,NSdy,rhoup,gy,gam)                  !!!y��������NSCBC�I��

							                             !!!NSCBC�֘A(z����)��call�J�n
!call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,rhoupRHSz,rhoup,a3,ad3)
!call forebackz(MX,MY,MZ,num,LNSz,UNSz,grhoupz,rhoupRHSz)
!call NSCBC_Lz(MX,MY,MZ,num,NSLz,grhoupz,rhoup,omegaNS,gam,c0z,cMXz,Ma0z,MaMXz,Pmugen,Lz)
!call NSCBC_dz(MX,MY,MZ,num,c0z,cMXz,NSLz,NSdz,rhoup)
!call NSCBC_dFz(MX,MY,MZ,num,NSdz,rhoup,gz,gam)                  !!!z��������NSCBC�I��

call Q_uT(MX,MY,MZ,num,gam,Ma,dQ2,Tuv)     !��������S���������߂� x����
call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,xmax,xmin,TuRHSx,Tuv,a3,ad3,LUTux,gTux)

call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,ymax,ymin,TuRHSy,Tuv,a3,ad3,LUTuy,gTuy)  !y�����S����
do o=1,num
	do i=0,MX
		do k=0,MZ
		gTuy(o,i,:,k)=gTuy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo

call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,zmax,zmin,TuRHSz,Tuv,a3,ad3)  !z�����S����
call forebackz(MX,MY,MZ,num,LTuz,UTuz,gTuz,TuRHSz)

call viscosityx(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvisx,gTux,gTuy,gTuz,mu)  !x���� �S������2�K����
call RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,xmax,xmin,TuRHS2x,Tuvisx,a3,ad3,LUTux,ggTux)

call viscosityy(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvisy,gTux,gTuy,gTuz,mu)  !y���� �S������2�K����
call RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,ymax,ymin,TuRHS2y,Tuvisy,a3,ad3,LUTuy,ggTuy)
do o=1,num
	do i=0,MX
		do k=0,MZ
		ggTuy(o,i,:,k)=ggTuy(o,i,:,k)*ysy(:)    !�i�q�L����K�p
		enddo
	enddo
enddo

call viscosityz(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvisz,gTux,gTuy,gTuz,mu)  !z���� �S������2�K����
call RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega0,zmax,zmin,TuRHS2z,Tuvisz,a3,ad3)
call forebackz(MX,MY,MZ,num,LTuz,UTuz,ggTuz,TuRHS2z)
call visout(MX,MY,MZ,num,ggTux,ggTuy,grhoupx,grhoupy,Tuvisx,Tuvisy,rhoup)

			Q(:,:,:,:)=(1d0/3d0)*Q(:,:,:,:)+(2d0/3d0)*dQ2(:,:,:,:)-(2d0/3d0)*dt*(gx(:,:,:,:)+gy(:,:,:,:)+gz(:,:,:,:)&
			           &-(ggTux(:,:,:,:)+ggTuy(:,:,:,:)+ggTuz(:,:,:,:)))	


call kakuran_Q(MX,MY,MZ,MXt,num,s,gam,Ma,Q,kakuran_u,kakuran_v,kakuran_w,times,u,v,w,T)


Q(2,:,0,:)=Q(1,:,0,:)*u(:,0,:)
Q(3,:,0,:)=Q(1,:,0,:)*v(:,0,:)    !!�f�B��������(�����ǖ�)
Q(4,:,0,:)=Q(1,:,0,:)*w(:,0,:) 
Q(5,:,0,:)=(Q(1,:,0,:)*T(:,0,:)/gam/(Ma**2.d0))/(gam-1.d0)+0.5d0*Q(1,:,0,:)&
             &*(u(:,0,:)*u(:,0,:)+v(:,0,:)*v(:,0,:)+w(:,0,:)*w(:,0,:))
 

!call NBound(MX,MY,Q)
call Q_rhoup(MX,MY,MZ,num,Q,rhoup,gam)
call Q_uT(MX,MY,MZ,num,gam,Ma,Q,Tuv)
			
!call NBound(MX,MY,Q)
!do k=0,MY									!!isnan ��ق�
!	do i=0,MX
!		do o=1,3
!		if(isnan(Q(o,i,k))) then
!		write(*,*) 't=', dt*(dble(j))
!		endif
!		enddo
!	enddo
!enddo


if(isnan(rhoup(1,3,My,1))) stop '�G���['
		

if(mod(times,1000)==0) then
filecount=filecount+1
write(Chara, '(i4.4)') filecount      
open(11,file='��'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(12,file='u'//trim(Chara)//'.txt',status='unknown',form='formatted')
!open(13,file='v'//trim(Chara)//'.txt',status='unknown',form='formatted')
!open(14,file='w'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(15,file='P'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(16,file='T'//trim(Chara)//'.txt',status='unknown',form='formatted')
open(17,file='udif'//trim(Chara)//'.txt',status='unknown',form='formatted')

!!u��z�����̕��ς����A���ۂ̑��x�Ƃ̍������---------------------------
sum_uz=0.d0

do k=0,MZ-1
sum_uz(:,:)=sum_uz(:,:)+rhoup(2,:,:,k)
enddo

u_ave=sum_uz/MZ

do k=0,MZ
u_margin(:,:,k)=rhoup(2,:,:,k)-u_ave(:,:)
enddo
!!------------------------------------------------------------------------

	do j=0,MY				!!t=0�̎��̃�,u,P���o�͂��邽��
		do i=0,MX		
		write(11,'(3f24.16)') x(i),y(j),rhoup(1,i,j,5)
		write(12,'(3f24.16)') x(i),y(j),rhoup(2,i,j,5)
!		write(13,'(3f24.16)') x(i),y(j),rhoup(3,i,j,5)
!		write(14,'(3f24.16)') x(i),y(j),rhoup(4,i,j,5)
		write(15,'(3f24.16)') x(i),y(j),rhoup(5,i,j,5)
		write(16,'(3f24.16)') x(i),y(j),Tuv(5,i,j,5)
		write(17,'(3f24.16)') x(i),y(j),u_margin(i,j,5)
		enddo
	write(11,*)
	write(12,*)
!	write(13,*)
!	write(14,*)
	write(15,*)
	write(16,*)
	write(17,*)
	enddo
	
	close(11)
	close(12)
!	close(13)
!	close(14)
	close(15)
	close(16)
	close(17)

	call cpu_time(CT2)
endif

if(mod(times,1000)==0) then
open(501,file='3D_��.txt',status='unknown')
open(502,file='3D_u.txt',status='unknown')
open(503,file='3D_v.txt',status='unknown')
open(504,file='3D_w.txt',status='unknown')
open(505,file='3D_P.txt',status='unknown')
open(506,file='3D_T.txt',status='unknown')

do k=0,MZ
	do j=0,MY
		do i=0,MX
		write(501,*) rhoup(1,i,j,k)
		write(502,*) rhoup(2,i,j,k)
		write(503,*) rhoup(3,i,j,k)
		write(504,*) rhoup(4,i,j,k)
		write(505,*) rhoup(5,i,j,k)
		write(506,*) Tuv(5,i,j,k)
		enddo
	enddo
enddo
close(501)
close(502)
close(503)
close(504)
close(505)
close(506)
endif		



!if(times >= 20000 .and. mod(times,200)==0) then
!call favre_ave(MX,MY,MZ,num,x,y,rhoup,Tuv,frhoup,count)
!
!call uzu(MX,MY,MZ,num,x,y,grhoupx,grhoupy,grhoupz,uzu_W,uzuW_ave,count2)
!
!call cf_tauw(MX,MY,MZ,x,Re,a,b,c,ad,bd,cd,omega,alpha,beta,alphaD,betaD,alpha3,alphaD3,ymax,ymin,ysy&
!			&,rho_read,u_read,frhoup,tauw,cf,count3)
!end if

		
		write(*,*) times
enddo

end program DCS

!!!�T�u���[�`���E�Ӓ�` x���� LU�O�i���------------------------------------------------------------------------
subroutine RHSx(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega,xmax,xmin,QR,F,a3,ad3,LUcom,dF)
implicit none
integer :: MX,MY,MZ,num
double precision :: ad,bd,cd,a,b,c,omega,a3,ad3
double precision :: f(1:num,0:MX,0:MY,0:MZ)
double precision :: QR(1:num,0:MX,0:MY,0:MZ)
double precision :: LUcom(-2:2,0:MX)
double precision :: Y(1:num,0:MX,0:MY,0:MZ)
double precision :: dF(1:num,0:MX,0:MY,0:MZ)
double precision :: xmax,xmin
double precision h
integer i

h=(xmax-xmin)/dble(MX)

	QR(:,0,:,:)=((-17.d0/6.d0)*f(:,0,:,:)+1.5d0*f(:,1,:,:)+1.5d0*f(:,2,:,:)-(1.0d0/6.d0)*f(:,3,:,:))/h
	QR(:,MX,:,:)=((1.d0/6.d0)*f(:,MX-3,:,:)-1.5d0*f(:,MX-2,:,:)-1.5d0*f(:,MX-1,:,:)+(17.d0/6.d0)*f(:,MX,:,:))/h


	QR(:,1,:,:)=a3*(-f(:,0,:,:)+f(:,2,:,:))/(2.d0*h)+omega*(ad3*(f(:,0,:,:)-2.d0*f(:,1,:,:)+f(:,2,:,:))/h)
	QR(:,MX-1,:,:)=a3*(-f(:,MX-2,:,:)+f(:,MX,:,:))/(2.d0*h)+omega*(ad3*(f(:,MX-2,:,:)-2.d0*f(:,MX-1,:,:)+f(:,MX,:,:))/h)

		do i=2,MX-2	
		QR(:,i,:,:)=a*(-f(:,i-1,:,:)+f(:,i+1,:,:))/(2.d0*h)+b*(-f(:,i-2,:,:)+f(:,i+2,:,:))/(4.d0*h)&
	     &+omega*(ad*(f(:,i-1,:,:)-2.d0*f(:,i,:,:)+f(:,i+1,:,:))/h+bd*(f(:,i-2,:,:)-2.d0*f(:,i,:,:)+f(:,i+2,:,:))/(4.d0*h))
		enddo 

!!  �O�i���
	Y(:,0,:,:)=QR(:,0,:,:)
	do i=1,MX
		Y(:,i,:,:)=QR(:,i,:,:)-LUcom(-1,i)*Y(:,i-1,:,:)
	enddo
!!  ��ޑ��
	dF(:,MX,:,:)=Y(:,MX,:,:)/LUcom(0,MX)
	dF(:,MX-1,:,:)=(Y(:,MX-1,:,:)-dF(:,MX,:,:)*LUcom(1,MX-1))/LUcom(0,MX-1)
	do i=MX-2,0,-1
		dF(:,i,:,:)=(Y(:,i,:,:)-LUcom(1,i)*dF(:,i+1,:,:)-LUcom(2,i)*dF(:,i+2,:,:))/LUcom(0,i)
	enddo

end subroutine RHSx

!!!�T�u���[�`���E�Ӓ�` y���� LU�̑O�i���------------------------------------------------------------------------
subroutine RHSy(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega,ymax,ymin,QR,F,a3,ad3,LUcom,dF)
integer :: MX,MY,MZ,num
double precision :: ad,bd,cd,a,b,c,omega,a3,ad3
double precision :: f(1:num,0:MX,0:MY,0:MZ)
double precision :: QR(1:num,0:MX,0:MY,0:MZ)
double precision :: LUcom(-2:2,0:MY)
double precision :: Y(1:num,0:MX,0:MY,0:MZ)
double precision :: dF(1:num,0:MX,0:MY,0:MZ)
double precision :: ymax,ymin
double precision h
integer j,i

h=(ymax-ymin)/dble(MY)

	QR(:,:,0,:)=((-17.d0/6.d0)*f(:,:,0,:)+1.5d0*f(:,:,1,:)+1.5d0*f(:,:,2,:)-(1.0d0/6.d0)*f(:,:,3,:))/h
	QR(:,:,MY,:)=((1.d0/6.d0)*f(:,:,MY-3,:)-1.5d0*f(:,:,MY-2,:)-1.5d0*f(:,:,MY-1,:)+(17.d0/6.d0)*f(:,:,MY,:))/h


	QR(:,:,1,:)=a3*(-f(:,:,0,:)+f(:,:,2,:))/(2.d0*h)+omega*(ad3*(f(:,:,0,:)-2.d0*f(:,:,1,:)+f(:,:,2,:))/h)
	QR(:,:,MY-1,:)=a3*(-f(:,:,MY-2,:)+f(:,:,MY,:))/(2.d0*h)+omega*(ad3*(f(:,:,MY-2,:)-2.d0*f(:,:,MY-1,:)+f(:,:,MY,:))/h)

		do i=2,MY-2	
		QR(:,:,i,:)=a*(-f(:,:,i-1,:)+f(:,:,i+1,:))/(2.d0*h)+b*(-f(:,:,i-2,:)+f(:,:,i+2,:))/(4.d0*h)&
	     &+omega*(ad*(f(:,:,i-1,:)-2.d0*f(:,:,i,:)+f(:,:,i+1,:))/h+bd*(f(:,:,i-2,:)-2.d0*f(:,:,i,:)+f(:,:,i+2,:))/(4.d0*h))
		enddo 

!!  �O�i���
	Y(:,:,0,:)=QR(:,:,0,:)
	do j=1,MY
		Y(:,:,j,:)=QR(:,:,j,:)-LUcom(-1,j) * Y(:,:,j-1,:)
	enddo
!!  ��ޑ��
	dF(:,:,MY,:)=Y(:,:,MY,:)/LUcom(0,MY)
	dF(:,:,MY-1,:)=(Y(:,:,MY-1,:)-dF(:,:,MY,:)*LUcom(1,MY-1))/LUcom(0,MY-1)
	do j=MY-2,0,-1
		dF(:,:,j,:)=(Y(:,:,j,:)-LUcom(1,j)*dF(:,:,j+1,:)-LUcom(2,j)*dF(:,:,j+2,:))/LUcom(0,j)
	enddo

end subroutine RHSy

!!!�T�u���[�`���E�Ӓ�` z����  LU�O�i���------------------------------------------------------------------------
subroutine RHSz(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega,zmax,zmin,QR,F,a3,ad3,LUcom,dF)
integer :: MX,MY,MZ,num
double precision :: ad,bd,cd,a,b,c,omega,a3,ad3
double precision :: f(1:num,0:MX,0:MY,0:MZ)
double precision :: QR(1:num,0:MX,0:MY,0:MZ)
double precision :: LUcom(-2:2,0:MZ)
double precision :: Y(1:num,0:MX,0:MY,0:MZ)
double precision :: dF(1:num,0:MX,0:MY,0:MZ)
double precision :: zmax,zmin
double precision h
integer k,i

h=(zmax-zmin)/dble(MZ)

	QR(:,:,:,0)=((-17.d0/6.d0)*f(:,:,:,0)+1.5d0*f(:,:,:,1)+1.5d0*f(:,:,:,2)-(1.0d0/6.d0)*f(:,:,:,3))/h
	QR(:,:,:,MZ)=((1.d0/6.d0)*f(:,:,:,MZ-3)-1.5d0*f(:,:,:,MZ-2)-1.5d0*f(:,:,:,MZ-1)+(17.d0/6.d0)*f(:,:,:,MZ))/h


	QR(:,:,:,1)=a3*(-f(:,:,:,0)+f(:,:,:,2))/(2.d0*h)+omega*(ad3*(f(:,:,:,0)-2.d0*f(:,:,:,1)+f(:,:,:,2))/h)
	QR(:,:,:,MZ-1)=a3*(-f(:,:,:,MZ-2)+f(:,:,:,MZ))/(2.d0*h)+omega*(ad3*(f(:,:,:,MZ-2)-2.d0*f(:,:,:,MZ-1)+f(:,:,:,MZ))/h)

		do i=2,MZ-2	
		QR(:,:,:,i)=a*(-f(:,:,:,i-1)+f(:,:,:,i+1))/(2.d0*h)+b*(-f(:,:,:,i-2)+f(:,:,:,i+2))/(4.d0*h)&
	     &+omega*(ad*(f(:,:,:,i-1)-2.d0*f(:,:,:,i)+f(:,:,:,i+1))/h+bd*(f(:,:,:,i-2)-2.d0*f(:,:,:,i)+f(:,:,:,i+2))/(4.d0*h))
		enddo 
		
		
!!  �O�i���
	Y(:,:,:,0)=QR(:,:,:,0)
	do k=1,MZ
		Y(:,:,:,k)=QR(:,:,:,k)-LUcom(-1,k) * Y(:,:,:,k-1)
	enddo
!!  ��ޑ��
	dF(:,:,:,MZ)=Y(:,:,:,MZ)/LUcom(0,MZ)
	dF(:,:,:,MZ-1)=(Y(:,:,:,MZ-1)-dF(:,:,:,MZ)*LUcom(1,MZ-1))/LUcom(0,MZ-1)
	do k=MZ-2,0,-1
		dF(:,:,:,k)=(Y(:,:,:,k)-LUcom(1,k)*dF(:,:,:,k+1)-LUcom(2,k)*dF(:,:,:,k+2))/LUcom(0,k)
	enddo

end subroutine RHSz

!!!�T�u���[�`���E�Ӓ�` z���� ����  LU�O�i���------------------------------------------------------------------------
subroutine RHSz_shuuki(MX,MY,MZ,num,a,b,c,ad,bd,cd,omega,zmax,zmin,QR,F,a3,ad3)
integer :: MX,MY,MZ,num
double precision :: ad,bd,cd,a,b,c,omega,a3,ad3
double precision :: f(1:num,0:MX,0:MY,0:MZ)
double precision :: QR(1:num,0:MX,0:MY,0:MZ)
double precision :: LUcom(-2:2,0:MZ)
double precision :: Y(1:num,0:MX,0:MY,0:MZ)
double precision :: dF(1:num,0:MX,0:MY,0:MZ)
double precision :: zmax,zmin
double precision h
integer im1,im2,im3,ip1,ip2,ip3,ip0
integer k,i

im1=0 ;im2=0 ;im3=0 ;ip1=0 ;ip2=0 ;ip3=0 ;ip0=0

h=(zmax-zmin)/dble(MZ)

		do i=0,MZ-1
		im1=mod((MZ-1)+i,MZ)       !i=0,N-1
		im2=mod((MZ-2)+i,MZ)       !i=0,N-2     
		im3=mod(MZ+i-3,MZ)         !i=0,N-3   
		ip1=mod(MZ+i+1,MZ)         !i=0,1 
		ip2=mod(MZ+i+2,MZ)         !i=0,2
		ip3=mod(MZ+i+3,MZ)         !i=0,3
		
		QR(:,:,:,i)=a*(-f(:,:,:,im1)+f(:,:,:,ip1))/(2.d0*h)+b*(-f(:,:,:,im2)+f(:,:,:,ip2))/(4.d0*h)&
	     &+omega*(ad*(f(:,:,:,im1)-2.d0*f(:,:,:,i)+f(:,:,:,ip1))/h+bd*(f(:,:,:,im2)-2.d0*f(:,:,:,i)+f(:,:,:,ip2))/(4.d0*h))
		enddo 
		
!!!  �O�i���
!	Y(:,:,:,0)=QR(:,:,:,0)
!	do k=1,MZ-1
!		Y(:,:,:,k)=QR(:,:,:,k)-LUcom(-1,k) * Y(:,:,:,k-1)
!	enddo
!!!  ��ޑ��
!	dF(:,:,:,MZ-1)=Y(:,:,:,MZ-1)/LUcom(0,MZ-1)
!	dF(:,:,:,MZ-2)=(Y(:,:,:,MZ-2)-dF(:,:,:,MZ-1)*LUcom(1,MZ-2))/LUcom(0,MZ-2)
!	do k=MZ-3,0,-1
!		dF(:,:,:,k)=(Y(:,:,:,k)-LUcom(1,k)*dF(:,:,:,k+1)-LUcom(2,k)*dF(:,:,:,k+2))/LUcom(0,k)
!	enddo

!dF(:,:,:,MZ)=dF(:,:,:,0)

end subroutine RHSz_shuuki


!!!�T�u���[�`�����Ӓ�`(A�̒�`)---------------------------------------------------------------
subroutine Ax(MX,AA,alpha,beta,alphaD,betaD,omega,alpha3,alphaD3)
integer :: MX
double precision :: alpha,beta,alphaD,betaD,omega,alpha3,alphaD3
double precision :: AA(0:MX,0:MX)
integer i

AA=0.0d0

do i=0,MX
	AA(i,i)=1.0d0
enddo

	AA(0,1)=3.0d0                !�Б�����
	AA(MX,MX-1)=3.0d0
	
	AA(1,0)=alpha3*(1.d0-alphaD3*omega)  !3��DCS
	AA(1,2)=alpha3*(1.d0+alphaD3*omega)
	AA(MX-1,MX)=alpha3*(1.d0+alphaD3*omega)
	AA(MX-1,MX-2)=alpha3*(1.d0-alphaD3*omega)

do i=2,MX-2								!5��DCS
	AA(i,i+1)=alpha*(1.d0+omega*alphaD)
	AA(i,i+2)=beta*(1.d0+omega*betaD)
	AA(i,i-1)=alpha*(1.d0-omega*alphaD)
	AA(i,i-2)=beta*(1.d0-omega*betaD)
enddo

end subroutine Ax

!!!�T�u���[�`�����Ӓ�` ����(CCS6�p) (A�̒�`)---------------------------------------------------------------
subroutine Ax_shuuki(MX,AA,alpha,beta,alphaD,betaD,omega)
integer, intent(in) :: MX
double precision, intent(in) :: alpha,beta,alphaD,betaD,omega
double precision, intent(out) :: AA(0:MX,0:MX)
integer i

AA=0.0d0
do i=0,MX-1
	AA(i,i)=1.0d0
enddo

do i=0,MX-2
	AA(i,i+1)=alpha*(1.d0+omega*alphaD)
enddo

do i=0,MX-3
	AA(i,i+2)=beta*(1.d0+omega*betaD)
enddo
do i=1,MX-1
	AA(i,i-1)=alpha*(1.d0-omega*alphaD)
enddo
do i=2,MX-1
	AA(i,i-2)=beta*(1.d0-omega*betaD)
enddo

AA(0,MX-1)=alpha*(1.d0-omega*alphaD)
AA(MX-1,0)=alpha*(1.d0+omega*alphaD)
AA(0,MX-2)=beta*(1.d0-omega*betaD)
AA(MX-2,0)=beta*(1.d0+omega*betaD)
AA(1,MX-1)=beta*(1.d0-omega*betaD)
AA(MX-1,1)=beta*(1.d0+omega*betaD)

end subroutine Ax_shuuki



!!!�T�u���[�`��LU���� x��y��z�������L �����āALU���k---------------------------------------------------------------------------
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
do i=0,MX              !���J�n
L(i,i)=1.0d0
enddo
	do i=0,MX		
		do j=0,MX
		s1=0.0d0
		s2=0.0d0
	   	if(i<=j) then        !��O�p�s������(U)
	   		do k=0,i-1 		
	   		s1=s1+L(i,k)*U(k,j)
	   		enddo
	   	U(i,j)=AA(i,j)-s1
	   	else if(i>=j) then   !���O�p�s������(L)
	   		do k=0,j-1
	   		s2=s2+L(i,k)*U(k,j)
	   		enddo
	   	L(i,j)=(AA(i,j)-s2)/U(j,j)
	   	endif
	   	enddo
	enddo
	
	do i = 0,MX                   !0�𔲂��Ĉ��k
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

!!!�T�u���[�`��LU���� x��y��z�������L(�����p)---------------------------------------------------------------------------
subroutine LU_shuuki(MX,AA,L,U)
implicit none
integer :: MX
double precision :: AA(0:MX,0:MX)
double precision  :: L(0:MX,0:MX),U(0:MX,0:MX)
double precision :: LUcom(-2:2,0:MX)
double precision s1,s2
integer i,j,k
L=0.0d0               
U=0.0d0
do i=0,MX-1              !���J�n
L(i,i)=1.0d0
enddo
	do i=0,MX-1		
		do j=0,MX-1
		s1=0.0d0
		s2=0.0d0
	   	if(i<=j) then        !��O�p�s������(U)
	   		do k=0,i-1 		
	   		s1=s1+L(i,k)*U(k,j)
	   		enddo
	   	U(i,j)=AA(i,j)-s1
	   	else if(i>=j) then   !���O�p�s������(L)
	   		do k=0,j-1
	   		s2=s2+L(i,k)*U(k,j)
	   		enddo
	   	L(i,j)=(AA(i,j)-s2)/U(j,j)
	   	endif
	   	enddo
	enddo
	
!	do i = 0,MX-1                   !0�𔲂��Ĉ��k
!		LUcom(0,i) = U(i,i)	
!	enddo	
!	do i = 1,MX-2
!		LUcom(-1,i) = L(i,i-1)	
!		LUcom(1,i) = U(i,i+1)
!	enddo
!	
!	LUcom(-1,0)= 0
!	LUcom(1,0)=U(0,1)
!	LUcom(-1,MX-1)=L(MX-1,MX-2)
!	LUcom(1,MX-1)=0

end subroutine LU_shuuki

!!!�T�u���[�`���O�i��� z����-------------------------------------------------------------------------
subroutine forebackz(MX,MY,MZ,num,L,U,g2,QR)
implicit none
integer :: MX,MY,MZ,num
double precision ::  L(0:MZ,0:MZ),U(0:MZ,0:MZ)
double precision :: g2(1:num,0:MX,0:MY,0:MZ)
double precision :: QR(1:num,0:MX,0:MY,0:MZ)
double precision s1,s2
integer i,j,k,o,m

!LU������L�ɒ���
s1=0.0d0
do m=0,MY
	do o=0,MX
		do k=1,num
		g2(k,o,m,0)=QR(k,o,m,0)
			do i=1,MZ-1
				do j=0,i-1
				s1=s1+L(i,j)*g2(k,o,m,j)
				enddo
				g2(k,o,m,i)=QR(k,o,m,i)-s1
				s1=0.0d0
			enddo
		enddo
	enddo
enddo
!Lg=f��g(�O)�����߁CUg(��)=g(�O)��g(��)�����߂�
s2=0.0d0
do m=0,MY
	do o=0,MX
		do k=1,num
			do i=MZ-1,0,-1   !!!!!!MZ-1�ɂ��ĂȂ�
				do j=i+1,MZ-1
				s2=s2+U(i,j)*g2(k,o,m,j)
				enddo
				g2(k,o,m,i)=(g2(k,o,m,i)-s2)/U(i,i)
				s2=0.0d0
			enddo
		enddo
	enddo
enddo

g2(:,:,:,MZ)=g2(:,:,:,0)

end subroutine forebackz


!!!x���� F(FP��FM)�̒�`--------------------------------------------------------------------
subroutine FPMx(zeta,MX,MY,MZ,num,gam,Q,F,FP,FM)
implicit none
integer MX,MY,MZ,num
double precision :: zeta,gam
double precision :: Q(1:num,0:MX,0:MY,0:MZ)
double precision :: F(1:num,0:MX,0:MY,0:MZ),FP(1:num,0:MX,0:MY,0:MZ),FM(1:num,0:MX,0:MY,0:MZ)
integer i,j

	F(1,:,:,:)=Q(2,:,:,:)
	F(2,:,:,:)=((Q(2,:,:,:)**2.d0)/Q(1,:,:,:))+(gam-1.0d0)&
			   &*(Q(5,:,:,:)-(Q(2,:,:,:)**2.d0+Q(3,:,:,:)**2.d0+Q(4,:,:,:)**2.d0)/(2.d0*Q(1,:,:,:)))
	F(3,:,:,:)=(Q(3,:,:,:)*Q(2,:,:,:)/Q(1,:,:,:))
	F(4,:,:,:)=(Q(4,:,:,:)*Q(2,:,:,:)/Q(1,:,:,:))
	F(5,:,:,:)=(Q(5,:,:,:)+(gam-1.0d0)*(Q(5,:,:,:)-(Q(2,:,:,:)**2.d0+Q(3,:,:,:)**2.d0+Q(4,:,:,:)**2.d0)/(2.d0*Q(1,:,:,:))))&
	           &*Q(2,:,:,:)/Q(1,:,:,:)


	FP=0.5d0*(F+zeta*Q)
	FM=0.5d0*(F-1.0d0*zeta*Q)
	
end subroutine FPMx

!!!y���� F(FP��FM)�̒�`--------------------------------------------------------------------
subroutine FPMy(zeta,MX,MY,MZ,num,gam,Q,F,FP,FM)
implicit none
integer MX,MY,MZ,num
double precision :: zeta,gam
double precision :: Q(1:num,0:MX,0:MY,0:MZ)
double precision :: F(1:num,0:MX,0:MY,0:MZ),FP(1:num,0:MX,0:MY,0:MZ),FM(1:num,0:MX,0:MY,0:MZ)
integer i,j

	F(1,:,:,:)=Q(3,:,:,:)
	F(2,:,:,:)=(Q(3,:,:,:)*Q(2,:,:,:)/Q(1,:,:,:))
	F(3,:,:,:)=((Q(3,:,:,:)**2.d0)/Q(1,:,:,:))+(gam-1.d0)&
			   &*(Q(5,:,:,:)-(Q(2,:,:,:)**2.d0+Q(3,:,:,:)**2.d0+Q(4,:,:,:)**2.d0)/(2.d0*Q(1,:,:,:)))
	F(4,:,:,:)=(Q(3,:,:,:)*Q(4,:,:,:)/Q(1,:,:,:))
	F(5,:,:,:)=(Q(5,:,:,:)+(gam-1.d0)*(Q(5,:,:,:)-(Q(2,:,:,:)**2.d0+Q(3,:,:,:)**2.d0+Q(4,:,:,:)**2.d0)/(2.d0*Q(1,:,:,:))))&
	           &*Q(3,:,:,:)/Q(1,:,:,:)


	FP=0.5d0*(F+zeta*Q)
	FM=0.5d0*(F-1.0d0*zeta*Q)
	
end subroutine FPMy

!!!z���� F(FP��FM)�̒�`--------------------------------------------------------------------
subroutine FPMz(zeta,MX,MY,MZ,num,gam,Q,F,FP,FM)
implicit none
integer MX,MY,MZ,num
double precision :: zeta,gam
double precision :: Q(1:num,0:MX,0:MY,0:MZ)
double precision :: F(1:num,0:MX,0:MY,0:MZ),FP(1:num,0:MX,0:MY,0:MZ),FM(1:num,0:MX,0:MY,0:MZ)
integer i,j

	F(1,:,:,:)=Q(4,:,:,:)
	F(2,:,:,:)=(Q(4,:,:,:)*Q(2,:,:,:)/Q(1,:,:,:))
	F(3,:,:,:)=(Q(3,:,:,:)*Q(4,:,:,:)/Q(1,:,:,:))
	F(4,:,:,:)=((Q(4,:,:,:)**2.d0)/Q(1,:,:,:))+(gam-1.d0)&
			   &*(Q(5,:,:,:)-(Q(2,:,:,:)**2.d0+Q(3,:,:,:)**2.d0+Q(4,:,:,:)**2.d0)/(2.d0*Q(1,:,:,:)))
	F(5,:,:,:)=(Q(5,:,:,:)+(gam-1.d0)*(Q(5,:,:,:)-(Q(2,:,:,:)**2.d0+Q(3,:,:,:)**2.d0+Q(4,:,:,:)**2.d0)/(2.d0*Q(1,:,:,:))))&
	           &*Q(4,:,:,:)/Q(1,:,:,:)


	FP=0.5d0*(F+zeta*Q)
	FM=0.5d0*(F-1.0d0*zeta*Q)
	
end subroutine FPMz

!!!NSCBC�֘A�̃T�u���[�`��----------------------------------------------------------------------------------------
!!!��,u,P���`����---------------------------------------------------------------------
subroutine Q_rhoup(MX,MY,MZ,num,Q,rhoup,gam)
implicit none
integer MX,MY,MZ,num
double precision :: gam
double precision :: Q(1:num,0:MX,0:MY,0:MZ)
double precision :: rhoup(1:num,0:MX,0:MY,0:MZ)
integer i,j

	rhoup(1,:,:,:)=Q(1,:,:,:)                                           !!��
	rhoup(2,:,:,:)=Q(2,:,:,:)/Q(1,:,:,:)                                  !!u(x�������x)
	rhoup(3,:,:,:)=Q(3,:,:,:)/Q(1,:,:,:)									!!v(y�������x)
	rhoup(4,:,:,:)=Q(4,:,:,:)/Q(1,:,:,:)									!!w(z�������x)
	rhoup(5,:,:,:)=(gam-1.0d0)*(Q(5,:,:,:)-(Q(2,:,:,:)**2.d0+Q(3,:,:,:)**2.d0+Q(4,:,:,:)**2.d0)/(2.d0*Q(1,:,:,:)))       !!P

end subroutine Q_rhoup

!!!x���� NSCBC��L�s��̒�`----------------------------------------------------------------
subroutine NSCBC_Lx(MX,MY,MZ,num,NSLx,drhoup,rhoup,omega,gam,c0,cMX,Ma0,MaMX,Pmugen,Lx)
integer MX,MY,MZ,num
double precision :: omega,gam
double precision :: drhoup(1:num,0:MX,0:MY,0:MZ),rhoup(1:num,0:MX,0:MY,0:MZ)
double precision :: NSLx(1:num,0:1,0:MY,0:MZ)
double precision :: c0(0:MY,0:MZ),cMX(0:MY,0:MZ),Ma0(0:MY,0:MZ),MaMX(0:MY,0:MZ)
double precision :: Pmugen,Lx
integer i

c0(:,:)=sqrt(gam*rhoup(5,0,:,:)/rhoup(1,0,:,:))                !!�n�_�ł̉���
cMX(:,:)=sqrt(gam*rhoup(5,MX,:,:)/rhoup(1,MX,:,:))             !!�I�_�ł̉���
Ma0(:,:)=rhoup(2,0,:,:)/c0(:,:)                                !!�n�_�̃}�b�n��
MaMX(:,:)=rhoup(2,MX,:,:)/cMX(:,:)                             !!�I�_�̃}�b�n��

!!�����˗��o����----------------------------------------------------------------------------------------
!NSLx(1,0,:,:)=(rhoup(2,0,:,:)-c0(:,:))*(-rhoup(1,0,:,:)*c0(:,:)*drhoup(2,0,:,:)+drhoup(5,0,:,:))   !!!�����͎n�_��
!NSLx(2,0,:,:)=rhoup(2,0,:,:)*((c0(:,:)**2.d0)*drhoup(1,0,:,:)-drhoup(5,0,:,:))
!NSLx(3,0,:,:)=rhoup(2,0,:,:)*drhoup(3,0,:,:)           
!NSLx(4,0,:,:)=rhoup(2,0,:,:)*drhoup(4,0,:,:)           
!NSLx(5,0,:,:)=omega*(1.d0-(Ma0(:,:)**2.d0))*c0(:,:)*(rhoup(5,0,:,:)-Pmugen)/Lx  

NSLx(1,1,:,:)=omega*(1.d0-(MaMX(:,:)**2.d0))*cMX(:,:)*(rhoup(5,MX,:,:)-Pmugen)/Lx  !!!�����͏I�_��
NSLx(2,1,:,:)=rhoup(2,MX,:,:)*((cMX(:,:)**2.d0)*drhoup(1,MX,:,:)-drhoup(5,MX,:,:))
NSLx(3,1,:,:)=rhoup(2,MX,:,:)*drhoup(3,MX,:,:) 
NSLx(4,1,:,:)=rhoup(2,MX,:,:)*drhoup(4,MX,:,:) 	 
NSLx(5,1,:,:)=(rhoup(2,MX,:,:)+cMX(:,:))*(rhoup(1,MX,:,:)*cMX(:,:)*drhoup(2,MX,:,:)+drhoup(5,MX,:,:))
!!--------------------------------------------------------------------------------------------------------

!!��������������------------------------------------------------------------------------------------------
NSLx(1,0,:,:)=(rhoup(2,0,:,:)-c0(:,:))*(-1.d0*rhoup(1,0,:,:)*c0(:,:)*drhoup(2,0,:,:)+drhoup(5,0,:,:))  !!!�����͎n�_��
NSLx(5,0,:,:)=NSLx(1,0,:,:)
NSLx(2,0,:,:)=0.5d0*(gam-1.d0)*(NSLx(1,0,:,:)+NSLx(5,0,:,:))
NSLx(3,0,:,:)=0.d0          
NSLx(4,0,:,:)=0.d0           
!!----------------------------------------------------------------------------------------------------------

end subroutine NSCBC_Lx

!!!x���� NSCBC��d�s��̒�`----------------------------------------------------------------
subroutine NSCBC_dx(MX,MY,MZ,num,c0,cMX,NSLx,NSdx,rhoup)
integer MX,MY,MZ,num
double precision :: c0(0:MY,0:MZ),cMX(0:MY,0:MZ)
double precision :: NSLx(1:num,0:1,0:MY,0:MZ),rhoup(1:num,0:MX,0:MY,0:MZ)
double precision :: NSdx(1:num,0:1,0:MY,0:MZ)
integer i

NSdx(1,0,:,:)=(1.d0/(c0(:,:)**2.d0))*(0.5d0*(NSLx(1,0,:,:)+NSLx(5,0,:,:))+NSLx(2,0,:,:))
NSdx(2,0,:,:)=0.5d0*(NSLx(1,0,:,:)+NSLx(5,0,:,:))
NSdx(3,0,:,:)=1.d0/(2.d0*rhoup(1,0,:,:)*c0(:,:))*(-NSLx(1,0,:,:)+NSLx(5,0,:,:))
NSdx(4,0,:,:)=NSLx(3,0,:,:)
NSdx(5,0,:,:)=NSLx(4,0,:,:)

NSdx(1,1,:,:)=(1.d0/(cMX(:,:)**2.d0))*(0.5d0*(NSLx(1,1,:,:)+NSLx(5,1,:,:))+NSLx(2,1,:,:))
NSdx(2,1,:,:)=0.5d0*(NSLx(1,1,:,:)+NSLx(5,1,:,:))
NSdx(3,1,:,:)=0.5d0/(rhoup(1,MX,:,:)*cMX(:,:))*(-NSLx(1,1,:,:)+NSLx(5,1,:,:))
NSdx(4,1,:,:)=NSLx(3,1,:,:)
NSdx(5,1,:,:)=NSLx(4,1,:,:)

end subroutine NSCBC_dx

!!!x����NSCBC x�����̔��������̍����ւ�------------------------------------------------
subroutine NSCBC_dFx(MX,MY,MZ,num,NSdx,rhoup,dFx,gam)
integer MX,MY,MZ,num
double precision :: NSdx(1:num,0:1,0:MY,0:MZ),rhoup(1:num,0:MX,0:MY,0:MZ)
double precision :: dFx(1:num,0:MX,0:MY,0:MZ)
double precision :: gam
integer i
 
dFx(1,0,:,:)=NSdx(1,0,:,:)
dFx(2,0,:,:)=rhoup(2,0,:,:)*NSdx(1,0,:,:)+rhoup(1,0,:,:)*NSdx(3,0,:,:)
dFx(3,0,:,:)=rhoup(3,0,:,:)*NSdx(1,0,:,:)+rhoup(1,0,:,:)*NSdx(4,0,:,:)
dFx(4,0,:,:)=rhoup(4,0,:,:)*NSdx(1,0,:,:)+rhoup(1,0,:,:)*NSdx(5,0,:,:)
dFx(5,0,:,:)=0.5d0*((rhoup(2,0,:,:)**2.d0)+(rhoup(3,0,:,:)**2.d0)+(rhoup(4,0,:,:)**2.d0))*NSdx(1,0,:,:)&
			 &+NSdx(2,0,:,:)/(gam-1.d0)&
		     &+rhoup(1,0,:,:)*(rhoup(2,0,:,:)*NSdx(3,0,:,:)+rhoup(3,0,:,:)*NSdx(4,0,:,:)+rhoup(4,0,:,:)*NSdx(5,0,:,:))

dFx(1,MX,:,:)=NSdx(1,1,:,:)
dFx(2,MX,:,:)=rhoup(2,MX,:,:)*NSdx(1,1,:,:)+rhoup(1,MX,:,:)*NSdx(3,1,:,:)
dFx(3,MX,:,:)=rhoup(3,MX,:,:)*NSdx(1,1,:,:)+rhoup(1,MX,:,:)*NSdx(4,1,:,:)
dFx(4,MX,:,:)=rhoup(4,MX,:,:)*NSdx(1,1,:,:)+rhoup(1,MX,:,:)*NSdx(5,1,:,:)
dFx(5,MX,:,:)=0.5d0*((rhoup(2,MX,:,:)**2.d0)+(rhoup(3,MX,:,:)**2.d0)+(rhoup(4,MX,:,:)**2.d0))*NSdx(1,1,:,:)&
			  &+NSdx(2,1,:,:)/(gam-1.d0)&
		      &+rhoup(1,MX,:,:)*(rhoup(2,MX,:,:)*NSdx(3,1,:,:)+rhoup(3,MX,:,:)*NSdx(4,1,:,:)+rhoup(4,MX,:,:)*NSdx(5,1,:,:))

!!!x��y������NSCBC�̏d�Ȃ镔���̒���-------------------------------------------------------		   
dFx(:,0,0,:)=0.5d0*dFx(:,0,0,:)     !(0.5����Ȃ��Ă������D�Ⴆ�ΕǏ����Ȃ�ǂ̉e�����������邽�߂����Ƒ傫���l�ŏ����Ă������Ddfx=2/3*dfx�Ƃ�)
dFx(:,MX,0,:)=0.5d0*dFx(:,MX,0,:)
dFx(:,0,MY,:)=0.5d0*dFx(:,0,MY,:)
dFx(:,MX,MY,:)=0.5d0*dFx(:,MX,MY,:)
!!!----------------------------------------------------------------------------------------	

end subroutine NSCBC_dFx

!!!y���� NSCBC��L�s��̒�`----------------------------------------------------------------
subroutine NSCBC_Ly(MX,MY,MZ,num,NSLy,drhoup,rhoup,omega,gam,c0,cMX,Ma0,MaMX,Pmugen,Ly)
integer MX,MY,MZ,num
double precision :: omega,gam
double precision :: drhoup(1:num,0:MX,0:MY,0:MZ),rhoup(1:num,0:MX,0:MY,0:MZ)
double precision :: NSLy(1:num,0:MX,0:1,0:MZ)
double precision :: c0(0:MX,0:MZ),cMX(0:MX,0:MZ),Ma0(0:MX,0:MZ),MaMX(0:MX,0:MZ)
double precision :: Pmugen,Ly
integer i

c0(:,:)=sqrt(gam*rhoup(5,:,0,:)/rhoup(1,:,0,:))                !!�n�_�ł̉���
cMX(:,:)=sqrt(gam*rhoup(5,:,MY,:)/rhoup(1,:,MY,:))             !!�I�_�ł̉���
Ma0(:,:)=rhoup(3,:,0,:)/c0(:,:)                                !!�n�_�̃}�b�n��
MaMX(:,:)=rhoup(3,:,MY,:)/cMX(:,:)                             !!�I�_�̃}�b�n��

!!�����˗��o����---------------------------------------------------------------------------------------- 
!NSLy(1,:,0,:)=(rhoup(3,:,0,:)-c0(:,:))*(-rhoup(1,:,0,:)*c0(:,:)*drhoup(3,:,0,:)+drhoup(5,:,0,:))   !!!�����͎n�_��
!NSLy(2,:,0,:)=rhoup(3,:,0,:)*((c0(:,:)**2.d0)*drhoup(1,:,0,:)-drhoup(5,:,0,:))
!NSLy(3,:,0,:)=rhoup(3,:,0,:)*drhoup(2,:,0,:)           
!NSLy(4,:,0,:)=rhoup(3,:,0,:)*drhoup(4,:,0,:)
!NSLy(5,:,0,:)=omega*(1.d0-(Ma0(:,:)**2.d0))*c0(:,:)*(rhoup(5,:,0,:)-Pmugen)/Ly   

NSLy(1,:,1,:)=omega*(1.d0-(MaMX(:,:)**2.d0))*cMX(:,:)*(rhoup(5,:,MY,:)-Pmugen)/Ly  !!!�����͏I�_��
NSLy(2,:,1,:)=rhoup(3,:,MY,:)*((cMX(:,:)**2.d0)*drhoup(1,:,MY,:)-drhoup(5,:,MY,:))
NSLy(3,:,1,:)=rhoup(3,:,MY,:)*drhoup(2,:,MY,:) 
NSLy(4,:,1,:)=rhoup(3,:,MY,:)*drhoup(4,:,MY,:)
NSLy(5,:,1,:)=(rhoup(3,:,MY,:)+cMX(:,:))*(rhoup(1,:,MY,:)*cMX(:,:)*drhoup(3,:,MY,:)+drhoup(5,:,MY,:))
!!---------------------------------------------------------------------------------------------------------

!!�������薳���ǖʏ���-----------------------------------------------------------------------------------
NSLy(1,:,0,:)=(rhoup(3,:,0,:)-c0(:,:))*(-1.d0*rhoup(1,:,0,:)*c0(:,:)*drhoup(3,:,0,:)+drhoup(5,:,0,:))  !!!�����͎n�_��
NSLy(2,:,0,:)=0.d0
NSLy(3,:,0,:)=0.d0           
NSLy(4,:,0,:)=0.d0           
NSLy(5,:,0,:)=NSLy(1,:,0,:)
!!-----------------------------------------------------------------------------------------------------------

end subroutine NSCBC_Ly

!!!y���� NSCBC��d�s��̒�`----------------------------------------------------------------
subroutine NSCBC_dy(MX,MY,MZ,num,c0,cMX,NSLy,NSdy,rhoup)
integer MX,MY,MZ,num
double precision :: c0(0:MX,0:MZ),cMX(0:MX,0:MZ)
double precision :: NSLy(1:num,0:MX,0:1,0:MZ),rhoup(1:num,0:MX,0:MY,0:MZ)
double precision :: NSdy(1:num,0:MX,0:1,0:MZ)
integer i

NSdy(1,:,0,:)=(1.d0/(c0(:,:)**2.d0))*(0.5d0*(NSLy(1,:,0,:)+NSLy(5,:,0,:))+NSLy(2,:,0,:))
NSdy(2,:,0,:)=0.5d0*(NSLy(1,:,0,:)+NSLy(5,:,0,:))
NSdy(3,:,0,:)=NSLy(3,:,0,:)
NSdy(4,:,0,:)=1.d0/(2.d0*rhoup(1,:,0,:)*c0(:,:))*(-NSLy(1,:,0,:)+NSLy(5,:,0,:))
NSdy(5,:,0,:)=NSLy(4,:,0,:)

NSdy(1,:,1,:)=(1.d0/(cMX(:,:)**2.d0))*(0.5d0*(NSLy(1,:,1,:)+NSLy(5,:,1,:))+NSLy(2,:,1,:))
NSdy(2,:,1,:)=0.5d0*(NSLy(1,:,1,:)+NSLy(5,:,1,:))
NSdy(3,:,1,:)=NSLy(3,:,1,:)
NSdy(4,:,1,:)=0.5d0/(rhoup(1,:,MY,:)*cMX(:,:))*(-NSLy(1,:,1,:)+NSLy(5,:,1,:))
NSdy(5,:,1,:)=NSLy(4,:,1,:)

end subroutine NSCBC_dy

!!!y����NSCBC y�����̔��������̍����ւ�------------------------------------------------
subroutine NSCBC_dFy(MX,MY,MZ,num,NSdy,rhoup,dFy,gam)
integer MX,MY,MZ,num
double precision :: NSdy(1:num,0:MX,0:1,0:MZ),rhoup(1:num,0:MX,0:MY,0:MZ)
double precision :: dFy(1:num,0:MX,0:MY,0:MZ)
double precision :: gam
integer i

dFy(1,:,0,:)=NSdy(1,:,0,:)
dFy(2,:,0,:)=rhoup(2,:,0,:)*NSdy(1,:,0,:)+rhoup(1,:,0,:)*NSdy(3,:,0,:)
dFy(3,:,0,:)=rhoup(3,:,0,:)*NSdy(1,:,0,:)+rhoup(1,:,0,:)*NSdy(4,:,0,:)
dFy(4,:,0,:)=rhoup(4,:,0,:)*NSdy(1,:,0,:)+rhoup(1,:,0,:)*NSdy(5,:,0,:)
dFy(5,:,0,:)=0.5d0*((rhoup(2,:,0,:)**2.d0)+(rhoup(3,:,0,:)**2.d0)+(rhoup(4,:,0,:)**2.d0))*NSdy(1,:,0,:)+NSdy(2,:,0,:)/(gam-1.d0)&
		  &+rhoup(1,:,0,:)*(rhoup(2,:,0,:)*NSdy(3,:,0,:)+rhoup(3,:,0,:)*NSdy(4,:,0,:)+rhoup(4,:,0,:)*NSdy(5,:,0,:))

do i=0,MZ
if(mod(i,2)==0) then 
dFy(1,0,0,i)=1.4d0*dFy(1,0,0,i)
endif
end do


dFy(1,:,MY,:)=NSdy(1,:,1,:)
dFy(2,:,MY,:)=rhoup(2,:,MY,:)*NSdy(1,:,1,:)+rhoup(1,:,MY,:)*NSdy(3,:,1,:)
dFy(3,:,MY,:)=rhoup(3,:,MY,:)*NSdy(1,:,1,:)+rhoup(1,:,MY,:)*NSdy(4,:,1,:)
dFy(4,:,MY,:)=rhoup(4,:,MY,:)*NSdy(1,:,1,:)+rhoup(1,:,MY,:)*NSdy(5,:,1,:)
dFy(5,:,MY,:)=0.5d0*((rhoup(2,:,MY,:)**2.d0)+(rhoup(3,:,MY,:)**2.d0)+(rhoup(4,:,MY,:)**2.d0))*NSdy(1,:,1,:)&
			  &+NSdy(2,:,1,:)/(gam-1.d0)&
		      &+rhoup(1,:,MY,:)*(rhoup(2,:,MY,:)*NSdy(3,:,1,:)+rhoup(3,:,MY,:)*NSdy(4,:,1,:)+rhoup(4,:,MY,:)*NSdy(5,:,1,:))

!!!x��y������NSCBC�̏d�Ȃ镔���̒���-------------------------------------------------------
dFy(:,0,0,:)=0.5d0*dFy(:,0,0,:)
dFy(:,0,MY,:)=0.5d0*dFy(:,0,MY,:)
dFy(:,MX,0,:)=0.5d0*dFy(:,MX,0,:)
dFy(:,MX,MY,:)=0.5d0*dFy(:,MX,MY,:)
!!!----------------------------------------------------------------------------------------

end subroutine NSCBC_dFy


!!!z���� NSCBC��L�s��̒�`----------------------------------------------------------------
subroutine NSCBC_Lz(MX,MY,MZ,num,NSLz,drhoup,rhoup,omega,gam,c0,cMX,Ma0,MaMX,Pmugen,Lz)
integer MX,MY,MZ,num
double precision :: omega,gam
double precision :: drhoup(1:num,0:MX,0:MY,0:MZ),rhoup(1:num,0:MX,0:MY,0:MZ)
double precision :: NSLz(1:num,0:MX,0:MY,0:1)
double precision :: c0(0:MX,0:MY),cMX(0:MX,0:MY),Ma0(0:MX,0:MY),MaMX(0:MX,0:MY)
double precision :: Pmugen,Lz
integer i

c0(:,:)=sqrt(gam*rhoup(5,:,:,0)/rhoup(1,:,:,0))                !!�n�_�ł̉���
cMX(:,:)=sqrt(gam*rhoup(5,:,:,MZ)/rhoup(1,:,:,MZ))             !!�I�_�ł̉���
Ma0(:,:)=rhoup(4,:,:,0)/c0(:,:)                                !!�n�_�̃}�b�n��
MaMX(:,:)=rhoup(4,:,:,MZ)/cMX(:,:)                             !!�I�_�̃}�b�n��

NSLz(1,:,:,0)=(rhoup(4,:,:,0)-c0(:,:))*(-rhoup(1,:,:,0)*c0(:,:)*drhoup(4,:,:,0)+drhoup(5,:,:,0))   !!!�����͎n�_��
NSLz(2,:,:,0)=rhoup(4,:,:,0)*((c0(:,:)**2.d0)*drhoup(1,:,:,0)-drhoup(5,:,:,0))
NSLz(3,:,:,0)=rhoup(4,:,:,0)*drhoup(2,:,:,0)           
NSLz(4,:,:,0)=rhoup(4,:,:,0)*drhoup(3,:,:,0)
NSLz(5,:,:,0)=omega*(1.d0-(Ma0(:,:)**2.d0))*c0(:,:)*(rhoup(5,:,:,0)-Pmugen)/Lz   

NSLz(1,:,:,1)=omega*(1.d0-(MaMX(:,:)**2.d0))*cMX(:,:)*(rhoup(5,:,:,MZ)-Pmugen)/Lz  !!!�����͏I�_��
NSLz(2,:,:,1)=rhoup(4,:,:,MZ)*((cMX(:,:)**2.d0)*drhoup(1,:,:,MZ)-drhoup(5,:,:,MZ))
NSLz(3,:,:,1)=rhoup(4,:,:,MZ)*drhoup(2,:,:,MZ) 
NSLz(4,:,:,1)=rhoup(4,:,:,MZ)*drhoup(3,:,:,MZ)
NSLz(5,:,:,1)=(rhoup(4,:,:,MZ)+cMX(:,:))*(rhoup(1,:,:,MZ)*cMX(:,:)*drhoup(4,:,:,MZ)+drhoup(5,:,:,MZ))

end subroutine NSCBC_Lz

!!!z���� NSCBC��d�s��̒�`----------------------------------------------------------------
subroutine NSCBC_dz(MX,MY,MZ,num,c0,cMX,NSLz,NSdz,rhoup)
integer MX,MY,MZ,num
double precision :: c0(0:MX,0:MY),cMX(0:MX,0:MY)
double precision :: NSLz(1:num,0:MX,0:MY,0:1),rhoup(1:num,0:MX,0:MY,0:MZ)
double precision :: NSdz(1:num,0:MX,0:MY,0:1)
integer i

NSdz(1,:,:,0)=(1.d0/(c0(:,:)**2.d0))*(0.5d0*(NSLz(1,:,:,0)+NSLz(5,:,:,0))+NSLz(2,:,:,0))
NSdz(2,:,:,0)=0.5d0*(NSLz(1,:,:,0)+NSLz(5,:,:,0))
NSdz(3,:,:,0)=NSLz(3,:,:,0)
NSdz(4,:,:,0)=NSLz(4,:,:,0)
NSdz(5,:,:,0)=1.d0/(2.d0*rhoup(1,:,:,0)*c0(:,:))*(-NSLz(1,:,:,0)+NSLz(5,:,:,0))

NSdz(1,:,:,1)=(1.d0/(cMX(:,:)**2.d0))*(0.5d0*(NSLz(1,:,:,1)+NSLz(5,:,:,1))+NSLz(2,:,:,1))
NSdz(2,:,:,1)=0.5d0*(NSLz(1,:,:,1)+NSLz(5,:,:,1))
NSdz(3,:,:,1)=NSLz(3,:,:,1)
NSdz(4,:,:,1)=NSLz(4,:,:,1)
NSdz(5,:,:,1)=0.5d0/(rhoup(1,:,:,1)*cMX(:,:))*(-NSLz(1,:,:,1)+NSLz(5,:,:,1))

end subroutine NSCBC_dz

!!!z����NSCBC z�����̔��������̍����ւ�------------------------------------------------
subroutine NSCBC_dFz(MX,MY,MZ,num,NSdz,rhoup,dFz,gam)
integer MX,MY,MZ,num
double precision :: NSdz(1:num,0:MX,0:MY,0:1),rhoup(1:num,0:MX,0:MY,0:MZ)
double precision :: dFz(1:num,0:MX,0:MY,0:MZ)
double precision :: gam
integer i

dFz(1,:,:,0)=NSdz(1,:,:,0)
dFz(2,:,:,0)=rhoup(2,:,:,0)*NSdz(1,:,:,0)+rhoup(1,:,:,0)*NSdz(3,:,:,0)
dFz(3,:,:,0)=rhoup(3,:,:,0)*NSdz(1,:,:,0)+rhoup(1,:,:,0)*NSdz(4,:,:,0)
dFz(4,:,:,0)=rhoup(4,:,:,0)*NSdz(1,:,:,0)+rhoup(1,:,:,0)*NSdz(5,:,:,0)
dFz(5,:,:,0)=0.5d0*((rhoup(2,:,:,0)**2.d0)+(rhoup(3,:,:,0)**2.d0)+(rhoup(4,:,:,0)**2.d0))*NSdz(1,:,:,0)+NSdz(2,:,:,0)/(gam-1.d0)&
		  &+rhoup(1,:,:,0)*(rhoup(2,:,:,0)*NSdz(3,:,:,0)+rhoup(3,:,:,0)*NSdz(4,:,:,0)+rhoup(4,:,:,0)*NSdz(5,:,:,0))

dFz(1,:,:,MZ)=NSdz(1,:,:,1)
dFz(2,:,:,MZ)=rhoup(2,:,:,MZ)*NSdz(1,:,:,1)+rhoup(1,:,:,MZ)*NSdz(3,:,:,1)
dFz(3,:,:,MZ)=rhoup(3,:,:,MZ)*NSdz(1,:,:,1)+rhoup(1,:,:,MZ)*NSdz(4,:,:,1)
dFz(4,:,:,MZ)=rhoup(4,:,:,MZ)*NSdz(1,:,:,1)+rhoup(1,:,:,MZ)*NSdz(5,:,:,1)
dFz(5,:,:,MZ)=0.5d0*((rhoup(2,:,:,MZ)**2.d0)+(rhoup(3,:,:,MZ)**2.d0)+(rhoup(4,:,:,MZ)**2.d0))*NSdz(1,:,:,1)&
			  &+NSdz(2,:,:,1)/(gam-1.d0)&
		      &+rhoup(1,:,:,MZ)*(rhoup(2,:,:,MZ)*NSdz(3,:,:,1)+rhoup(3,:,:,MZ)*NSdz(4,:,:,1)+rhoup(4,:,:,MZ)*NSdz(5,:,:,1))

end subroutine NSCBC_dFz
!!!!!!!!!!!!!!!!!!!NSCBC�֘A �����܂�-----------------------------------------------------------

!!!Dirichlet�̋��E����--------------------------------------------------------------------
subroutine DBound(MX,Q,Bound1,BoundMX)
implicit none
integer MX
double precision :: Q(-2:MX+2,-2:MX+2)
double precision :: Bound1(-2:MX+2),BoundMX(-2:MX+2)
integer i
do i=1,4
Q(i,0)=Bound1(i)
Q(i,MX)=BoundMX(i)
enddo

end subroutine DBound

!!!Neumann�̋��E����--------------------------------------------------------------------
subroutine NBound(MX,MY,Q)
implicit none
integer MX,MY
double precision :: Q(-2:MX+2,-2:MX+2,-2:MX+2)
integer i,j
do j=1,MY-1
	do i=1,4
	Q(i,0,j)=Q(i,1,j)
	Q(i,MX,j)=Q(i,MX-1,j)
	enddo
enddo

do j=1,MX-1
	do i=1,4
	Q(i,j,0)=Q(i,j,1)
	Q(i,j,MY)=Q(i,j,MY-1)
	enddo
enddo

end subroutine NBound


!!!Q����u,v,T�̃T�u���[�`��---------------------------------------------------------------------
subroutine Q_uT(MX,MY,MZ,num,gam,Ma,Q,Tuv)
implicit none
integer :: MX,MY,MZ,num
double precision :: gam,Ma
double precision :: Tuv(1:num,0:MX,0:MY,0:MZ),Q(1:num,0:MX,0:MY,0:MZ)

	Tuv(2,:,:,:)=Q(2,:,:,:)/Q(1,:,:,:)                          !u�̓��o
	Tuv(3,:,:,:)=Q(3,:,:,:)/Q(1,:,:,:)                          !v�̓��o
	Tuv(4,:,:,:)=Q(4,:,:,:)/Q(1,:,:,:)                          !w�̓��o
	Tuv(5,:,:,:)=gam*(Ma**2.d0)*(gam-1.d0)*(Q(5,:,:,:)&
	             &-((Q(4,:,:,:)**2.d0)+(Q(3,:,:,:)**2.d0)+(Q(2,:,:,:)**2.d0))/(2.d0*Q(1,:,:,:)))/Q(1,:,:,:)  !T�̓��o

end subroutine Q_uT

!!!du��dT����S�����̓��o����єS�x�ʂ̓��o x����------------------------------------------------------------
subroutine viscosityx(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvis,dTux,dTuy,dTuz,mu)
implicit none
integer :: MX,MY,MZ,num
double precision :: gam,Re,Pr,Ma
double precision :: mu(0:MX,0:MY,0:MZ),Tuvis(1:num,0:MX,0:MY,0:MZ),Sc(0:MX,0:MY,0:MZ)
double precision :: Tuv(1:num,0:MX,0:MY,0:MZ),dTux(1:num,0:MX,0:MY,0:MZ),dTuy(1:num,0:MX,0:MY,0:MZ),dTuz(1:num,0:MX,0:MY,0:MZ)

	mu(:,:,:)=Tuv(5,:,:,:)**(2.d0/3.d0)

	Tuvis(2,:,:,:)=(2.d0/3.d0)*mu(:,:,:)*(2.d0*dTux(2,:,:,:)-dTuy(3,:,:,:)-dTuz(4,:,:,:))/Re
	Tuvis(3,:,:,:)=mu(:,:,:)*(dTuy(2,:,:,:)+dTux(3,:,:,:))/Re
	Tuvis(4,:,:,:)=mu(:,:,:)*(dTuz(2,:,:,:)+dTux(4,:,:,:))/Re
	Tuvis(5,:,:,:)=Tuvis(2,:,:,:)*Tuv(2,:,:,:)+Tuvis(3,:,:,:)*Tuv(3,:,:,:)+Tuvis(4,:,:,:)*Tuv(4,:,:,:)&
				   &+mu(:,:,:)*dTux(5,:,:,:)/((gam-1.d0)*Re*Pr*(Ma**2.d0))

end subroutine viscosityx

!!!du��dT����S�����̓��o����єS�x�ʂ̓��o y����------------------------------------------------------------
subroutine viscosityy(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvis,dTux,dTuy,dTuz,mu)
implicit none
integer :: MX,MY,MZ,num
double precision :: gam,Re,Pr,Ma
double precision :: mu(0:MX,0:MY,0:MZ),Tuvis(1:num,0:MX,0:MY,0:MZ),Sc(0:MX,0:MY,0:MZ)
double precision :: Tuv(1:num,0:MX,0:MY,0:MZ),dTux(1:num,0:MX,0:MY,0:MZ),dTuy(1:num,0:MX,0:MY,0:MZ),dTuz(1:num,0:MX,0:MY,0:MZ)

	mu(:,:,:)=Tuv(5,:,:,:)**(2.d0/3.d0)

	Tuvis(2,:,:,:)=mu(:,:,:)*(dTux(3,:,:,:)+dTuy(2,:,:,:))/Re
	Tuvis(3,:,:,:)=(2.d0/3.d0)*mu(:,:,:)*(-1.0d0*dTux(2,:,:,:)+2.d0*dTuy(3,:,:,:)-1.0d0*dTuz(4,:,:,:))/Re
	Tuvis(4,:,:,:)=mu(:,:,:)*(dTuy(4,:,:,:)+dTuz(3,:,:,:))/Re
	Tuvis(5,:,:,:)=Tuvis(2,:,:,:)*Tuv(2,:,:,:)+Tuvis(3,:,:,:)*Tuv(3,:,:,:)+Tuvis(4,:,:,:)*Tuv(4,:,:,:)&
				   &+mu(:,:,:)*dTuy(5,:,:,:)/((gam-1.d0)*Re*Pr*(Ma**2.d0))
				   
end subroutine viscosityy

!!!�����˗��o�����̃T�u���[�`��----------------------------------------------------------------------------------
subroutine visout(MX,MY,MZ,num,ggTux,ggTuy,grhoupx,grhoupy,Tuvisx,Tuvisy,rhoup)
integer :: MX,MY,MZ,num
double precision :: ggTux(1:num,0:MX,0:MY,0:MZ),ggTuy(1:num,0:MX,0:MY,0:MZ),grhoupx(1:num,0:MX,0:MY,0:MZ)
double precision :: grhoupy(1:num,0:MX,0:MY,0:MZ)
double precision :: Tuvisx(1:num,0:MX,0:MY,0:MZ),Tuvisy(1:num,0:MX,0:MY,0:MZ),rhoup(1:num,0:MX,0:MY,0:MZ)

ggTux(3,MX,:,:)=0.d0
ggTux(4,MX,:,:)=0.d0
ggTux(5,MX,:,:)=ggTux(2,MX,:,:)*rhoup(2,MX,:,:)+grhoupx(2,MX,:,:)*Tuvisx(2,MX,:,:)+grhoupx(3,MX,:,:)*Tuvisx(3,MX,:,:)&
				&+grhoupx(4,MX,:,:)*Tuvisx(4,MX,:,:)
				
ggTuy(2,:,MY,:)=0.d0
ggTuy(4,:,MY,:)=0.d0
ggTuy(5,:,MY,:)=ggTuy(3,:,MY,:)*rhoup(3,:,MY,:)+grhoupy(2,:,MY,:)*Tuvisy(2,:,MY,:)+grhoupy(3,:,MY,:)*Tuvisy(3,:,MY,:)&
				&+grhoupy(4,:,MY,:)*Tuvisy(4,:,MY,:)

end subroutine visout

!!!du��dT����S�����̓��o����єS�x�ʂ̓��o z����------------------------------------------------------------
subroutine viscosityz(MX,MY,MZ,num,gam,Sc,Re,Pr,Ma,Tuv,Tuvis,dTux,dTuy,dTuz,mu)
implicit none
integer :: MX,MY,MZ,num
double precision :: gam,Re,Pr,Ma
double precision :: mu(0:MX,0:MY,0:MZ),Tuvis(1:num,0:MX,0:MY,0:MZ),Sc(0:MX,0:MY,0:MZ)
double precision :: Tuv(1:num,0:MX,0:MY,0:MZ),dTux(1:num,0:MX,0:MY,0:MZ),dTuy(1:num,0:MX,0:MY,0:MZ),dTuz(1:num,0:MX,0:MY,0:MZ)

	mu(:,:,:)=Tuv(5,:,:,:)**(2.d0/3.d0)

	Tuvis(2,:,:,:)=mu(:,:,:)*(dTuz(2,:,:,:)+dTux(4,:,:,:))/Re
	Tuvis(3,:,:,:)=mu(:,:,:)*(dTuz(3,:,:,:)+dTuy(4,:,:,:))/Re
	Tuvis(4,:,:,:)=(2.d0/3.d0)*mu(:,:,:)*(-1.0d0*dTux(2,:,:,:)-1.d0*dTuy(3,:,:,:)+2.0d0*dTuz(4,:,:,:))/Re
	Tuvis(5,:,:,:)=Tuvis(2,:,:,:)*Tuv(2,:,:,:)+Tuvis(3,:,:,:)*Tuv(3,:,:,:)+Tuvis(4,:,:,:)*Tuv(4,:,:,:)&
				   &+mu(:,:,:)*dTuz(5,:,:,:)/((gam-1.d0)*Re*Pr*(Ma**2.d0))

!!�����˗��o����----------------------------------
!	Tuvis(2,:,:,0)=0.d0
!	Tuvis(2,:,:,MZ)=0.d0
!	
!	Tuvis(3,:,:,0)=0.d0
!	Tuvis(3,:,:,MZ)=0.d0
!	
!	dTuz(5,:,:,0)=0.d0
!	dTuz(5,:,:,MZ)=0.d0	
!	Tuvis(5,:,:,:)=Tuvis(2,:,:,:)*Tuv(2,:,:,:)+Tuvis(3,:,:,:)*Tuv(3,:,:,:)+Tuvis(4,:,:,:)*Tuv(4,:,:,:)&
!				   &+mu(:,:,:)*dTuz(5,:,:,:)/((gam-1.d0)*Re*Pr*(Ma**2.d0))


end subroutine viscosityz


!!!�����_�������̒�`-----------------------------------------------------------------------------------------------
subroutine kakuran_Q(MX,MY,MZ,MXt,num,s,gam,Ma,Q,kakuran_u,kakuran_v,kakuran_w,times,u,v,w,T)
implicit none
integer :: MX,MY,MZ,MXt,num,s,times
double precision :: gam,Ma
double precision :: Q(1:num,0:MX,0:MY,0:MZ),kari
double precision :: kakuran_u(0:MXt,0:MY,0:MZ),kakuran_v(0:MXt,0:MY,0:MZ),kakuran_w(0:MXt,0:MY,0:MZ)
double precision :: u(0:MX,0:MY,0:MZ),v(0:MX,0:MY,0:MZ),w(0:MX,0:MY,0:MZ),T(0:MX,0:MY,0:MZ)
integer i,j,k

kari=0.065d0*dble(times)/450d0
if(times > 450) kari=0.065d0
i=0
do k=0,MZ                                              !�����ő��֐��̖������ʂ���(0<Ly<5��10<Ly<20�́A�������Ȃ��B������)
		do j=0,MY
!		if(j <= 79 .or. j > 112 ) then 
!			Q(1,i,j,k)=rhoup(1,i,j,k)
!			Q(2,i,j,k)=rhoup(1,i,j,k)*(rhoup(2,i,j,k))
!			Q(3,i,j,k)=rhoup(1,i,j,k)*(rhoup(3,i,j,k))
!			Q(4,i,j,k)=rhoup(1,i,j,k)*(rhoup(4,i,j,k))
!			Q(5,i,j,k)=rhoup(5,i,j,k)/(gam-1.d0)+0.5d0*rhoup(1,i,j,k)&
!			           &*(rhoup(2,i,j,k)**2.d0+rhoup(3,i,j,k)**2.d0+rhoup(4,i,j,k)**2.d0)	                  
!		else if(j > 79 .and. j <= 112) then


!			Q(2,i,j,k)=rhoup(1,i,j,k)*(rhoup(2,i,j,k)+kakuran_u(s,j,k))
!			Q(3,i,j,k)=rhoup(1,i,j,k)*(rhoup(3,i,j,k)+kakuran_v(s,j,k))
!			Q(4,i,j,k)=rhoup(1,i,j,k)*(rhoup(4,i,j,k)+kakuran_w(s,j,k))
!			Q(5,i,j,k)=rhoup(5,i,j,k)/(gam-1.d0)+0.5d0*rhoup(1,i,j,k)&
!			           &*((rhoup(2,i,j,k)+kakuran_u(s,j,k))**2.d0+(rhoup(3,i,j,k)+kakuran_v(s,j,k))**2.d0&
!			           &+(rhoup(4,i,j,k)+kakuran_w(s,j,k))**2.d0)	
!

			Q(2,i,j,k)=Q(1,i,j,k)*(u(i,j,k)+kakuran_u(s,j,k)*kari)
			Q(3,i,j,k)=Q(1,i,j,k)*(v(i,j,k)+kakuran_v(s,j,k)*kari)
			Q(4,i,j,k)=Q(1,i,j,k)*(w(i,j,k)+kakuran_w(s,j,k)*kari)
			Q(5,i,j,k)=(Q(1,i,j,k)/(gam*(Ma**2.d0)*(gam-1.d0))*T(i,j,k)&
                       &+0.5d0*Q(1,i,j,k)&
		               &*((u(i,j,k)+kakuran_u(s,j,k)*kari)**2.d0+(v(i,j,k)&
			           &+kakuran_v(s,j,k)*kari)**2.d0+(w(i,j,k)+kakuran_w(s,j,k)*kari)**2.d0)	)

!			Q(2,i,j,k)=Q(1,i,j,k)*(Q(2,i,j,k)/Q(1,i,j,k)+kakuran_u(s,j,k))
!			Q(3,i,j,k)=Q(1,i,j,k)*(Q(3,i,j,k)/Q(1,i,j,k)+kakuran_v(s,j,k))
!			Q(4,i,j,k)=Q(1,i,j,k)*(Q(4,i,j,k)/Q(1,i,j,k)+kakuran_w(s,j,k))
!			Q(5,i,j,k)=(Q(1,i,j,k)*(gam*(Ma**2.d0)*(gam-1.d0)*(Q(5,i,j,k)&
!	                   &-((Q(4,i,j,k)**2.d0)+(Q(3,i,j,k)**2.d0)+(Q(2,i,j,k)**2.d0))/(2.d0*Q(1,i,j,k)))/Q(1,i,j,k))& 
!			           &/gam/(Ma**2.d0))/(gam-1.d0)&
!                       &+0.5d0*Q(1,i,j,k)&
!		               &*((Q(2,i,j,k)/Q(1,i,j,k)+kakuran_u(s,j,k))**2.d0+(Q(3,i,j,k)/Q(1,i,j,k)&
!			           &+kakuran_v(s,j,k))**2.d0+(Q(4,i,j,k)/Q(1,i,j,k)+kakuran_w(s,j,k))**2.d0)	                   
!		endif
		enddo
enddo

end subroutine kakuran_Q

!!!�t�@�[�u������-----------------------------------------------------------------------------------------------------
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
open(41,file='��_heikin'//trim(Chara)//'.txt',status='unknown',form='formatted')
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

!!�Q�\���̉����p--------------------------------------------------------------------------------
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
!!z�����̕��ς����---------------------------
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

!!!�T�u���[�`���E�Ӓ�` y����(�s��x,y�����̂�) LU�̑O�i���------------------------------------------------------------------------
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

!!  �O�i���
	Y(:,:,0)=QR(:,:,0)
	do j=1,MY
		Y(:,:,j)=QR(:,:,j)-LUcom(-1,j) * Y(:,:,j-1)
	enddo
!!  ��ޑ��
	dF(:,:,MY)=Y(:,:,MY)/LUcom(0,MY)
	dF(:,:,MY-1)=(Y(:,:,MY-1)-dF(:,:,MY)*LUcom(1,MY-1))/LUcom(0,MY-1)
	do j=MY-2,0,-1
		dF(:,:,j)=(Y(:,:,j)-LUcom(1,j)*dF(:,:,j+1)-LUcom(2,j)*dF(:,:,j+2))/LUcom(0,j)
	enddo

end subroutine RHS2y


!!����f���͋y�і��C�W��Cf----------------------------------------------------------------------------------
subroutine cf_tauw(MX,MY,MZ,x,Re,a,b,c,ad,bd,cd,omega,alpha,beta,alphaD,betaD,alpha3,alphaD3,ymax,ymin,ysy&
                  &,rho_read,u_read,frhoup,tauw,cf,count3)
implicit none
integer :: MX,MY,MZ
double precision :: ad,bd,cd,a,b,c,omega,a3,ad3
double precision :: alpha,beta,alphaD,betaD,alpha3,alphaD3
double precision :: ymax,ymin
double precision :: ysy(0:MY),rho_read(0:MX,0:MY),u_read(0:MX,0:MY),frhoup(1:6,0:MX,0:MY),x(0:MX)
double precision :: tauw(0:MX),cf(0:MX)
double precision :: AAuave_y(0:MY,0:MY),LUuave_y(0:MY,0:MY)
double precision :: QRuave_y(1:6,0:MX,0:MY),guave_y(1:6,0:MX,0:MY)
double precision :: Re
double precision :: mu1
double precision :: Rex(0:MX),Cftheo(0:MX)
character(len=16) :: Chara
integer :: i,j,o,count3

mu1=1.d0/Re

call Ax(MY,AAuave_y,alpha,beta,alphaD,betaD,0d0,alpha3,alphaD3)
call LU(MY,AAuave_y,LUuave_y)
call RHS2y(MX,MY,MZ,6,a,b,c,ad,bd,cd,0d0,ymax,ymin,QRuave_y,frhoup,a3,ad3,LUuave_y,guave_y)
do o=1,6
	do i=0,MX
	guave_y(o,i,:)=guave_y(o,i,:)*ysy(:)    !�i�q�L����K�p
	enddo
enddo

tauw=mu1*guave_y(2,:,0)       !!��*du/dy


cf(:)=tauw(:)/(0.5d0*rho_read(:,0)*(u_read(:,0)**2.d0))

do i=0,MX
Rex(i)=(Re**2.d0)/(1.721**2.d0)+Re*x(i)
enddo

Cftheo(:)=0.455d0/((dlog10(Rex(:)))**2.58d0)    !!���_��

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

