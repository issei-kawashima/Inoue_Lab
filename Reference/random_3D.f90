!!!kakusan
!!!森山史孝　1024
program kakusan
  !$use omp_lib
  implicit none
  integer,parameter::NX=100,NY=70,NZ=50
  integer,parameter::Kmx=10,Kmy=10,Kmz=10
  double precision,parameter::PI=dacos(-1d0)
  double precision,parameter::Lx=150d0,Ly=8d0,Lz=12d0
  double precision,parameter::dx=Lx/dble(NX)
  double precision,parameter::dy=Ly/dble(NY)
  double precision,parameter::dz=Lz/dble(NZ)
  double precision,parameter::Kmax=1d0

  integer::i,j,k
  integer::Kx,Ky,Kz
  double precision::x(0:NX),y(0:NY),ys(0:NY),z(0:NZ-1)
  double precision,dimension(0:NY,0:NZ-1)::WF3d1,WF3d2,WF3d
  double precision,dimension(0:NX,0:NY,0:NZ-1)::u_d3,v_d3,w_d3
  double precision::ES3d(Kmx,Kmy,Kmz)
  double precision::Rand_x,Rand_y,Rand_z
  double precision::theta_x,theta_y,theta_z
  double precision::max_E
  double precision::max_u,max_v,max_w
  double precision::abs_k,NAN=1,ran=1,a1=0.13d0
  double precision,dimension(0:NX,0:NY,0:NZ-1)::kakuran_u,kakuran_v,kakuran_w

!!!初期座標および格子伸長
  do i=0,NX
     x(i)=dx*dble(i)
  end do
  do i=0,NY
     y(i)=dy*dble(i)
  end do
  ys=(Ly)*dexp(-a1*(Ly-y)) - (Ly - y)*dexp(-a1*(Ly))
  do i=0,NZ-1
     z(i)=dz*dble(i)
  end do
  WF3d2=0d0
  !ここって必要ないループなんじゃ？
 !  do k=0,NZ-1
 !     do j=0,NY
 ! !       WF3d1(j,k)=dtanh(19d0*ys(j))*dexp(-10.8d0*(ys(j)**(6d0)))
 ! !       wf3d2(j,k)=dtanh(19d0*(-ys(j)+2d0))*dexp(-10.8d0*((-ys(j)+2d0)**(6d0)))
 !     end do
 !  end do
 !  wf3d=wf3d1+wf3d2

  WF3d(0,:)=0d0
  WF3d(NY,:)=0d0
wf3d=1d0
!ここって必要ないループなんじゃ？
!   do j=0,NY
! !     WF3d(j,:)=(dcos(pi*ys(j)*0.5d0)**2)*dexp(-2d-5*ys(j)**6)
! !     wf3d(j,:)=1d0-ys(j)**10d0
!   end do

  open(11,file='WF3D.csv')
  do j=0,NY
     write(11,*) ys(j),',',WF3d(j,0),',',WF3d(j,NZ-1)
  end do
    do kx=1,Kmx
       do ky=1,Kmy
          do kz=1,Kmz
             abs_k=dsqrt(dble(Kx*Kx)+dble(Ky*Ky)+dble(Kz*Kz))
             ES3d(Kx,Ky,Kz)=(abs_k/Kmax)**4*dexp(-2d0*(abs_k/Kmax)**2)
          end do
       end do
    end do

  max_E=maxval(ES3d)
if(ran==1)then
  u_d3=0d0
  v_d3=0d0
  w_d3=0d0
if (Nan==1)then
  call random_seed()
       !$omp parallel sections
     !$omp section
    do j=0,NY
       do i=0,NX
          do k=0,NZ-1
             do Kx=1,Kmx
                do ky=1,Kmy
                   do kz=1,Kmz
                      call random_number(Rand_x)
                      call random_number(Rand_y)
                      call random_number(Rand_z)
                      theta_x=Rand_x*Lx
                      theta_y=Rand_y*Ly
                      theta_z=Rand_z*Lz
                      u_d3(i,j,k)=u_d3(i,j,k)&
                       +dsqrt(ES3D(Kx,Ky,Kz))*dsin(2d0*pi*Kx*(x(i)+theta_x)/Lx)&
                           *dsin(2d0*pi*Ky*(ys(j)+theta_y)/Ly)&
                           *dsin(2d0*pi*Kz*(z(k)+theta_z)/Lz)
                   end do
                end do
             end do
          end do
       end do
    end do
    u_d3(NX,:,:)=u_d3(0,:,:)
    u_d3(:,NY,:)=u_d3(:,0,:)
    u_d3(:,:,NZ-1)=u_d3(:,:,0)
    write(*,*)'step-1'
     !$omp section
    !!!            v
    do j=0,NY
       do i=0,NX
          do k=0,NZ-1
             do Kx=1,Kmx
                do ky=1,Kmy
                   do kz=1,Kmz
                      call random_number(Rand_x)
                      call random_number(Rand_y)
                      call random_number(Rand_z)
                      theta_x=Rand_x*Lx
                      theta_y=Rand_y*Ly
                      theta_z=Rand_z*Lz
                      v_d3(i,j,k)=v_d3(i,j,k)&
                       +dsqrt(ES3D(Kx,Ky,Kz))*dsin(2d0*pi*Kx*(x(i)+theta_x)/Lx)&
                           *dsin(2d0*pi*Ky*(ys(j)+theta_y)/Ly)&
                           *dsin(2d0*pi*Kz*(z(k)+theta_z)/Lz)
                   end do
                end do
             end do
          end do
       end do
    end do
    v_d3(NX,:,:)=v_d3(0,:,:)
    v_d3(:,NY,:)=v_d3(:,0,:)
    v_d3(:,:,NZ-1)=v_d3(:,:,0)
    write(*,*)'step-2'
!$omp section
!!!!!!!!!!!!!!!!!!!!!!                  w
    do j=0,NY
       do i=0,NX
          do k=0,NZ-1
             do Kx=1,Kmx
                do ky=1,Kmy
                   do kz=1,Kmz
                      call random_number(Rand_x)
                      call random_number(Rand_y)
                      call random_number(Rand_z)
                      theta_x=Rand_x*Lx
                      theta_y=Rand_y*Ly
                      theta_z=Rand_z*Lz
                      w_d3(i,j,k)=w_d3(i,j,k)&
                       +dsqrt(ES3D(Kx,Ky,Kz))*dsin(2d0*pi*Kx*(x(i)+theta_x)/Lx)&
                           *dsin(2d0*pi*Ky*(ys(j)+theta_y)/Ly)&
                           *dsin(2d0*pi*Kz*(z(k)+theta_z)/Lz)
                   end do
                end do
             end do
          end do
       end do
    end do
    w_d3(NX,:,:)=w_d3(0,:,:)
    w_d3(:,NY,:)=w_d3(:,0,:)
    w_d3(:,:,NZ-1)=w_d3(:,:,0)
    write(*,*)'step-3'
     !$omp end parallel sections
    open(120,file='randam_tasikame.csv')
    do i=0,NY
       do j=0,NZ-1
          write(120,*)i,j,u_d3(1,i,j)
       end do
    end do

 open(21,file='kakuran3D_u.txt',status='replace')
 open(22,file='kakuran3D_v.txt',status='replace')
 open(23,file='kakuran3D_w.txt',status='replace')


 do i=0,NX
    max_u=maxval(dabs(u_d3(i,:,:)*WF3D))
    max_v=maxval(dabs(v_d3(i,:,:)*WF3D))
    max_w=maxval(dabs(w_d3(i,:,:)*WF3D))
    kakuran_u(i,:,:)=u_d3(i,:,:)/max_u
    kakuran_v(i,:,:)=v_d3(i,:,:)/max_v
    kakuran_w(i,:,:)=w_d3(i,:,:)/max_w
 end do
      do i=0,NX
            do j=0,NY
                 do k=0,NZ-1
 write(21,*) kakuran_u(i,j,k)
 write(22,*) kakuran_v(i,j,k)
 write(23,*) kakuran_w(i,j,k)
end do
end do
end do
 close(21)
 close(22)
 close(23)

 open(100,file='tasikame.csv')
 do i=5,5
    do j=0,NY
       do k=0,NZ-1
          write(100,*)j,k,kakuran_u(i,j,k)!*WF3d(j,k)
       end do
    end do
 end do
end if
end if
end program kakusan
