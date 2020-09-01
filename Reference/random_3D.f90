!森山が作成したランダム撹乱作成コード(境界層用)
!を河島がジェット用に改変(2020/08/31)
program kakusan
  !$use omp_lib
  implicit none
  integer,parameter :: NX = 360
  integer,parameter :: NY = 200
  integer,parameter :: NZ = 20
  !x方向格子伸長
  double precision,parameter :: x_width=10.8d0
  double precision,parameter :: a1=1d0/14d0
  double precision,parameter :: a2=7.d0
  double precision,parameter :: b1=1.d0/1.4d0
  !y方向格子伸長
  double precision,parameter :: y_width=3.d0!!!以下はxと同じなのでコメントアウト!!!a1=1d0/14d0;a2=7d0;b1=1.d0/1.4d0
  double precision,parameter :: b = 1.d0!Jet半径は1で固定してしまう
  double precision,parameter :: Cx = 24.d0*b !x軸の幅の設定
  double precision,parameter :: Cy = 8.d0*b !y軸の幅の設定
  !Buffer領域の幅は常にWx<=Cxと計算領域よりも小さくてはならない
  double precision,parameter :: Wrx = 12.d0*b!Buffer領域x方向右側の幅
  double precision,parameter :: Wlx = Wrx!Buffer領域x方向左側の幅
  double precision,parameter :: Wry = 2.d0*b!Buffer領域y方向右側の幅
  double precision,parameter :: Wly = Wry!Buffer領域y方向左側の幅
  double precision,parameter :: Lx =  Cx+Wrx!+Wlx x方向の長さを定義.x軸左側にもbufferをかけるなら変更が必要
  double precision,parameter :: Ly = 2.d0*Cy+Wry+Wly!y方向の長さを定義 計算領域がy軸対称なのでCyは*2にしている
  double precision,parameter :: Lz = 1.d0
  double precision :: Ymin

  integer,parameter::Kmx=10,Kmy=10,Kmz=10!kx,ky,kz(打ち切り波数)
  double precision,parameter::PI=dacos(-1d0)
  double precision,parameter::dx=Lx/dble(NX)
  double precision,parameter::dy=Ly/dble(NY)
  double precision,parameter::dz=Lz/dble(NZ)
  double precision,parameter::Kmax=1d0

  integer::i,j,k
  integer::Kx,Ky,Kz
  double precision::x(0:NX),xs(0:NX),y(0:NY),ys(0:NY),z(0:NZ-1)
  double precision,dimension(0:NX,0:NY,0:NZ-1)::u_d3,v_d3,w_d3
  double precision::ES3d(Kmx,Kmy,Kmz)
  double precision::Rand_x,Rand_y,Rand_z
  double precision::theta_x,theta_y,theta_z
  double precision::max_E
  double precision::max_u,max_v,max_w
  double precision::abs_k
  double precision,dimension(0:NX,0:NY,0:NZ-1)::kakuran_u,kakuran_v,kakuran_w

!!!初期座標および格子伸長
  !x座標設定
  do i=0,NX
     x(i)=dx*dble(i)
  end do
  !x方向の格子伸長
  do i= 0,NX
   xs(i) = b1 * (1.7d0*x(i)-a1*&
    (-dlog(dcosh(a2*(x(i) - x_width))) + dlog(dcosh(a2*(x(i) + x_width)))))
  enddo

  !y座標設定
  Ymin = -(Ly/2.d0)
  do i=0,NY
     y(i)=Ymin + dy*dble(i)
  end do
  !y方向の格子伸長
  do i= 0,NY
    ys(i) = b1 * (1.7d0*y(i)-a1*&
  (-dlog(dcosh(a2*(y(i) - y_width))) + dlog(dcosh(a2*(y(i)+ y_width)))))
  enddo

  !z座標設定
  do i=0,NZ-1
     z(i)=dz*dble(i)
  end do

  !撹乱エネルギースペクトルE(k)の計算
  do kx=1,Kmx
     do ky=1,Kmy
        do kz=1,Kmz
          !k=(kx,ky,kz)T(転置)
          !|k|=(kx^2+ky^2+kz^2)
           abs_k=dsqrt(dble(Kx*Kx)+dble(Ky*Ky)+dble(Kz*Kz))
          !E(k)=(|k|/Kmax)^4 * exp(-2*(|k|/Kmax)^2)
           ES3d(Kx,Ky,Kz)=(abs_k/Kmax)**4*dexp(-2d0*(abs_k/Kmax)**2)
        end do
     end do
  end do

  max_E=maxval(ES3d)
  u_d3=0d0
  v_d3=0d0
  w_d3=0d0
  call random_seed()
       !$omp parallel sections
     !$omp section
    !ランダム撹乱のメイン部分
    do j=0,NY
       do i=0,NX
          do k=0,NZ-1
             do Kx=1,Kmx
                do ky=1,Kmy
                   do kz=1,Kmz
                     !0<=Rand_x<=1の範囲でランダムに値を与える
                      call random_number(Rand_x)
                      call random_number(Rand_y)
                      call random_number(Rand_z)
                      theta_x=Rand_x*Lx
                      theta_y=Rand_y*Ly
                      theta_z=Rand_z*Lz
                      !撹乱の関数f(x,y,z)
                      !xとyの座標はxs,ysと格子伸長適用済みのものに変更
                      u_d3(i,j,k)=u_d3(i,j,k)&
                       +dsqrt(ES3D(Kx,Ky,Kz))*dsin(2d0*pi*Kx*(xs(i)+theta_x)/Lx)&
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
    write(*,*)'u方向撹乱計算完了'
     !$omp section
    !!!v方向
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
                       +dsqrt(ES3D(Kx,Ky,Kz))*dsin(2d0*pi*Kx*(xs(i)+theta_x)/Lx)&
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
    write(*,*)'v方向撹乱計算完了'
!$omp section
!w方向
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
                       +dsqrt(ES3D(Kx,Ky,Kz))*dsin(2d0*pi*Kx*(xs(i)+theta_x)/Lx)&
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
    write(*,*)'w方向撹乱計算完了'
     !$omp end parallel sections
    open(120,file='random_check.csv')
    !u方向の撹乱関数を確認(x=1)
    do i=0,NY
       do j=0,NZ-1
          write(120,*)i,j,u_d3(1,i,j)
       end do
    end do

 open(21,file='kakuran3D_u.txt',status='replace')
 open(22,file='kakuran3D_v.txt',status='replace')
 open(23,file='kakuran3D_w.txt',status='replace')


 do i=0,NX
    max_u=maxval(dabs(u_d3(i,:,:)))
    max_v=maxval(dabs(v_d3(i,:,:)))
    max_w=maxval(dabs(w_d3(i,:,:)))
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

 open(100,file='check.csv')
 !u方向の撹乱関数を確認
 do i=5,5
    do j=0,NY
       do k=0,NZ-1
          write(100,*)j,k,kakuran_u(i,j,k)
       end do
    end do
 end do
end program kakusan
