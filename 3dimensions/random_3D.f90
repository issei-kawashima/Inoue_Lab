!森山が作成したランダム撹乱作成コード(境界層用)
!を河島がジェット用に改変(2020/08/31)
program kakusan
  !$use omp_lib
  implicit none
  ! integer,parameter :: NX = 360!X=0で流入させるので、不要
  integer,parameter :: NY = 200
  integer,parameter :: NZ = 20
  !y方向格子伸長
  double precision,parameter :: y_width=3.d0
  double precision,parameter :: a1=1d0/14d0
  double precision,parameter :: a2=7.d0
  double precision,parameter :: b1=1.d0/1.4d0
  double precision,parameter :: b = 1.d0!Jet半径は1で固定してしまう
  double precision,parameter :: Cx = 24.d0*b !x軸の幅の設定
  double precision,parameter :: Cy = 8.d0*b !y軸の幅の設定
  !Buffer領域の幅は常にWx<=Cxと計算領域よりも小さくてはならない
  double precision,parameter :: Wrx = 12.d0*b!Buffer領域x方向右側の幅
  double precision,parameter :: Wlx = Wrx!Buffer領域x方向左側の幅
  double precision,parameter :: Wry = 2.d0*b!Buffer領域y方向右側の幅
  double precision,parameter :: Wly = Wry!Buffer領域y方向左側の幅
  double precision,parameter :: Lx =  Cx+Wrx!+Wlx x方向の長さを定義.x軸左側にもbufferをかけるなら変更が必要
  ! double precision,parameter :: Ly = 2.d0*Cy+Wry+Wly!y方向の長さを定義 計算領域がy軸対称なのでCyは*2にしている
  !撹乱を全て+にするために、y座標を絶対値で計算するようにした。
  !したがって、幅は半分になるので、Lyも半分になる
  double precision,parameter :: Ly = Cy+Wry
  double precision,parameter :: Lz = 1.d0
  double precision :: Ymin

  integer,parameter::Kmx=10,Kmy=10,Kmz=10!kx,ky,kz(打ち切り波数)
  double precision,parameter::pi= acos(-1d0)
  !Lyが半分なので、ここは2Lyになる
  double precision,parameter::dy= 2.d0*Ly/dble(NY)
  double precision,parameter::dz=Lz/dble(NZ)
  double precision,parameter::Kmax = dsqrt(3.d0)

  integer::i,j,k
  integer::Kx,Ky,Kz
  double precision,dimension(0:NY):: y,ys
  double precision,dimension(0:NZ-1) :: z
  double precision,dimension(0:NY,0:NZ-1)::u_d3,v_d3,w_d3
  double precision::ES3d(Kmx,Kmy,Kmz)
  double precision::Rand_x,Rand_y,Rand_z
  double precision::theta_x,theta_y,theta_z
  double precision::max_E
  double precision::max_u,max_v,max_w
  double precision::abs_k
  double precision,dimension(0:NY,0:NZ-1)::kakuran_u,kakuran_v,kakuran_w

!!!初期座標および格子伸長
  !x座標は流入部にしか撹乱を導入しないので、不要

  !y座標設定
  ! Ymin = -(Ly/2.d0)
  !Lyは半分なので、Yminも半分
  Ymin = -Ly
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
  open(120,file='dirturbance_conditions/E(k)_check.csv')
  !エネルギースペクトルを確認
  do Kz=1,Kmz
    do Ky=1,Kmy
       do Kx=1,Kmx
         abs_k=dsqrt(dble(Kx**2)+dble(Ky**2)+dble(Kz**2))
         write(120,*)abs_k,",",ES3d(Kx,Ky,Kz)/max_E
      enddo
     end do
  end do
  close(120)

  u_d3=0d0
  v_d3=0d0
  w_d3=0d0
  call random_seed()
!$omp parallel sections
  !$omp section
    !ランダム撹乱のメイン部分
    do k=0,NZ-1
      do j=0,NY
         do Kz=1,Kmz
            do Ky=1,Kmy
               do Kx=1,Kmx
                 !0<=Rand_x<=1の範囲でランダムに値を与える
                  call random_number(Rand_x)!波数が変わるたびに位相が変わるようにする
                  call random_number(Rand_y)
                  call random_number(Rand_z)
                  theta_x=Rand_x*Lx
                  theta_y=Rand_y*Ly
                  theta_z=Rand_z*Lz
                  !撹乱の関数f(x,y,z)
                  !xについては関数にはxの変動も入るが、使用するのはx=0の位置の撹乱だけなので、x=0の箇所しか計算しない
                  !yの座標はysと格子伸長適用済みのものに変更
                  u_d3(j,k)=u_d3(j,k)&
                   +dsqrt(ES3D(Kx,Ky,Kz))*dsin(2d0*pi*Kx*(0.d0+theta_x)/Lx)&
                       *dsin(2d0*pi*Ky*(ys(j)+theta_y)/Ly)&
                       *dsin(2d0*pi*Kz*(z(k)+theta_z)/Lz)
               end do
            end do
         end do
      end do
    end do
    u_d3(NY,:)=u_d3(0,:)
    u_d3(:,NZ-1)=u_d3(:,0)
    write(*,*)'u方向撹乱計算完了'
    !$omp section
    !!!v方向
    do k=0,NZ-1
      do j=0,NY
         do Kz=1,Kmz
            do Ky=1,Kmy
               do Kx=1,Kmx
                  call random_number(Rand_x)
                  call random_number(Rand_y)
                  call random_number(Rand_z)
                  theta_x=Rand_x*Lx
                  theta_y=Rand_y*Ly
                  theta_z=Rand_z*Lz
                  v_d3(j,k)=v_d3(j,k)&
                   +dsqrt(ES3D(Kx,Ky,Kz))*dsin(2d0*pi*Kx*(0.d0+theta_x)/Lx)&
                       *dsin(2d0*pi*Ky*(ys(j)+theta_y)/Ly)&
                       *dsin(2d0*pi*Kz*(z(k)+theta_z)/Lz)
               end do
            end do
         end do
      end do
    end do
    v_d3(NY,:)=v_d3(0,:)
    v_d3(:,NZ-1)=v_d3(:,0)
    write(*,*)'v方向撹乱計算完了'
    !$omp section
!w方向
    do k=0,NZ-1
      do j=0,NY
        do Kz=1,Kmz
          do Ky=1,Kmy
             do Kx=1,Kmx
              call random_number(Rand_x)
              call random_number(Rand_y)
              call random_number(Rand_z)
              theta_x=Rand_x*Lx
              theta_y=Rand_y*Ly
              theta_z=Rand_z*Lz
              w_d3(j,k)=w_d3(j,k)&
               +dsqrt(ES3D(Kx,Ky,Kz))*dsin(2d0*pi*Kx*(0.d0+theta_x)/Lx)&
                   *dsin(2d0*pi*Ky*(ys(j)+theta_y)/Ly)&
                   *dsin(2d0*pi*Kz*(z(k)+theta_z)/Lz)
             end do
          end do
        end do
      end do
    end do
    w_d3(NY,:)=w_d3(0,:)
    w_d3(:,NZ-1)=w_d3(:,0)
    write(*,*)'w方向撹乱計算完了'
!$omp end parallel sections

 open(21,file='dirturbance_conditions/kakuran3D_u.txt',status='replace')
 open(22,file='dirturbance_conditions/kakuran3D_v.txt',status='replace')
 open(23,file='dirturbance_conditions/kakuran3D_w.txt',status='replace')


 ! do k=0,NZ-1
   do j=0,Ny
    max_u=maxval(dabs(u_d3(j,:)))
    max_v=maxval(dabs(v_d3(j,:)))
    max_w=maxval(dabs(w_d3(j,:)))
    kakuran_u(j,:)=u_d3(j,:)/max_u
    kakuran_v(j,:)=v_d3(j,:)/max_v
    kakuran_w(j,:)=w_d3(j,:)/max_w
  enddo
 ! end do

 do k=0,NZ-1
   do j=0,NY
     write(21,*) kakuran_u(j,k)
     write(22,*) kakuran_v(j,k)
     write(23,*) kakuran_w(j,k)
   end do
 end do
 close(21)
 close(22)
 close(23)

 open(100,file='dirturbance_conditions/check.csv')
 !u方向の撹乱関数を確認
do k=0,NZ-1
   do j=0,NY
      write(100,*)j,",",k,",",kakuran_v(j,k)
   end do
end do
close(100)
end program kakusan
