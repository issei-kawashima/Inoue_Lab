!dx,dyの格子伸長の倍率は0.5倍/1.2倍に設定してある
!2019.05.18 Nx=180 Ny=100 dt=5.d-3で実行すれば卒論の格子伸長適用時の条件になる
!M=75000(T=150)まで計算をすることが亜音速では少なくともできた。(3dim.f90で)
!2020.06.03 ファイル出力形式を.dから.txtにした。これによってpara viewで可視化できるようになるし、gnuplotでも可視化できる。
!2020.06.04 Tjet=1.12Tempでは適性膨張ジェットなのでTjet=1.4*Tempへ変更した
!計算時間短縮のために、NSCBCのx=Nxの部分は計算しなくて済むようにNSCBC_xのsubroutineを_0と_Nxに分割した。
!Pr=0.71>Pr=1へ変更。
!Nx=180,Ny=100,Nz=20,dt=2.d-3で計算。
!2020.06.11 超音速のFirst Schockを捉えるために、グリッド数をあげる必要がある。
!しかし現状ではメモリ制限(コンパイラーのせい)で容量が超えてしまうので、allocateに配列を書き換える。
!Nx=360,Ny=299,Nz=20でも計算開始できた。したがって、allocateで本当にグリッド数の限界を突破した
!2020.06.13 Pr=1の理由を探す事に。今回は森山が現実的という0.71で計算する
!計算高速化のために以下のことをした。
!1.dif_zのsubrouineをmodを使用しない形に変更。
!2. M=の出力をif文とmodを使用しないで行うように変更
!3.計算時間の計測をやめた。(実用的な意味がないし、三日間とかになると現状では桁不足だから)
!2020.06.20 dFy,dFzの1/3の角処理も廃止。
!variable_settingのUVWT,dUVWTの(0,:,:,:)は不要だが、微分の際に形式があっていないと同じsubroutineを
!使用できないので、仕方なく今回は廃止を見逃す。将来的にdif_x,y,zを0:4ごとなどに縮小できたらUVWTの(0)は廃止可能
!2020.08.25 音響成分のdiv_uを実装
!2020.08.26 速度勾配テンソルの第二不変量Qを実装
!2020.09.01 ランダム撹乱を読み込み窓関数を適用し、Jetの流入条件と足し合わせるコードの追加を開始
!2020.09.06 窓関数には、y方向にランダム撹乱の窓関数はtop-hat型ジェットの関数をそのまま使用。
!これにより、ランダム撹乱は完全にジェットの中にのみ、存在することにした。
!2020.09.07 ランダム撹乱の定義を変更し、配列を縮小したので、それに合わせてコードを書き換えた
!ujetがLxに達するまでの時間=36秒になるまで徐々にランダム撹乱を強くするようにした
!inflow subrouitneを改変して、top-hat Jet+撹乱を流入させるようにした
!撹乱の強さはジェット中心速度ujetの5%とdis_strengthを設定
!2020.09.08 NSCBC_x_0を適用することにした。これは矩形ジェットのため、ジェットが流入しないz面での流入境界条件は引き続き必要であるかである。
!しかし矩形ジェット流入部にはdirichlet条件で、top-hatジェットと撹乱を入れているので、そこのQのみはinflow subroutineで上書きする方式になっている
!矩形ジェットの計算条件は、矩形ジェットをある一部分だけ切り取った平面ジェットで計算することを指すので、0/toは0/19になる
!=>つまり今まで考えていた矩形ジェットとは計算条件が違う
!dt=1.d-2にした。　これで計算時間が1/5になった=>計算回れば計算精度OK
!2020.09.09NSCBCの流入/流出条件を見直し、それに伴い撹乱uと流入条件も見直した
!超音速なので、x=Nxの境界条件では逆流が起きないとして、NSCBCの無反射流出条件で逆流を示すL1を=0にする
!また、流入条件も領域の中から外に逆流するものがないつまり、L1=0、またuに関してはランダム撹乱を入れない
!そうすることで、du/dtの計算がなくなるので、結果的に、x=0の超音速流入条件はFx(0)=0だけになる、
!y方向に関しては亜音速流出条件を適用し、NSCBCの角・縁処理は一旦やらないでおく。
!2020.09.10 格子数を変更する際にはNUxも変更しなくてはいけない。
!Nx=360でNUx=213とする。Nx=180ではNUx=90で良い。
!加えて、ランダム撹乱の生成の際に、格子点と格子伸長の関数でメインプログラムと同一のものを使用しているので、そちらも直さなければいけない


module three_omp
  !連続の式、Eulerの運動方程式、エネルギー方程式を並列に並べた行列Q,Fの設定等をする
  !これらの式をまとめて基礎式と呼ぶ
  implicit none
  !計算条件値の設定
  double precision,parameter :: gamma = 1.4d0
  integer,parameter :: t_end = 150 !時刻tの設定
  integer,parameter :: p_output = 10 !時間毎の局所圧力を出力させる際のステップ間隔
  integer,parameter :: Nx = 360
  integer,parameter :: Ny = 200
  integer,parameter :: Nz = 20
  double precision,parameter :: dt = 1.d-2
  integer,parameter :: NUx = 213!buffer_xのUxで流入側のUxを0にする座標(格子点番号)Nx=180ならNUx=90,Nx=360ならNUx=213
  integer,parameter :: Mmax = t_end / dt
  integer,parameter :: output_count = int(1.d0/dt)!出力ファイルを1sec間隔で出力するように設定
  double precision,parameter :: b = 1.d0!Jet半径は1で固定してしまう
  double precision,parameter :: Cx = 24.d0*b !x軸の幅の設定
  double precision,parameter :: Cy = 8.d0*b !y軸の幅の設定
  !Buffer領域の幅は常にWx<=Cxと計算領域よりも小さくてはならない
  double precision,parameter :: Wrx = 12.d0*b!Buffer領域x方向右側の幅
  double precision,parameter :: Wlx = Wrx!Buffer領域x方向左側の幅
  double precision,parameter :: Wry = 2.d0*b!Buffer領域y方向右側の幅
  double precision,parameter :: Wly = Wry!Buffer領域y方向左側の幅
  double precision,parameter :: Lx =  Cx+Wrx!x方向の長さを定義=>流入部にはBufferを入れてはいけないので、Wlxは不要
  double precision,parameter :: Ly = 2.d0*Cy+Wry+Wly!y方向の長さを定義 計算領域がy軸対称なのでCyは*2にしている
  double precision,parameter :: Lz = 1.d0

  double precision,parameter :: psigma = -0.25d0
  double precision,parameter :: msigma = 0.25d0
  double precision,parameter :: NS_sigma = 0.25d0
  double precision,parameter :: ccs_sigma = 0.d0
  double precision,parameter :: c = 1.d0
  double precision,parameter :: Pr = 0.71d0
  double precision,parameter :: Ma = 1.6d0
  double precision,parameter :: Temp = 1.d0
  double precision,parameter :: Tjet = 1.4d0*Temp
  double precision,parameter :: ujet = 1.d0
  double precision,parameter :: dis_strength = 5.d-2*ujet!ジェット中心速度の5%撹乱
  integer,parameter :: times = int((Lx/ujet)/dt)!流入撹乱の時間変動基準(timesを超えたらフルパワー)
  integer,parameter :: observe_start_time = int(50.d0/dt)!ランダム撹乱で乱流化したかどうかを時間変動で、集計する開始時刻
  integer,parameter :: observe_end_time = int(55.d0/dt)!ランダム撹乱で乱流化したかどうかを時間変動で、集計する終了時刻
  double precision,parameter :: Sc = 120.d0 / (273.15d0 + 18.d0)
  double precision,parameter :: zeta = 1.d0
  double precision,parameter :: pi = acos(-1.d0)
  double precision,parameter :: Re = 1.d3
  double precision,parameter :: A2 = ujet*5.d-2!v方向の撹乱の振幅
  double precision,parameter :: T2 = 1.d0!v方向の撹乱の周期 T2=1で1秒に1周期
contains
  !初期条件G(rho,u,p)を用いてQ行列の設定
  subroutine Q_matrix(G,Q)
    integer i,j,k
    double precision,allocatable,dimension(:,:,:,:) :: Q,G
    !$omp parallel do
      do k=0,Nz-1
         do i=0,Ny
           do j=0,Nx
            Q(0,j,i,k) = G(0,j,i,k)
            Q(1,j,i,k) = G(0,j,i,k) * G(1,j,i,k)
            Q(2,j,i,k) = G(0,j,i,k) * G(2,j,i,k)
            Q(3,j,i,k) = G(0,j,i,k) * G(3,j,i,k)
            Q(4,j,i,k) = G(4,j,i,k) / (gamma -1.d0) + G(0,j,i,k) * (G(1,j,i,k)**2.d0+&
                        G(2,j,i,k)**2.d0 +G(3,j,i,k)**2.d0) * 0.5d0
          end do
        end do
      enddo
    !$omp end parallel do
  endsubroutine Q_matrix
  !求めQを用いてF行列の設定
  subroutine F_matrix(Q,Fpx,Fmx,Fpy,Fmy,Fpz,Fmz)
    double precision,allocatable,dimension(:,:,:,:) :: Q,Fpx,Fmx,Fpy,Fmy,Fpz,Fmz
    double precision,allocatable,dimension(:,:,:,:)::Fx,Fy,Fz
    integer i,j,k,l
    allocate(Fx(0:4,0:Nx,0:Ny,0:Nz-1),Fy(0:4,0:Nx,0:Ny,0:Nz-1),Fz(0:4,0:Nx,0:Ny,0:Nz-1))
    Fx=0.d0
    Fy=0.d0
    Fz=0.d0
    !$omp parallel do
      do k=0,Nz-1
         do i=0,Ny
           do j=0,Nx
       !F行列の設定(x方向)
        Fx(0,j,i,k) = Q(1,j,i,k)
        Fx(1,j,i,k) = 1.d0/(2.d0*Q(0,j,i,k))*((3.d0-gamma)*(Q(1,j,i,k)**2.d0)+&
                    (1.d0-gamma)*(Q(2,j,i,k)**2.d0 + Q(3,j,i,k)**2.d0))+(gamma-1.d0)*Q(4,j,i,k)
        Fx(2,j,i,k) = Q(1,j,i,k)*Q(2,j,i,k)/Q(0,j,i,k)
        Fx(3,j,i,k) = Q(1,j,i,k)*Q(3,j,i,k)/Q(0,j,i,k)
        Fx(4,j,i,k) = gamma * Q(1,j,i,k) * Q(4,j,i,k) / Q(0,j,i,k) &
                   &+ (1.d0- gamma)*Q(1,j,i,k)*(Q(1,j,i,k)**2.d0+Q(2,j,i,k)**2.d0&
                   +Q(3,j,i,k)**2.d0)/(2.d0 * (Q(0,j,i,k)**2.d0))
         !F行列の設定(y方向)
         Fy(0,j,i,k) = Q(2,j,i,k)
         Fy(1,j,i,k) = Q(2,j,i,k)*Q(1,j,i,k)/Q(0,j,i,k)
         Fy(2,j,i,k) = 1.d0/(2.d0*Q(0,j,i,k))*((3.d0 - gamma)*(Q(2,j,i,k)**2.d0) + &
                       &(1.d0 - gamma)*(Q(1,j,i,k)**2.d0+Q(3,j,i,k)**2.d0))+(gamma - 1.d0)*Q(4,j,i,k)
         Fy(3,j,i,k) = Q(2,j,i,k)*Q(3,j,i,k)/Q(0,j,i,k)
         Fy(4,j,i,k) = gamma * Q(2,j,i,k) * Q(4,j,i,k) / Q(0,j,i,k) &
                   &+ (1.d0- gamma)*Q(2,j,i,k)*(Q(1,j,i,k)**2.d0+Q(2,j,i,k)**2.d0&
                   +Q(3,j,i,k)**2.d0)/(2.d0 * (Q(0,j,i,k)**2.d0))
         !F行列の設定(z方向)
         Fz(0,j,i,k) = Q(3,j,i,k)
         Fz(1,j,i,k) = Q(3,j,i,k)*Q(1,j,i,k)/Q(0,j,i,k)
         Fz(2,j,i,k) = Q(3,j,i,k)*Q(2,j,i,k)/Q(0,j,i,k)
         Fz(3,j,i,k) = 1.d0/(2.d0*Q(0,j,i,k))*((3.d0 - gamma)*(Q(3,j,i,k)**2.d0) + &
         &(1.d0 - gamma)*(Q(1,j,i,k)**2.d0+Q(2,j,i,k)**2.d0))+(gamma - 1.d0)*Q(4,j,i,k)
         Fz(4,j,i,k) = gamma * Q(3,j,i,k) * Q(4,j,i,k) / Q(0,j,i,k) &
                   &+ (1.d0- gamma)*Q(3,j,i,k)*(Q(1,j,i,k)**2.d0+Q(2,j,i,k)**2.d0&
                   +Q(3,j,i,k)**2.d0)/(2.d0 * (Q(0,j,i,k)**2.d0))
           end do
         end do
       enddo
   !$omp end parallel do

   !=======Fp,FmはFを使用するので、分割して並列化======================================
    !求めたFを特製速度の正負によって分割する
    !Lax-Friedrichの流速分割を用いる
    !$omp parallel do
      do k=0,Nz-1
           do i=0,Ny
             do j=0,Nx
               do l=0,4
                Fpx(l,j,i,k) = (Fx(l,j,i,k) + zeta * Q(l,j,i,k)) * 0.5d0
                Fmx(l,j,i,k) = (Fx(l,j,i,k) - zeta * Q(l,j,i,k)) * 0.5d0
                Fpy(l,j,i,k) = (Fy(l,j,i,k) + zeta * Q(l,j,i,k)) * 0.5d0
                Fmy(l,j,i,k) = (Fy(l,j,i,k) - zeta * Q(l,j,i,k)) * 0.5d0
                Fpz(l,j,i,k) = (Fz(l,j,i,k) + zeta * Q(l,j,i,k)) * 0.5d0
                Fmz(l,j,i,k) = (Fz(l,j,i,k) - zeta * Q(l,j,i,k)) * 0.5d0
            end do
          end do
        enddo
      enddo
    !$omp end parallel do
        deallocate(Fx,Fy,Fz)
      endsubroutine F_matrix
!du/dx,dT.dxを求めるためにまずはuとTの設定。
!UVWT,myuの計算はrho,u,pをそのまま代入しても良いがその場合は求めたtQから毎回rho_u_pのsubroutineを
!呼び出して計算しなければいけないので今回はQから直接計算できるようにプログラムを組んだ
    subroutine variable_setting(UVWT,Q,myu)
      !V行列を設定する際にμの計算が複雑になっているのでそれを簡略に示すために別でμを計算するsubroutine
      double precision,allocatable,dimension(:,:,:,:) :: UVWT,Q
      double precision,allocatable,dimension(:,:,:) :: myu
      integer i,j,k
      !$omp parallel do
        do k=0,Nz-1
           do i=0,Ny
             do j=0,Nx
               !UVWTはcall前に毎回0クリアされているのでUVWT0=0の代入は二度手間
        ! UVWT(0,j,i,k) = 0.d0
        UVWT(1,j,i,k) = Q(1,j,i,k) / Q(0,j,i,k)!u
        UVWT(2,j,i,k) = Q(2,j,i,k) / Q(0,j,i,k)!v
        UVWT(3,j,i,k) = Q(3,j,i,k) / Q(0,j,i,k)!w
        UVWT(4,j,i,k) = (Ma**2.d0)*(gamma * (gamma - 1.d0) * (Q(4,j,i,k) - (((Q(1,j,i,k) **2.d0)+ &
        &(Q(2,j,i,k) **2.d0)+(Q(3,j,i,k)**2.d0)) / (2.d0 * Q(0,j,i,k))))) / Q(0,j,i,k)!T
        !Tの値はMa^2*gamma*p/rhoこれは速度で無次元化したもの
            end do
          end do
        enddo
      !$omp end parallel do
!!==UVWT(4)を使用しているし、Doループが大きいので、並列化にあたり分けた。 UVWT(4)を埋め込むと計算量が増えるので、行わない===
      !$omp parallel do
        do k=0,Nz-1
           do i=0,Ny
             do j=0,Nx
            !UVWTからTの値を代入することで計算を簡略化している
            myu(j,i,k) = (UVWT(4,j,i,k) ** 1.5d0) * (1.d0 + Sc) / (UVWT(4,j,i,k) + Sc)
            end do
          end do
        enddo
      !$omp end parallel do
    end subroutine variable_setting

    subroutine V_matrix(Vx,Vy,Vz,myu,UVWT,dUVWTx,dUVWTy,dUVWTz)
      !基礎式右辺をV行列として設定
      double precision,allocatable,dimension(:,:,:,:) :: Vx,Vy,Vz
      double precision,allocatable,dimension(:,:,:,:) :: dUVWTx,dUVWTy,dUVWTz,UVWT
      double precision,allocatable,dimension(:,:,:) :: myu
      integer i,j,k
      !$omp parallel do
        do k=0,Nz-1
           do i=0,Ny
             do j=0,Nx
               !粘性項の設定(x方向)
               !Vx,y,zはcallする前に0クリアされているので、V0に0を代入する作業は二度手間
               !したがってコメントアウトしてしまう(削除するとわからなくなるから消さない)
        ! Vx(0,j,i,k) = 0.d0
        Vx(1,j,i,k) = (2.d0*myu(j,i,k)/(3.d0*Re)) * (2.d0 * dUVWTx(1,j,i,k) - dUVWTy(2,j,i,k)-dUVWTz(3,j,i,k))
        Vx(2,j,i,k) = (myu(j,i,k) / Re) * (dUVWTx(2,j,i,k) + dUVWTy(1,j,i,k))
        Vx(3,j,i,k) = (myu(j,i,k) / Re) * (dUVWTx(3,j,i,k) + dUVWTz(1,j,i,k))

                !粘性項の設定(y方向)
        ! Vy(0,j,i,k) = 0.d0
        Vy(1,j,i,k) = (myu(j,i,k) / Re) * (dUVWTx(2,j,i,k) + dUVWTy(1,j,i,k))
        Vy(2,j,i,k) = (2.d0*myu(j,i,k)/(3.d0*Re))*(-dUVWTx(1,j,i,k) + 2.d0*dUVWTy(2,j,i,k)-dUVWTz(3,j,i,k))
        Vy(3,j,i,k) = (myu(j,i,k) / Re) * (dUVWTy(3,j,i,k) + dUVWTz(2,j,i,k))

                !粘性項の設定(z方向)
        ! Vz(0,j,i,k) = 0.d0
        Vz(1,j,i,k) = (myu(j,i,k) / Re) * (dUVWTx(3,j,i,k) + dUVWTz(1,j,i,k))
        Vz(2,j,i,k) = (myu(j,i,k) / Re) * (dUVWTy(3,j,i,k) + dUVWTz(2,j,i,k))
        Vz(3,j,i,k) = (2.d0*myu(j,i,k)/(3.d0*Re))*(-dUVWTx(1,j,i,k) - dUVWTy(2,j,i,k)+2.d0*dUVWTz(3,j,i,k))
            end do
          end do
        enddo
      !$omp end parallel do
      !====V4だけが唯一V1~V3を必要とし、並列化不可能&Doループが大きいので分割して並列化する====
      !$omp parallel do
        do k=0,Nz-1
           do i=0,Ny
             do j=0,Nx
               !V4_x
               Vx(4,j,i,k) =Vx(1,j,i,k)*UVWT(1,j,i,k) + Vx(2,j,i,k)* UVWT(2,j,i,k)+ &
                             Vx(3,j,i,k)*UVWT(3,j,i,k)+((myu(j,i,k) * dUVWTx(4,j,i,k))&
                             / ((gamma - 1.d0)*Re*Pr*(Ma ** 2.d0)))
              !V4_y
               Vy(4,j,i,k) = Vy(1,j,i,k)*UVWT(1,j,i,k) + Vy(2,j,i,k)* UVWT(2,j,i,k)+ &
                             Vy(3,j,i,k)*UVWT(3,j,i,k)+((myu(j,i,k) * dUVWTy(4,j,i,k))&
                             / ((gamma - 1.d0)*Re*Pr*(Ma ** 2.d0)))
              !V4_z
               Vz(4,j,i,k) = Vz(1,j,i,k)*UVWT(1,j,i,k) + Vz(2,j,i,k)* UVWT(2,j,i,k)+ &
                             Vz(3,j,i,k)*UVWT(3,j,i,k)+((myu(j,i,k) * dUVWTz(4,j,i,k))&
                             / ((gamma - 1.d0)*Re*Pr*(Ma ** 2.d0)))
            end do
          end do
        enddo
      !$omp end parallel do
    endsubroutine V_matrix
    !求めたQnから時間毎のrho,u,v,w,pを求めるサブルーチン
    !ある指定した時間の時の計算結果のみを取り出せばグラフが作成できるので毎回は使用しない
    subroutine rho_u_p(G,Qn)
      double precision,allocatable,dimension(:,:,:,:) :: Qn,G
      integer i,j,k
        G=0.d0
      !$omp parallel do
        do k=0,Nz-1
           do i=0,Ny
             do j=0,Nx
        G(0,j,i,k) = Qn(0,j,i,k)!ρ
        G(1,j,i,k) = Qn(1,j,i,k) / Qn(0,j,i,k)!u
        G(2,j,i,k) = Qn(2,j,i,k) / Qn(0,j,i,k)!v
        G(3,j,i,k) = Qn(3,j,i,k) / Qn(0,j,i,k)!w
        G(4,j,i,k) = (gamma-1.d0)*(Qn(4,j,i,k)-(Qn(1,j,i,k)**2.d0+&
                      Qn(2,j,i,k)**2.d0+Qn(3,j,i,k)**2.d0)/(2.d0*Qn(0,j,i,k)))!p
            end do
          end do
        enddo
      !$omp end parallel do
    endsubroutine rho_u_p
    !DCS用の行列Aの設定(左辺の設定)
    !５次精度DCSと6次精度CCSの両方に対応
    subroutine LU_DecompoNonP(N,sigma,LU)
      integer i,N,j,k
      double precision sigma
      double precision Usum,Lsum
      double precision,allocatable,dimension(:,:) :: A,L,U,LU
      double precision,parameter :: alpha3 = 0.25d0
      double precision,parameter :: alpha5 = 1.d0/3.d0
      double precision,parameter :: Dalpha = 1.d0  !5次精度のDCSとなるための係数設定
      allocate(A(0:N,0:N),L(0:N,0:N),U(0:N,0:N))
      A=0.d0;L=0.d0;U=0.d0!Nx,Nyで同じsubroutineを用いているため0クリアが必要
      LU=0.d0
      !非周期条件よりi=0,Nでは片側DCS、i=1,N-1では3次精度DCSを用いる
      !その他の箇所は5次精度DCSを使う(sigmaを0にすれば6次CCSとなる)
      !0行目 片側DCS
      A(0,0) = 1.d0
      A(0,1) = 3.d0
      !1行目 3次精度DCS
      A(1,0) = alpha3 * (1.d0 - Dalpha * sigma)
      A(1,1) = 1.d0
      A(1,2) = alpha3 * (1.d0 + Dalpha * sigma)

      !2からN-2行目まで
      !$omp parallel do
      do i = 2,N-2
        A(i,i-1) = (1.d0 - Dalpha * sigma) * alpha5
        A(i,i) = 1.d0
        A(i,i+1) = (1.d0 + Dalpha * sigma) * alpha5
      enddo
      !$omp end parallel do

      !N-1行目 3次精度DCS
      A(N-1,N-2) = alpha3 * (1.d0 - Dalpha * sigma)
      A(N-1,N-1) = 1.d0
      A(N-1,N) = alpha3 * (1.d0 + Dalpha * sigma)
      !N行目 片側DCS
      A(N,N-1) = 3.d0
      A(N,N) = 1.d0
    !LU分解
      !まずL,Uの初期値を設定
      !===========並列化しない=====================================================
      !L(i,0)を求める際のU(0,0)だが、doループ内で求める値なので、並列化するとNaNになってしまう
      !U(0,i)とL(i,0),L(i,i)でDoループを分けることも可能だが、並列化の利点よりもループを増やす欠点の方が
      !大きいと判断し、ここでは並列化は行わない。
      do i = 0,N
        U(0,i) = A(0,i)
        L(i,0) = A(i,0) / U(0,0)
        L(i,i) = 1.d0
      enddo
      !===========並列化しない=====================================================

      !========並列化不可能======================================================
      !Uは行ごとに、Lは列ごとに求めていく。
      !ただしUの２行目、 Lの２列目、Uの3行目、Lの３列目といった順番
      do i = 1,N !初期条件の結果を利用してUのi列、Lのi行の順に求めていく
        do j = i,N !Uのi行の列要素を求めていく。i<=jを考慮する
          Usum = 0.d0
          do k = 0, i-1
            Usum = Usum + L(i,k) * U(k,j)
          enddo
          U(i,j) = A(i,j) - Usum
        enddo
        !U(i,j)はi=0,Nでj=i,Nで回しながら行要素を求めるが
        !Lはj=0,Nでi=j+1,Nで回しながら列要素を求めなければならない
        !しかしそれではdoループをi=0,Nとj=0,Nで分けなければならない
        !LはL(j,i)として考えるとi=0,Nの一つのdo文でU,L両方を定義できる
        !その結果Uの２行目、 Lの２列目、Uの3行目、Lの３列目といった順に求めるコードが書ける
        !そうでなければNの数だけコードを書かなければならず非現実的になってしまう
        !それがCroutのアルゴリズム
        do j = i+1,N !Lのi列の行要素を求めていく。j<iを考慮
          Lsum = 0.d0
          do k = 0,i-1
            Lsum = Lsum + L(j,k) * U(k,i)
          enddo
          L(j,i) = (A(j,i) - Lsum) / U(i,i)
        enddo
      enddo
      !========並列化不可能=======================================================
!LU行列の圧縮　L,U行列の対角成分のみをLU行列に保存する
!LU(-1,i)にL行列の対角成分L(i,i-1)を保存
!L(i,i)=1.d0は全身代入では使用しないので圧縮した行列には含めない
!LU(0,i)にU_i,iの成分を保存
!LU(1,i)にU_i,i+1の成分を保存
!これは6次制度CCSと5次制度DCSにのみ対応しているのでそれ以上の精度で計算する場合には変更が必要
      LU(-1,0) = 0.d0
      LU(0,N) = U(N,N)
      LU(1,N) = 0.d0
      !$omp parallel do
      do i = 1,N
        LU(-1,i) = L(i,i-1)
      enddo
      !$omp end parallel do

      !$omp parallel do
      do j = 0,1
        do i = 0,N-1
          LU(j,i) = U(i,i+j)
        enddo
      enddo
      !$omp end parallel do
      deallocate(A,L,U)
    endsubroutine LU_DecompoNonP

    !z方向用の周期条件で計算するA,L,Uの設定subroutine
    !DCS用の行列Aの設定(左辺の設定)
    !５次精度DCSと6次精度CCSの両方に対応
    subroutine LU_DecompoPiriodic(N,sigma,LU)
      integer i,N,j,k
      double precision sigma
      double precision Usum,Lsum
      double precision,allocatable,dimension(:,:) :: A,L,U,LU
      double precision,parameter :: alpha5 = 1.d0/3.d0
      double precision,parameter :: Dalpha = 1.d0  !5次精度のDCSとなるための係数設定
      allocate(A(0:N-1,0:N-1),L(0:N-1,0:N-1),U(0:N-1,0:N-1))
      A=0.d0;L=0.d0;U=0.d0!Nx,Nyで同じsubroutineを用いているため0クリアが必要
      LU=0.d0
      !周期条件なので5次精度DCSを使う(sigmaを0にすれば6次CCSとなる)
      !0行目
      A(0,0) = 1.d0
      A(0,1) = alpha5*(1.d0+sigma*Dalpha)
      A(0,N-1) = alpha5*(1.d0-sigma*Dalpha)

      !$omp parallel do
      !1からN-2行目まで
      do i = 1,N-2
        A(i,i-1) = (1.d0 - Dalpha * sigma) * alpha5
        A(i,i) = 1.d0
        A(i,i+1) = (1.d0 + Dalpha * sigma) * alpha5
      enddo
      !$omp end parallel do

      !N行目
      A(N-1,0) = alpha5 * (1.d0 + Dalpha * sigma)
      A(N-1,N-2) = alpha5*(1.d0-sigma*Dalpha)
      A(N-1,N-1) = 1.d0
    !LU分解
      !まずL,Uの初期値を設定
      !===============L(i,0)はU(0,0)を必要とするので、並列化しない(Doループが小さいから)====
      do i = 0,N-1
        U(0,i) = A(0,i)
        L(i,0) = A(i,0) / U(0,0)
        L(i,i) = 1.d0
      enddo
      !===============並列化しない================================================

      !========並列化不可能=======================================================
      !Uは行ごとに、Lは列ごとに求めていく。
      !ただしUの２行目、 Lの２列目、Uの3行目、Lの３列目といった順番
      do i = 1,N-1 !初期条件の結果を利用してUのi列、Lのi行の順に求めていく
        do j = i,N-1 !Uのi行の列要素を求めていく。i<=jを考慮する
          Usum = 0.d0
          do k = 0, i-1
            Usum = Usum + L(i,k) * U(k,j)
          enddo
          U(i,j) = A(i,j) - Usum
        enddo
        !U(i,j)はi=0,Nでj=i,Nで回しながら行要素を求めるが
        !Lはj=0,Nでi=j+1,Nで回しながら列要素を求めなければならない
        !しかしそれではdoループをi=0,Nとj=0,Nで分けなければならない
        !LはL(j,i)として考えるとi=0,Nの一つのdo文でU,L両方を定義できる
        !その結果Uの２行目、 Lの２列目、Uの3行目、Lの３列目といった順に求めるコードが書ける
        !そうでなければNの数だけコードを書かなければならず非現実的になってしまう
        !それがCroutのアルゴリズム
        do j = i+1,N-1 !Lのi列の行要素を求めていく。j<iを考慮
          Lsum = 0.d0
          do k = 0,i-1
            Lsum = Lsum + L(j,k) * U(k,i)
          enddo
          L(j,i) = (A(j,i) - Lsum) / U(i,i)
        enddo
      enddo
      !========並列化不可能=======================================================
!===============================================
!LU圧縮は周期なので-(2:2,0:N-1)の幅でLUをとればできる
!===============================================
!LU行列の圧縮　L,U行列の対角成分のみをLU行列に保存する
!LU(-1,i)にL行列の対角成分L(i,i-1)を保存
!L(i,i)=1.d0は全身代入では使用しないので圧縮した行列には含めない
!LU(0,i)にU_i,iの成分を保存
!LU(1,i)にU_i,i+1の成分を保存
!これは6次制度CCSと5次制度DCSにのみ対応しているのでそれ以上の精度で計算する場合には変更が必要
      LU(-2,N-2:N-1) = 0.d0
      LU(-1,0) = 0.d0
      LU(0,N-1) = U(N-1,N-1)
      !LU(-2),LU(2)は対角成分と被りがあるのでその部分は圧縮行列には入れずに0とする
      LU(1,N-1) = 0.d0
      LU(2,N-2:N-1) = 0.d0

      !$omp parallel do
      do i= 0,N-3
        LU(-2,i) = U(i,N-1)
        LU(2,i)  = L(N-1,i)
      enddo
      !$omp end parallel do

      !$omp parallel do
      do i = 1,N-1
        LU(-1,i) = L(i,i-1)
      enddo
      !$omp end parallel do

      !$omp parallel do
      do j = 0,1
        do i = 0,N-2
          LU(j,i) = U(i,i+j)
        enddo
      enddo
      !$omp end parallel do
      deallocate(A,L,U)
    endsubroutine LU_DecompoPiriodic

  !DCS右辺の計算(RHS)サブルーチン
  !x方向
    subroutine dif_x(sigma,dx,Fx,dFzeta,LU,dzeta_inx)
      integer i,j,k,l
      double precision,allocatable,dimension(:,:,:,:) :: Fx,dFzeta
      double precision,allocatable,dimension(:) :: dzeta_inx
      double precision,allocatable,dimension(:,:,:,:):: D2,D4,D6,D8
      double precision,allocatable,dimension(:,:,:,:):: x,y,RHS_x
      double precision dx,sigma,dxinv
      double  precision,parameter :: ra = 14.d0/9.d0, rb = 1.d0/9.d0&
      &,da = 4.d0 / 9.d0,db = 2.d0 / 9.d0 !5次精度のDCSとなるための係数設定
      double precision,allocatable,dimension(:,:) :: LU

      allocate(D2(0:4,2:Nx-2,0:Ny,0:Nz-1),D4(0:4,2:Nx-2,0:Ny,0:Nz-1),&
      D6(0:4,2:Nx-2,0:Ny,0:Nz-1),D8(0:4,2:Nx-2,0:Ny,0:Nz-1))
      allocate(x(0:4,0:Nx,0:Ny,0:Nz-1),y(0:4,0:Nx,0:Ny,0:Nz-1),&
      RHS_x(0:4,0:Nx,0:Ny,0:Nz-1))

      D2=0.d0;D4=0.d0;D6=0.d0;D8=0.d0;RHS_x=0.d0;y=0.d0;x=0.d0
      !片側DCS,3次精度DCSも入れた非周期条件の際のbの設定
      dxinv = 1.d0/dx

      !$omp parallel do
        do k=0,Nz-1
           do i=0,Ny
             do l=0,4
       !片側DCSの右辺設定
        RHS_x(l,0,i,k)=((-17.d0/6.d0)*Fx(l,0,i,k)+1.5d0*(Fx(l,1,i,k)+&
                        Fx(l,2,i,k))-Fx(l,3,i,k)/6.d0)*dxinv
        RHS_x(l,Nx,i,k)=((1.d0/6.d0)*Fx(l,Nx-3,i,k)-1.5d0*(Fx(l,Nx-2,i,k)+&
                        Fx(l,Nx-1,i,k))+(17.d0/6.d0)*Fx(l,Nx,i,k))*dxinv

        !3次精度DCSの右辺設定
        RHS_x(l,1,i,k)=((1.5d0)*(-Fx(l,0,i,k)+Fx(l,2,i,k))*0.5d0*dxinv)+sigma*&
                      ((Fx(l,0,i,k)-2.d0*Fx(l,1,i,k)+Fx(l,2,i,k))*(0.5d0*dxinv))
        RHS_x(l,Nx-1,i,k)=((1.5d0)*(-Fx(l,Nx-2,i,k)+Fx(l,Nx,i,k))*(0.5d0*dxinv))+&
                      sigma*((Fx(l,Nx-2,i,k)-2.d0*Fx(l,Nx-1,i,k)+Fx(l,Nx,i,k))*(0.5d0*dxinv))
              end do
          end do
        enddo
      !$omp end parallel do

         !5次精度DCSの右辺設定
      !$omp parallel do
      do k=0,Nz-1
       do i=0,Ny
        do j=2,Nx-2
          do l=0,4
         D2(l,j,i,k) = (-Fx(l,j-1,i,k)+Fx(l,j+1,i,k)) * (0.5d0*dxinv)
         D4(l,j,i,k) = (-Fx(l,j-2,i,k)+Fx(l,j+2,i,k)) * (0.25d0*dxinv)
         !D6(l,j,i,k) = (-Fx(l,j-3,i,k)+Fx(l,j+3,i,k)) / (6.d0*dx)
         !7次精度DCS用のため不要
         D6(l,j,i,k) = (Fx(l,j-1,i,k)+Fx(l,j+1,i,k)- 2.d0* Fx(l,j,i,k)) * dxinv
         D8(l,j,i,k) = (Fx(l,j-2,i,k)+Fx(l,j+2,i,k)- 2.d0* Fx(l,j,i,k)) * (0.25d0*dxinv)
         ! D12(l,j,i,k) = (Fl(:,j-i,k,:)+Fl(:,j+i,k,:)- 2.d0* Fx(l,j,i,k)) / (9.d0*dx)
         !7次精度DCS用のため不要
           end do
         end do
        end do
       enddo
     !$omp end parallel do
     !==========RHSにはD2~D8が必要で、Doループが大きいので、分割して並列化する===============
     !$omp parallel do
     do k=0,Nz-1
      do i=0,Ny
       do j=2,Nx-2
         do l=0,4
        RHS_x(l,j,i,k)=ra*D2(l,j,i,k)+rb*D4(l,j,i,k)+sigma*(da*D6(l,j,i,k)+db*D8(l,j,i,k))
         end do
       end do
      end do
     enddo
    !$omp end parallel do
     !前進代入法、後退代入法の計算サブルーチン(x方向)
          !前進代入
          !$omp parallel do
          do k=0,Nz-1
            do i=0,Ny
              do l=0,4
                y(l,0,i,k) = RHS_x(l,0,i,k)!例外の境界値
              end do
            end do
          end do
          !$omp end parallel do
    !=============前後のyを使用するので、並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝===＝＝＝======
            do j = 1,Nx!y(:,k)でk=0は上で定義したので残りの1〜Nxを定義する
              y(:,j,:,:) = RHS_x(:,j,:,:) - LU(-1,j)*y(:,j-1,:,:)
            enddo
    !=============並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝===＝＝＝======
            !後退代入
            !$omp parallel do
            do k=0,Nz-1
              do i=0,Ny
                do l=0,4
                  x(l,Nx,i,k) = y(l,Nx,i,k) / LU(0,Nx)!例外の境界値(Nx)
                end do
              end do
            end do
            !$omp end parallel do
            !x(:,j)でj=Nxは定義したので残りのj=0~Nx-1を定義する
    !=============前後のxを使用するので、並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝===＝＝＝======
            do j = Nx-1, 0, -1!後退するので-1ずつ進む
              x(:,j,:,:) = (y(:,j,:,:) - LU(1,j)*x(:,j+1,:,:)) / LU(0,j)
            enddo
     !============並列化不可能＝＝＝＝＝=＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝==＝＝＝=======

          call combine_x(dzeta_inx,x,dFzeta)
          deallocate(D2,D4,D6,D8,x,y,RHS_x)
      end subroutine dif_x
   !DCS右辺の計算(RHS)サブルーチン
   !y方向
     subroutine dif_y(sigma,dy,Fy,dFzeta,LU,dzeta_iny)
       integer i,j,k,l
       double precision,allocatable,dimension(:,:,:,:):: Fy,dFzeta
       double precision,allocatable,dimension(:):: dzeta_iny
       double precision,allocatable,dimension(:,:,:,:):: x,y,RHS_y
       double precision,allocatable,dimension(:,:,:,:):: D2,D4,D6,D8
       double precision dy,sigma,dyinv
       double  precision,parameter :: ra = 14.d0/9.d0, rb = 1.d0/9.d0&
       &,da = 4.d0 / 9.d0,db = 2.d0 / 9.d0 !5次精度のDCSとなるための係数設定
       double precision,allocatable,dimension(:,:) :: LU
       allocate(x(0:4,0:Nx,0:Ny,0:Nz-1),y(0:4,0:Nx,0:Ny,0:Nz-1),&
       RHS_y(0:4,0:Nx,0:Ny,0:Nz-1))
       allocate(D2(0:4,0:Nx,2:Ny-2,0:Nz-1),D4(0:4,0:Nx,2:Ny-2,0:Nz-1),&
       D6(0:4,0:Nx,2:Ny-2,0:Nz-1),D8(0:4,0:Nx,2:Ny-2,0:Nz-1))
       D2=0.d0;D4=0.d0;D6=0.d0;D8=0.d0;RHS_y=0.d0;y=0.d0;x=0.d0
       dyinv = 1.d0 / dy
       !片側DCS,3次精度DCSも入れた非周期条件の際のbの設定
       !$omp parallel do
         do k=0,Nz-1
          do j=0,Nx
            do l=0,4
              !片側DCSの右辺設定
         RHS_y(l,j,0,k) = ((-17.d0/6.d0)*Fy(l,j,0,k)+(1.5d0)*Fy(l,j,1,k)+&
         (1.5d0)*Fy(l,j,2,k)-Fy(l,j,3,k)/6.d0)*dyinv
         RHS_y(l,j,Ny,k)=((1.d0/6.d0)*Fy(l,j,Ny-3,k)-(1.5d0)&
         &*Fy(l,j,Ny-2,k)-(1.5d0)*Fy(l,j,Ny-1,k)+(17.d0/6.d0)*Fy(l,j,Ny,k))*dyinv

         !3次精度DCSの右辺設定
         RHS_y(l,j,1,k)=((1.5d0)*(-Fy(l,j,0,k)+Fy(l,j,2,k))*(0.5d0*dyinv))+sigma*&
         ((Fy(l,j,0,k)-2.d0*Fy(l,j,1,k)+Fy(l,j,2,k))*(0.5d0*dyinv))
         RHS_y(l,j,Ny-1,k)=((1.5d0)*(-Fy(l,j,Ny-2,k)+Fy(l,j,Ny,k))*(0.5d0*dyinv))&
            +sigma*((Fy(l,j,Ny-2,k)-2.d0*Fy(l,j,Ny-1,k)+Fy(l,j,Ny,k))*(0.5d0*dyinv))
              end do
            enddo
          enddo
        !$omp end parallel do

    !$omp parallel do
    do k=0,Nz-1
      do i = 2,Ny-2
        do j=0,Nx
          do l=0,4
        !5次精度DCSの右辺設定
        D2(l,j,i,k) = (-Fy(l,j,i-1,k)+Fy(l,j,i+1,k)) * (0.5d0*dyinv)
        D4(l,j,i,k) = (-Fy(l,j,i-2,k)+Fy(l,j,i+2,k)) * (0.25d0*dyinv)
        !D6(l,j,i,k) = (-Fy(l,j,i-3,k)+Fy(l,j,i+3,k)) / (6.d0*dy)
        !7次精度DCS用のため不要
        D6(l,j,i,k) = (Fy(l,j,i-1,k)+Fy(l,j,i+1,k)- 2.d0* Fy(l,j,i,k)) * dyinv
        D8(l,j,i,k) = (Fy(l,j,i-2,k)+Fy(l,j,i+2,k)- 2.d0* Fy(l,j,i,k)) * (0.25d0*dyinv)
        ! D12(l,j,i,k) = (Fy(l,j,i-3,k)+Fy(l,j,i+3,k)- 2.d0* Fy(l,j,i,k)) / (9.d0*dy)
        !7次精度DCS用のため不要
          end do
        end do
      enddo
    enddo
    !$omp end parallel do
    !==========RHSにはD2~D8が必要で、Doループが大きいので、分割して並列化する===============
    !$omp parallel do
      do k=0,Nz-1
        do i = 2,Ny-2
          do j=0,Nx
            do l=0,4
        RHS_y(l,j,i,k)=ra*D2(l,j,i,k)+rb*D4(l,j,i,k)+sigma*(da*D6(l,j,i,k)+db*D8(l,j,i,k))
            end do
          end do
        enddo
      enddo
    !$omp end parallel do

      !前進代入法、後退代入法の計算サブルーチン(y方向)
        !前進代入
      !$omp parallel do
        do k=0,Nz-1
          do j=0,Nx
            do l=0,4
              y(l,j,0,k) = RHS_y(l,j,0,k)
            end do
          end do
        end do
      !$omp end parallel do
    !=============前後のyを使うので、並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝=====
        do i = 1,Ny!y(:,:,i,:)でi=0は上で定義したので残りの1〜Nyを定義する
          y(:,:,i,:) = RHS_y(:,:,i,:) - LU(-1,i)*y(:,:,i-1,:)
        enddo
    !=============並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝====＝＝＝======
        !後退代入
        !$omp parallel do
          do k=0,Nz-1
            do j=0,Nx
              do l=0,4
                x(l,j,Ny,k) = y(l,j,Ny,k) / LU(0,Ny)!境界値(Ny)
              end do
            end do
          end do
        !$omp end parallel do
    !=============前後のxを使うので、並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
        !x(:,:,i)でi=Nyは定義したので残りのi=0~Ny-1を定義する
        do i = Ny-1, 0, -1!後退するので-1ずつ進む
          x(:,:,i,:) = (y(:,:,i,:) - LU(1,i)*x(:,:,i+1,:)) / LU(0,i)
        enddo
    !=============並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝====＝＝＝======

        call combine_y(dzeta_iny,x,dFzeta)!y方向の格子伸長を適用　微分変換をしている
        deallocate(D2,D4,D6,D8,x,y,RHS_y)
       end subroutine dif_y

       subroutine dif_z(sigma,dz,Fz,x,LU)
         integer i,j,k,l
         double precision,allocatable,dimension(:,:,:,:):: Fz,x
         double precision,allocatable,dimension(:,:,:,:):: D2,D4,D6,D8
         double precision,allocatable,dimension(:,:,:,:):: y,RHS_z
         double precision,allocatable,dimension(:,:) :: LU
         double precision dz,sigma,dzinv,Lsum
         double  precision,parameter :: ra = 14.d0/9.d0, rb = 1.d0/9.d0&
         &,da = 4.d0 / 9.d0,db = 2.d0 / 9.d0 !5次精度のDCSとなるための係数設定

         allocate(D2(0:4,0:Nx,0:Ny,2:Nz-3),D4(0:4,0:Nx,0:Ny,2:Nz-3),&
         D6(0:4,0:Nx,0:Ny,2:Nz-3),D8(0:4,0:Nx,0:Ny,2:Nz-3))
         allocate(y(0:4,0:Nx,0:Ny,0:Nz-1),RHS_z(0:4,0:Nx,0:Ny,0:Nz-1))

         D2=0.d0;D4=0.d0;D6=0.d0;D8=0.d0;RHS_z=0.d0;y=0.d0;x=0.d0;Lsum=0.d0
         dzinv = 1.d0 / dz

         !$omp parallel do
          do i=0,Ny
            do j=0,Nx
              do l=0,4
          !5次精度DCSの右辺設定
          RHS_z(l,j,i,0) = ra*dzinv*((-Fz(l,j,i,Nz-1)+Fz(l,j,i,1))*0.5d0)+&
                          rb*dzinv*((-Fz(l,j,i,Nz-2)+Fz(l,j,i,2))*0.25d0)+&
                          sigma*dzinv*(da*(Fz(l,j,i,Nz-1)+Fz(l,j,i,1)-2.d0*Fz(l,j,i,0))+&
                          db*((Fz(l,j,i,Nz-2)+Fz(l,j,i,2)-2.d0*Fz(l,j,i,0))*0.25d0))

          RHS_z(l,j,i,1) = ra*dzinv*((-Fz(l,j,i,0)+Fz(l,j,i,2))*0.5d0)+&
                          rb*dzinv*((-Fz(l,j,i,Nz-1)+Fz(l,j,i,3))*0.25d0)+&
                          sigma*dzinv*(da*(Fz(l,j,i,0)+Fz(l,j,i,2)-2.d0*Fz(l,j,i,1))+&
                          db*((Fz(l,j,i,Nz-1)+Fz(l,j,i,3)-2.d0*Fz(l,j,i,1))*0.25d0))

          RHS_z(l,j,i,Nz-2) =  ra*dzinv*((-Fz(l,j,i,Nz-3)+Fz(l,j,i,Nz-1))*0.5d0)+&
                          rb*dzinv*((-Fz(l,j,i,Nz-4)+Fz(l,j,i,0))*0.25d0)+&
                          sigma*dzinv*(da*(Fz(l,j,i,Nz-3)+Fz(l,j,i,Nz-1)-2.d0*Fz(l,j,i,Nz-2))+&
                          db*((Fz(l,j,i,Nz-4)+Fz(l,j,i,0)-2.d0*Fz(l,j,i,Nz-2))*0.25d0))

          RHS_z(l,j,i,Nz-1) =  ra*dzinv*((-Fz(l,j,i,Nz-2)+Fz(l,j,i,0))*0.5d0)+&
                          rb*dzinv*((-Fz(l,j,i,Nz-3)+Fz(l,j,i,1))*0.25d0)+&
                          sigma*dzinv*(da*(Fz(l,j,i,Nz-2)+Fz(l,j,i,0)-2.d0*Fz(l,j,i,Nz-1))+&
                          db*((Fz(l,j,i,Nz-3)+Fz(l,j,i,1)-2.d0*Fz(l,j,i,Nz-1))*0.25d0))
                end do
              enddo
            enddo
          !$omp end parallel do

          !$omp parallel do
           do k = 2,Nz-3
            do i=0,Ny
              do j=0,Nx
                do l=0,4
             D2(l,j,i,k) = (-Fz(l,j,i,k-1)+Fz(l,j,i,k+1)) * (0.5d0*dzinv)
             D4(l,j,i,k) = (-Fz(l,j,i,k-2)+Fz(l,j,i,k+2)) * (0.25d0*dzinv)

             D6(l,j,i,k) = (Fz(l,j,i,k-1)+Fz(l,j,i,k+1)- 2.d0* Fz(l,j,i,k)) * dzinv
             D8(l,j,i,k) = (Fz(l,j,i,k-2)+Fz(l,j,i,k+2)- 2.d0* Fz(l,j,i,k)) * (0.25d0*dzinv)
                end do
              end do
            end do
           enddo
         !$omp end parallel do
         !==========RHSにはD2~D8が必要で、Doループが大きいので、分割して並列化する
         !$omp parallel do
         do k = 2,Nz-3
           do i=0,Ny
            do j=0,Nx
              do l=0,4
            RHS_z(l,j,i,k)=ra*D2(l,j,i,k)+rb*D4(l,j,i,k)+sigma*(da*D6(l,j,i,k)+db*D8(l,j,i,k))
              end do
            end do
           end do
         enddo
         !$omp end parallel do

        !前進代入法、後退代入法の計算サブルーチン(z方向)
          !前進代入
          !$omp parallel do
            do i=0,Ny
             do j=0,Nx
               do l=0,4
                 y(l,j,i,0) = RHS_z(l,j,i,0)!例外の境界値
               end do
             end do
            end do
          !$omp end parallel do
    !=============前後のyを使うので、並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
              do k = 1,Nz-2!y(:,k)でk=0は上で定義したので残りの1〜Nz-2を定義する
                !これはfill-inのない通常部分
                y(:,:,:,k) = RHS_z(:,:,:,k) - LU(-1,k)*y(:,:,:,k-1)
              enddo
              !Fill-in部はΣの計算が必要になるので追加 Lの最後の行i=Nz-1のみ別で計算
              do i= 0,Ny
                do j = 0,Nx
                  do l= 0,4
                    Lsum = 0.d0
                    do k = 0,Nz-3
                      Lsum = Lsum +LU(2,k)*y(l,j,i,k)
                    enddo
                      y(l,j,i,Nz-1) = RHS_z(l,j,i,Nz-1) - LU(-1,Nz-1)*y(l,j,i,Nz-2)-Lsum
                  enddo
                enddo
              enddo
    !=============並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝====＝＝＝======
              !後退代入
              !fill-inの計算の範囲外なので別で計算
              !$omp parallel do
                do i=0,Ny
                 do j=0,Nx
                   do l=0,4
                     x(l,j,i,Nz-1) = y(l,j,i,Nz-1) / LU(0,Nz-1)!例外の境界値
                     x(l,j,i,Nz-2) = (y(l,j,i,Nz-2) -LU(1,Nz-2)*(y(l,j,i,Nz-1)/LU(0,Nz-1)))/ LU(0,Nz-2)
                     !本来は以下のようだったが、x(Nz-1)が計算前で並列化できないので、埋め込んでDoループ1つにまとめた
                     ! x(:,:,:,Nz-2) = (y(:,:,:,Nz-2) -LU(1,Nz-2)*x(:,:,:,Nz-1))/ LU(0,Nz-2)
                   end do
                 end do
                end do
              !$omp end parallel do
    !=============前後のxを使うので、並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
              !x(:,i)でi=Nxは定義したので残りのi=0~Nx-1を定義する
              do k = Nz-3, 0, -1!後退するので-1ずつ進む
                x(:,:,:,k) = (y(:,:,:,k) - LU(1,k)*x(:,:,:,k+1)-LU(-2,k)*x(:,:,:,Nz-1)) / LU(0,k)
              enddo
    !=============並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝====＝＝＝======

          deallocate(D2,D4,D6,D8,y,RHS_z)
       endsubroutine dif_z
    !境界条件をNSCBCで設定
    !L,d行列を設定することでdfx(0:2,0)とdFx(0:2,Nx)の値を定める
    !定めたdFxをxm,xpをそれぞれ求めた時点でi=0とNxで置き換える
    !これが境界条件となる。neumannはNx-1とNx,1と0を同じにしたがNSCBC_xはdFxの境界を定めることで
    !境界条件を設定している
    !そしてそのままQ1,Q2,Qnを求める
    !まずはx方向用のNSCBC　subrouitineを作成
    !(:,:,:)等で計算を行わせる際は:の配列数が対応していることが要確認
    subroutine NSCBC_x_0_super(dFx)
      !超音速流入条件
      !u,v,w,Tはtop-hat,Crocce-Busemannとランダム撹乱により流入条件として固定してる(imposed, 課されている)ので、
      !このNSCBCでは密度ρのみを求めるものである。
      !密度ρはQ(0)である。Q(0)はFx,y,z(0)とVx,y,z(0)から求められる
      ! NSCBCではFxの書き換えを行う。その中で、必要なのは、Fx(0)のみである。
      !Q(1:4)は上の流入条件で最終的に上書きしてしまうので、Fx(1:4)を求めたとしてもQになってから全て上書き消去されるので
      !わざわざNSCBCで計算=>上書きしても無駄。=>そのためこのsubroutineでは計算しない!!!!!!!!!!!!!!!!
      double precision,allocatable,dimension(:,:,:,:):: dFx
      integer i,k
      !超音速流入では、領域内部から、外部に逆流するものがないので、L1=0したがって、L5=0
      !*uが時間変動しないtop-hatのみなので、du/dt=0のため
      !よって、L2=0そのため、d1=0
      !*ただし、流入条件のTはCrocce-Busemannで定常なので、時間変動しないdT/dt=0のため
      !Fx(0)=d1なので、Fx(0)=0。したがって、今回はそれだけを書き換えている
    !$omp parallel do
      do k = 0,Nz-1
        do i = 0,Ny
      !設定したdからNSCBCで置き換える境界地点のdFxを定義する
      !i=0の時の差し替えdFx
      dFx(0,0,i,k) = 0.d0
        end do
      end do
    !$omp end parallel do
  endsubroutine NSCBC_x_0_super
  subroutine NSCBC_x_0_sub(G,dGx,dFx)
    !亜音速流入条件
    !u,v,w,Tはtop-hat,Crocce-Busemannとランダム撹乱により流入条件として固定してる(imposed, 課されている)ので、
    !このNSCBCでは密度ρのみを求めるものである。
    !密度ρはQ(0)である。Q(0)はFx,y,z(0)とVx,y,z(0)から求められる
    ! NSCBCではFxの書き換えを行う。その中で、必要なのは、Fx(0)のみである。
    !Q(1:4)は上の流入条件で最終的に上書きしてしまうので、Fx(1:4)を求めたとしてもQになってから全て上書き消去されるので
    !わざわざNSCBCで計算=>上書きしても無駄。=>そのためこのsubroutineでは計算しない!!!!!!!!!!!!!!!!
    double precision,allocatable,dimension(:,:,:,:):: dFx
    double precision,allocatable,dimension(:,:,:,:):: G,dGx
    double precision,allocatable,dimension(:,:,:):: L0
    double precision,allocatable,dimension(:,:):: c_NS,Ma_NS
    integer i,k!,l
    ! double precision p0x_infty
    allocate(L0(1:5,0:Ny,0:Nz-1))
    allocate(c_NS(0:Ny,0:Nz-1),Ma_NS(0:Ny,0:Nz-1))

    L0=0.d0;c_NS=0.d0;Ma_NS=0.d0

  !========並列化しない(Maにはcが必要なため、Doループを分割しないといけないから)＝＝＝＝＝＝＝＝＝
    do k = 0,Nz-1
      do i = 0,Ny
        !音速cはi=0,Nxの両点においてそれぞれ定義しなければならない
        c_NS(i,k) = sqrt(gamma * G(4,0,i,k) / G(0,0,i,k))
        !マッハ数Ma_NSはi=0,Nxで使うので別々に定義する
        Ma_NS(i,k) = G(1,0,i,k) / c_NS(i,k)!uを使う
      end do
    end do
  !========並列化しない==========================================================

  !============L0(2)は同時並列化できないので、並列化をしない===========================
    do k = 0,Nz-1
      do i = 0,Ny
  !   x方向右側つまりi=0の点において亜音速流入条件でL行列を設定する
    L0(1,i,k)=(G(1,0,i,k)-c_NS(i,k))*(-G(0,0,i,k)*c_NS(i,k)*dGx(1,0,i,k)+dGx(4,0,i,k))
  !   論文によるとL3,L4は不要
  !   L1=0.d0&uはtop-hatジェットで程上流なので、du/dt=0
    L0(5,i,k)=L0(1,i,k)!-2.d0*c_NS(i,k)*du/dt!
  !   流入速度uを時間変動させないので今回はdu/dt=0となるため省略
    L0(2,i,k)=(0.5d0)*(gamma-1.d0)*(L0(5,i,k)+L0(1,i,k))!+G(0,0,i,k)*c_NS(i,k)**2.d0/T*dT/dtが本来はあるが
  !   流入条件のTは時間変動させずに、Ceocco-Busemannのやつで固定なので、dT/dt＝0となり計算不要
      end do
    end do
  !============L0(2)は同時並列化できないので、並列化をしない===========================

    !$omp parallel do
      do k = 0,Nz-1
        do i = 0,Ny
          !設定したL行列からd1をi=0において設定する
          !設定したd1からNSCBCで置き換える境界地点のdFxを定義する
          !i=0の時の差し替えdFx
          ! d1(i,k) = (1.d0 / (c_NS(i,k) **2.d0)) * ((L0(1,i,k)+L0(5,i,k))*0.5d0 + L0(2,i,k))
          ! dFx(0,0,i,k) = d1(i,k)
          !本来は一度d1に格納するが、dFx(0)=d1なので、d1を省略してしまう
          dFx(0,0,i,k) = (1.d0 / (c_NS(i,k) **2.d0)) * ((L0(1,i,k)+L0(5,i,k))*0.5d0 + L0(2,i,k))
        end do
      end do
    !$omp end parallel do

    ! !$omp parallel do
    !   do k = 0,Nz-1
    !     do l = 0,4
    !   !NSCBCの角処理(x方向,y方向で設定した境界値が重複するため1/2ずつ加える)
    !   !z方向に関しyては、周期で今回はz方向の縁の重複の話なので、関係ない。
    !   dFx(l,0,0,k) = dFx(l,0,0,k) / 2.d0
    !   dFx(l,0,Ny,k) = dFx(l,0,Ny,k) / 2.d0
    !     end do
    !   end do
    ! !$omp end parallel do

      deallocate(L0,c_NS,Ma_NS)
    endsubroutine NSCBC_x_0_sub

      subroutine NSCBC_x_Nx_super(G,dGx,dFx)
        !超音速無反射流出条件
        double precision,allocatable,dimension(:,:,:,:):: G,dGx,dFx
        double precision,allocatable,dimension(:,:,:):: LNx,dNx
        double precision,allocatable,dimension(:,:):: c_NS,Ma_NS
        integer i,k!,l
        allocate(LNx(1:5,0:Ny,0:Nz-1),dNx(1:5,0:Ny,0:Nz-1))
        allocate(c_NS(0:Ny,0:Nz-1),Ma_NS(0:Ny,0:Nz-1))
        LNx=0.d0;dNx=0.d0;c_NS=0.d0;Ma_NS=0.d0

      !========並列化しない(Maにはcが必要なため、Doループを分割しないといけないから)＝＝＝＝＝＝＝＝＝
        do k = 0,Nz-1
          do i = 0,Ny
        !音速cはi=0,Nxの両点においてそれぞれ定義しなければならない
        c_NS(i,k) = sqrt(gamma * G(4,Nx,i,k) / G(0,Nx,i,k))
        !マッハ数Ma_NSはi=0,Nxで使うので別々に定義する
        Ma_NS(i,k) = G(1,Nx,i,k) / c_NS(i,k)
          end do
        end do
      !========並列化しない＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝=＝＝＝＝＝＝＝＝＝

      !$omp parallel do
        do k = 0,Nz-1
          do i = 0,Ny
        !x方向左側つまりi=Nxの点において無反射流出条件でL行列を設定する
        ! LNx(1,i,k)=NS_sigma*c_NS(i,k)*(1.d0-(Ma_NS(i,k)**2.d0))*(G(4,Nx,i,k)-&
        ! &pNx_infty)/Lx
        !超音速無反射流出条件にするので、領域外から内部に流入する波が無いとする。したがって、L１=0とする
        LNx(1,i,k)=0.d0
        LNx(2,i,k)=G(1,Nx,i,k)*((c_NS(i,k)**2.d0)*dGx(0,Nx,i,k)-dGx(4,Nx,i,k))
        LNx(3,i,k)=G(1,Nx,i,k)*dGx(2,Nx,i,k)
        LNx(4,i,k)=G(1,Nx,i,k)*dGx(3,Nx,i,k)
        LNx(5,i,k)=(G(1,Nx,i,k)+c_NS(i,k))*(G(0,Nx,i,k)*c_NS(i,k)*dGx(1,Nx,i,k)+dGx(4,Nx,i,k))
          end do
        end do
      !$omp end parallel do

        !$omp parallel do
          do k = 0,Nz-1
            do i = 0,Ny
        !設定したL行列からd1~5をi=Nxにおいて設定する
          dNx(1,i,k) = (1.d0 / (c_NS(i,k) **2.d0)) * ((LNx(1,i,k)+LNx(5,i,k))*0.5d0 + LNx(2,i,k))
          dNx(2,i,k) = (LNx(1,i,k)+LNx(5,i,k))*0.5d0
          dNx(3,i,k) = 0.5d0/(G(0,Nx,i,k) * c_NS(i,k)) * (-LNx(1,i,k) + LNx(5,i,k))
          dNx(4,i,k) = LNx(3,i,k)
          dNx(5,i,k) = LNx(4,i,k)
            end do
          end do
        !$omp end parallel do

        !$omp parallel do
          do k = 0,Nz-1
            do i = 0,Ny
        !設定したdからNxSCBCで置き換える境界地点のdFxを定義する
        !i=Nxの時の差し替えF
        dFx(0,Nx,i,k) = dNx(1,i,k)
        dFx(1,Nx,i,k) = (G(1,Nx,i,k)*dNx(1,i,k)) + (G(0,Nx,i,k)*dNx(3,i,k))
        dFx(2,Nx,i,k) = (G(2,Nx,i,k)*dNx(1,i,k)) + (G(0,Nx,i,k)*dNx(4,i,k))
        dFx(3,Nx,i,k) = (G(3,Nx,i,k)*dNx(1,i,k)) + (G(0,Nx,i,k)*dNx(5,i,k))
        dFx(4,Nx,i,k) =(0.5d0)*((G(1,Nx,i,k)**2.d0)+(G(2,Nx,i,k)**2.d0)+&
                      (G(3,Nx,i,k)**2.d0))*dNx(1,i,k)+dNx(2,i,k)/(gamma-1.d0)+ &
        (G(0,Nx,i,k)*G(1,Nx,i,k)*dNx(3,i,k))+(G(0,Nx,i,k)*G(2,Nx,i,k)*dNx(4,i,k))+&
                      (G(0,Nx,i,k)*G(3,Nx,i,k)*dNx(5,i,k))
            end do
          end do
      !$omp end parallel do

    !   !$omp parallel do
    !     do k = 0,Nz-1
    !       do l = 0,4
    !         !NSCBCの角処理(x方向,y方向で設定した境界値が重複するため1/2ずつ加える)
    !         dFx(l,Nx,0,k) = dFx(l,Nx,0,k) / 2.d0
    !         dFx(l,Nx,Ny,k) = dFx(l,Nx,Ny,k) / 2.d0
    !       end do
    !     end do
    ! !$omp end parallel do
      deallocate(LNx,dNx,c_NS,Ma_NS)
    endsubroutine NSCBC_x_Nx_super

    subroutine NSCBC_x_Nx_sub(G,dGx,dFx,pNx_infty)
      !亜音速無反射流出条件
      double precision,allocatable,dimension(:,:,:,:):: G,dGx,dFx
      double precision,allocatable,dimension(:,:,:):: LNx,dNx
      double precision,allocatable,dimension(:,:):: c_NS,Ma_NS
      double precision pNx_infty
      integer i,k!,l
      allocate(LNx(1:5,0:Ny,0:Nz-1),dNx(1:5,0:Ny,0:Nz-1))
      allocate(c_NS(0:Ny,0:Nz-1),Ma_NS(0:Ny,0:Nz-1))
      LNx=0.d0;dNx=0.d0;c_NS=0.d0;Ma_NS=0.d0

    !========並列化しない(Maにはcが必要なため、Doループを分割しないといけないから)＝＝＝＝＝＝＝＝＝
      do k = 0,Nz-1
        do i = 0,Ny
      !音速cはi=0,Nxの両点においてそれぞれ定義しなければならない
      c_NS(i,k) = sqrt(gamma * G(4,Nx,i,k) / G(0,Nx,i,k))
      !マッハ数Ma_NSはi=0,Nxで使うので別々に定義する
      Ma_NS(i,k) = G(1,Nx,i,k) / c_NS(i,k)
        end do
      end do
    !========並列化しない＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝=＝＝＝＝＝＝＝＝＝

    !$omp parallel do
      do k = 0,Nz-1
        do i = 0,Ny
      !x方向左側つまりi=Nxの点において無反射流出条件でL行列を設定する
      ! LNx(1,i,k)=NS_sigma*c_NS(i,k)*(1.d0-(Ma_NS(i,k)**2.d0))*(G(4,Nx,i,k)-&
      ! &pNx_infty)/Lx
      !超音速無反射流出条件にするので、領域外から内部に流入する波が無いとする。したがって、L１=0とする
      LNx(1,i,k)=0.d0
      LNx(2,i,k)=G(1,Nx,i,k)*((c_NS(i,k)**2.d0)*dGx(0,Nx,i,k)-dGx(4,Nx,i,k))
      LNx(3,i,k)=G(1,Nx,i,k)*dGx(2,Nx,i,k)
      LNx(4,i,k)=G(1,Nx,i,k)*dGx(3,Nx,i,k)
      LNx(5,i,k)=(G(1,Nx,i,k)+c_NS(i,k))*(G(0,Nx,i,k)*c_NS(i,k)*dGx(1,Nx,i,k)+dGx(4,Nx,i,k))
        end do
      end do
    !$omp end parallel do

      !$omp parallel do
        do k = 0,Nz-1
          do i = 0,Ny
      !設定したL行列からd1~5をi=Nxにおいて設定する
        dNx(1,i,k) = (1.d0 / (c_NS(i,k) **2.d0)) * ((LNx(1,i,k)+LNx(5,i,k))*0.5d0 + LNx(2,i,k))
        dNx(2,i,k) = (LNx(1,i,k)+LNx(5,i,k))*0.5d0
        dNx(3,i,k) = 0.5d0/(G(0,Nx,i,k) * c_NS(i,k)) * (-LNx(1,i,k) + LNx(5,i,k))
        dNx(4,i,k) = LNx(3,i,k)
        dNx(5,i,k) = LNx(4,i,k)
          end do
        end do
      !$omp end parallel do

      !$omp parallel do
        do k = 0,Nz-1
          do i = 0,Ny
      !設定したdからNxSCBCで置き換える境界地点のdFxを定義する
      !i=Nxの時の差し替えF
      dFx(0,Nx,i,k) = dNx(1,i,k)
      dFx(1,Nx,i,k) = (G(1,Nx,i,k)*dNx(1,i,k)) + (G(0,Nx,i,k)*dNx(3,i,k))
      dFx(2,Nx,i,k) = (G(2,Nx,i,k)*dNx(1,i,k)) + (G(0,Nx,i,k)*dNx(4,i,k))
      dFx(3,Nx,i,k) = (G(3,Nx,i,k)*dNx(1,i,k)) + (G(0,Nx,i,k)*dNx(5,i,k))
      dFx(4,Nx,i,k) =(0.5d0)*((G(1,Nx,i,k)**2.d0)+(G(2,Nx,i,k)**2.d0)+&
                    (G(3,Nx,i,k)**2.d0))*dNx(1,i,k)+dNx(2,i,k)/(gamma-1.d0)+ &
      (G(0,Nx,i,k)*G(1,Nx,i,k)*dNx(3,i,k))+(G(0,Nx,i,k)*G(2,Nx,i,k)*dNx(4,i,k))+&
                    (G(0,Nx,i,k)*G(3,Nx,i,k)*dNx(5,i,k))
          end do
        end do
    !$omp end parallel do

    !   !$omp parallel do
    !     do k = 0,Nz-1
    !       do l = 0,4
    !         !NSCBCの角処理(x方向,y方向で設定した境界値が重複するため1/2ずつ加える)
    !         dFx(l,Nx,0,k) = dFx(l,Nx,0,k) / 2.d0
    !         dFx(l,Nx,Ny,k) = dFx(l,Nx,Ny,k) / 2.d0
    !       end do
    !     end do
    ! !$omp end parallel do
    deallocate(LNx,dNx,c_NS,Ma_NS)
    endsubroutine NSCBC_x_Nx_sub

    !次にy方向のNSCBC　sunrouineを作成
    subroutine NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
      double precision,allocatable,dimension(:,:,:,:):: G,dGy,dFy
      double precision,allocatable,dimension(:,:,:):: L0,d0,L1,d1
      double precision,allocatable,dimension(:,:):: c_NS0,Ma_NS0,c_NS1,Ma_NS1
      double precision pNy_infty,p0y_infty
      integer j,k

      allocate(L0(1:5,0:Nx,0:Nz-1),d0(1:5,0:Nx,0:Nz-1),L1(1:5,0:Nx,0:Nz-1),d1(1:5,0:Nx,0:Nz-1))
      allocate(c_NS0(0:Nx,0:Nz-1),Ma_NS0(0:Nx,0:Nz-1),c_NS1(0:Nx,0:Nz-1),Ma_NS1(0:Nx,0:Nz-1))

      L0=0.d0;d0=0.d0;c_NS0=0.d0;Ma_NS0=0.d0
      L1=0.d0;d1=0.d0;c_NS1=0.d0;Ma_NS1=0.d0

    !========並列化しない(Maにはcが必要なため、Doループを分割しないといけないから)＝＝＝＝＝＝＝＝＝
      do k=0,Nz-1
        do j=0,Nx
      !音速cはi=0,Nxの両点においてそれぞれ定義しなければならない
      c_NS0(j,k) = sqrt(gamma * G(4,j,0,k) / G(0,j,0,k))
      c_NS1(j,k) = sqrt(gamma * G(4,j,Ny,k) / G(0,j,Ny,k))
      !マッハ数Maはi=0,Nxで使うので別々に定義する
      Ma_NS0(j,k) = G(2,j,0,k) / c_NS0(j,k)!vを使う
      Ma_NS1(j,k) = G(2,j,Ny,k) / c_NS1(j,k)
        end do
      end do
    !========並列化しない＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
!亜音速流出条件
!主流垂直方向なので、主流方向と違い、流速が超音速ではないと推測し、亜音速流出条件にする
    !$omp parallel do
      do k=0,Nz-1
        do j=0,Nx
      !y方向右側つまりi=0の点において無反射流出条件でL行列を設定する
      L0(1,j,k) = (G(2,j,0,k) - c_NS0(j,k)) * (-G(0,j,0,k)*c_NS0(j,k)*dGy(2,j,0,k)+dGy(4,j,0,k))
      L0(2,j,k) = G(2,j,0,k) * ((c_NS0(j,k) ** 2.d0)*dGy(0,j,0,k) - dGy(4,j,0,k))
      L0(3,j,k) = G(2,j,0,k) * dGy(1,j,0,k)
      L0(4,j,k) = G(2,j,0,k) * dGy(3,j,0,k)
      L0(5,j,k) = NS_sigma * c_NS0(j,k) * (1.d0 - (Ma_NS0(j,k) ** 2.d0))*(G(4,j,0,k) - &
                    &p0y_infty)/Ly
        end do
      end do
    !$omp end parallel do

    !$omp parallel do
      do k=0,Nz-1
        do j=0,Nx
      !y方向左側つまりi=Nyの点において無反射流出条件でL行列を設定する
      L1(1,j,k) = NS_sigma * c_NS1(j,k) * (1.d0 - (Ma_NS1(j,k) ** 2.d0))*(G(4,j,Ny,k) - &
    &  pNy_infty)/Ly
      L1(2,j,k) = G(2,j,Ny,k) * ((c_NS1(j,k) ** 2.d0)*dGy(0,j,Ny,k) - dGy(4,j,Ny,k))
      L1(3,j,k) = G(2,j,Ny,k) * dGy(1,j,Ny,k)
      L1(4,j,k) = G(2,j,Ny,k) * dGy(3,j,Ny,k)
      L1(5,j,k) = (G(2,j,Ny,k) + c_NS1(j,k)) * (G(0,j,Ny,k)*c_NS1(j,k)*dGy(2,j,Ny,k) + dGy(4,j,Ny,k))
        end do
      end do
    !$omp end parallel do

      !$omp parallel do
        do k=0,Nz-1
          do j=0,Nx
      !設定したL行列からd1~5をi=0,Nyの両点においてそれぞれ設定する
        d0(1,j,k) = (1.d0 / (c_NS0(j,k) **2.d0)) * ((L0(1,j,k)+L0(5,j,k))*0.5d0 + L0(2,j,k))
        d0(2,j,k) = (0.5d0) * (L0(1,j,k)+L0(5,j,k))
        d0(3,j,k) = L0(3,j,k)

        d1(1,j,k) = (1.d0 / (c_NS1(j,k) **2.d0)) * ((L1(1,j,k)+L1(5,j,k))*0.5d0 + L1(2,j,k))
        d1(2,j,k) = (0.5d0) * (L1(1,j,k)+L1(5,j,k))
        d1(3,j,k) = L1(3,j,k)
        !d4のみrhoを含むので個別で設定しなければいけない
        d0(4,j,k) = 0.5d0/(G(0,j,0,k) * c_NS0(j,k)) * (-L0(1,j,k) + L0(5,j,k))
        d1(4,j,k) = 0.5d0/(G(0,j,Ny,k) * c_NS1(j,k)) * (-L1(1,j,k) + L1(5,j,k))

        d0(5,j,k) = L0(4,j,k)
        d1(5,j,k) = L1(4,j,k)
      !設定したdからNySCBCで置き換える境界地点のdFyを定義する
          end do
        end do
      !$omp end parallel do

    !$omp parallel do
      do k=0,Nz-1
        do j=0,Nx
      !i=0の時の差し替えdFy
      dFy(0,j,0,k) = d0(1,j,k)
      dFy(1,j,0,k) = G(1,j,0,k)*d0(1,j,k)+G(0,j,0,k)*d0(3,j,k)
      dFy(2,j,0,k) = G(2,j,0,k)*d0(1,j,k)+G(0,j,0,k)*d0(4,j,k)
      dFy(3,j,0,k) = G(3,j,0,k)*d0(1,j,k)+G(0,j,0,k)*d0(5,j,k)
      dFy(4,j,0,k) =(0.5d0)*((G(1,j,0,k)**2.d0)+(G(2,j,0,k)**2.d0)+(G(3,j,0,k)**2.d0))*d0(1,j,k)+&
                    d0(2,j,k)/(gamma-1.d0)+(G(0,j,0,k)*G(1,j,0,k)*d0(3,j,k))+&
          (G(0,j,0,k)*G(2,j,0,k)*d0(4,j,k))+(G(0,j,0,k)*G(3,j,0,k)*d0(5,j,k))
        end do
      end do
    !$omp end parallel do

    !$omp parallel do
      do k=0,Nz-1
        do j=0,Nx
      !i=Nyの時の差し替えdFy
      dFy(0,j,Ny,k) = d1(1,j,k)
      dFy(1,j,Ny,k) = (G(1,j,Ny,k)*d1(1,j,k)) + (G(0,j,Ny,k)*d1(3,j,k))
      dFy(2,j,Ny,k) = (G(2,j,Ny,k)*d1(1,j,k)) + (G(0,j,Ny,k)*d1(4,j,k))
      dFy(3,j,Ny,k) = (G(3,j,Ny,k)*d1(1,j,k)) + (G(0,j,Ny,k)*d1(5,j,k))
      dFy(4,j,Ny,k) =(0.5d0)*((G(1,j,Ny,k)**2.d0)+(G(2,j,Ny,k)**2.d0)+(G(3,j,Ny,k)**2.d0))*d1(1,j,k)+&
                  d1(2,j,k)/(gamma-1.d0)+(G(0,j,Ny,k)*G(1,j,Ny,k)*d1(3,j,k))+&
      (G(0,j,Ny,k)*G(2,j,Ny,k)*d1(4,j,k))+(G(0,j,Ny,k)*G(3,j,Ny,k)*d1(5,j,k))
        end do
      end do
    !$omp end parallel do

      !NSCBCの縁処理(x方向,y方向で設定した境界値が重複するため1/2ずつ加える)
      ! !$omp parallel do
      ! do k=0,Nz-1
      !   do l=0,4
      !     dFy(l,0,0,k) = dFy(l,0,0,k) / 2.d0
      !     dFy(l,0,Ny,k) = dFy(l,0,Ny,k) / 2.d0
      !     dFy(l,Nx,0,k) = dFy(l,Nx,0,k) / 2.d0
      !     dFy(l,Nx,Ny,k) = dFy(l,Nx,Ny,k) / 2.d0
      !   end do
      ! end do
      ! !$omp end parallel do
      deallocate(L0,d0,c_NS0,Ma_NS0)
      deallocate(L1,d1,c_NS1,Ma_NS1)
    endsubroutine NSCBC_y
    !超音速・亜音速に関係なく、全体にNeumann条件を設定したい時に使うsubroutine
    subroutine Q_boundary(Q)
      double precision,allocatable,dimension(:,:,:,:):: Q
      integer i,j,k,l
      !$omp parallel do
        do k=0,Nz-1
          do i=0,Ny
            do l=0,4
              Q(l,0,i,k) = Q(l,1,i,k)
              Q(l,Nx,i,k) = Q(l,Nx-1,i,k)
            end do
          enddo
        end do
      !$omp end parallel do

      !$omp parallel do
        do k=0,Nz-1
          do j=0,Nx
            do l=0,4
              Q(l,j,0,k) = Q(l,j,1,k)
              Q(l,j,Ny,k) = Q(l,j,Ny-1,k)
            end do
          enddo
        end do
      !$omp end parallel do
    endsubroutine Q_boundary

    subroutine inflow(M,Q,in_G0,in_G1_top,in_G2,in_G3)
    ! subroutine inflow(M,Q,in_G0,in_G1_top,in_G1_du,in_G2,in_G3)
      double precision,allocatable,dimension(:,:,:,:):: Q
      double precision,allocatable,dimension(:,:):: in_G0,in_G1_top,in_G2,in_G3
      ! double precision,allocatable,dimension(:,:):: in_G1_du
      double precision :: fluct_dis_strength
      integer i,k,M
      if (M < times) then
        fluct_dis_strength = dis_strength*dble(M)/dble(times)
      else
        fluct_dis_strength = dis_strength
      endif
      !$omp parallel do
        do k =0,Nz-1
          do i=0,Ny
            !Q(0)に関しては、NSCBCを使用して求めたF,Vから求めたQ(0)の密度を使用する
        !uに撹乱を入れないパターン=>これで、流入条件の計算で、時間変動を気にしなくて良くなる
       Q(1,0,i,k) = in_G0(i,k)*in_G1_top(i,k)!rho*u_in
       Q(2,0,i,k) = in_G0(i,k)*fluct_dis_strength*in_G2(i,k)!rho*kakuran_v
       Q(3,0,i,k) = in_G0(i,k)*fluct_dis_strength*in_G3(i,k)!rho*kakuran_w
       Q(4,0,i,k) = 1.d0/((Ma**2.d0)*gamma*(gamma-1.d0))&
                   +in_G0(i,k)*((in_G1_top(i,k))**2.d0&
                   +(fluct_dis_strength*in_G2(i,k))**2.d0&
                   +(fluct_dis_strength*in_G3(i,k))**2.d0)*0.5d0!Et
       !uに撹乱を入れるパターン=>そのため、NSCBC流入条件の計算で、時間変動を計算しなくてはいけない
       ! Q(1,0,i,k) = in_G0(i,k)*(in_G1_top(i,k)+fluct_dis_strength*in_G1_du(i,k))!rho*(u_in+kakuran_u)
       ! Q(4,0,i,k) = 1.d0/((Ma**2.d0)*gamma*(gamma-1.d0))&
       !             +in_G0(i,k)*((in_G1_top(i,k)+fluct_dis_strength*in_G1_du(i,k))**2.d0&
       !             +(fluct_dis_strength*in_G2(i,k))**2.d0&
       !             +(fluct_dis_strength*in_G3(i,k))**2.d0)*0.5d0!Et
          enddo
        end do
      !$omp end parallel do
    endsubroutine inflow
    !NSCBC_x_Nxが不要ならoutflowも不要なので、outflowをxとyに分割する
    !矩型Jetなどを流入させるようになったら、部分的に必要なので、修正して適用する
    subroutine outflow_x(UVWT,dUVWTx,Vx,dVx)
      double precision,allocatable,dimension(:,:,:,:):: Vx,dUVWTx,UVWT,dVx
      integer i,k
      !無反射流出条件の時の条件を設定するsubroutine
      !境界のdVx,dVyに条件を設定するのでdVの計算ができた後に境界値のみ上書きをする
      !τ11=Vx(1),τ12=Vx(2),τ13=Vx(3)
      !τ21=Vy(1),τ22=Vy(2),τ23=Vy(3)
      !dτ12/dx,dτ13/dx,dτ21/dy,dτ23/dy,dq/dx,dq/dyの全てが0(境界のみ)
      !$omp parallel do
        do k=0,Nz-1
          do i=0,Ny
        !x方向右側の条件設定
        dVx(2,Nx,i,k) = 0.d0!dτ12/dx
        dVx(3,Nx,i,k) = 0.d0!dτ13/dx
        dVx(4,Nx,i,k) = dVx(1,Nx,i,k)*UVWT(1,Nx,i,k)+dUVWTx(1,Nx,i,k)*Vx(1,Nx,i,k)&
                      +dUVWTx(2,Nx,i,k)*Vx(2,Nx,i,k)+dUVWTx(3,Nx,i,k)*Vx(3,Nx,i,k)
        !dτ11/dx*u+τ11*du/dx+dτ12/dx(=0)*v+τ12*dv/dx+dτ13/dx(=0)*w+τ13*dw/dx+dq/dx(=0)
        !=dτ11/dx*u+τ11*du/dx+τ12*dv/dx+τ13*dw/dx
          end do
        end do
      !$omp end parallel do
    endsubroutine outflow_x
    subroutine outflow_y(UVWT,dUVWTy,Vy,dVy)
      double precision,allocatable,dimension(:,:,:,:):: Vy,dUVWTy,UVWT,dVy
      integer j,k
      !無反射流出条件の時の条件を設定するsubroutine
      !$omp parallel do
        do k=0,Nz-1
            do j=0,Nx
        !y方向左側の条件設定
        dVy(1,j,0,k) = 0.d0!dτ21/dy
        dVy(3,j,0,k) = 0.d0!dτ23/dy
        dVy(4,j,0,k) = dVy(2,j,0,k)*UVWT(2,j,0,k)+dUVWTy(1,j,0,k)*Vy(1,j,0,k)&
                      +dUVWTy(2,j,0,k)*Vy(2,j,0,k)+dUVWTy(3,j,0,k)*Vy(3,j,0,k)
        !dτ21/dy(=0)*u+τ21*du/dy+dτ22/dy*v+τ22*dv/dy+dτ23/dy(=0)*w+τ23*dw/dy+dq/dy(=0)
        !τ21*du/dy+dτ22/dy*v+τ22*dv/dy+τ23*dw/dy
        !y方向右側の条件設定
        dVy(1,j,Ny,k) = 0.d0
        dVy(3,j,Ny,k) = 0.d0
        dVy(4,j,Ny,k) = dVy(2,j,Ny,k)*UVWT(2,j,Ny,k)+dUVWTy(1,j,Ny,k)*Vy(1,j,Ny,k)&
                      +dUVWTy(2,j,Ny,k)*Vy(2,j,Ny,k)+dUVWTy(3,j,Ny,k)*Vy(3,j,Ny,k)
                      !dq/dy
          end do
        end do
      !$omp end parallel do
    endsubroutine outflow_y
    !buffer領域の設定subroutine
    !x方向
    !ここでまず計算に必要なU(x)とσ(x)を定義している
    !係数等はタサイさんの資料を参考にしている
    subroutine buffer_x(c_infty,Ux,sigma_x,zeta_fx)
      integer i
      double precision c_infty,Xmax,Xmin
      double precision,allocatable,dimension(:):: Ux,sigma_x
      double precision,allocatable,dimension(:):: zeta_fx
      double precision,parameter ::alpha_u=1.5d0,alpha_sigma=1.125d0,beta_r=0.01d0,beta_l=0.01d0
      Xmax = Lx;Xmin = 0.d0
      !$omp parallel do
        do i = 0,Nx
          Ux(i) = alpha_u*c_infty*(dtanh(dble(atanh(beta_r/alpha_u-1.d0))*(zeta_fx(i)-Xmax)/(-Wrx))&
          -dtanh(dble(atanh(beta_l/alpha_u-1.d0))*(zeta_fx(i)-Xmin)/Wlx))

          if(zeta_fx(i)<(Xmax-Wrx)) then
            sigma_x(i) = 0.d0!流出部にのみBufferをつけるのでx左側もσは0となる
          else
            sigma_x(i) = alpha_sigma*c_infty*((zeta_fx(i)-(Xmax-Wrx))/Wrx)**3.d0
          endif
        enddo
      !$omp end parallel do
       !dx不要かも？
     Ux(0:NUx) = 0.d0!x左側のBufferを取るためにWlxの範囲のUxを確実に0に設定している
     ! open(100,file="ux-check.csv")
     ! do i=0,Nx
     ! write(100,*) zeta_fx(i),",",Ux(i)
     ! enddo
     ! close(100)
      !格子伸長が入っているためほぼ計算領域の真ん中になる値を調べて代入した
    endsubroutine buffer_x
    !y方向
    subroutine buffer_y(c_infty,Uy,sigma_y,zeta_fy)
      integer i
      double precision Ymax,Ymin,c_infty
      double precision,allocatable,dimension(:):: Uy,sigma_y
      double precision,allocatable,dimension(:):: zeta_fy
      double precision,parameter ::alpha_u=1.15d0,alpha_sigma=1.125d0,beta_r=0.01d0,beta_l=0.01d0
      Ymax = (Ly/2.d0);Ymin = -(Ly/2d0)
      !格子伸長を行うので新しい座標ζ_yを用いてUyとsigma_yを設定する
      !$omp parallel do
        do i = 0,Ny
          Uy(i) = alpha_u*c_infty*(dtanh(dble(atanh(beta_r/alpha_u-1.d0))*(zeta_fy(i)-Ymax)/(-Wry))&
          -dtanh(dble(atanh(beta_l/alpha_u-1.d0))*(zeta_fy(i)-Ymin)/Wly))

          if(zeta_fy(i)<(Wly+Ymin)) then
            sigma_y(i) = alpha_sigma*c_infty*((-zeta_fy(i)+Ymin+Wly)/Wly)**3.d0
          elseif((zeta_fy(i)>=(Wly+Ymin)).and.(zeta_fy(i)<(Ymax-Wry))) then
            sigma_y(i) = 0.d0
          elseif(zeta_fy(i)>=(Ymax-Wry)) then
            sigma_y(i) = alpha_sigma*c_infty*((zeta_fy(i)-(Ymax-Wry))/Wry)**3.d0
          endif
        enddo
      !$omp end parallel do
    endsubroutine buffer_y
    !ζ,dζ/dxの定義subroutine
    subroutine lattice_x(dx,zeta_fx,dzeta_inx)
      integer i
      double precision,allocatable,dimension(:):: zeta_fx
      double precision,allocatable,dimension(:):: dzeta,dzeta_inx
      double precision dx,width,a1,a2,b1
      allocate(dzeta(0:Nx))
      dzeta=0.d0;width=10.8d0;a1=1d0/14d0;a2=7.d0;b1=1.d0/1.4d0
      !widthは格子間隔を細かくする範囲。この式では-width<=x<=widthの範囲で適用される
      !a2は粗い所と細かい所の境界の傾きの大きさを設定している
      !a1はどの程度の格子数の差をつけるかを設定する係数。このパラメータは条件によって適宜調整する

      !並列化するために、xの座標をzeta_fxやdzetaの式の中に組み込んだ
      !$omp parallel do
        do i= 0,Nx
         zeta_fx(i) = b1 * ((1.7d0*(dx*dble(i))) - a1 * &
          (-dlog(dcosh(a2*((dx*dble(i)) - width))) + dlog(dcosh(a2*((dx*dble(i)) + width)))))
          dzeta(i) = b1 * (1.7d0 - (a1*a2) * &
          (-dtanh(a2*((dx*dble(i)) - width)) + dtanh(a2*((dx*dble(i)) + width))))
        enddo
      !$omp end parallel do
      dzeta_inx = 1.d0/dzeta
      deallocate(dzeta)
    endsubroutine lattice_x
    !ζ,dζ/dyの定義subroutine
    subroutine lattice_y(dy,zeta_fy,dzeta_iny)
      integer i
      double precision,allocatable,dimension(:):: zeta_fy
      double precision,allocatable,dimension(:):: dzeta,dzeta_iny
      double precision dy,width,a1,a2,b1,Ymin
      allocate(dzeta(0:Ny))
      dzeta=0.d0;width=3.d0;a1=1d0/14d0;a2=7d0;b1=1.d0/1.4d0
      !widthは格子間隔を細かくする範囲。この式では-width<=y<=widthの範囲で適用される
      !a2は粗い所と細かい所の境界の傾きの大きさを設定している
      !a1はどの程度の格子数の差をつけるかを設定する係数
      Ymin = -(Ly/2.d0)

      !並列化するために、yの座標をzeta_fyやdzetaの式の中に組み込んだ
      !$omp parallel do
        do i= 0,Ny
          zeta_fy(i) = b1 * ((1.7d0*(Ymin + dy*dble(i))) - a1 * &
        (-dlog(dcosh(a2*((Ymin + dy*dble(i)) - width))) + dlog(dcosh(a2*((Ymin + dy*dble(i)) + width)))))
        dzeta(i) = b1 * (1.7d0 - (a1*a2) * &
        (-dtanh(a2*((Ymin + dy*dble(i)) - width)) + dtanh(a2*((Ymin + dy*dble(i)) + width))))
        enddo
      !$omp end parallel do

      dzeta_iny = 1.d0/dzeta
      deallocate(dzeta)
    endsubroutine lattice_y
    !作成したdx/dζをdF/dyなどに掛けて微分変換を行うsubroutine
    subroutine combine_x(dzeta_in,dF,dFzeta)
      integer i,j,k,l
      double precision,allocatable,dimension(:):: dzeta_in
      double precision,allocatable,dimension(:,:,:,:):: dFzeta,dF
      !$omp parallel do
      do k=0,Nz-1
       do i=0,Ny
         do j=0,Nx
           do l=0,4
             dFzeta(l,j,i,k) = dF(l,j,i,k) * dzeta_in(j)
           enddo
         end do
       end do
      end do
      !$omp end parallel do
    endsubroutine combine_x
    subroutine combine_y(dzeta_in,dF,dFzeta)
      integer i,j,k,l
      double precision,allocatable,dimension(:):: dzeta_in
      double precision,allocatable,dimension(:,:,:,:):: dFzeta,dF
      !$omp parallel do
      do k=0,Nz-1
       do i=0,Ny
         do j=0,Nx
           do l=0,4
             dFzeta(l,j,i,k) = dF(l,j,i,k) * dzeta_in(i)
           enddo
         end do
       end do
      end do
      !$omp end parallel do
    endsubroutine combine_y
end module three_omp

    program main
      use three_omp
      implicit none
      character(len = 16) filename
      character(len = 16) z_name
      !時間更新毎に出力ファイルを変更するためのファイル名設定
      !NSCBCでdrho,du,dv,dw,dpを計算する必要があるので計算しやすいようにG行列を作成しG(rho,u,v,w,p)
      !となるように設定した。すなわちdGは密度、速度、圧力のそれぞれの微分項を含む行列となる
      !x方向
      !Runge-Kutta法のためQ1,Q2用の値を設定
      double precision,allocatable,dimension(:,:,:,:) :: G,Q,Q0,Q1,Q2,Qn,Fpx,Fmx,xp,xm,oldG
      !y方向
      double precision,allocatable,dimension(:,:,:,:) :: Fpy,Fmy,yp,ym
      !z方向
      double precision,allocatable,dimension(:,:,:,:) :: Fpz,Fmz,zp,zm
      double precision,allocatable,dimension(:,:,:) :: myu
      !A,L,U,x,y,bはF+とF-用に分けるため別々で定義
      !x方向
      double precision,allocatable,dimension(:,:) :: LUmx,LUpx
      !y方向
      double precision,allocatable,dimension(:,:) :: LUmy,LUpy
      !z方向
      double precision,allocatable,dimension(:,:) :: LUmz,LUpz
      !粘性項の計算に使う行列(x方向)
      double precision,allocatable,dimension(:,:,:,:) :: Vx,dVx,UVWT,dUVWTx
      double precision,allocatable,dimension(:,:) :: LUccsx
      !y方向
      double precision,allocatable,dimension(:,:,:,:) :: Vy,dVy,dUVWTy
      double precision,allocatable,dimension(:,:) :: LUccsy
      !z方向
      double precision,allocatable,dimension(:,:,:,:) :: Vz,dVz,dUVWTz
      double precision,allocatable,dimension(:,:) :: LUccsz
      ! NSCBC用
      double precision,allocatable,dimension(:,:) :: in_G0,in_G1_top,in_G2,in_G3
      ! double precision,allocatable,dimension(:,:) :: in_G1_du
      !x方向
      double precision,allocatable,dimension(:,:,:,:) :: dGx,dFx
      double precision  dx,pNx_infty
      !y方向
      double precision,allocatable,dimension(:,:,:,:) :: dGy,dFy
      double precision dy,p0y_infty,pNy_infty
      !z方向
      double precision,allocatable,dimension(:,:,:,:) :: dFz,dGz
      double precision dz,z
      integer i,j,k,l,M,ii,jj,kk
      ! double precision theta !時間での周期的撹乱用の変数
      double precision c_infty
      double precision,allocatable,dimension(:) :: ur,Tu
      double precision,allocatable,dimension(:) :: Ux,sigma_x,Uy,sigma_y
      double precision,allocatable,dimension(:,:,:,:) :: dQx,dQy
      double precision,allocatable,dimension(:) :: dzeta_iny,dzeta_inx
      double precision,allocatable,dimension(:) :: zeta_fx,zeta_fy
      double precision,allocatable,dimension(:,:,:) :: omega_1,omega_2,omega_3,dp!渦度と圧力変動差を入れる配列
      double precision,allocatable,dimension(:,:,:) :: div_u,Invariant_2 !音響成分と渦構造(第二不変量)を入れる配列
      double precision,allocatable,dimension(:,:) :: kakuran_u,kakuran_v,kakuran_w!ランダム撹乱を入れる配列

      allocate(G(0:4,0:Nx,0:Ny,0:Nz-1),Q(0:4,0:Nx,0:Ny,0:Nz-1),Q0(0:4,0:Nx,0:Ny,0:Nz-1)&
      ,Q1(0:4,0:Nx,0:Ny,0:Nz-1),Q2(0:4,0:Nx,0:Ny,0:Nz-1),Qn(0:4,0:Nx,0:Ny,0:Nz-1)&
      ,Fpx(0:4,0:Nx,0:Ny,0:Nz-1),Fmx(0:4,0:Nx,0:Ny,0:Nz-1),xp(0:4,0:Nx,0:Ny,0:Nz-1)&
      ,xm(0:4,0:Nx,0:Ny,0:Nz-1),oldG(0:4,0:Nx,0:Ny,0:Nz-1))

      allocate(Fpy(0:4,0:Nx,0:Ny,0:Nz-1),Fmy(0:4,0:Nx,0:Ny,0:Nz-1),&
      yp(0:4,0:Nx,0:Ny,0:Nz-1),ym(0:4,0:Nx,0:Ny,0:Nz-1))

      allocate(Fpz(0:4,0:Nx,0:Ny,0:Nz-1),Fmz(0:4,0:Nx,0:Ny,0:Nz-1),&
      zp(0:4,0:Nx,0:Ny,0:Nz-1),zm(0:4,0:Nx,0:Ny,0:Nz-1))

      allocate(myu(0:Nx,0:Ny,0:Nz-1))
      allocate(Vx(0:4,0:Nx,0:Ny,0:Nz-1),dVx(0:4,0:Nx,0:Ny,0:Nz-1),&
      UVWT(0:4,0:Nx,0:Ny,0:Nz-1),dUVWTx(0:4,0:Nx,0:Ny,0:Nz-1))

      allocate(Vy(0:4,0:Nx,0:Ny,0:Nz-1),dVy(0:4,0:Nx,0:Ny,0:Nz-1)&
      ,dUVWTy(0:4,0:Nx,0:Ny,0:Nz-1))

      allocate(Vz(0:4,0:Nx,0:Ny,0:Nz-1),dVz(0:4,0:Nx,0:Ny,0:Nz-1),&
      dUVWTz(0:4,0:Nx,0:Ny,0:Nz-1))

      allocate(in_G0(0:Ny,0:Nz-1),in_G1_top(0:Ny,0:Nz-1),in_G2(0:Ny,0:Nz-1),in_G3(0:Ny,0:Nz-1))
      ! allocate(in_G1_du(0:Ny,0:Nz-1))
      allocate(dGx(0:4,0:Nx,0:Ny,0:Nz-1),dFx(0:4,0:Nx,0:Ny,0:Nz-1))
      allocate(dGy(0:4,0:Nx,0:Ny,0:Nz-1),dFy(0:4,0:Nx,0:Ny,0:Nz-1))
      allocate(dFz(0:4,0:Nx,0:Ny,0:Nz-1),dGz(0:4,0:Nx,0:Ny,0:Nz-1))
      allocate(Ux(0:Nx),sigma_x(0:Nx),Uy(0:Ny),sigma_y(0:Ny),&
      dQx(0:4,0:Nx,0:Ny,0:Nz-1),dQy(0:4,0:Nx,0:Ny,0:Nz-1))

      allocate(dzeta_inx(0:Nx),dzeta_iny(0:Ny))
      allocate(omega_1(0:Nx,0:Ny,0:Nz-1),omega_2(0:Nx,0:Ny,0:Nz-1),&
      omega_3(0:Nx,0:Ny,0:Nz-1),dp(0:Nx,0:Ny,0:Nz-1),div_u(0:Nx,0:Ny,0:Nz-1),&
      Invariant_2(0:Nx,0:Ny,0:Nz-1))
      allocate(kakuran_u(0:Ny,0:Nz-1),kakuran_v(0:Ny,0:Nz-1),&
      kakuran_w(0:Ny,0:Nz-1))

      allocate(zeta_fx(0:Nx),zeta_fy(0:Ny))
      allocate(ur(0:Ny),Tu(0:Ny))
      !x_axis
      allocate(LUmx(-1:1,0:Nx),LUpx(-1:1,0:Nx))
      allocate(LUccsx(-1:1,0:Nx))
      !y_axis
      allocate(LUmy(-1:1,0:Ny),LUpy(-1:1,0:Ny))
      allocate(LUccsy(-1:1,0:Ny))
      !z_axis
      allocate(LUmz(-2:2,0:Nz-1),LUpz(-2:2,0:Nz-1))
      allocate(LUccsz(-2:2,0:Nz-1))
      dx = Lx /dble(Nx)
      dy = Ly /dble(Ny)
      dz = Lz / dble(Nz)
!一応ゼロクリア
      G=0.d0;Q=0.d0;Qn=0.d0;Q0=0.d0;Q1=0.d0;Q2=0.d0
      pNx_infty=0.d0;p0y_infty=0.d0;pNy_infty=0.d0;ur=0.d0;Tu=0.d0
      in_G0=0.d0;in_G1_top=0.d0;in_G2=0.d0;in_G3=0.d0
      ! in_G1_du=0.d0
      Ux=0.d0;sigma_x=0.d0;Uy=0.d0;sigma_y=0.d0;zeta_fy=0.d0;dzeta_iny=0.d0
      zeta_fx=0.d0;dzeta_inx=0.d0;omega_1=0.d0;omega_2=0.d0;omega_3=0.d0;dp=0.d0;oldG=0.d0
      div_u=0.d0;Invariant_2=0.d0;kakuran_u=0.d0;kakuran_v=0.d0;kakuran_w=0.d0

      !============座標設定======================================================
      !y方向の格子伸長のための座標設定
      call lattice_y(dy,zeta_fy,dzeta_iny)
      !x方向も
      call lattice_x(dx,zeta_fx,dzeta_inx)
      !=========================================================================
  !!!!!!============流入条件設定==================================================
    !top-hat型ジェットの導出・計算
    !初期条件もζ_yの座標系で設定する
      ur(Ny/2) = ujet
      Tu(Ny/2) = Tjet

    !========並列化しない=========================================================
    !Tuを求めるのにurが必要なため1つのDoループで並列化不可能。ここはDoループが小さいので、並列化しない
      do i = (Ny/2)+1,Ny
        !Top-hat型の分布になるような式を設定
        ur(i) = ujet/2.d0*(1.d0 - dtanh((12.5d0/4.d0)*((zeta_fy(i)/b)- (b/zeta_fy(i)))))
        !Crocco-Busemanの関係式と主流速度分布を用いて温度分布Tuを設定
        Tu(i) = Ma**2.d0*(gamma-1.d0)/2.d0*(ur(i)*ujet-ur(i)**2.d0)/ujet+&
                Tjet*ur(i)/ujet+Temp*(ujet-ur(i))/ujet
      enddo
      !========並列化しない========================================================
    !初期分布をx軸対象になるようにする。
    !そのためにy軸正の範囲の値を負の範囲に軸対象になるようにコピーする
    !$omp parallel do
      do i = 0,(Ny/2)-1
        ur(i) = ur(Ny-i)
        Tu(i) = Tu(Ny-i)
      enddo
    !$omp end parallel do
    !===============================流入ジェットの計算終了===========================
    !==============================流入ランダム撹乱の読み込み========================
      ! open(31,file='dirturbance_conditions/kakuran3D_u.txt',status='old')
      open(32,file='dirturbance_conditions/kakuran3D_v.txt',status='old')
      open(33,file='dirturbance_conditions/kakuran3D_w.txt',status='old')
      !=====読み込みは順番が大切だろうから、並列化しない==================================
      do k=0,Nz-1
        do i=0,Ny
          ! read(31,*) kakuran_u(i,k)
          read(32,*) kakuran_v(i,k)
          read(33,*) kakuran_w(i,k)
        enddo
      enddo
      !=====読み込みは順番が大切だろうから、並列化しない==================================
      ! close(31)
      close(32)
      close(33)

      ! open(34,file='result_omp/kakkuran_kakunin.txt',status='replace')
      ! do k=0,Nz-1
      !   z = dz*dble(k)
      !   do i=0,Ny
      !     write(34,'(3f24.16)') zeta_fy(i),z,kakuran_v(i,k)
      !   enddo
      !   write(34,*)
      ! enddo
      ! close(34)
   !============================================================================
   !Bufferの計算のための初期値を用いて無限遠方での音速を定義
   c_infty = sqrt(Temp/Ma**2.d0)
   !Buffer領域の計算に使うUx,Uy,sigma_x,sigma_yの計算
   call buffer_x(c_infty,Ux,sigma_x,zeta_fx)
   call buffer_y(c_infty,Uy,sigma_y,zeta_fy)
 !==========初期値の設定==========================================================
 !流入条件
 !x=0の軸上にのみ流入条件を適用することでここからどんどん流入が起こる
 !矩型JetをZ方向 k=8~10のみに流入させる
 !ランダム撹乱の窓関数はtop-hat型ジェットの関数をそのまま使用
 !これにより、ランダム撹乱は完全にジェットの中にのみ、存在する
!$omp parallel do
 do k =0,Nz-1
   do i=0,Ny
     in_G0(i,k) = 1.d0/Tu(i)!密度ρは理想気体状態方程式に従うから
     in_G1_top(i,k) = ur(i)
     !本当は窓関数を別で適用するが、今回は窓関数=ur(i)なので、計算量削減のためそのまま適用
     ! in_G1_du = dis_strength*kakuran_u*窓関数
     ! in_G1_du(i,k) = dis_strength*kakuran_u(i,k)*ur(i)
     in_G2(i,k) = dis_strength*kakuran_v(i,k)*ur(i)
     in_G3(i,k) = dis_strength*kakuran_w(i,k)*ur(i)
   end do
 enddo
!$omp end parallel do
 !まず最初に、流入条件をGのx=0の場所にのみ適用する
 !G(1,2,3)に撹乱入れない。時間進行によって徐々に強くするから初期条件で撹乱は0
 !$omp parallel do
  do k =0,Nz-1
    do i=0,Ny
      !Crocco-Busemannの関係式より
      G(0,0,i,k) = in_G0(i,k)!ρ
      !top-hat Jetのみ
      G(1,0,i,k) = in_G1_top(i,k)!u
      G(4,0,i,k) = 1.d0*Temp/((Ma**2.d0)*gamma)!p
    enddo
  enddo
 !$omp end parallel do

 !次にx=1~Nxの残りの部分に一括で、初期値を入れる
!=====大きいDoループなので並列化する=================================================
  !$omp parallel do
   do k=0,Nz-1
     do i=0,Ny
       do j=1,Nx
         G(0,j,i,k) = 1.d0!ρ
         !以下は上で0クリアしているので、再度やるのは無駄
         ! G(1,j,i,k) = 0.d0!u
         ! G(2,j,i,k) = 0.d0!v
         ! G(3,j,i,k) = 0.d0!w
         !G(0)は単なる定数なので、分けて並列化するよりも、数値で代入してしまった方が計算量は少ない
         G(4,j,i,k) = 1.d0*Temp/((Ma**2.d0)*gamma)!p
         ! G(4,j,i,k) = G(0,j,i,k)*Temp/((Ma**2.d0)*gamma)!p
       end do
     enddo
   enddo
  !$omp end parallel do
  !=====大きいDoループなので並列化する==============================================
  !初期値の出力
  !まずt=0はループ外で個別に作成
  !もちろん出力もζ_y座標系とζ_x座標系で行う
  !=======ファイルへの書き出しはもちろん順番が大切なので、並列化不可能====================
       do k=0,Nz-1
         z = dz*dble(k)
         write(z_name, '(i2.2)') k
         open(10, file = "result_omp/parameter000000_"//trim(z_name)//".txt")
          ! z = dz*dble(Nz/2)
          do i = 0,Ny
            do j = 0,Nx
              write(10,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
              &f24.16,",",f24.16)') zeta_fx(j),zeta_fy(i),z,&
              G(0,j,i,k),div_u(j,i,k),Invariant_2(j,i,k),dp(j,i,k)/dt
            enddo
            write(10,*)
          enddo
          close(10)
       enddo
   !=======ファイルへの書き出しはもちろん順番が大切なので、並列化不可能=======================
   open(41, file = "result_omp/turbulent_check_1.csv")
   open(42, file = "result_omp/turbulent_check_2.csv")
   open(43, file = "result_omp/turbulent_check_3.csv")
   open(44, file = "result_omp/turbulent_check_4.csv")


      !p_inftyの定義
      pNx_infty = G(4,Nx,0,0)
      p0y_infty = G(4,0,0,0)
      pNy_infty = G(4,0,Ny,0)
        !粘性項の計算はCCSを用いるためA,L,U行列がsigma=0となる
        !そのためAp,Amなどとはまた別に設定する
        !CCS用のA,L,U行列はNの値がNx,Nyで異なるので別々に設定する
        call LU_DecompoNonP(Nx,ccs_sigma,LUccsx)
        call LU_DecompoNonP(Ny,ccs_sigma,LUccsy)
        call LU_DecompoPiriodic(Nz,ccs_sigma,LUccsz)
        !A行列の設定をF+とF-のそれぞれで行う
        !LU分解の過程は同じだが用いる値が異なるので別々にcall
        call LU_DecompoNonP(Nx,psigma,LUpx)
        call LU_DecompoNonP(Ny,psigma,LUpy)
        call LU_DecompoPiriodic(Nz,psigma,LUpz)
        call LU_DecompoNonP(Nx,msigma,LUmx)
        call LU_DecompoNonP(Ny,msigma,LUmy)
        call LU_DecompoPiriodic(Nz,msigma,LUmz)
        call Q_matrix(G,Q)
        Q0 = Q!Bufferの計算で使う初期値のQを保存
      !NSCBCのdGを求める微分も粘性項と同じCCSなのでA,L,U行列の使い回し可能
      do M =1,Mmax
        Fpx=0.d0;Fmx=0.d0;xp=0.d0;xm=0.d0;Fpy=0.d0;Fmy=0.d0;yp=0.d0;ym=0.d0;Fpz=0.d0;Fmz=0.d0;zp=0.d0;zm=0.d0
        UVWT=0.d0;dUVWTx=0.d0;dUVWTy=0.d0;dUVWTz=0.d0;dVx=0.d0;dVy=0.d0;dVz=0.d0;dQx=0.d0;dQy=0.d0;myu=0.d0
        Vx=0.d0;Vy=0.d0;Vz=0.d0;dGx=0.d0;dGy=0.d0;dGz=0.d0;dFx=0.d0;dFy=0.d0;dFz=0.d0
        !x(i)を導出するのにf(i)を用いるがそれが時間ごとに変化するため幾つかのsubroutineは
        !このDoループ内に入れる必要がある。
            !移流方程式の導出
            !3次精度Runge-Kutta法での導出！
        !Q1
        !F行列のdfx/dxの計算
        !x_axis
        call F_matrix(Q,Fpx,Fmx,Fpy,Fmy,Fpz,Fmz)
        call dif_x(psigma,dx,Fpx,xp,LUpx,dzeta_inx)
        call dif_x(msigma,dx,Fmx,xm,LUmx,dzeta_inx)
        !y_axis dFy/dyの計算
        call dif_y(psigma,dy,Fpy,yp,LUpy,dzeta_iny)
        call dif_y(msigma,dy,Fmy,ym,LUmy,dzeta_iny)
        !z_axis
        call dif_z(psigma,dz,Fpz,zp,LUpz)
        call dif_z(msigma,dz,Fmz,zm,LUmz)
        !NSCBCの境界条件を適用させるためにQ1を求める前にdFを定義してその両端に境界条件を
        !適用しなければならない
        !$omp parallel do
          do k=0,Nz-1
            do i=0,Ny
              do j=0,Nx
                do l=0,4
                dFx(l,j,i,k) = xm(l,j,i,k) + xp(l,j,i,k)
                dFy(l,j,i,k) = ym(l,j,i,k) + yp(l,j,i,k)
                dFz(l,j,i,k) = zm(l,j,i,k) + zp(l,j,i,k)
                end do
              end do
            enddo
          enddo
        !$omp end parallel do
        !粘性項V行列のdv/dxの計算
        !まずはdu/dx,dT/dxの導出とμの設定
        call variable_setting(UVWT,Q,myu)
        !x_axis dUVWT/dxの計算
        call dif_x(ccs_sigma,dx,UVWT,dUVWTx,LUccsx,dzeta_inx)
        !y_axis dUVWT/dyの計算
        call dif_y(ccs_sigma,dy,UVWT,dUVWTy,LUccsy,dzeta_iny)
        !z_axis dUVWT/dzの計算
        call dif_z(ccs_sigma,dz,UVWT,dUVWTz,LUccsz)
        !V行列を構成する値が揃ったのでV行列の設定とdV/dx,y,zの計算
        call V_matrix(Vx,Vy,Vz,myu,UVWT,dUVWTx,dUVWTy,dUVWTz)
        !x_axis
        call dif_x(ccs_sigma,dx,Vx,dVx,LUccsx,dzeta_inx)
        !y_axis
        call dif_y(ccs_sigma,dy,Vy,dVy,LUccsy,dzeta_iny)
        !z_axis
        call dif_z(ccs_sigma,dz,Vz,dVz,LUccsz)
        !NSCBCの計算開始
        !x方向のNSCBCの計算
        call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
        call NSCBC_x_0_super(dFx)
        ! call NSCBC_x_0_sub(G,dGx,dFx)
        call NSCBC_x_Nx_super(G,dGx,dFx)
        ! call NSCBC_x_Nx_sub(G,dGx,dFx,pNx_infty)
        call outflow_x(UVWT,dUVWTx,Vx,dVx)
        !y方向
        call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
        call NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
        !無反射流出条件の際の境界での粘性項の条件を設定
        call outflow_y(UVWT,dUVWTy,Vy,dVy)
        !Buffer領域の計算
        !計算に必要なdQ/dx,dQ/dyをCCSで導出。今まではdQ/dtしか求めていなかった
        !dQ/dx,dQ/dy自体はdG/dx,dG/dyを組み合わせて作ることができるのでそうして作成すると微分せずに済み計算の短縮に繋がる
        !しかしEtの微分だけは足し算が多くなってしまい遅くなるのでEtはそれで新規に微分、Qのそれ以外の値は組み合わせて作るのが最善
        !だがこの場合Etは単体で定義していないので今回はそれができない。そのためQを直接微分する
        call dif_x(ccs_sigma,dx,Q,dQx,LUccsx,dzeta_inx)
        call dif_y(ccs_sigma,dy,Q,dQy,LUccsy,dzeta_iny)
       !計算して求めたdF,dVそしてBuffer領域の計算のための値などを組み合わせ、代入してdQ/dtを求める
     !$omp parallel do
       do k=0,Nz-1
        do i=0,Ny
          do j=0,Nx
            do l=0,4
            Q1(l,j,i,k) = Q(l,j,i,k) + c*dt*(dVx(l,j,i,k)+dVy(l,j,i,k)+dVz(l,j,i,k)-dFx(l,j,i,k)&
              &-dFy(l,j,i,k)-dFz(l,j,i,k)-sigma_x(j)*(Q(l,j,i,k)-Q0(l,j,i,k))&
              -sigma_y(i)*(Q(l,j,i,k)-Q0(l,j,i,k))-Ux(j)*dQx(l,j,i,k)-Uy(i)*dQy(l,j,i,k))
            end do
          end do
        enddo
      enddo
    !$omp end parallel do
      !call Q_boundary(Q1)
      !i=0で流入条件させるのでその部分のQ1を上書きして流入させ続ける
      ! call inflow(M,Q1,in_G0,in_G1_top,in_G1_du,in_G2,in_G3)!dirichlet条件で流入部を固定
      call inflow(M,Q1,in_G0,in_G1_top,in_G2,in_G3)!dirichlet条件で流入部を固定
      !Q2(Q,F,x+-,y+-,f+-はそれぞれの計算過程において分ける必要がある。
      !またL,Uなどは DCSという方法が変わらないので同じものを使用できる)
      !dF/dxの計算
      Fpx=0.d0;Fmx=0.d0;xp=0.d0;xm=0.d0;Fpy=0.d0;Fmy=0.d0;yp=0.d0;ym=0.d0;Fpz=0.d0;Fmz=0.d0;zp=0.d0;zm=0.d0
      UVWT=0.d0;dUVWTx=0.d0;dUVWTy=0.d0;dUVWTz=0.d0;dVx=0.d0;dVy=0.d0;dVz=0.d0;dQx=0.d0;dQy=0.d0;myu=0.d0
      Vx=0.d0;Vy=0.d0;Vz=0.d0;dGx=0.d0;dGy=0.d0;dGz=0.d0;dFx=0.d0;dFy=0.d0;dFz=0.d0
      call F_matrix(Q1,Fpx,Fmx,Fpy,Fmy,Fpz,Fmz)
      call dif_x(psigma,dx,Fpx,xp,LUpx,dzeta_inx)
      call dif_x(msigma,dx,Fmx,xm,LUmx,dzeta_inx)
      !y_axis dFy/dyの計算
      call dif_y(psigma,dy,Fpy,yp,LUpy,dzeta_iny)
      call dif_y(msigma,dy,Fmy,ym,LUmy,dzeta_iny)
      !z_axis dFz/dzの計算
      call dif_z(psigma,dz,Fpz,zp,LUpz)
      call dif_z(msigma,dz,Fmz,zm,LUmz)
      !NSCBCの境界条件を適用させるためにQ1を求める前にdFを定義してその両端に境界条件を
      !適用しなければならない
      !$omp parallel do
        do k=0,Nz-1
          do i=0,Ny
            do j=0,Nx
              do l=0,4
              dFx(l,j,i,k) = xm(l,j,i,k) + xp(l,j,i,k)
              dFy(l,j,i,k) = ym(l,j,i,k) + yp(l,j,i,k)
              dFz(l,j,i,k) = zm(l,j,i,k) + zp(l,j,i,k)
              end do
            end do
          enddo
        enddo
      !$omp end parallel do

        !粘性項V行列のdv/dxの計算
        !まずはdu/dx,dT/dxの導出とμの設定
        call variable_setting(UVWT,Q1,myu)
        !x_axis dUVWT/dxの計算
        call dif_x(ccs_sigma,dx,UVWT,dUVWTx,LUccsx,dzeta_inx)
        !y_axis dUVWT/dyの計算
        call dif_y(ccs_sigma,dy,UVWT,dUVWTy,LUccsy,dzeta_iny)
        !z_axis dUVWT/dzの計算
        call dif_z(ccs_sigma,dz,UVWT,dUVWTz,LUccsz)
        !V行列を構成する値が揃ったのでV行列の設定とdV/dxの計算
        call V_matrix(Vx,Vy,Vz,myu,UVWT,dUVWTx,dUVWTy,dUVWTz)
        !x_axis
        call dif_x(ccs_sigma,dx,Vx,dVx,LUccsx,dzeta_inx)
        !y_axis
        call dif_y(ccs_sigma,dy,Vy,dVy,LUccsy,dzeta_iny)
        !z_axis
        call dif_z(ccs_sigma,dz,Vz,dVz,LUccsz)
        !NSCBCの計算開始
        call rho_u_p(G,Q1)
        !x方向のNSCBCの計算
        call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
        call NSCBC_x_0_super(dFx)!超音速流入
        ! call NSCBC_x_0_sub(G,dGx,dFx)!亜音速流入
        call NSCBC_x_Nx_super(G,dGx,dFx)
        ! call NSCBC_x_Nx_sub(G,dGx,dFx,pNx_infty)
        call outflow_x(UVWT,dUVWTx,Vx,dVx)
        !y方向
        call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
        call NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
        call outflow_y(UVWT,dUVWTy,Vy,dVy)
        !Buffer領域の計算
        call dif_x(ccs_sigma,dx,Q1,dQx,LUccsx,dzeta_inx)
        call dif_y(ccs_sigma,dy,Q1,dQy,LUccsy,dzeta_iny)
      !$omp parallel do
        do k=0,Nz-1
         do i=0,Ny
           do j=0,Nx
             do l=0,4
              Q2(l,j,i,k) = (0.75d0)*Q(l,j,i,k) +(0.25d0) * Q1(l,j,i,k) + c* (dt*0.25d0)&
                  &*(dVx(l,j,i,k)+dVy(l,j,i,k)+dVz(l,j,i,k)-dFx(l,j,i,k)-dFy(l,j,i,k)&
                  -dFz(l,j,i,k)-sigma_x(j)*(Q1(l,j,i,k)-Q0(l,j,i,k))&
                  -sigma_y(i)*(Q1(l,j,i,k)-Q0(l,j,i,k))-Ux(j)*dQx(l,j,i,k)-Uy(i)*dQy(l,j,i,k))
               end do
             end do
           enddo
         enddo
       !$omp end parallel do

!        call Q_boundary(Q2)
        ! call inflow(M,Q2,in_G0,in_G1_top,in_G1_du,in_G2,in_G3)
        call inflow(M,Q2,in_G0,in_G1_top,in_G2,in_G3)
      !Qn
      !dF/dxの計算
      Fpx=0.d0;Fmx=0.d0;xp=0.d0;xm=0.d0;Fpy=0.d0;Fmy=0.d0;yp=0.d0;ym=0.d0;Fpz=0.d0;Fmz=0.d0;zp=0.d0;zm=0.d0
      UVWT=0.d0;dUVWTx=0.d0;dUVWTy=0.d0;dUVWTz=0.d0;dVx=0.d0;dVy=0.d0;dVz=0.d0;dQx=0.d0;dQy=0.d0;myu=0.d0
      Vx=0.d0;Vy=0.d0;Vz=0.d0;dGx=0.d0;dGy=0.d0;dGz=0.d0;dFx=0.d0;dFy=0.d0;dFz=0.d0
      call F_matrix(Q2,Fpx,Fmx,Fpy,Fmy,Fpz,Fmz)
      call dif_x(psigma,dx,Fpx,xp,LUpx,dzeta_inx)
      call dif_x(msigma,dx,Fmx,xm,LUmx,dzeta_inx)
      !y_axis dFy/dyの計算
      call dif_y(psigma,dy,Fpy,yp,LUpy,dzeta_iny)
      call dif_y(msigma,dy,Fmy,ym,LUmy,dzeta_iny)
      !z_axis dFz/dzの計算
      call dif_z(psigma,dz,Fpz,zp,LUpz)
      call dif_z(msigma,dz,Fmz,zm,LUmz)
      !NSCBCの境界条件を適用させるためにQ2を求める前にdFを定義してその両端に境界条件を
      !適用しなければならない
      !$omp parallel do
        do k=0,Nz-1
          do i=0,Ny
            do j=0,Nx
              do l=0,4
              dFx(l,j,i,k) = xm(l,j,i,k) + xp(l,j,i,k)
              dFy(l,j,i,k) = ym(l,j,i,k) + yp(l,j,i,k)
              dFz(l,j,i,k) = zm(l,j,i,k) + zp(l,j,i,k)
              end do
            end do
          enddo
        enddo
      !$omp end parallel do

        !粘性項V行列のdv/dxの計算
        !まずはdu/dx,dT/dxの導出とμの設定
        call variable_setting(UVWT,Q2,myu)
        !x_axis dUVWT/dxの計算
        call dif_x(ccs_sigma,dx,UVWT,dUVWTx,LUccsx,dzeta_inx)
        !y_axis dUVWT/dyの計算
        call dif_y(ccs_sigma,dy,UVWT,dUVWTy,LUccsy,dzeta_iny)
        !z_axis dUVWT/dzの計算
        call dif_z(ccs_sigma,dz,UVWT,dUVWTz,LUccsz)
        !V行列を構成する値が揃ったのでV行列の設定とdV/dxの計算
        call V_matrix(Vx,Vy,Vz,myu,UVWT,dUVWTx,dUVWTy,dUVWTz)
        !x_axis
        call dif_x(ccs_sigma,dx,Vx,dVx,LUccsx,dzeta_inx)
        !y_axis
        call dif_y(ccs_sigma,dy,Vy,dVy,LUccsy,dzeta_iny)
        !z_axis
        call dif_z(ccs_sigma,dz,Vz,dVz,LUccsz)
        !NSCBCの計算開始
        call rho_u_p(G,Q2)
        !=====================================================
        !流入を全てのZ座標から行う場合はinflowのsubroutineで、
        !密度すらもDirichlet条件で固定しているのでx方向のNSCBCは全て不要
        !=====================================================
        !x方向のNSCBCの計算
        call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
        call NSCBC_x_0_super(dFx)
        ! call NSCBC_x_0_sub(G,dGx,dFx)
        call NSCBC_x_Nx_super(G,dGx,dFx)
        ! call NSCBC_x_Nx_sub(G,dGx,dFx,pNx_infty)
        call outflow_x(UVWT,dUVWTx,Vx,dVx)
        !y方向
        call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
        call NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
        call outflow_y(UVWT,dUVWTy,Vy,dVy)
        !Buffer領域の計算
        call dif_x(ccs_sigma,dx,Q2,dQx,LUccsx,dzeta_inx)
        call dif_y(ccs_sigma,dy,Q2,dQy,LUccsy,dzeta_iny)

      !$omp parallel do
        do k=0,Nz-1
         do i=0,Ny
           do j=0,Nx
             do l=0,4
               Qn(l,j,i,k)=Q(l,j,i,k)/3.d0+(2.d0/3.d0)*Q2(l,j,i,k)+c*&
               &((2.d0*dt)/3.d0)*(dVx(l,j,i,k)+dVy(l,j,i,k)+dVz(l,j,i,k)-dFx(l,j,i,k)&
               -dFy(l,j,i,k)-dFz(l,j,i,k)-sigma_x(j)*(Q2(l,j,i,k)-Q0(l,j,i,k))&
               -sigma_y(i)*(Q2(l,j,i,k)-Q0(l,j,i,k))-Ux(j)*dQx(l,j,i,k)-Uy(i)*dQy(l,j,i,k))
               end do
             end do
           enddo
         enddo
       !$omp end parallel do

!        call Q_boundary(Qn)
        ! call inflow(M,Qn,in_G0,in_G1_top,in_G1_du,in_G2,in_G3)
        call inflow(M,Qn,in_G0,in_G1_top,in_G2,in_G3)
        call rho_u_p(G,Qn)
        if((M >= observe_start_time).and.(observe_end_time >= M)) then
          call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
          write(41,'(f24.16)') dGx(1,5*Nx/6,Ny/2,Nz/2)
          write(42,'(f24.16)') dGx(1,5*Nx/6,Ny/4,Nz/4)
          write(43,'(f24.16)') dGx(1,5*Nx/6,3*Ny/4,3*Nz/4)
          write(44,'(f24.16)') dGx(1,3*Nx/4,3*Ny/4,3*Nz/4)
        endif

        if(mod(M,output_count) == 0) then!dt=1.d-4で0.01秒刻みで出力するためにMの条件を設定
         !渦度用のdGの計算
         dGx=0.d0;dGy=0.d0;dGz=0.d0
         call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
         call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
         call dif_z(ccs_sigma,dz,G,dGz,LUccsz)

         !$omp parallel do
           do k=0,Nz-1
             do i=0,Ny
               do j=0,Nx
                 !音響成分div u
                 div_u(j,i,k) = dGx(1,j,i,k)+dGy(2,j,i,k)+dGz(3,j,i,k)
                 !速度勾配テンソルの第二不変量Q
                 Invariant_2(j,i,k) = dGx(1,j,i,k)*dGy(2,j,i,k)+&
                     dGy(2,j,i,k)*dGz(3,j,i,k)+dGz(3,j,i,k)*dGx(1,j,i,k)-&
                     (dGz(1,j,i,k)*dGx(3,j,i,k)+dGy(3,j,i,k)*dGz(2,j,i,k)+&
                     dGx(2,j,i,k)*dGy(1,j,i,k))
                 ! omega_1(j,i,k) = dGy(3,j,i,k) - dGz(2,j,i,k)!dw/dy-dv/dz
                 ! omega_2(j,i,k) = dGz(1,j,i,k) - dGx(3,j,i,k)!du/dz-dw/dx
                 ! omega_3(j,i,k) = dGx(2,j,i,k) - dGy(1,j,i,k)!dv/dx-du/dy
               end do
             enddo
           enddo
         !$omp end parallel do

         call rho_u_p(oldG,Q)
         !$omp parallel do
           do k=0,Nz-1
             do i=0,Ny
               do j=0,Nx
                 dp(j,i,k) = (G(4,j,i,k) - oldG(4,j,i,k))/dt
               end do
             enddo
           enddo
         !$omp end parallel do

           write(filename, '(i6.6)') M
           !Mの計算毎に出力ファイル名を変更して出力する
           !i5.5で5桁分の数字を表示できるのでdt=1.d-5以下で計算するならここも変更が必要
      !=======ファイルへの書き出しはもちろん順番が大切なので、並列化不可能====================
           do kk= 0,Nz-1
             write(z_name, '(i2.2)') kk
             open(10, file = "result_omp/parameter"//trim(filename)//"_"//trim(z_name)//".txt")
             z=dz*dble(kk)
             do ii = 0,Ny
               do jj = 0,Nx
                 write(10,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
                 &f24.16,",",f24.16)') zeta_fx(jj),zeta_fy(ii),z,&
                 G(0,jj,ii,kk),div_u(jj,ii,kk),Invariant_2(jj,ii,kk),dp(jj,ii,kk)
               enddo
               write(10,*)
               !一度に全てを出力する際にはデータの切れ目として空白を一行挿入しなくてはいけない
             enddo
             close(10)
           enddo
     !=======ファイルへの書き出しはもちろん順番が大切なので、並列化不可能====================
         endif
        !計算が破綻している場合に計算を終了させるプログラム
      !=======最初にNaNでなっている場所を特定したいので、並列化不可能======================
        do k = 0,Nz-1
            do j = 0,Ny
              do i = 0,Nx
                  if(isnan(Qn(0,i,j,k))) then
                    !渦度用のdGの計算
                    oldG=0;dGx=0.d0;dGy=0.d0;dGz=0.d0
                    call rho_u_p(oldG,Q)
                    call dif_x(ccs_sigma,dx,oldG,dGx,LUccsx,dzeta_inx)
                    call dif_y(ccs_sigma,dy,oldG,dGy,LUccsy,dzeta_iny)
                    call dif_z(ccs_sigma,dz,oldG,dGz,LUccsz)

                    !$omp parallel do
                    do kk=0,Nz-1
                      do ii=0,Ny
                        do jj=0,Nx
                          !音響成分du/dx+dv/dy+dw/dz
                          div_u(jj,ii,kk) =dGx(1,jj,ii,kk)+dGy(2,jj,ii,kk)+dGz(3,jj,ii,kk)
                          Invariant_2(jj,ii,kk) = dGx(1,jj,ii,kk)*dGy(2,jj,ii,kk)+&
                          dGy(2,jj,ii,kk)*dGz(3,jj,ii,kk)+dGz(3,jj,ii,kk)*dGx(1,jj,ii,kk)-&
                          (dGz(1,jj,ii,kk)*dGx(3,jj,ii,kk)+dGy(3,jj,ii,kk)*dGz(2,jj,ii,kk)+&
                          dGx(2,jj,ii,kk)*dGy(1,jj,ii,kk))
                          ! !dw/dy-dv/dz
                          ! omega_1(jj,ii,kk) = dGy(3,jj,ii,kk) - dGz(2,jj,ii,kk)
                          ! !du/dz-dw/dx
                          ! omega_2(jj,ii,kk) = dGz(1,jj,ii,kk) - dGx(3,jj,ii,kk)
                          ! !dv/dx-du/dy
                          ! omega_3(jj,ii,kk) = dGx(2,jj,ii,kk) - dGy(1,jj,ii,kk)
                        end do
                      end do
                    end do
                    !$omp end parallel do

                    !Gは計算破綻時間には常にNaNなのでdpは計算不能

                    !i5.5で5桁分の数字を表示できるのでdt=1.d-5以下で計算するならここも変更が必要
                    write(filename, '(i6.6)') M-1
                    !Mの計算毎に出力ファイル名を変更して出力する
                    !計算破綻直前の値を出力するので1step前の結果になる
                    do kk= 0,Nz-1
                      z=dz*dble(kk)
                      write(z_name, '(i2.2)') kk
                      open(10, file = "result_omp/parameter"//trim(filename)//"_"//trim(z_name)//".txt")
                      do ii = 0,Ny
                        do jj = 0,Nx
                          write(10,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16,",",&
                          &f24.16)') zeta_fx(jj),zeta_fy(ii),z,&
                          oldG(0,jj,ii,kk),div_u(jj,ii,kk),Invariant_2(jj,ii,kk)
                        enddo
                        write(10,*)
                        !一度に全てを出力する際にはデータの切れ目として空白を一行挿入しなくてはいけない
                      enddo
                      close(10)
                    enddo
                    write(*,*) "x=",i,"y=",j,"z=",k,"M=",M
                    stop "rho becomes NAN"
                  endif
                enddo
              enddo
          enddo
      !=======最初にNaNでなっている場所を特定したいので、並列化不可能======================

!              Q = Q1!オイラー法の時間の更新
        !RK法の時間の更新
              Q = Qn
        write(*,*) "M=",M!計算に時間がかかるので進行状況の確認用に出力
      enddo
! ===========メイン計算終了========================================================
     close(41)
     close(42)
     close(43)
     close(44)
      deallocate(G,Q,Q0,Q1,Q2,Qn,Fpx,Fmx,xp,xm,oldG)
      deallocate(Fpy,Fmy,yp,ym,Fpz,Fmz,zp,zm,myu)
      deallocate(LUmx,LUpx,LUmy,LUpy,LUmz,LUpz,LUccsx,LUccsy,LUccsz)
      deallocate(Vx,dVx,UVWT,dUVWTx,Vy,dVy,dUVWTy,Vz,dVz,dUVWTz)
      deallocate(in_G0,in_G1_top,in_G2,in_G3,dGx,dFx,dGy,dFy,dGz,dFz)
      ! deallocate(in_G1_du)
      deallocate(Ux,sigma_x,Uy,sigma_y,dQx,dQy,dzeta_iny,dzeta_inx)
      deallocate(omega_1,omega_2,omega_3,dp,div_u,Invariant_2)
    end program main
