!dx,dyの格子伸長の倍率は0.5倍/1.2倍に設定してある
!2019.05.18 Nx=180 Ny=100 dt=5.d-3で実行すれば卒論の格子伸長適用時の条件になる
!2019.05.20 Nx,Ny,dtをプログラム内で定義する仕様に変更開始>終了!!
!>> 496行目付近でUxを確実に0にしている部分の座標をNxを変えるなら確認しないといけない
!また定置観測している箇所の座標もNx,Nyを変えるなら調節しないといけないwrite文二つ使用
!2019.05.27 粘性項の係数μの設定でTがMa数を反映していない状況になっていた。修正済
!2019.05.27 NSCBC_xのi=0の流入部でのL3の設定がv方向に撹乱を入れない前提で0にしていたが今回は入れているので修正した
!2019.06.20　計算結果の出力に関して渦度を出していたが誤ったG(Q2から導出したG)でdGx,dGyを計算していたのでQnから計算したdGx,dGyで渦度を求めるように変更した
!2020.06.04 計算結果出力形式を.dから.txtへ変更した
!2020.06.13 Pr=1の理由を探す事に。今回は森山が現実的という0.71で計算する
module dimension
  !連続の式、Eulerの運動方程式、エネルギー方程式を並列に並べた行列Q,Fの設定等をする
  !これらの式をまとめて基礎式と呼ぶ
  implicit none
  !計算条件値の設定
  double precision,parameter :: gamma = 1.4d0
  integer,parameter :: t_end = 100 !時刻tの設定
  integer,parameter :: p_output = 10 !時間毎の局所圧力を出力させる際のステップ間隔
  integer,parameter :: Nx = 180
  integer,parameter :: Ny = 100
  double precision,parameter :: dt = 2.d-3
  double precision,parameter :: b = 1.d0!Jet半径は1で固定してしまう
  double precision,parameter :: Cx = 24.d0*b !x軸の幅の設定
  double precision,parameter :: Cy = 8.d0*b !y軸の幅の設定
  !Buffer領域の幅は常にWx<=Cxと計算領域よりも小さくてはならない
  double precision,parameter :: Wrx = 12.d0*b!Buffer領域x方向右側の幅
  double precision,parameter :: Wlx = Wrx!Buffer領域x方向左側の幅
  double precision,parameter :: Wry = 2.d0*b!Buffer領域y方向右側の幅
  double precision,parameter :: Wly = Wry!Buffer領域y方向左側の幅
  double precision,parameter :: Lx =  Cx+Wrx!+Wlx x方向の長さを定義.x軸左側にもbufferをかけるなら変更が必要
  double precision,parameter :: Ly = Cy+Wry!y方向の長さを定義 計算領域がy軸対称なのでWyも片方だけ

  double precision,parameter :: psigma = -0.25d0
  double precision,parameter :: msigma = 0.25d0
  double precision,parameter :: NS_sigma = 0.25d0
  double precision,parameter :: ccs_sigma = 0.d0
  double precision,parameter :: c = 1.d0
  !亜音速流入のためRe数は小さめに
  double precision,parameter :: Pr = 0.71d0
  double precision,parameter :: Ma = 5.d-1
  !Ma数も同様に小さめに
  double precision,parameter :: Temp = 1.d0
  double precision,parameter :: Tjet = 1.12d0*Temp
  double precision,parameter :: ujet = 1.d0
  double precision,parameter :: Sc = 120.d0 / (273.15d0 + 18.d0)
  double precision,parameter :: zeta = 1.d0
  double precision,parameter :: pi = acos(-1.d0)
  double precision,parameter :: Re = 1.d3
  double precision,parameter :: A2 = ujet*5.d-2!v方向の撹乱の振幅
  double precision,parameter :: T2 = 1.d0!v方向の撹乱の周期 T2=1で1秒に1周期
contains
  !初期条件G(rho,u,p)を用いてQ行列の設定
  subroutine Q_matrix(G,Q)
    double precision,dimension(0:3,0:Nx,0:Ny) :: Q,G
      Q(0,:,:) = G(0,:,:)
      Q(1,:,:) = G(0,:,:) * G(1,:,:)
      Q(2,:,:) = G(0,:,:) * G(2,:,:)
      Q(3,:,:) = G(3,:,:) / (gamma -1.d0) + G(0,:,:) * (G(1,:,:)**2.d0 + &
      &G(2,:,:)**2.d0) * 0.5d0
  endsubroutine Q_matrix
  !求めQを用いてF行列の設定(x方向)
  subroutine Fx_matrix(Q,Fpx,Fmx)
    double precision,dimension(0:3,0:Nx,0:Ny) :: Q,Fx,Fpx,Fmx
    Fx=0.d0
        Fx(0,:,:) = Q(1,:,:)
       Fx(1,:,:) = 1.d0/(2.d0*Q(0,:,:))*((3.d0-gamma)*(Q(1,:,:)**2.d0)+&
                  (1.d0-gamma)*(Q(2,:,:)**2.d0))+(gamma-1.d0)*Q(3,:,:)
        Fx(2,:,:) = Q(1,:,:)*Q(2,:,:)/Q(0,:,:)
       Fx(3,:,:) = gamma * Q(1,:,:) * Q(3,:,:) / Q(0,:,:) &
                 &- (gamma - 1.d0) * (Q(1,:,:)**3.d0 + Q(1,:,:)*(Q(2,:,:)**2.d0))&
                 & / (2.d0 * (Q(0,:,:)**2.d0))
    !求めたFを特製速度の正負によって分割する
    !Lax-Friedrichの流速分割を用いる
        Fpx(:,:,:) = (Fx(:,:,:) + zeta * Q(:,:,:)) * 0.5d0
        Fmx(:,:,:) = (Fx(:,:,:) - zeta * Q(:,:,:)) * 0.5d0
  endsubroutine Fx_matrix
    !求めQを用いてF行列の設定(y方向)
    subroutine Fy_matrix(Q,Fpy,Fmy)
      double precision,dimension(0:3,0:Nx,0:Ny) :: Q,Fy,Fpy,Fmy
      Fy=0.d0
          Fy(0,:,:) = Q(2,:,:)
          Fy(1,:,:) = Q(1,:,:)*Q(2,:,:)/Q(0,:,:)
          Fy(2,:,:) = 1.d0/(2.d0*Q(0,:,:))*((3.d0 - gamma)*(Q(2,:,:)**2.d0) + &
          &(1.d0 - gamma)*(Q(1,:,:)**2.d0))+(gamma - 1.d0)*Q(3,:,:)
          Fy(3,:,:) = gamma * Q(2,:,:) * Q(3,:,:) / Q(0,:,:) &
                    &- (gamma - 1.d0) * (Q(2,:,:)**3.d0 + Q(2,:,:)*(Q(1,:,:)**2.d0))&
                     &/ (2.d0 * (Q(0,:,:)**2.d0))
      !求めたFを特製速度の正負によって分割する
      !Lax-Friedrichの流速分割を用いる
          Fpy(:,:,:) = (Fy(:,:,:) + zeta * Q(:,:,:)) * 0.5d0
          Fmy(:,:,:) = (Fy(:,:,:) - zeta * Q(:,:,:)) * 0.5d0
      endsubroutine Fy_matrix
!du/dx,dT.dxを求めるためにまずはuとTの設定。一行目は他と係数を合わせ易くするために0だが含めている
!UVT,myuの計算はrho,u,pをそのまま代入しても良いがその場合は求めたtQから毎回rho_u_pのsubroutineを
!呼び出して計算しなければいけないので今回はQから直接計算できるようにプログラムを組んだ
    subroutine variable_setting(UVT,Q,myu)
      double precision,dimension(0:3,0:Nx,0:Ny) :: UVT, Q
      double precision,dimension(0:Nx,0:Ny) :: myu
        UVT(0,:,:) = 0.d0
        UVT(1,:,:) = Q(1,:,:) / Q(0,:,:)
        UVT(2,:,:) = Q(2,:,:) / Q(0,:,:)
        UVT(3,:,:) = (Ma**2.d0)*(gamma * (gamma - 1.d0) * (Q(3,:,:) - (((Q(1,:,:) **2.d0)+ &
        &(Q(2,:,:) **2.d0)) / (2.d0 * Q(0,:,:))))) / Q(0,:,:)
        !Tの値はMa^2*gamma*p/rhoこれは速度で無次元化したもの
!V行列を設定する際にμの計算が複雑になっているのでそれを簡略に示すために別でμを計算するsubroutine
!UVTからTの値を代入することで計算を簡略化している
        myu(:,:) = (UVT(3,:,:) ** 1.5d0) * (1.d0 + Sc) / (UVT(3,:,:) + Sc)
    end subroutine variable_setting

    subroutine Vx_matrix(Vx,myu,UVT,dUVTx,dUVTy)
      !粘性項の設定(x方向)
      !基礎式右辺をV行列として設定
      double precision,dimension(0:3,0:Nx,0:Ny) :: Vx,dUVTx,dUVTy,UVT
      double precision,dimension(0:Nx,0:Ny) :: myu
        Vx(0,:,:) = 0.d0
        Vx(1,:,:) = (2.d0*myu(:,:) / (3.d0*Re)) * (2.d0 * dUVTx(1,:,:) - dUVTy(2,:,:))
        Vx(2,:,:) = (myu(:,:) / Re) * (dUVTx(2,:,:) + dUVTy(1,:,:))
        Vx(3,:,:) = Vx(1,:,:)*UVT(1,:,:) + Vx(2,:,:)* UVT(2,:,:)+ &
        ((myu(:,:) * dUVTx(3,:,:)) / ((gamma - 1.d0)*Re*Pr*(Ma ** 2.d0)))
    endsubroutine Vx_matrix

    subroutine Vy_matrix(Vy,myu,UVT,dUVTx,dUVTy)
      !粘性項の設定(y方向)
      !基礎式右辺をV行列として設定
      double precision,dimension(0:3,0:Nx,0:Ny) :: Vy,dUVTx,dUVTy,UVT
      double precision,dimension(0:Nx,0:Ny) :: myu
        Vy(0,:,:) = 0.d0
        Vy(1,:,:) = (myu(:,:) / Re) * (dUVTx(2,:,:) + dUVTy(1,:,:))
        Vy(2,:,:) = (myu(:,:) / Re) * (2.d0/3.d0)*(-dUVTx(1,:,:) + 2.d0*dUVTy(2,:,:))
        Vy(3,:,:) = Vy(1,:,:)*UVT(1,:,:)+ Vy(2,:,:) * UVT(2,:,:)+&
        ((myu(:,:) * dUVTy(3,:,:)) / ((gamma - 1.d0)*Re*Pr*(Ma ** 2.d0)))
    endsubroutine Vy_matrix
    !求めたQnから時間毎のrho,u,pを求めるサブルーチン
    !ある指定した時間の時の計算結果のみを取り出せばグラフが作成できるので毎回は使用しない
    subroutine rho_u_p(cG,Qn)
      double precision,dimension(0:3,0:Nx,0:Ny) :: Qn,cG
        cG(:,:,:) = 0.d0
        cG(0,:,:) = Qn(0,:,:)
        cG(1,:,:) = Qn(1,:,:) / Qn(0,:,:)
        cG(2,:,:) = Qn(2,:,:) / Qn(0,:,:)
        cG(3,:,:) = (gamma - 1.d0) * (Qn(3,:,:) - (Qn(1,:,:) ** 2.d0 + &
        &Qn(2,:,:) ** 2.d0) / (2.d0 * Qn(0,:,:)))
    endsubroutine rho_u_p

    !DCS用の行列Aの設定(左辺の設定)
    !５次精度DCSと6次精度CCSの両方に対応
    subroutine LU_Decompo(N,sigma,LU)
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
      !0行目
      A(0,0) = 1.d0
      A(0,1) = 3.d0
      !1行目
      A(1,0) = alpha3 * (1.d0 - Dalpha * sigma)
      A(1,1) = 1.d0
      A(1,2) = alpha3 * (1.d0 + Dalpha * sigma)
      !2からN-2行目まで
      do i = 2,N-2
        A(i,i-1) = (1.d0 - Dalpha * sigma) * alpha5
        A(i,i) = 1.d0
        A(i,i+1) = (1.d0 + Dalpha * sigma) * alpha5
      enddo
      !N-1行目
      A(N-1,N-2) = alpha3 * (1.d0 - Dalpha * sigma)
      A(N-1,N-1) = 1.d0
      A(N-1,N) = alpha3 * (1.d0 + Dalpha * sigma)
      !N行目
      A(N,N-1) = 3.d0
      A(N,N) = 1.d0
    !LU分解
      !まずL,Uの初期値を設定
      do i = 0,N
        U(0,i) = A(0,i)
        L(i,0) = A(i,0) / U(0,0)
        L(i,i) = 1.d0
      enddo
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
!LU行列の圧縮　L,U行列の対角成分のみをLU行列に保存する
!LU(-1,i)にL行列の対角成分L(i,i-1)を保存
!L(i,i)=1.d0は全身代入では使用しないので圧縮した行列には含めない
!LU(0,i)にU_i,iの成分を保存
!LU(1,i)にU_i,i+1の成分を保存
!これは6次制度CCSと5次制度DCSにのみ対応しているのでそれ以上の精度で計算する場合には変更が必要
      LU(-1,0) = 0.d0
      LU(0,N) = U(N,N)
      LU(1,N) = 0.d0
      do i = 1,N
        LU(-1,i) = L(i,i-1)
      enddo
      do j = 0,1
        do i = 0,N-1
          LU(j,i) = U(i,i+j)
        enddo
      enddo
      deallocate(A,L,U)
    endsubroutine LU_Decompo

  !DCS右辺の計算(RHS)サブルーチン
  !x方向
    subroutine dif_x(sigma,dx,Fx,dFzeta,LU,dzeta_inx)
      integer i
      double  precision,dimension(0:3,0:Nx,0:Ny) :: RHS_x,Fx
      double  precision,dimension(0:3,2:Nx-2,0:Ny) :: D2,D4,D6,D8
      double precision dx,sigma
      double  precision,parameter :: ra = 14.d0/9.d0, rb = 1.d0/9.d0&
      &,da = 4.d0 / 9.d0,db = 2.d0 / 9.d0 !5次精度のDCSとなるための係数設定
      double precision,allocatable,dimension(:,:) :: LU
      double precision,dimension(0:3,0:Nx,0:Ny) :: x,y
      double precision,dimension(0:3,0:Nx,0:Ny) :: dzeta_inx,dFzeta
      D2=0.d0;D4=0.d0;D6=0.d0;D8=0.d0;RHS_x=0.d0;y=0.d0;x=0.d0
      !片側DCS,3次精度DCSも入れた非周期条件の際のbの設定
        !片側DCSの右辺設定
        RHS_x(:,0,:)=((-17.d0/6.d0)*Fx(:,0,:)+(1.5d0)*Fx(:,1,:)+(1.5d0)&
                      *Fx(:,2,:)-Fx(:,3,:)/6.d0)/dx
        RHS_x(:,Nx,:)=((1.d0/6.d0)*Fx(:,Nx-3,:)-(1.5d0)* Fx(:,Nx-2,:)-&
                      (1.5d0)*Fx(:,Nx-1,:)+(17.d0/6.d0)*Fx(:,Nx,:))/dx
        !3次精度DCSの右辺設定
        RHS_x(:,1,:)=((1.5d0)*(-Fx(:,0,:)+Fx(:,2,:))/(2.d0*dx))+sigma*&
                      ((Fx(:,0,:)-2.d0*Fx(:,1,:)+Fx(:,2,:))/(2.d0*dx))
        RHS_x(:,Nx-1,:)=((1.5d0)*(-Fx(:,Nx-2,:)+Fx(:,Nx,:))/(2.d0*dx))+&
                      sigma*((Fx(:,Nx-2,:)-2.d0*Fx(:,Nx-1,:)+Fx(:,Nx,:))/(2.d0*dx))
         !5次精度DCSの右辺設定
       do i = 2,Nx-2
         D2(:,i,:) = (-Fx(:,i-1,:)+Fx(:,i+1,:)) / (2.d0*dx)
         D4(:,i,:) = (-Fx(:,i-2,:)+Fx(:,i+2,:)) / (4.d0*dx)
         !D6(:,i,:) = (-Fx(:,i-3,:)+Fx(:,i+3,:)) / (6.d0*dx)
         !7次精度DCS用のため不要
         D6(:,i,:) = (Fx(:,i-1,:)+Fx(:,i+1,:)- 2.d0* Fx(:,i,:)) / dx
         D8(:,i,:) = (Fx(:,i-2,:)+Fx(:,i+2,:)- 2.d0* Fx(:,i,:)) / (4.d0*dx)
         ! D12(:,i,:) = (Fx(:,i-3,:)+Fx(:,i+3,:)- 2.d0* Fx(:,i,:)) / (9.d0*dx)
         !7次精度DCS用のため不要
         RHS_x(:,i,:)=ra*D2(:,i,:)+rb*D4(:,i,:)+sigma*(da*D6(:,i,:)+db*D8(:,i,:))
       enddo
     !前進代入法、後退代入法の計算サブルーチン(x方向)
            !前進代入
            y(:,0,:) = RHS_x(:,0,:)!例外の境界値
            do i = 1,Nx!y(:,k)でk=0は上で定義したので残りの1〜Nxを定義する
              y(:,i,:) = RHS_x(:,i,:) - LU(-1,i)*y(:,i-1,:)
            enddo
            !後退代入
            x(:,Nx,:) = y(:,Nx,:) / LU(0,Nx)!例外の境界値
            !x(:,i)でi=Nxは定義したので残りのi=0~Nx-1を定義する
            do i = Nx-1, 0, -1!後退するので-1ずつ進む
              x(:,i,:) = (y(:,i,:) - LU(1,i)*x(:,i+1,:)) / LU(0,i)
            enddo
          call combine(dzeta_inx,x,dFzeta)
      end subroutine dif_x
   !DCS右辺の計算(RHS)サブルーチン
   !y方向
     subroutine dif_y(sigma,dy,Fy,dFzeta,LU,dzeta_iny)
       integer i
       double  precision,dimension(0:3,0:Nx,0:Ny) :: RHS_y,Fy
       double  precision,dimension(0:3,0:Nx,2:Ny-2) :: D2,D4,D6,D8
       double precision dy,sigma
       double  precision,parameter :: ra = 14.d0/9.d0, rb = 1.d0/9.d0&
       &,da = 4.d0 / 9.d0,db = 2.d0 / 9.d0 !5次精度のDCSとなるための係数設定
       double precision,allocatable,dimension(:,:) :: LU
       double precision,dimension(0:3,0:Nx,0:Ny) :: x,y
       double precision,dimension(0:3,0:Nx,0:Ny) :: dzeta_iny,dFzeta
       D2=0.d0;D4=0.d0;D6=0.d0;D8=0.d0;RHS_y=0.d0;y=0.d0;x=0.d0
       !片側DCS,3次精度DCSも入れた非周期条件の際のbの設定
         !片側DCSの右辺設定
         RHS_y(:,:,0) = ((-17.d0/6.d0)*Fy(:,:,0)+(1.5d0)*Fy(:,:,1)+&
         (1.5d0)*Fy(:,:,2)-Fy(:,:,3)/6.d0)/dy
         RHS_y(:,:,Ny)=((1.d0/6.d0)*Fy(:,:,Ny-3)-(1.5d0)&
         &*Fy(:,:,Ny-2)-(1.5d0)*Fy(:,:,Ny-1)+(17.d0/6.d0)*Fy(:,:,Ny))/dy
         !3次精度DCSの右辺設定
         RHS_y(:,:,1)=((1.5d0)*(-Fy(:,:,0)+Fy(:,:,2))/(2.d0*dy))+sigma*&
         ((Fy(:,:,0)-2.d0*Fy(:,:,1)+Fy(:,:,2))/(2.d0*dy))
         RHS_y(:,:,Ny-1)=((1.5d0)*(-Fy(:,:,Ny-2)+Fy(:,:,Ny))/&
         &(2.d0*dy))+sigma*((Fy(:,:,Ny-2)-2.d0*Fy(:,:,Ny-1)+Fy(:,:,Ny))/(2.d0*dy))
          !5次精度DCSの右辺設定
      do i = 2,Ny-2
        D2(:,:,i) = (-Fy(:,:,i-1)+Fy(:,:,i+1)) / (2.d0*dy)
        D4(:,:,i) = (-Fy(:,:,i-2)+Fy(:,:,i+2)) / (4.d0*dy)
        !D6(:,:,i) = (-Fy(:,:,i-3)+Fy(:,:,i+3)) / (6.d0*dy)
        !7次精度DCS用のため不要
        D6(:,:,i) = (Fy(:,:,i-1)+Fy(:,:,i+1)- 2.d0* Fy(:,:,i)) / dy
        D8(:,:,i) = (Fy(:,:,i-2)+Fy(:,:,i+2)- 2.d0* Fy(:,:,i)) / (4.d0*dy)
        ! D12(:,:,i) = (Fy(:,:,i-3)+Fy(:,:,i+3)- 2.d0* Fy(:,:,i)) / (9.d0*dy)
        !7次精度DCS用のため不要
        RHS_y(:,:,i)=ra*D2(:,:,i)+rb*D4(:,:,i)+sigma*(da*D6(:,:,i)+db*D8(:,:,i))
      enddo
      !前進代入法、後退代入法の計算サブルーチン(y方向)
        !前進代入
        y(:,:,0) = RHS_y(:,:,0)
        do i = 1,Ny!y(:,:,i)でi=0は上で定義したので残りの1〜Nyを定義する
          y(:,:,i) = RHS_y(:,:,i) - LU(-1,i)*y(:,:,i-1)
        enddo
        !後退代入
        x(:,:,Ny) = y(:,:,Ny) / LU(0,Ny)
        !x(:,:,i)でi=Nyは定義したので残りのi=0~Ny-1を定義する
        do i = Ny-1, 0, -1!後退するので-1ずつ進む
          x(:,:,i) = (y(:,:,i) - LU(1,i)*x(:,:,i+1)) / LU(0,i)
        enddo
        call combine(dzeta_iny,x,dFzeta)!y方向の格子伸長を適用　微分変換をしている
       end subroutine dif_y

    !境界条件をNSCBCで設定
    !L,d行列を設定することでdfx(0:2,0)とdFx(0:2,Nx)の値を定める
    !定めたdFxをxm,xpをそれぞれ求めた時点でi=0とNxで置き換える
    !これが境界条件となる。neumannはNx-1とNx,1と0を同じにしたがNSCBC_xはdFxの境界を定めることで
    !境界条件を設定している
    !そしてそのままQ1,Q2,Qnを求める
    !まずはx方向用のNSCBC　subrouitineを作成
    !(:,:,:)等で計算を行わせる際は:の配列数が対応していることが要確認
    subroutine NSCBC_x(G,dGx,dFx,pNx_infty)
      double precision,dimension(0:3,0:Nx,0:Ny) :: G,dGx,dFx
      double precision,dimension(1:5,0:Ny,0:1) :: L,d
      double precision,dimension(0:Ny,0:1) :: c_NS,Ma_NS
      double precision pNx_infty
      L=0.d0;d=0.d0;c_NS=0.d0;Ma_NS=0.d0
      !音速cはi=0,Nxの両点においてそれぞれ定義しなければならない
      c_NS(:,0) = sqrt(gamma * G(3,0,:) / G(0,0,:))
      c_NS(:,1) = sqrt(gamma * G(3,Nx,:) / G(0,Nx,:))
      !マッハ数Ma_NSはi=0,Nxで使うので別々に定義する
      Ma_NS(:,0) = G(1,0,:) / c_NS(:,0)!uを使う
      Ma_NS(:,1) = G(1,Nx,:) / c_NS(:,1)
      !x方向右側つまりi=0の点において亜音速流入条件でL行列を設定する
      L(1,:,0)=(G(1,0,:)-c_NS(:,0))*(-G(0,0,:)*c_NS(:,0)*dGx(1,0,:)+dGx(3,0,:))
!      L(2,:,0)=G(1,0,:)*((c_NS(:,0)**2.d0)*dGx(0,0,:)-dGx(3,0,:))
!      L(3,:,0)=0.d0!L3=u*dv/dxで流入速度
      L(3,:,0)=G(1,0,:)*dGx(2,0,:)
      L(4,:,0)=0.d0 !!今回は2次元なのでz方向成分は0
      L(5,:,0)=L(1,:,0) !-2.d0*c_NS(:,0)*du/dtが本来はあるが
      !流入速度uを時間変動させないので今回はdu/dt=0となるため省略
!      L(5,:,0)=msigma*c_NS(:,0)*(1.d0-(Ma_NS(:,0)**2.d0))*(G(3,0,:)-p0x_infty)/Lx
      L(2,:,0)=(5.d-1)*(gamma-1.d0)*(L(5,:,0)+L(1,:,0))!+G(0,0,:)*c_NS(:,0)**2.d0/T*dT/dtが本来はあるが
      !境界温度Tを時間変動させないので今回はdT/dt=0となるため省略
      !x方向左側つまりi=Nxの点において無反射流入条件でL行列を設定する
      L(1,:,1)=NS_sigma*c_NS(:,1)*(1.d0-(Ma_NS(:,1)**2.d0))*(G(3,Nx,:)-&
      &pNx_infty)/Lx
      L(2,:,1)=G(1,Nx,:)*((c_NS(:,1)**2.d0)*dGx(0,Nx,:)-dGx(3,Nx,:))
      L(3,:,1)=G(1,Nx,:)*dGx(2,Nx,:)
      L(4,:,1)=0.d0!!今回は2次元なのでz方向成分は0
      L(5,:,1)=(G(1,Nx,:)+c_NS(:,1))*(G(0,Nx,:)*c_NS(:,1)*dGx(1,Nx,:)+dGx(3,Nx,:))
      !設定したL行列からd1~5をi=0,Nxの両点においてそれぞれ設定する
        d(1,:,:) = (1.d0 / (c_NS(:,:) **2.d0)) * ((L(1,:,:)+L(5,:,:))/2.d0 + L(2,:,:))
        d(2,:,:) = (L(1,:,:)+L(5,:,:))*0.5d0
        !d3のみrhoを含むので個別で設定しなければいけない
        d(3,:,0) = (1.d0/ (2.d0 * G(0,0,:) * c_NS(:,0))) * (-L(1,:,0) + L(5,:,0))
        d(3,:,1) = (1.d0/ (2.d0 * G(0,Nx,:) * c_NS(:,1))) * (-L(1,:,1) + L(5,:,1))
        d(4,:,:) = L(3,:,:)
        d(5,:,:) = L(4,:,:)
      !設定したdからNxSCBCで置き換える境界地点のdFxを定義する
      !i=0の時の差し替えdFx
      dFx(0,0,:) = d(1,:,0)
      dFx(1,0,:) = G(1,0,:)*d(1,:,0)+G(0,0,:)*d(3,:,0)
      dFx(2,0,:) = G(2,0,:)*d(1,:,0)+G(0,0,:)*d(4,:,0)
      dFx(3,0,:) =(5.d-1)*((G(1,0,:)**2.d0)+(G(2,0,:)**2.d0))*d(1,:,0)+d(2,:,0)/(gamma-1.d0)&
              &  + (G(0,0,:)*G(1,0,:)*d(3,:,0)) + (G(0,0,:)*G(2,0,:)*d(4,:,0))

      !i=Nxの時の差し替えF
      dFx(0,Nx,:) = d(1,:,1)
      dFx(1,Nx,:) = (G(1,Nx,:)*d(1,:,1)) + (G(0,Nx,:)*d(3,:,1))
      dFx(2,Nx,:) = (G(2,Nx,:)*d(1,:,1)) + (G(0,Nx,:)*d(4,:,1))
      dFx(3,Nx,:) =(5.d-1)*((G(1,Nx,:)**2.d0)+(G(2,Nx,:)**2.d0))*d(1,:,1)+d(2,:,1)/(gamma-1.d0)&
              &  + (G(0,Nx,:)*G(1,Nx,:)*d(3,:,1)) + (G(0,Nx,:)*G(2,Nx,:)*d(4,:,1))
    !NSCBCの角処理(x方向,y方向で設定した境界値が重複するため半分ずつ加える)
                  dFx(:,0,0) = dFx(:,0,0) * 0.5d0
                  dFx(:,0,Ny) = dFx(:,0,Ny) * 0.5d0
                  dFx(:,Nx,0) = dFx(:,Nx,0) * 0.5d0
                  dFx(:,Nx,Ny) = dFx(:,Nx,Ny) * 0.5d0
    endsubroutine NSCBC_x
    !次にy方向のNSCBC　sunrouineを作成
    subroutine NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
      double precision,dimension(0:3,0:Nx,0:Ny) :: G,dGy,dFy
      double precision,dimension(1:5,0:Nx,0:1) :: L,d
      double precision,dimension(0:Nx,0:1) :: c_NS,Ma_NS
      double precision pNy_infty,p0y_infty
      L=0.d0;d=0.d0;c_NS=0.d0;Ma_NS=0.d0
      !音速cはi=0,Nxの両点においてそれぞれ定義しなければならない
      c_NS(:,0) = sqrt(gamma * G(3,:,0) / G(0,:,0))
      c_NS(:,1) = sqrt(gamma * G(3,:,Ny) / G(0,:,Ny))
      !マッハ数Maはi=0,Nxで使うので別々に定義する
      Ma_NS(:,0) = G(2,:,0) / c_NS(:,0)!vを使う
      Ma_NS(:,1) = G(2,:,Ny) / c_NS(:,1)
      !y方向右側つまりi=0の点において無反射流出条件でL行列を設定する
      L(1,:,0) = (G(2,:,0) - c_NS(:,0)) * (-G(0,:,0)*c_NS(:,0)*dGy(2,:,0)+dGy(3,:,0))
      L(2,:,0) = G(2,:,0) * ((c_NS(:,0) ** 2.d0)*dGy(0,:,0) - dGy(3,:,0))
      L(3,:,0) = G(2,:,0) * dGy(1,:,0)
      L(4,:,0) = 0.d0 !!今回は2次元なのでz方向成分は0
      L(5,:,0) = NS_sigma * c_NS(:,0) * (1.d0 - (Ma_NS(:,0) ** 2.d0))*(G(3,:,0) - &
      &p0y_infty) / (2.d0*Ly)
      !y方向左側つまりi=Nyの点において無反射流出条件でL行列を設定する
      L(1,:,1) = NS_sigma * c_NS(:,1) * (1.d0 - (Ma_NS(:,1) ** 2.d0))*(G(3,:,Ny) - &
    &  pNy_infty) / (2.d0*Ly)
      L(2,:,1) = G(2,:,Ny) * ((c_NS(:,1) ** 2.d0)*dGy(0,:,Ny) - dGy(3,:,Ny))
      L(3,:,1) = G(2,:,Ny) * dGy(1,:,Ny)
      L(4,:,1) = 0.d0 !!今回は2次元なのでz方向成分は0
      L(5,:,1) = (G(2,:,Ny) + c_NS(:,1)) * (G(0,:,Ny)*c_NS(:,1)*dGy(2,:,Ny) + dGy(3,:,Ny))
      !設定したL行列からd1~5をi=0,Nyの両点においてそれぞれ設定する
        d(1,:,:) = (1.d0 / (c_NS(:,:) **2.d0)) * ((L(1,:,:)+L(5,:,:))/2.d0 + L(2,:,:))
        d(2,:,:) = (5.d-1) * (L(1,:,:)+L(5,:,:))
        d(3,:,:) = L(3,:,:)
        !d4のみrhoを含むので個別で設定しなければいけない
        d(4,:,0) = (1.d0/ (2.d0 * G(0,:,0) * c_NS(:,0))) * (-L(1,:,0) + L(5,:,0))
        d(4,:,1) = (1.d0/ (2.d0 * G(0,:,Ny) * c_NS(:,1))) * (-L(1,:,1) + L(5,:,1))
        d(5,:,:) = L(4,:,:)
      !設定したdからNySCBCで置き換える境界地点のdFyを定義する
      !i=0の時の差し替えdFy
      dFy(0,:,0) = d(1,:,0)
      dFy(1,:,0) = G(1,:,0)*d(1,:,0)+G(0,:,0)*d(3,:,0)
      dFy(2,:,0) = G(2,:,0)*d(1,:,0)+G(0,:,0)*d(4,:,0)
      dFy(3,:,0) =(5.d-1)*((G(1,:,0)**2.d0)+(G(2,:,0)**2.d0))*d(1,:,0)+&
      d(2,:,0)/(gamma-1.d0)+(G(0,:,0)*G(1,:,0)*d(3,:,0))+(G(0,:,0)*G(2,:,0)*d(4,:,0))

      !i=Nyの時の差し替えdFy
      dFy(0,:,Ny) = d(1,:,1)
      dFy(1,:,Ny) = (G(1,:,Ny)*d(1,:,1)) + (G(0,:,Ny)*d(3,:,1))
      dFy(2,:,Ny) = (G(2,:,Ny)*d(1,:,1)) + (G(0,:,Ny)*d(4,:,1))
      dFy(3,:,Ny) =(5.d-1)*((G(1,:,Ny)**2.d0)+(G(2,:,Ny)**2.d0))*d(1,:,1)+&
      d(2,:,1)/(gamma-1.d0)+(G(0,:,Ny)*G(1,:,Ny)*d(3,:,1))+(G(0,:,Ny)*G(2,:,Ny)*d(4,:,1))
      !NSCBCの角処理(x方向,y方向で設定した境界値が重複するため半分ずつ加える)
                    dFy(:,0,0) = dFy(:,0,0) * 0.5d0
                    dFy(:,0,Ny) = dFy(:,0,Ny) * 0.5d0
                    dFy(:,Nx,0) = dFy(:,Nx,0) * 0.5d0
                    dFy(:,Nx,Ny) = dFy(:,Nx,Ny) * 0.5d0
    endsubroutine NSCBC_y
    subroutine inflow(Q,in_G)
      double precision,dimension(0:3,0:Nx,0:Ny) :: Q
      double precision,dimension(0:2,0:Ny) :: in_G
      Q(0,0,:) = in_G(0,:)!今までと違いrhoをNSCBCで求めずにdirichlet条件で固定してしまう
      Q(1,0,:) = in_G(0,:)*in_G(1,:)!rho*u
      Q(2,0,:) = in_G(0,:)*in_G(2,:)!rho*v
      Q(3,0,:) = 1.d0/((Ma**2.d0)*gamma*(gamma-1.d0))&
                  +in_G(0,:)*(in_G(1,:)**2.d0 + in_G(2,:)**2.d0)/2.d0!Et
    endsubroutine inflow
    subroutine outflow(UVT,dUVTx,Vx,dVx,dUVTy,Vy,dVy)
      double precision,dimension(0:3,0:Nx,0:Ny) :: Vx,Vy,dUVTx,dUVTy,UVT,dVx,dVy
      !無反射流出条件の時の条件を設定するsubroutine
      !境界のdVx,dVyに条件を設定するのでdVの計算ができた後に境界値のみ上書きをする
      !x方向右側の条件設定
      dVx(2,Nx,:) = 0.d0
      dVx(3,Nx,:) = dVx(1,Nx,:)*UVT(1,Nx,:)+dUVTx(1,Nx,:)*Vx(1,Nx,:)+dUVTx(2,Nx,:)*Vx(2,Nx,:)
      !y方向左側の条件設定
      dVy(1,:,0) = 0.d0
      dVy(3,:,0) = dVy(2,:,0)*UVT(2,:,0)+dUVTy(1,:,0)*Vy(1,:,0)+dUVTy(2,:,0)*Vy(2,:,0)
      !y方向右側の条件設定
      dVy(1,:,Ny) = 0.d0
      dVy(3,:,Ny) = dVy(2,:,Ny)*UVT(2,:,Ny)+dUVTy(1,:,Ny)*Vy(1,:,Ny)+dUVTy(2,:,Ny)*Vy(2,:,Ny)
    endsubroutine outflow
    !buffer領域の設定subroutine
    !x方向
    !ここでまず計算に必要なU(x)とσ(x)を定義している
    !係数等はタサイさんの資料を参考にしている
    subroutine buffer_x(c_infty,Ux,sigma_x,zeta_fx)
      integer i
      double precision c_infty,Xmax,Xmin,x1
      double precision,dimension(0:3,0:Nx,0:Ny) :: Ux,sigma_x
      double precision,dimension(0:Nx) :: zeta_fx
      double  precision,parameter ::alpha_u=1.5d0,alpha_sigma=1.125d0,beta_r=0.01d0,beta_l=0.01d0
      Xmax = Lx;Xmin = 0.d0
      do i = 0,Nx
        x1 = zeta_fx(i)
        !Uxは本来0:Nxの一次の配列で十分だが座標変換の際に他の配列と計算する際に配列が揃っていないと:で省略して計算できないので
        !無駄に同じ値を入れて4次の配列にしている
        Ux(:,i,:) = alpha_u*c_infty*(dtanh(dble(atanh(beta_r/alpha_u-1.d0))*(x1-Xmax)/(-Wrx))&
        -dtanh(dble(atanh(beta_l/alpha_u-1.d0))*(x1-Xmin)/Wlx))
        if(x1<(Xmax-Wrx)) then
          sigma_x(:,i,:) = 0.d0!流出部にのみBufferをつけるのでx左側もσは0となる
        else
          sigma_x(:,i,:) = alpha_sigma*c_infty*((x1-(Xmax-Wrx))/Wrx)**3.d0
        endif
      enddo
       !dx不要かも？
     Ux(:,0:90,:) = 0.d0!x左側のBufferを取るためにWlxの範囲のUxを確実に0に設定している
      !open(100,file="ux-check.csv")
      !do i=0,Nx
      !write(100,*) zeta_fx(i),",",Ux(0,i,0)
      !enddo
      !close(100)
      !格子伸長が入っているためほぼ計算領域の真ん中になる値を調べて代入した
    endsubroutine buffer_x
    !y方向
    subroutine buffer_y(c_infty,Uy,sigma_y,zeta_fy)
      integer i
      double precision Ymax,Ymin,c_infty,y1
      double precision,dimension(0:3,0:Nx,0:Ny) :: Uy,sigma_y
      double precision,dimension(0:Ny) :: zeta_fy
      double  precision,parameter ::alpha_u=1.15d0,alpha_sigma=1.125d0,beta_r=0.01d0,beta_l=0.01d0
      Ymax = Ly;Ymin = -Ly
      !格子伸長を行うので新しい座標ζ_yを用いてUyとsigma_yを設定する
      do i = 0,Ny
        y1 = zeta_fy(i)
        Uy(:,:,i) = alpha_u*c_infty*(dtanh(dble(atanh(beta_r/alpha_u-1.d0))*(y1-Ymax)/(-Wry))&
        -dtanh(dble(atanh(beta_l/alpha_u-1.d0))*(y1-Ymin)/Wly))
        if(y1<(Wly+Ymin)) then
          sigma_y(:,:,i) = alpha_sigma*c_infty*((-y1+Ymin+Wly)/Wly)**3.d0
        elseif((y1>=(Wly+Ymin)).and.(y1<(Ymax-Wry))) then
          sigma_y(:,:,i) = 0.d0
        elseif(y1>=(Ymax-Wry)) then
          sigma_y(:,:,i) = alpha_sigma*c_infty*((y1-(Ymax-Wry))/Wry)**3.d0
        endif
      enddo
    endsubroutine buffer_y
    !ζ,dζ/dxの定義subroutine
    subroutine lattice_x(dx,zeta_fx,dzeta_inx)
      integer i
      double precision,dimension(0:Nx) :: zeta_fx
      double precision,dimension(0:3,0:Nx,0:Ny) :: dzeta,dzeta_inx
      double precision dx,x,width,a1,a2
      dzeta=0.d0;width=10.8d0;a1=1d0/14d0;a2=7.d0
      !widthは格子間隔を細かくする範囲。この式では-width<=x<=widthの範囲で適用される
      !a2は粗い所と細かい所の境界の傾きの大きさを設定している
      !a1はどの程度の格子数の差をつけるかを設定する係数。このパラメータは条件によって適宜調整する
      do i= 0,Nx
        x = dx*dble(i)
       zeta_fx(i) = (1.d0/1.4d0) * ((1.7d0*x) - a1 * &
        (-dlog(dcosh(a2*(x - width))) + dlog(dcosh(a2*(x + width)))))
        dzeta(:,i,:) = (1.d0/1.4d0) * (1.7d0 - (a1*a2) * &
        (-dtanh(a2*(x - width)) + dtanh(a2*(x + width))))
      enddo
      dzeta_inx(:,:,:) = 1.d0/dzeta(:,:,:)
    endsubroutine lattice_x
    !ζ,dζ/dyの定義subroutine
    subroutine lattice_y(dy,zeta_fy,dzeta_iny)
      integer i
      double precision,dimension(0:3,0:Nx,0:Ny) :: dzeta,dzeta_iny
      double precision,dimension(0:Ny) :: zeta_fy
      double precision dy,y,width,a1,a2
      dzeta=0.d0;width=3.d0;a1=1d0/14d0;a2=7d0
      !widthは格子間隔を細かくする範囲。この式では-width<=y<=widthの範囲で適用される
      !a2は粗い所と細かい所の境界の傾きの大きさを設定している
      !a1はどの程度の格子数の差をつけるかを設定する係数
      do i= 0,Ny
        y = -Ly + dy*dble(i)
        zeta_fy(i) = (1.d0/1.4d0) * ((1.7d0*y) - a1 * &
        (-dlog(dcosh(a2*(y - width))) + dlog(dcosh(a2*(y + width)))))
        dzeta(:,:,i) = (1.d0/1.4d0) * (1.7d0 - (a1*a2) * &
        (-dtanh(a2*(y - width)) + dtanh(a2*(y + width))))
      enddo
      dzeta_iny(:,:,:) = 1.d0/dzeta(:,:,:)
    endsubroutine lattice_y
    !作成したdx/dζをdF/dyなどに掛けて微分変換を行うsubroutine
    subroutine combine(dzeta_in,dF,dFzeta)
      double precision,dimension(0:3,0:Nx,0:Ny) :: dzeta_in,dFzeta,dF
      dFzeta(:,:,:) = dF(:,:,:) * dzeta_in(:,:,:)
    endsubroutine combine
end module dimension

    program main
      use dimension
      implicit none
      character(len = 16) filename
      !時間更新毎に出力ファイルを変更するためのファイル名設定
      !NSCBCでdrho,du,dpを計算する必要があるので計算しやすいようにG行列を作成しG(rho,u,p)
      !となるように設定した。すなわちdGは密度、速度、圧力のそれぞれの微分項を含む行列となる
      !x方向
      !Runge-Kutta法のためQ1,Q2用の値を設定
      double precision,dimension(0:3,0:Nx,0:Ny) :: G,Q,Q0,Q1,Q2,Qn,Fpx,Fmx,xp,xm,oldG
      !y方向
      double precision,dimension(0:3,0:Nx,0:Ny) :: Fpy,Fmy,yp,ym
      double precision,dimension(0:Nx,0:Ny) :: myu
      !A,L,U,x,y,bはF+とF-用に分けるため別々で定義
      !x方向
      double precision,allocatable,dimension(:,:) :: LUmx,LUpx
      !y方向
      double precision,allocatable,dimension(:,:) :: LUmy,LUpy
      !粘性項の計算に使う行列(x方向)
      double precision,dimension(0:3,0:Nx,0:Ny) :: Vx,dVx,UVT,dUVTx
      double precision,allocatable,dimension(:,:) :: LUccsx
      !y方向
      double precision,dimension(0:3,0:Nx,0:Ny) :: Vy,dVy,dUVTy
      double precision,allocatable,dimension(:,:) :: LUccsy
      ! NSCBC用
      double precision,dimension(0:2,0:Ny) :: in_G
      !x方向
      double  precision,dimension(0:3,0:Nx,0:Ny) :: dGx,dFx
      double precision  dx,x,pNx_infty
      !y方向
      double precision,dimension(0:3,0:Nx,0:Ny) :: dGy,dFy
      double precision dy,y,p0y_infty,pNy_infty
      integer i,j,M,Mmax,output_count
      double precision t0,t1,theta
      double precision c_infty
      double precision,dimension(0:Ny) :: ur,Tu
      double precision,dimension(0:3,0:Nx,0:Ny) :: Ux,sigma_x,Uy,sigma_y,dQx,dQy
      double precision,dimension(0:3,0:Nx,0:Ny) :: dzeta_iny,dzeta_inx
      double precision,dimension(0:Nx) :: zeta_fx
      double precision,dimension(0:Ny) :: zeta_fy
      double precision,dimension(0:Nx,0:Ny) :: omega_z,dp!渦度と圧力変動差を入れる配列
      !計算にかかる時間をCPU時間で計測する
      call cpu_time(t0)

      !x_axis
      allocate(LUmx(-1:1,0:Nx),LUpx(-1:1,0:Nx))
      allocate(LUccsx(-1:1,0:Nx))
      !y_axis
      allocate(LUmy(-1:1,0:Ny),LUpy(-1:1,0:Ny))
      allocate(LUccsy(-1:1,0:Ny))
      dx = Lx /dble(Nx)
      dy = 2.d0*Ly /dble(Ny)
      Mmax = t_end / dt
      output_count = int(0.1d0/dt)!出力ファイルを0.1sec間隔で出力するように設定
!一応ゼロクリア
      G=0.d0;Q=0.d0;Qn=0.d0;Fpx=0.d0;Fmx=0.d0;xp=0.d0;xm=0.d0;Q0=0.d0
      Q1=0.d0;Q2=0.d0;Fpy=0.d0;Fmy=0.d0;yp=0.d0;ym=0.d0
      myu=0.d0;pNx_infty=0.d0;p0y_infty=0.d0;pNy_infty=0.d0
      Vx=0.d0;dVx=0.d0;UVT=0.d0;dUVTx=0.d0;dFx=0.d0;dGx=0.d0
      Vy=0.d0;dVy=0.d0;dUVTy=0.d0;dFy=0.d0;dGy=0.d0;in_G=0.d0;ur=0.d0;Tu=0.d0
      Ux=0.d0;sigma_x=0.d0;Uy=0.d0;sigma_y=0.d0;dQx=0.d0;dQy=0.d0;zeta_fy=0.d0;dzeta_iny=0.d0
      zeta_fx=0.d0;dzeta_inx=0.d0;omega_z=0.d0;dp=0.d0;oldG=0.d0
      !y方向の格子伸長のための座標設定
      call lattice_y(dy,zeta_fy,dzeta_iny)
      !x方向も
      call lattice_x(dx,zeta_fx,dzeta_inx)
    !流入条件の導出
    !初期条件もζ_yの座標系で設定する
    ur(Ny/2) = ujet
    Tu(Ny/2) = Tjet
    do i = (Ny/2)+1,Ny
      !Top-hat型の分布になるような式を設定
      ur(i) = ujet/2.d0*(1.d0 - dtanh((12.5d0/4.d0)*((zeta_fy(i)/b)- (b/zeta_fy(i)))))
      !Crocco-Busemanの関係式と主流速度分布を用いて温度分布Tuを設定
      Tu(i) = Ma**2.d0*(gamma-1.d0)/2.d0*(ur(i)*ujet-ur(i)**2.d0)/ujet+&
              Tjet*ur(i)/ujet+Temp*(ujet-ur(i))/ujet
    enddo
    !初期分布をx軸対象になるようにする。
    !そのためにy軸正の範囲の値を負の範囲に軸対象になるようにコピーする
    do i = 0,(Ny/2)-1
      ur(i) = ur(Ny-i)
      Tu(i) = Tu(Ny-i)
    enddo
    !初期値
    !全体にまず初期値を入れてしまう
      G(0,:,:) = 1.d0!ρ
      G(1,:,:) = 0.d0!u
      G(2,:,:) = 0.d0!v
      G(3,:,:) = G(0,:,:)*Temp/((Ma**2.d0)*gamma)!p
      !Bufferの計算のための初期値を用いて無限遠方での音速を定義
      c_infty = sqrt(Temp/Ma**2.d0)
      !Buffer領域の計算に使うUx,Uy,sigma_x,sigma_yの計算
      call buffer_x(c_infty,Ux,sigma_x,zeta_fx)
      call buffer_y(c_infty,Uy,sigma_y,zeta_fy)
    !流入条件
    !x=0の軸上にのみ流入条件を適用することでここからどんどん流入が起こる
      in_G(0,:) = 1.d0/Tu(:)!密度ρは理想気体状態方程式に従うから
      in_G(1,:) = ur(:)! NSCBC流入条件で使う流入値のみを保存する配列
      !in_G(2,:) = 0.d0! NSCBC流入条件で使う流入値のみを保存する配列
      !Doループ内でin_G(2,:)をDirichlet条件で設定して撹乱を導入しているのでここでは設定しない
      !初期値の出力
      !まずt=0はループ外で個別に作成
      open(10, file = "result_2D/parameter000000.txt")
      !もちろん出力もζ_y座標系とζ_x座標系で行う
      do i = 0,Ny
        do j = 0,Nx
          write(10,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16)') &
          zeta_fx(j),zeta_fy(i),G(0,j,i),omega_z(j,i),dp(j,i)/dt
        enddo
        write(10,*)
      enddo
      !コメント行はファイルの最後に置くことで読み込みの際にエラーが起きる可能性をなくすようにする
      write(10,'(2A1,1I1)') "#","M",0!#を入れているのはgnuplotでコメントアウトするため
      write(10,'(6A10)')"#","x","y","rho","vorticity","dp/dt"
      close(10)
      open(20,file = "result_2D/1pressure.txt")
      write(20,'(1I1,",",1f24.16)') 0,G(3,162,Ny/2)!(23,7)を指定しているが実際は(22.89,6.97)にずれてしまう
      !p_inftyの定義
      pNx_infty = G(3,Nx,0)
      p0y_infty = G(3,0,0)
      pNy_infty = G(3,0,Ny)
        !粘性項の計算はCCSを用いるためA,L,U行列がsigma=0となる
        !そのためAp,Amなどとはまた別に設定する
        !CCS用のA,L,U行列はNの値がNx,Nyで異なるので別々に設定する
        call LU_Decompo(Nx,ccs_sigma,LUccsx)
        call LU_Decompo(Ny,ccs_sigma,LUccsy)
        !A行列の設定をF+とF-のそれぞれで行う
        !LU分解の過程は同じだが用いる値が異なるので別々にcall
        call LU_Decompo(Nx,psigma,LUpx)
        call LU_Decompo(Ny,psigma,LUpy)
        call LU_Decompo(Nx,msigma,LUmx)
        call LU_Decompo(Ny,msigma,LUmy)
        call Q_matrix(G,Q)
        Q0 = Q!Bufferの計算で使う初期値のQを保存
      !NSCBCのdGを求める微分も粘性項と同じCCSなのでA,L,U行列の使い回し可能
      do M =1,Mmax
        !x(i)を導出するのにf(i)を用いるがそれが時間ごとに変化するため幾つかのsubroutineは
        !このDoループ内に入れる必要がある。
            !移流方程式の導出
            !3次精度Runge-Kutta法での導出！
        !v方向から流入させる撹乱のためのsin波の流入条件を設定。時間変動させている
        theta = 2.d0*pi*dble(M)*dt
        in_G(2,:)=0.d0
        do i = 0,Ny
          if(zeta_fy(i) < -b) then
            in_G(2,i) = 0.d0
          elseif((zeta_fy(i) >= -b).and.(zeta_fy(i) < b)) then
            in_G(2,i) = A2*sin(T2*theta)!1秒で1周期になる。1/4*θなら4秒で1周期
          elseif(zeta_fy(i) >= b) then
            in_G(2,i) = 0.d0
          endif
        enddo
        !Q1
        !F行列のdfx/dxの計算
        !x_axis
        call Fx_matrix(Q,Fpx,Fmx)
        call dif_x(psigma,dx,Fpx,xp,LUpx,dzeta_inx)
        call dif_x(msigma,dx,Fmx,xm,LUmx,dzeta_inx)
        !y_axis dFy/dyの計算
        call Fy_matrix(Q,Fpy,Fmy)
        call dif_y(psigma,dy,Fpy,yp,LUpy,dzeta_iny)
        call dif_y(msigma,dy,Fmy,ym,LUmy,dzeta_iny)
        !NSCBCの境界条件を適用させるためにQ1を求める前にdFを定義してその両端に境界条件を
        !適用しなければならない
                dFx(:,:,:) = xm(:,:,:) + xp(:,:,:)
                dFy(:,:,:) = ym(:,:,:) + yp(:,:,:)
        !粘性項V行列のdv/dxの計算
        !まずはdu/dx,dT/dxの導出とμの設定
        call variable_setting(UVT,Q,myu)
        !x_axis dUVT/dxの計算
        call dif_x(ccs_sigma,dx,UVT,dUVTx,LUccsx,dzeta_inx)
        !y_axis dUVT/dyの計算
        call dif_y(ccs_sigma,dy,UVT,dUVTy,LUccsy,dzeta_iny)
        !V行列を構成する値が揃ったのでV行列の設定とdV/dxの計算
        !x_axis
        call Vx_matrix(Vx,myu,UVT,dUVTx,dUVTy)
        call dif_x(ccs_sigma,dx,Vx,dVx,LUccsx,dzeta_inx)
        !y_axis
        call Vy_matrix(Vy,myu,UVT,dUVTx,dUVTy)
        call dif_y(ccs_sigma,dy,Vy,dVy,LUccsy,dzeta_iny)
        !NSCBCの計算開始
        !x方向のNSCBCの計算
        call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
        call NSCBC_x(G,dGx,dFx,pNx_infty)
        !y方向
        call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
        call NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
        !無反射流出条件の際の境界での粘性項の条件を設定
        call outflow(UVT,dUVTx,Vx,dVx,dUVTy,Vy,dVy)
        !Buffer領域の計算
        !計算に必要なdQ/dx,dQ/dyをCCSで導出。今まではdQ/dtしか求めていなかった
        !dQ/dx,dQ/dy自体はdG/dx,dG/dyを組み合わせて作ることができるのでそうして作成すると微分せずに済み計算の短縮に繋がる
        !しかしEtの微分だけは足し算が多くなってしまい遅くなるのでEtはそれで新規に微分、Qのそれ以外の値は組み合わせて作るのが最善
        !だがこの場合Etは単体で定義していないので今回はそれができない。そのためQを直接微分する
        call dif_x(ccs_sigma,dx,Q,dQx,LUccsx,dzeta_inx)
        call dif_y(ccs_sigma,dy,Q,dQy,LUccsy,dzeta_iny)
       !計算して求めたdF,dVそしてBuffer領域の計算のための値などを組み合わせ、代入してdQ/dtを求める
          Q1(:,:,:) = Q(:,:,:) + c*dt*(dVx(:,:,:) + dVy(:,:,:) - dFx(:,:,:)&
              & - dFy(:,:,:)-sigma_x(:,:,:)*(Q(:,:,:)-Q0(:,:,:))&
              -sigma_y(:,:,:)*(Q(:,:,:)-Q0(:,:,:))-Ux(:,:,:)*dQx(:,:,:)-Uy(:,:,:)*dQy(:,:,:))
      !i=0で流入条件させるのでその部分のQ1を上書きして流入させ続ける
      call inflow(Q1,in_G)!dirichlet条件で流入部を固定
      !Q2(Q,F,x+-,y+-,f+-はそれぞれの計算過程において分ける必要がある。
      !またL,Uなどは DCSという方法が変わらないので同じものを使用できる)
      !dF/dxの計算
      Fpx=0.d0;Fmx=0.d0;xp=0.d0;xm=0.d0;Fpy=0.d0;Fmy=0.d0;yp=0.d0;ym=0.d0
      call Fx_matrix(Q1,Fpx,Fmx)
      call dif_x(psigma,dx,Fpx,xp,LUpx,dzeta_inx)
      call dif_x(msigma,dx,Fmx,xm,LUmx,dzeta_inx)
      !y_axis dFy/dyの計算
      call Fy_matrix(Q1,Fpy,Fmy)
      call dif_y(psigma,dy,Fpy,yp,LUpy,dzeta_iny)
      call dif_y(msigma,dy,Fmy,ym,LUmy,dzeta_iny)
      !NSCBCの境界条件を適用させるためにQ1を求める前にdFを定義してその両端に境界条件を
      !適用しなければならない
              dFx(:,:,:) = xm(:,:,:) + xp(:,:,:)
              dFy(:,:,:) = ym(:,:,:) + yp(:,:,:)
        !粘性項V行列のdv/dxの計算
        !まずはdu/dx,dT/dxの導出とμの設定
        call variable_setting(UVT,Q1,myu)
        !x_axis dUVT/dxの計算
        call dif_x(ccs_sigma,dx,UVT,dUVTx,LUccsx,dzeta_inx)
        !y_axis dUVT/dyの計算
        call dif_y(ccs_sigma,dy,UVT,dUVTy,LUccsy,dzeta_iny)
        !V行列を構成する値が揃ったのでV行列の設定とdV/dxの計算
        !x_axis
        call Vx_matrix(Vx,myu,UVT,dUVTx,dUVTy)
        call dif_x(ccs_sigma,dx,Vx,dVx,LUccsx,dzeta_inx)
        !y_axis
        call Vy_matrix(Vy,myu,UVT,dUVTx,dUVTy)
        call dif_y(ccs_sigma,dy,Vy,dVy,LUccsy,dzeta_iny)
        !NSCBCの計算開始
        call rho_u_p(G,Q1)
        !x方向のNSCBCの計算
        call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
        call NSCBC_x(G,dGx,dFx,pNx_infty)
        !y方向
        call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
        call NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
        call outflow(UVT,dUVTx,Vx,dVx,dUVTy,Vy,dVy)
        !Buffer領域の計算
        call dif_x(ccs_sigma,dx,Q1,dQx,LUccsx,dzeta_inx)
        call dif_y(ccs_sigma,dy,Q1,dQy,LUccsy,dzeta_iny)
        Q2(:,:,:) = (0.75d0)*Q(:,:,:) +(0.25d0) * Q1(:,:,:) + c* (dt*0.25d0)&
            &*(dVx(:,:,:) + dVy(:,:,:) - dFx(:,:,:) - dFy(:,:,:)-sigma_x(:,:,:)*(Q1(:,:,:)-Q0(:,:,:))&
            -sigma_y(:,:,:)*(Q1(:,:,:)-Q0(:,:,:))-Ux(:,:,:)*dQx(:,:,:)-Uy(:,:,:)*dQy(:,:,:))
        call inflow(Q2,in_G)
      !Qn
      !dF/dxの計算
      Fpx=0.d0;Fmx=0.d0;xp=0.d0;xm=0.d0;Fpy=0.d0;Fmy=0.d0;yp=0.d0;ym=0.d0
      call Fx_matrix(Q2,Fpx,Fmx)
      call dif_x(psigma,dx,Fpx,xp,LUpx,dzeta_inx)
      call dif_x(msigma,dx,Fmx,xm,LUmx,dzeta_inx)
      !y_axis dFy/dyの計算
      call Fy_matrix(Q2,Fpy,Fmy)
      call dif_y(psigma,dy,Fpy,yp,LUpy,dzeta_iny)
      call dif_y(msigma,dy,Fmy,ym,LUmy,dzeta_iny)
      !NSCBCの境界条件を適用させるためにQ2を求める前にdFを定義してその両端に境界条件を
      !適用しなければならない
              dFx(:,:,:) = xm(:,:,:) + xp(:,:,:)
              dFy(:,:,:) = ym(:,:,:) + yp(:,:,:)
        !粘性項V行列のdv/dxの計算
        !まずはdu/dx,dT/dxの導出とμの設定
        call variable_setting(UVT,Q2,myu)
        !x_axis dUVT/dxの計算
        call dif_x(ccs_sigma,dx,UVT,dUVTx,LUccsx,dzeta_inx)
        !y_axis dUVT/dyの計算
        call dif_y(ccs_sigma,dy,UVT,dUVTy,LUccsy,dzeta_iny)
        !V行列を構成する値が揃ったのでV行列の設定とdV/dxの計算
        !x_axis
        call Vx_matrix(Vx,myu,UVT,dUVTx,dUVTy)
        call dif_x(ccs_sigma,dx,Vx,dVx,LUccsx,dzeta_inx)
        !y_axis
        call Vy_matrix(Vy,myu,UVT,dUVTx,dUVTy)
        call dif_y(ccs_sigma,dy,Vy,dVy,LUccsy,dzeta_iny)
        !NSCBCの計算開始
        call rho_u_p(G,Q2)
        !x方向のNSCBCの計算
        call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
        call NSCBC_x(G,dGx,dFx,pNx_infty)
        !y方向
        call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
        call NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
        call outflow(UVT,dUVTx,Vx,dVx,dUVTy,Vy,dVy)
        !Buffer領域の計算
        call dif_x(ccs_sigma,dx,Q2,dQx,LUccsx,dzeta_inx)
        call dif_y(ccs_sigma,dy,Q2,dQy,LUccsy,dzeta_iny)
            Qn(:,:,:) = Q(:,:,:) /3.d0 + (2.d0/3.d0) * Q2(:,:,:) + c* &
            &((2.d0*dt)/3.d0) *(dVx(:,:,:) + dVy(:,:,:) - dFx(:,:,:) - dFy(:,:,:)-sigma_x(:,:,:)*(Q2(:,:,:)-Q0(:,:,:))&
            -sigma_y(:,:,:)*(Q2(:,:,:)-Q0(:,:,:))-Ux(:,:,:)*dQx(:,:,:)-Uy(:,:,:)*dQy(:,:,:))
        call inflow(Qn,in_G)
        call rho_u_p(G,Qn)
        if (mod(M,p_output) == 0) then
          write(20,'(1I7,",",1f24.16)') M,G(3,162,Ny/2)
        endif
          if(mod(M,output_count) == 0) then!dt=1.d-4で0.01秒刻みで出力するためにMの条件を設定
            !渦度用のdGの計算
            dGx=0.d0;dGy=0.d0
            call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
            call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
            omega_z(:,:) = dGx(2,:,:) - dGy(1,:,:)
            call rho_u_p(oldG,Q)
            dp(:,:) = (G(3,:,:) - oldG(3,:,:))/dt
            write(filename, '(i6.6)') M
            !Mの計算毎に出力ファイル名を変更して出力する
            !i5.5で5桁分の数字を表示できるのでdt=1.d-5以下で計算するならここも変更が必要
            open(10, file = "result_2D/parameter"//trim(filename)//".txt")
            do i = 0,Ny
              do j = 0,Nx
                write(10,'(f24.16,",",f24.16,",",f24.16,",",f24.16,",",f24.16)') &
                zeta_fx(j),zeta_fy(i),G(0,j,i),omega_z(j,i),dp(j,i)
              enddo
              write(10,*)
              !一度に全てを出力する際にはデータの切れ目として空白を一行挿入しなくてはいけない
            enddo
            write(10,'(2A1,1I7)') "#","M",M
            write(10,'(6A10)')"#","x","y","rho","vorticity","dp/dt"
            close(10)
          endif
        !計算が破綻している場合に計算を終了させるプログラム
          do j = 0,Ny
            do i = 0,Nx
              if(isnan(Qn(0,i,j))) then
                write(*,*) "x=",i,"y=",j,"M=",M
                call cpu_time(t1)
                write(*,'("Time required = ",i3,"min",f4.1,"sec")') &
                &int((((t1-t0) - mod(t1-t0,60.d0)) /60.d0)), mod(t1-t0,60.d0)
                stop "rho becomes NAN"
              endif
            enddo
          enddo
          !RK法の時間の更新
              Q(:,:,:) = Qn(:,:,:)
        if(mod(M,1) == 0) then
          write(*,*) "M=",M!計算に時間がかかるので進行状況の確認用に出力
        endif
      enddo
      call cpu_time(t1)
      write(*,'("Time required = ",i3,"min",f4.1,"sec")') &
      &int((((t1-t0) - mod(t1-t0,60.d0)) /60.d0)), mod(t1-t0,60.d0)
      write(20,'(2A1,1f24.16,1A1,1f24.16)') "#","x",zeta_fx(Nx/2),"y",zeta_fy(Ny/2)
      write(20,'(3A15)')"#","M","Local Pressure"
      close(20)
    end program main
