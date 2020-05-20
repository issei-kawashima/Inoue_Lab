!dx,dyの格子伸長の倍率は0.5倍/1.2倍に設定してある
!2019.05.18 Nx=180 Ny=100 dt=5.d-3で実行すれば卒論の格子伸長適用時の条件になる
!2019.05.20 Nx,Ny,dtをプログラム内で定義する仕様に変更開始>終了!!
!>> 496行目付近でUxを確実に0にしている部分の座標をNxを変えるなら確認しないといけない
!また定置観測している箇所の座標もNx,Nyを変えるなら調節しないといけないwrite文二つ使用
!2019.05.21　３次元対応のコードの作成を開始
!2019.05.23 UVWT subrouitineでmyuの設定をUVWT(4)を使用した簡単な形式に変更した。後日それでいいか確認する>OK
!2019.05.24 Q,Fxyz,Vxyz,UVWT,rho_u_pのsubroutineは三次元用に変更ずみ
!2019.05.25 dif_x,yの3次元化完了. NSCBC_x,yのL,d,dFも変更完了
!2019.05.25 inflow,outflowのsubroutineは行列の次元は変更したが正しいかどうかの確証がない
!2019.05.25 dif_zは作成したが周期条件でのsubroutineとなるため行列の次元の設定等が未完了
!2019.06.12 z方向周期用のLU subroutineとdif_zが完成　別プログラムで計算させて確認済み。また周期方向のLU分解も圧縮済み
!2019.06.20 計算結果の出力に関して渦度を出していたが誤ったG(Q2から導出したG)でdGx,dGyを計算していたのでQnから計算したdGx,dGyで渦度を求めるように変更した
!2019.06.20 渦度を三次元のものに変更した。ただ成分が三つなのでomega_1~3の３つの配列に分けた。
!2019.06.20 NSCBCではなくneumannにしたらM=101まで回るようになった。しかしM=102(dt=5.d-3)で毎回Nanとなるdt=2.d-3に変えるとM=109でNaNとなった。
!2019.06.22 NSCBCで光速の設定でpをwで計算していた。またdpの計算もwになっていたので修正した。
!2019.06.22 Runge-Kutta+NSCBCで試行。M=120でnan
!2019.06.24 Euler+NSCBCで試行 M=319でNan => Nan直前の M=318の密度と渦度の出力の仕方を見つけて可視化したが、x=Nxの流出境界で密度が振動し、
!0になってしまいNanになったと見られる。その直前までは徐々に流入・移流している様子がz=Nz/2で観測できているのでz方向の微分自体はできているように思われる。
!2019.07.29 Euler+NeumannでBufferなし、Nx=180,Ny=100で施行。M=116でNan
!2019.12.10 Runge-Kutta+NSCBC+Bufer Nx=180 Ny=100 Nz=20,dt =2.d-3で計算してM=307でNaN
!2020.04.23 Fzの定義式が正しいことを確認
!2020.04.30 outflowのsubroutineが抜けていたことを発見。τとqの要素が=0となる条件で3次元で計算し直し、再定義。
 !>>抜けていたわけではなく、NSCBCのsubroutineに組み込んでいた。わかりにくいのでoutflow subroutineとして分けたままにした。(05.10)
!2020.04.30 Runge-Kutta+NSCBC+Bufer Nx=180, Ny=120, Nz=20, dt=2.d-3で計算開始。M=307でNaN
!2020.05.01 可視化するとz=0でのみrhoの値がx=0で1,x=>1で6から徐々に下がるという挙動をしていたので流入値の設定の問題のような気がした。
!2020.05.01 4/30と同条件&in_G(0)は1/Tu,　in_G(1)はurのジェットの条件に変更して計算。(それ以外に信頼できる流入条件を知らなかった)
!2020.05.01 M=307でNaN ただし、rhoの流入に関しては正常になったように見える。またどのz面で可視化しても結果は同じだった。
!2020.05.10 Tjet=1.12Tempでは適性膨張ジェットなのでTjet=1.4*Tempへ変更した。渡辺さん・姫野さんはMa=1.4で実験と計算を行なっているのでMa=1.4を採用
!2020.05.10 dt=5.d-3,NSCBC_x>Neumannにして計算。超音速なので、x流出方向の大体は超音速になり、逆流が起こらないと仮定した。
!2020.05.10 この条件で2次元で計算したら、上手くできるので、3次元ならどうか挑戦。M=89でNan。Z方向は全て同じ&少しづつ流入しているので何が問題かわからない。
!2020.05.13 格子伸長を外したら計算がより早く死んだので、格子幅の問題かと思ってNxを細くしてみた。dt=2.d-3, Nx=360,Ny=100, Ma=0.5,Tjet=1.12Tempで計算。
!2020.05.13 Nxを変更したのでBuffer_xのUxも90から180まで0にする範囲を変更した。
!2020.05.13 Ma数を亜音速にしたので、Neumannを廃止、NSCBC_xを復活させた。Ma=135でNan.
!2020.05.13 Ma=0.5,Nx=180,Ny100,Tjet=1.12Temp,dt=2.d-3,NSCBC_x&y,Ux(0~90)=0で計算。5/1と同条件なので同じ場所M=307でNan.
!2020.05.13 格子伸長を外したらより早くNanになったので、格子幅をより細かくしたら行けるかも。。
!ただそのためにはメモリが足りないのでallocatableに配列を変更しないといけない。
!2020.05.14 一旦allocateは変更しないで、Nx=360,その他は5/13(最後)と同様にして計算。M=135でNan.精度あげたのに...
!2020.05.14 格子幅をx方向は細かくしたのに、より早くNanになったので孔子幅の問題ではないと仮定して、配列の動的設定(allocate)に関しては一旦辞める。
!2020.05.20 周期条件のコードが間違ってるかもしれないので、その部分を見直す。
!2020.05.20 ああああああああああ
module threedim
  !連続の式、Eulerの運動方程式、エネルギー方程式を並列に並べた行列Q,Fの設定等をする
  !これらの式をまとめて基礎式と呼ぶ
  implicit none
  !計算条件値の設定
  double precision,parameter :: gamma = 1.4d0
  integer,parameter :: t_end = 150 !時刻tの設定
  integer,parameter :: p_output = 10 !時間毎の局所圧力を出力させる際のステップ間隔
  integer,parameter :: Nx = 180
  integer,parameter :: Ny = 100
  integer,parameter :: Nz = 20
  double precision,parameter :: dt = 2.d-3
  integer,parameter :: NUx = 90!buffer_xのUxで流入側のUxを0にする座標(格子点番号)Nx=180ならNUx=90
  integer,parameter :: Mmax = t_end / dt
  integer,parameter :: output_count = int(0.1d0/dt)!出力ファイルを1sec間隔で出力するように設定
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
  double precision,parameter :: Lz = 1.d0

  double precision,parameter :: psigma = -0.25d0
  double precision,parameter :: msigma = 0.25d0
  double precision,parameter :: NS_sigma = 0.25d0
  double precision,parameter :: ccs_sigma = 0.d0
  double precision,parameter :: c = 1.d0
  !亜音速流入のためRe数は小さめに
  double precision,parameter :: Pr = 0.71d0
  double precision,parameter :: Ma = 0.5d0
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
    double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Q,G
      Q(0,:,:,:) = G(0,:,:,:)
      Q(1,:,:,:) = G(0,:,:,:) * G(1,:,:,:)
      Q(2,:,:,:) = G(0,:,:,:) * G(2,:,:,:)
      Q(3,:,:,:) = G(0,:,:,:) * G(3,:,:,:)
      Q(4,:,:,:) = G(4,:,:,:) / (gamma -1.d0) + G(0,:,:,:) * (G(1,:,:,:)**2.d0+&
                  G(2,:,:,:)**2.d0 +G(3,:,:,:)**2.d0) * 0.5d0
  endsubroutine Q_matrix
  !求めQを用いてF行列の設定
  subroutine F_matrix(Q,Fpx,Fmx,Fpy,Fmy,Fpz,Fmz)
    double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Q,Fx,Fpx,Fmx,Fy,Fpy,Fmy,Fz,Fpz,Fmz
    Fx=0.d0
    Fy=0.d0
    Fz=0.d0
    !F行列の設定(x方向)
      Fx(0,:,:,:) = Q(1,:,:,:)
      Fx(1,:,:,:) = 1.d0/(2.d0*Q(0,:,:,:))*((3.d0-gamma)*(Q(1,:,:,:)**2.d0)+&
                  (1.d0-gamma)*(Q(2,:,:,:)**2.d0 + Q(3,:,:,:)**2.d0))+(gamma-1.d0)*Q(4,:,:,:)
      Fx(2,:,:,:) = Q(1,:,:,:)*Q(2,:,:,:)/Q(0,:,:,:)
      Fx(3,:,:,:) = Q(1,:,:,:)*Q(3,:,:,:)/Q(0,:,:,:)
      Fx(4,:,:,:) = gamma * Q(1,:,:,:) * Q(4,:,:,:) / Q(0,:,:,:) &
                 &+ (1.d0- gamma)*Q(1,:,:,:)*(Q(1,:,:,:)**2.d0+Q(2,:,:,:)**2.d0&
                 +Q(3,:,:,:)**2.d0)/(2.d0 * (Q(0,:,:,:)**2.d0))
    !求めたFを特製速度の正負によって分割する
    !Lax-Friedrichの流速分割を用いる
        Fpx(:,:,:,:) = (Fx(:,:,:,:) + zeta * Q(:,:,:,:)) * 0.5d0
        Fmx(:,:,:,:) = (Fx(:,:,:,:) - zeta * Q(:,:,:,:)) * 0.5d0

    !F行列の設定(y方向)
          Fy(0,:,:,:) = Q(2,:,:,:)
          Fy(1,:,:,:) = Q(2,:,:,:)*Q(1,:,:,:)/Q(0,:,:,:)
          Fy(2,:,:,:) = 1.d0/(2.d0*Q(0,:,:,:))*((3.d0 - gamma)*(Q(2,:,:,:)**2.d0) + &
                        &(1.d0 - gamma)*(Q(1,:,:,:)**2.d0+Q(3,:,:,:)**2.d0))+(gamma - 1.d0)*Q(4,:,:,:)
          Fy(3,:,:,:) = Q(2,:,:,:)*Q(3,:,:,:)/Q(0,:,:,:)
          Fy(4,:,:,:) = gamma * Q(2,:,:,:) * Q(4,:,:,:) / Q(0,:,:,:) &
                    &+ (1.d0- gamma)*Q(2,:,:,:)*(Q(1,:,:,:)**2.d0+Q(2,:,:,:)**2.d0&
                    +Q(3,:,:,:)**2.d0)/(2.d0 * (Q(0,:,:,:)**2.d0))
      !求めたFを特製速度の正負によって分割する
      !Lax-Friedrichの流速分割を用いる
          Fpy(:,:,:,:) = (Fy(:,:,:,:) + zeta * Q(:,:,:,:)) * 0.5d0
          Fmy(:,:,:,:) = (Fy(:,:,:,:) - zeta * Q(:,:,:,:)) * 0.5d0

      !F行列の設定(z方向)
            Fz(0,:,:,:) = Q(3,:,:,:)
            Fz(1,:,:,:) = Q(3,:,:,:)*Q(1,:,:,:)/Q(0,:,:,:)
            Fz(2,:,:,:) = Q(3,:,:,:)*Q(2,:,:,:)/Q(0,:,:,:)
            Fz(3,:,:,:) = 1.d0/(2.d0*Q(0,:,:,:))*((3.d0 - gamma)*(Q(3,:,:,:)**2.d0) + &
            &(1.d0 - gamma)*(Q(1,:,:,:)**2.d0+Q(2,:,:,:)**2.d0))+(gamma - 1.d0)*Q(4,:,:,:)
            Fz(4,:,:,:) = gamma * Q(3,:,:,:) * Q(4,:,:,:) / Q(0,:,:,:) &
                      &+ (1.d0- gamma)*Q(3,:,:,:)*(Q(1,:,:,:)**2.d0+Q(2,:,:,:)**2.d0&
                      +Q(3,:,:,:)**2.d0)/(2.d0 * (Q(0,:,:,:)**2.d0))
        !求めたFを特製速度の正負によって分割する
        !Lax-Friedrichの流速分割を用いる
            Fpz(:,:,:,:) = (Fz(:,:,:,:) + zeta * Q(:,:,:,:)) * 0.5d0
            Fmz(:,:,:,:) = (Fz(:,:,:,:) - zeta * Q(:,:,:,:)) * 0.5d0
      endsubroutine F_matrix
!du/dx,dT.dxを求めるためにまずはuとTの設定。一行目は他と係数を合わせ易くするために0だが含めている
!UVWT,myuの計算はrho,u,pをそのまま代入しても良いがその場合は求めたtQから毎回rho_u_pのsubroutineを
!呼び出して計算しなければいけないので今回はQから直接計算できるようにプログラムを組んだ
    subroutine variable_setting(UVWT,Q,myu)
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: UVWT,Q
      double precision,dimension(0:Nx,0:Ny,0:Nz-1) :: myu
        UVWT(0,:,:,:) = 0.d0
        UVWT(1,:,:,:) = Q(1,:,:,:) / Q(0,:,:,:)!u
        UVWT(2,:,:,:) = Q(2,:,:,:) / Q(0,:,:,:)!v
        UVWT(3,:,:,:) = Q(3,:,:,:) / Q(0,:,:,:)!w
        UVWT(4,:,:,:) = (Ma**2.d0)*(gamma * (gamma - 1.d0) * (Q(4,:,:,:) - (((Q(1,:,:,:) **2.d0)+ &
        &(Q(2,:,:,:) **2.d0)+(Q(3,:,:,:)**2.d0)) / (2.d0 * Q(0,:,:,:))))) / Q(0,:,:,:)!T
        !Tの値はMa^2*gamma*p/rhoこれは速度で無次元化したもの
!V行列を設定する際にμの計算が複雑になっているのでそれを簡略に示すために別でμを計算するsubroutine
!UVWTからTの値を代入することで計算を簡略化している
        myu(:,:,:) = (UVWT(4,:,:,:) ** 1.5d0) * (1.d0 + Sc) / (UVWT(4,:,:,:) + Sc)
    end subroutine variable_setting

    subroutine V_matrix(Vx,Vy,Vz,myu,UVWT,dUVWTx,dUVWTy,dUVWTz)
      !粘性項の設定(x方向)
      !基礎式右辺をV行列として設定
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Vx,Vy,Vz,dUVWTx,dUVWTy,dUVWTz,UVWT
      double precision,dimension(0:Nx,0:Ny,0:Nz-1) :: myu
        Vx(0,:,:,:) = 0.d0
        Vx(1,:,:,:) = (2.d0*myu(:,:,:)/(3.d0*Re)) * (2.d0 * dUVWTx(1,:,:,:) - dUVWTy(2,:,:,:)-dUVWTz(3,:,:,:))
        Vx(2,:,:,:) = (myu(:,:,:) / Re) * (dUVWTx(2,:,:,:) + dUVWTy(1,:,:,:))
        Vx(3,:,:,:) = (myu(:,:,:) / Re) * (dUVWTx(3,:,:,:) + dUVWTz(1,:,:,:))
        Vx(4,:,:,:) =Vx(1,:,:,:)*UVWT(1,:,:,:) + Vx(2,:,:,:)* UVWT(2,:,:,:)+ &
                      Vx(3,:,:,:)*UVWT(3,:,:,:)+((myu(:,:,:) * dUVWTx(4,:,:,:))&
                      / ((gamma - 1.d0)*Re*Pr*(Ma ** 2.d0)))
      !粘性項の設定(y方向)
        Vy(0,:,:,:) = 0.d0
        Vy(1,:,:,:) = (myu(:,:,:) / Re) * (dUVWTx(2,:,:,:) + dUVWTy(1,:,:,:))
        Vy(2,:,:,:) = (2.d0*myu(:,:,:)/(3.d0*Re))*(-dUVWTx(1,:,:,:) + 2.d0*dUVWTy(2,:,:,:)-dUVWTz(3,:,:,:))
        Vy(3,:,:,:) = (myu(:,:,:) / Re) * (dUVWTy(3,:,:,:) + dUVWTz(2,:,:,:))
        Vy(4,:,:,:) = Vy(1,:,:,:)*UVWT(1,:,:,:) + Vy(2,:,:,:)* UVWT(2,:,:,:)+ &
                      Vy(3,:,:,:)*UVWT(3,:,:,:)+((myu(:,:,:) * dUVWTy(4,:,:,:))&
                      / ((gamma - 1.d0)*Re*Pr*(Ma ** 2.d0)))
      !粘性項の設定(z方向)
        Vz(0,:,:,:) = 0.d0
        Vz(1,:,:,:) = (myu(:,:,:) / Re) * (dUVWTx(3,:,:,:) + dUVWTz(1,:,:,:))
        Vz(2,:,:,:) = (myu(:,:,:) / Re) * (dUVWTy(3,:,:,:) + dUVWTz(2,:,:,:))
        Vz(3,:,:,:) = (2.d0*myu(:,:,:)/(3.d0*Re))*(-dUVWTx(1,:,:,:) - dUVWTy(2,:,:,:)+2.d0*dUVWTz(3,:,:,:))
        Vz(4,:,:,:) = Vz(1,:,:,:)*UVWT(1,:,:,:) + Vz(2,:,:,:)* UVWT(2,:,:,:)+ &
                      Vz(3,:,:,:)*UVWT(3,:,:,:)+((myu(:,:,:) * dUVWTz(4,:,:,:))&
                      / ((gamma - 1.d0)*Re*Pr*(Ma ** 2.d0)))
    endsubroutine V_matrix
    !求めたQnから時間毎のrho,u,v,w,pを求めるサブルーチン
    !ある指定した時間の時の計算結果のみを取り出せばグラフが作成できるので毎回は使用しない
    subroutine rho_u_p(G,Qn)
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Qn,G
        G=0.d0
        G(0,:,:,:) = Qn(0,:,:,:)
        G(1,:,:,:) = Qn(1,:,:,:) / Qn(0,:,:,:)
        G(2,:,:,:) = Qn(2,:,:,:) / Qn(0,:,:,:)
        G(3,:,:,:) = Qn(3,:,:,:) / Qn(0,:,:,:)
        G(4,:,:,:) = (gamma-1.d0)*(Qn(4,:,:,:)-(Qn(1,:,:,:)**2.d0+&
                      Qn(2,:,:,:)**2.d0+Qn(3,:,:,:)**2.d0)/(2.d0*Qn(0,:,:,:)))!p
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
      do i = 2,N-2
        A(i,i-1) = (1.d0 - Dalpha * sigma) * alpha5
        A(i,i) = 1.d0
        A(i,i+1) = (1.d0 + Dalpha * sigma) * alpha5
      enddo
      !N-1行目 3次精度DCS
      A(N-1,N-2) = alpha3 * (1.d0 - Dalpha * sigma)
      A(N-1,N-1) = 1.d0
      A(N-1,N) = alpha3 * (1.d0 + Dalpha * sigma)
      !N行目 片側DCS
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
      !1からN-2行目まで
      do i = 1,N-2
        A(i,i-1) = (1.d0 - Dalpha * sigma) * alpha5
        A(i,i) = 1.d0
        A(i,i+1) = (1.d0 + Dalpha * sigma) * alpha5
      enddo
      !N行目
      A(N-1,0) = alpha5 * (1.d0 + Dalpha * sigma)
      A(N-1,N-2) = alpha5*(1.d0-sigma*Dalpha)
      A(N-1,N-1) = 1.d0
    !LU分解
      !まずL,Uの初期値を設定
      do i = 0,N-1
        U(0,i) = A(0,i)
        L(i,0) = A(i,0) / U(0,0)
        L(i,i) = 1.d0
      enddo
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
      do i= 0,N-3
        LU(-2,i) = U(i,N-1)
        LU(2,i)  = L(N-1,i)
      enddo
      do i = 1,N-1
        LU(-1,i) = L(i,i-1)
      enddo
      do j = 0,1
        do i = 0,N-2
          LU(j,i) = U(i,i+j)
        enddo
      enddo
      deallocate(A,L,U)
    endsubroutine LU_DecompoPiriodic

  !DCS右辺の計算(RHS)サブルーチン
  !x方向
    subroutine dif_x(sigma,dx,Fx,dFzeta,LU,dzeta_inx)
      integer i
      double  precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: RHS_x,Fx
      double  precision,dimension(0:4,2:Nx-2,0:Ny,0:Nz-1) :: D2,D4,D6,D8
      double precision dx,sigma,dxinv
      double  precision,parameter :: ra = 14.d0/9.d0, rb = 1.d0/9.d0&
      &,da = 4.d0 / 9.d0,db = 2.d0 / 9.d0 !5次精度のDCSとなるための係数設定
      double precision,allocatable,dimension(:,:) :: LU
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: x,y
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: dzeta_inx,dFzeta
      D2=0.d0;D4=0.d0;D6=0.d0;D8=0.d0;RHS_x=0.d0;y=0.d0;x=0.d0
      !片側DCS,3次精度DCSも入れた非周期条件の際のbの設定
      dxinv = 1.d0/dx
        !片側DCSの右辺設定
        RHS_x(:,0,:,:)=((-17.d0/6.d0)*Fx(:,0,:,:)+(1.5d0)*Fx(:,1,:,:)+(1.5d0)&
                      *Fx(:,2,:,:)-Fx(:,3,:,:)/6.d0)*dxinv
        RHS_x(:,Nx,:,:)=((1.d0/6.d0)*Fx(:,Nx-3,:,:)-(1.5d0)* Fx(:,Nx-2,:,:)-&
                      (1.5d0)*Fx(:,Nx-1,:,:)+(17.d0/6.d0)*Fx(:,Nx,:,:))*dxinv
        !3次精度DCSの右辺設定
        RHS_x(:,1,:,:)=((1.5d0)*(-Fx(:,0,:,:)+Fx(:,2,:,:))*0.5d0*dxinv)+sigma*&
                      ((Fx(:,0,:,:)-2.d0*Fx(:,1,:,:)+Fx(:,2,:,:))*(0.5d0*dxinv))
        RHS_x(:,Nx-1,:,:)=((1.5d0)*(-Fx(:,Nx-2,:,:)+Fx(:,Nx,:,:))/(0.5d0*dxinv))+&
                      sigma*((Fx(:,Nx-2,:,:)-2.d0*Fx(:,Nx-1,:,:)+Fx(:,Nx,:,:))*(0.5d0*dxinv))
         !5次精度DCSの右辺設定
       do i = 2,Nx-2
         D2(:,i,:,:) = (-Fx(:,i-1,:,:)+Fx(:,i+1,:,:)) * (0.5d0*dxinv)
         D4(:,i,:,:) = (-Fx(:,i-2,:,:)+Fx(:,i+2,:,:)) * (0.25d0*dxinv)
         !D6(:,i,:,:) = (-Fx(:,i-3,:,:)+Fx(:,i+3,:,:)) / (6.d0*dx)
         !7次精度DCS用のため不要
         D6(:,i,:,:) = (Fx(:,i-1,:,:)+Fx(:,i+1,:,:)- 2.d0* Fx(:,i,:,:)) * dxinv
         D8(:,i,:,:) = (Fx(:,i-2,:,:)+Fx(:,i+2,:,:)- 2.d0* Fx(:,i,:,:)) * (0.25d0*dxinv)
         ! D12(:,i,:,:) = (Fx(:,i-3,:,:)+Fx(:,i+3,:,:)- 2.d0* Fx(:,i,:,:)) / (9.d0*dx)
         !7次精度DCS用のため不要
         RHS_x(:,i,:,:)=ra*D2(:,i,:,:)+rb*D4(:,i,:,:)+sigma*(da*D6(:,i,:,:)+db*D8(:,i,:,:))
       enddo
     !前進代入法、後退代入法の計算サブルーチン(x方向)
            !前進代入
            y(:,0,:,:) = RHS_x(:,0,:,:)!例外の境界値
            do i = 1,Nx!y(:,k)でk=0は上で定義したので残りの1〜Nxを定義する
              y(:,i,:,:) = RHS_x(:,i,:,:) - LU(-1,i)*y(:,i-1,:,:)
            enddo
            !後退代入
            x(:,Nx,:,:) = y(:,Nx,:,:) / LU(0,Nx)!例外の境界値
            !x(:,i)でi=Nxは定義したので残りのi=0~Nx-1を定義する
            do i = Nx-1, 0, -1!後退するので-1ずつ進む
              x(:,i,:,:) = (y(:,i,:,:) - LU(1,i)*x(:,i+1,:,:)) / LU(0,i)
            enddo
          call combine(dzeta_inx,x,dFzeta)
      end subroutine dif_x
   !DCS右辺の計算(RHS)サブルーチン
   !y方向
     subroutine dif_y(sigma,dy,Fy,dFzeta,LU,dzeta_iny)
       integer i
       double  precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: RHS_y,Fy
       double  precision,dimension(0:4,0:Nx,2:Ny-2,0:Nz-1) :: D2,D4,D6,D8
       double precision dy,sigma,dyinv
       double  precision,parameter :: ra = 14.d0/9.d0, rb = 1.d0/9.d0&
       &,da = 4.d0 / 9.d0,db = 2.d0 / 9.d0 !5次精度のDCSとなるための係数設定
       double precision,allocatable,dimension(:,:) :: LU
       double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: x,y
       double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: dzeta_iny,dFzeta
       D2=0.d0;D4=0.d0;D6=0.d0;D8=0.d0;RHS_y=0.d0;y=0.d0;x=0.d0
       dyinv = 1.d0 / dy
       !片側DCS,3次精度DCSも入れた非周期条件の際のbの設定
         !片側DCSの右辺設定
         RHS_y(:,:,0,:) = ((-17.d0/6.d0)*Fy(:,:,0,:)+(1.5d0)*Fy(:,:,1,:)+&
         (1.5d0)*Fy(:,:,2,:)-Fy(:,:,3,:)/6.d0)*dyinv
         RHS_y(:,:,Ny,:)=((1.d0/6.d0)*Fy(:,:,Ny-3,:)-(1.5d0)&
         &*Fy(:,:,Ny-2,:)-(1.5d0)*Fy(:,:,Ny-1,:)+(17.d0/6.d0)*Fy(:,:,Ny,:))*dyinv
         !3次精度DCSの右辺設定
         RHS_y(:,:,1,:)=((1.5d0)*(-Fy(:,:,0,:)+Fy(:,:,2,:))*(0.5d0*dyinv))+sigma*&
         ((Fy(:,:,0,:)-2.d0*Fy(:,:,1,:)+Fy(:,:,2,:))*(0.5d0*dyinv))
         RHS_y(:,:,Ny-1,:)=((1.5d0)*(-Fy(:,:,Ny-2,:)+Fy(:,:,Ny,:))*(0.5d0*dyinv))&
            +sigma*((Fy(:,:,Ny-2,:)-2.d0*Fy(:,:,Ny-1,:)+Fy(:,:,Ny,:))*(0.5d0*dyinv))
          !5次精度DCSの右辺設定
      do i = 2,Ny-2
        D2(:,:,i,:) = (-Fy(:,:,i-1,:)+Fy(:,:,i+1,:)) * (0.5d0*dyinv)
        D4(:,:,i,:) = (-Fy(:,:,i-2,:)+Fy(:,:,i+2,:)) * (0.25d0*dyinv)
        !D6(:,:,i,:) = (-Fy(:,:,i-3,:)+Fy(:,:,i+3,:)) / (6.d0*dy)
        !7次精度DCS用のため不要
        D6(:,:,i,:) = (Fy(:,:,i-1,:)+Fy(:,:,i+1,:)- 2.d0* Fy(:,:,i,:)) * dyinv
        D8(:,:,i,:) = (Fy(:,:,i-2,:)+Fy(:,:,i+2,:)- 2.d0* Fy(:,:,i,:)) * (0.25d0*dyinv)
        ! D12(:,:,i,:) = (Fy(:,:,i-3,:)+Fy(:,:,i+3,:)- 2.d0* Fy(:,:,i,:)) / (9.d0*dy)
        !7次精度DCS用のため不要
        RHS_y(:,:,i,:)=ra*D2(:,:,i,:)+rb*D4(:,:,i,:)+sigma*(da*D6(:,:,i,:)+db*D8(:,:,i,:))
      enddo
      !前進代入法、後退代入法の計算サブルーチン(y方向)
        !前進代入
        y(:,:,0,:) = RHS_y(:,:,0,:)
        do i = 1,Ny!y(:,:,i)でi=0は上で定義したので残りの1〜Nyを定義する
          y(:,:,i,:) = RHS_y(:,:,i,:) - LU(-1,i)*y(:,:,i-1,:)
        enddo
        !後退代入
        x(:,:,Ny,:) = y(:,:,Ny,:) / LU(0,Ny)
        !x(:,:,i)でi=Nyは定義したので残りのi=0~Ny-1を定義する
        do i = Ny-1, 0, -1!後退するので-1ずつ進む
          x(:,:,i,:) = (y(:,:,i,:) - LU(1,i)*x(:,:,i+1,:)) / LU(0,i)
        enddo
        call combine(dzeta_iny,x,dFzeta)!y方向の格子伸長を適用　微分変換をしている
       end subroutine dif_y

       subroutine dif_z(sigma,dz,Fz,x,LU)
         integer i,j,r,k,im1,im2,ip1,ip2
         double  precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: RHS_z,Fz
         double  precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: D2,D4,D6,D8
         double precision dz,sigma,dzinv,Lsum
         double  precision,parameter :: ra = 14.d0/9.d0, rb = 1.d0/9.d0&
         &,da = 4.d0 / 9.d0,db = 2.d0 / 9.d0 !5次精度のDCSとなるための係数設定
         double precision,allocatable,dimension(:,:) :: LU
         double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: x,y
         D2=0.d0;D4=0.d0;D6=0.d0;D8=0.d0;RHS_z=0.d0;y=0.d0;x=0.d0;Lsum=0.d0
         dzinv = 1.d0 / dz
           !5次精度DCSの右辺設定
           do i = 0,Nz-1
             im1 = mod(Nz+i-1, Nz)
             im2 = mod(Nz+i-2, Nz)
             ip1 = mod(Nz+i+1, Nz)
             ip2 = mod(Nz+i+2, Nz)
             D2(:,:,:,i) = (-Fz(:,:,:,im1)+Fz(:,:,:,ip1)) * (0.5d0*dzinv)
             D4(:,:,:,i) = (-Fz(:,:,:,im2)+Fz(:,:,:,ip2)) * (0.25d0*dzinv)

             D6(:,:,:,i) = (Fz(:,:,:,im1)+Fz(:,:,:,ip1)- 2.d0* Fz(:,:,:,i)) * dzinv
             D8(:,:,:,i) = (Fz(:,:,:,im2)+Fz(:,:,:,ip2)- 2.d0* Fz(:,:,:,i)) * (0.25d0*dzinv)
           RHS_z(:,:,:,i)=ra*D2(:,:,:,i)+rb*D4(:,:,:,i)+sigma*(da*D6(:,:,:,i)+db*D8(:,:,:,i))
           enddo
           !前進代入法、後退代入法の計算サブルーチン(x方向)
              !前進代入
              y(:,:,:,0) = RHS_z(:,:,:,0)!例外の境界値
              do i = 1,Nz-2!y(:,k)でk=0は上で定義したので残りの1〜Nz-2を定義する
                !これはfill-inのない通常部分
                y(:,:,:,i) = RHS_z(:,:,:,i) - LU(-1,i)*y(:,:,:,i-1)
              enddo
              !Fill-in部はΣの計算が必要になるので追加 Lの最後の行i=Nz-1のみ別で計算
              do k= 0,Ny
                do r = 0,Nx
                  do j= 0,4
                    Lsum = 0.d0
                    do i = 0,Nz-3
                      Lsum = Lsum +LU(2,i)*y(j,r,k,i)
                    enddo
                      y(j,r,k,Nz-1) = RHS_z(j,r,k,Nz-1) - LU(-1,Nz-1)*y(j,r,k,Nz-2)-Lsum
                  enddo
                enddo
              enddo
              !後退代入
              x(:,:,:,Nz-1) = y(:,:,:,Nz-1) / LU(0,Nz-1)!例外の境界値
              !fill-inの計算の範囲外なので別で計算
              x(:,:,:,Nz-2) = (y(:,:,:,Nz-2) -LU(1,Nz-2)*x(:,:,:,Nz-1))/ LU(0,Nz-2)
              !x(:,i)でi=Nxは定義したので残りのi=0~Nx-1を定義する
              do i = Nz-3, 0, -1!後退するので-1ずつ進む
                x(:,:,:,i) = (y(:,:,:,i) - LU(1,i)*x(:,:,:,i+1)-LU(-2,i)*x(:,:,:,Nz-1)) / LU(0,i)
              enddo
       endsubroutine dif_z
    !境界条件をNSCBCで設定
    !L,d行列を設定することでdfx(0:2,0)とdFx(0:2,Nx)の値を定める
    !定めたdFxをxm,xpをそれぞれ求めた時点でi=0とNxで置き換える
    !これが境界条件となる。neumannはNx-1とNx,1と0を同じにしたがNSCBC_xはdFxの境界を定めることで
    !境界条件を設定している
    !そしてそのままQ1,Q2,Qnを求める
    !まずはx方向用のNSCBC　subrouitineを作成
    !(:,:,:)等で計算を行わせる際は:の配列数が対応していることが要確認
    subroutine NSCBC_x(G,dGx,dFx,pNx_infty)
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: G,dGx,dFx
      double precision,dimension(1:5,0:Ny,0:Nz-1,0:1) :: L,d
      double precision,dimension(0:Ny,0:Nz-1,0:1) :: c_NS,Ma_NS
      double precision pNx_infty
      L=0.d0;d=0.d0;c_NS=0.d0;Ma_NS=0.d0
      !音速cはi=0,Nxの両点においてそれぞれ定義しなければならない
      c_NS(:,:,0) = sqrt(gamma * G(4,0,:,:) / G(0,0,:,:))
      c_NS(:,:,1) = sqrt(gamma * G(4,Nx,:,:) / G(0,Nx,:,:))
      !マッハ数Ma_NSはi=0,Nxで使うので別々に定義する
      Ma_NS(:,:,0) = G(1,0,:,:) / c_NS(:,:,0)!uを使う
      Ma_NS(:,:,1) = G(1,Nx,:,:) / c_NS(:,:,1)
      !x方向右側つまりi=0の点において亜音速流入条件でL行列を設定する
      L(1,:,:,0)=(G(1,0,:,:)-c_NS(:,:,0))*(-G(0,0,:,:)*c_NS(:,:,0)*dGx(1,0,:,:)+dGx(4,0,:,:))
!      L(2,:,:,0)=G(1,0,:,:)*((c_NS(:,:,0)**2.d0)*dGx(0,0,:,:)-dGx(3,0,:,:))
!      L(3,:,:,0)=0.d0!L3=u*dv/dxで流入速度
      L(3,:,:,0)=0.d0!G(1,0,:,:)*dGx(2,0,:,:)
      L(4,:,:,0)=0.d0!G(1,0,:,:)*dGx(3,0,:,:)
      L(5,:,:,0)=L(1,:,:,0) !-2.d0*c_NS(:,:,0)*du/dtが本来はあるが
      !流入速度uを時間変動させないので今回はdu/dt=0となるため省略
!      L(5,:,:,0)=msigma*c_NS(:,:,0)*(1.d0-(Ma_NS(:,:,0)**2.d0))*(G(3,0,:,:)-p0x_infty)/Lx
      L(2,:,:,0)=(0.5d0)*(gamma-1.d0)*(L(5,:,:,0)+L(1,:,:,0))!+G(0,0,:,:)*c_NS(:,:,0)**2.d0/T*dT/dtが本来はあるが
      !境界温度Tを時間変動させないので今回はdT/dt=0となるため省略
      !x方向左側つまりi=Nxの点において無反射流出条件でL行列を設定する
      L(1,:,:,1)=NS_sigma*c_NS(:,:,1)*(1.d0-(Ma_NS(:,:,1)**2.d0))*(G(4,Nx,:,:)-&
      &pNx_infty)/Lx
      L(2,:,:,1)=G(1,Nx,:,:)*((c_NS(:,:,1)**2.d0)*dGx(0,Nx,:,:)-dGx(4,Nx,:,:))
      L(3,:,:,1)=G(1,Nx,:,:)*dGx(2,Nx,:,:)
      L(4,:,:,1)=G(1,Nx,:,:)*dGx(3,Nx,:,:)
      L(5,:,:,1)=(G(1,Nx,:,:)+c_NS(:,:,1))*(G(0,Nx,:,:)*c_NS(:,:,1)*dGx(1,Nx,:,:)+dGx(4,Nx,:,:))
      !設定したL行列からd1~5をi=0,Nxの両点においてそれぞれ設定する
        d(1,:,:,:) = (1.d0 / (c_NS(:,:,:) **2.d0)) * ((L(1,:,:,:)+L(5,:,:,:))*0.5d0 + L(2,:,:,:))
        d(2,:,:,:) = (L(1,:,:,:)+L(5,:,:,:))*0.5d0
        !d3のみrhoを含むので個別で設定しなければいけない
        d(3,:,:,0) = (1.d0/ (2.d0 * G(0,0,:,:) * c_NS(:,:,0))) * (-L(1,:,:,0) + L(5,:,:,0))
        d(3,:,:,1) = (1.d0/ (2.d0 * G(0,Nx,:,:) * c_NS(:,:,1))) * (-L(1,:,:,1) + L(5,:,:,1))
        d(4,:,:,:) = L(3,:,:,:)
        d(5,:,:,:) = L(4,:,:,:)
      !設定したdからNxSCBCで置き換える境界地点のdFxを定義する
      !i=0の時の差し替えdFx
      dFx(0,0,:,:) = d(1,:,:,0)
      dFx(1,0,:,:) = G(1,0,:,:)*d(1,:,:,0)+G(0,0,:,:)*d(3,:,:,0)
      dFx(2,0,:,:) = G(2,0,:,:)*d(1,:,:,0)+G(0,0,:,:)*d(4,:,:,0)
      dFx(3,0,:,:) = G(3,0,:,:)*d(1,:,:,0)+G(0,0,:,:)*d(5,:,:,0)
      dFx(4,0,:,:) =(0.5d0)*((G(1,0,:,:)**2.d0)+(G(2,0,:,:)**2.d0)+&
                    (G(3,0,:,:)**2.d0))*d(1,:,:,0)+d(2,:,:,0)/(gamma-1.d0)+&
      (G(0,0,:,:)*G(1,0,:,:)*d(3,:,:,0)) + (G(0,0,:,:)*G(2,0,:,:)*d(4,:,:,0))+&
                    (G(0,0,:,:)*G(3,0,:,:)*d(5,:,:,0))

      !i=Nxの時の差し替えF
      dFx(0,Nx,:,:) = d(1,:,:,1)
      dFx(1,Nx,:,:) = (G(1,Nx,:,:)*d(1,:,:,1)) + (G(0,Nx,:,:)*d(3,:,:,1))
      dFx(2,Nx,:,:) = (G(2,Nx,:,:)*d(1,:,:,1)) + (G(0,Nx,:,:)*d(4,:,:,1))
      dFx(3,Nx,:,:) = (G(3,Nx,:,:)*d(1,:,:,1)) + (G(0,Nx,:,:)*d(5,:,:,1))
      dFx(4,Nx,:,:) =(0.5d0)*((G(1,Nx,:,:)**2.d0)+(G(2,Nx,:,:)**2.d0)+&
                    (G(3,Nx,:,:)**2.d0))*d(1,:,:,1)+d(2,:,:,1)/(gamma-1.d0)+ &
      (G(0,Nx,:,:)*G(1,Nx,:,:)*d(3,:,:,1))+(G(0,Nx,:,:)*G(2,Nx,:,:)*d(4,:,:,1))+&
                    (G(0,Nx,:,:)*G(3,Nx,:,:)*d(5,:,:,1))
    !NSCBCの角処理(x方向,y方向で設定した境界値が重複するため半分ずつ加える)
                  dFx(:,0,0,:) = dFx(:,0,0,:) * 0.5d0
                  dFx(:,0,Ny,:) = dFx(:,0,Ny,:) * 0.5d0
                  dFx(:,Nx,0,:) = dFx(:,Nx,0,:) * 0.5d0
                  dFx(:,Nx,Ny,:) = dFx(:,Nx,Ny,:) * 0.5d0
    endsubroutine NSCBC_x
    !次にy方向のNSCBC　sunrouineを作成
    subroutine NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: G,dGy,dFy
      double precision,dimension(1:5,0:Nx,0:Nz-1,0:1) :: L,d
      double precision,dimension(0:Nx,0:Nz-1,0:1) :: c_NS,Ma_NS
      double precision pNy_infty,p0y_infty
      L=0.d0;d=0.d0;c_NS=0.d0;Ma_NS=0.d0
      !音速cはi=0,Nxの両点においてそれぞれ定義しなければならない
      c_NS(:,:,0) = sqrt(gamma * G(4,:,0,:) / G(0,:,0,:))
      c_NS(:,:,1) = sqrt(gamma * G(4,:,Ny,:) / G(0,:,Ny,:))
      !マッハ数Maはi=0,Nxで使うので別々に定義する
      Ma_NS(:,:,0) = G(2,:,0,:) / c_NS(:,:,0)!vを使う
      Ma_NS(:,:,1) = G(2,:,Ny,:) / c_NS(:,:,1)
      !y方向右側つまりi=0の点において無反射流出条件でL行列を設定する
      L(1,:,:,0) = (G(2,:,0,:) - c_NS(:,:,0)) * (-G(0,:,0,:)*c_NS(:,:,0)*dGy(2,:,0,:)+dGy(4,:,0,:))
      L(2,:,:,0) = G(2,:,0,:) * ((c_NS(:,:,0) ** 2.d0)*dGy(0,:,0,:) - dGy(4,:,0,:))
      L(3,:,:,0) = G(2,:,0,:) * dGy(1,:,0,:)
      L(4,:,:,0) = G(2,:,0,:) * dGy(3,:,0,:)
      L(5,:,:,0) = NS_sigma * c_NS(:,:,0) * (1.d0 - (Ma_NS(:,:,0) ** 2.d0))*(G(4,:,0,:) - &
      &p0y_infty) / (2.d0*Ly)
      !y方向左側つまりi=Nyの点において無反射流出条件でL行列を設定する
      L(1,:,:,1) = NS_sigma * c_NS(:,:,1) * (1.d0 - (Ma_NS(:,:,1) ** 2.d0))*(G(4,:,Ny,:) - &
    &  pNy_infty) / (2.d0*Ly)
      L(2,:,:,1) = G(2,:,Ny,:) * ((c_NS(:,:,1) ** 2.d0)*dGy(0,:,Ny,:) - dGy(4,:,Ny,:))
      L(3,:,:,1) = G(2,:,Ny,:) * dGy(1,:,Ny,:)
      L(4,:,:,1) = G(2,:,Ny,:) * dGy(3,:,Ny,:)
      L(5,:,:,1) = (G(2,:,Ny,:) + c_NS(:,:,1)) * (G(0,:,Ny,:)*c_NS(:,:,1)*dGy(2,:,Ny,:) + dGy(4,:,Ny,:))
      !設定したL行列からd1~5をi=0,Nyの両点においてそれぞれ設定する
        d(1,:,:,:) = (1.d0 / (c_NS(:,:,:) **2.d0)) * ((L(1,:,:,:)+L(5,:,:,:))*0.5d0 + L(2,:,:,:))
        d(2,:,:,:) = (0.5d0) * (L(1,:,:,:)+L(5,:,:,:))
        d(3,:,:,:) = L(3,:,:,:)
        !d4のみrhoを含むので個別で設定しなければいけない
        d(4,:,:,0) = (1.d0/ (2.d0 * G(0,:,0,:) * c_NS(:,:,0))) * (-L(1,:,:,0) + L(5,:,:,0))
        d(4,:,:,1) = (1.d0/ (2.d0 * G(0,:,Ny,:) * c_NS(:,:,1))) * (-L(1,:,:,1) + L(5,:,:,1))
        d(5,:,:,:) = L(4,:,:,:)
      !設定したdからNySCBCで置き換える境界地点のdFyを定義する
      !i=0の時の差し替えdFy
      dFy(0,:,0,:) = d(1,:,:,0)
      dFy(1,:,0,:) = G(1,:,0,:)*d(1,:,:,0)+G(0,:,0,:)*d(3,:,:,0)
      dFy(2,:,0,:) = G(2,:,0,:)*d(1,:,:,0)+G(0,:,0,:)*d(4,:,:,0)
      dFy(3,:,0,:) = G(3,:,0,:)*d(1,:,:,0)+G(0,:,0,:)*d(5,:,:,0)
      dFy(4,:,0,:) =(0.5d0)*((G(1,:,0,:)**2.d0)+(G(2,:,0,:)**2.d0)+(G(3,:,0,:)**2.d0))*d(1,:,:,0)+&
                    d(2,:,:,0)/(gamma-1.d0)+(G(0,:,0,:)*G(1,:,0,:)*d(3,:,:,0))+&
          (G(0,:,0,:)*G(2,:,0,:)*d(4,:,:,0))+(G(0,:,0,:)*G(3,:,0,:)*d(5,:,:,0))

      !i=Nyの時の差し替えdFy
      dFy(0,:,Ny,:) = d(1,:,:,1)
      dFy(1,:,Ny,:) = (G(1,:,Ny,:)*d(1,:,:,1)) + (G(0,:,Ny,:)*d(3,:,:,1))
      dFy(2,:,Ny,:) = (G(2,:,Ny,:)*d(1,:,:,1)) + (G(0,:,Ny,:)*d(4,:,:,1))
      dFy(3,:,Ny,:) = (G(3,:,Ny,:)*d(1,:,:,1)) + (G(0,:,Ny,:)*d(5,:,:,1))
      dFy(4,:,Ny,:) =(0.5d0)*((G(1,:,Ny,:)**2.d0)+(G(2,:,Ny,:)**2.d0)+(G(3,:,Ny,:)**2.d0))*d(1,:,:,1)+&
                  d(2,:,:,1)/(gamma-1.d0)+(G(0,:,Ny,:)*G(1,:,Ny,:)*d(3,:,:,1))+&
      (G(0,:,Ny,:)*G(2,:,Ny,:)*d(4,:,:,1))+(G(0,:,Ny,:)*G(3,:,Ny,:)*d(5,:,:,1))
      !NSCBCの角処理(x方向,y方向で設定した境界値が重複するため半分ずつ加える)
                    dFy(:,0,0,:) = dFy(:,0,0,:) * 0.5d0
                    dFy(:,0,Ny,:) = dFy(:,0,Ny,:) * 0.5d0
                    dFy(:,Nx,0,:) = dFy(:,Nx,0,:) * 0.5d0
                    dFy(:,Nx,Ny,:) = dFy(:,Nx,Ny,:) * 0.5d0
    endsubroutine NSCBC_y
    !x方向のi=0の流入部はdirichlet条件で固定。i=Nxの流出条件はNeumann条件を設定する。
    !なぜなら超音速のため流入部ではLが全て0になり、dFxは全て0になり、計算の意味そのものがなくなってしまうから。
    subroutine  Neumann(Q)
      integer i
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Q
        Q(:,Nx,:,:) = Q(:,Nx-1,:,:)
    endsubroutine Neumann
    subroutine inflow(Q,in_G)
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Q
      double precision,dimension(0:3,0:Ny,0:Nz-1) :: in_G
     Q(0,0,:,:) = in_G(0,:,:)!今までと違いrhoをNSCBCで求めずにdirichlet条件で固定してしまう
     Q(1,0,:,:) = in_G(0,:,:)*in_G(1,:,:)!rho*u
     Q(2,0,:,:) = in_G(0,:,:)*in_G(2,:,:)!rho*v
     Q(3,0,:,:) = in_G(0,:,:)*in_G(3,:,:)!rho*w
     Q(4,0,:,:) = 1.d0/((Ma**2.d0)*gamma*(gamma-1.d0))&
                 +in_G(0,:,:)*(in_G(1,:,:)**2.d0+in_G(2,:,:)**2.d0+in_G(3,:,:)**2.d0)*0.5d0!Et
    !まずは簡単な流入条件で試すために密度ρはNSCBCで求めたものを使うようにする
      ! Q(1,0,:,:) = Q(0,0,:,:)*in_G(1,:,:)!rho*u
      ! Q(2,0,:,:) = Q(0,0,:,:)*in_G(2,:,:)!rho*v
      ! Q(3,0,:,:) = Q(0,0,:,:)*in_G(3,:,:)!rho*w
      ! Q(4,0,:,:) = 1.d0/((Ma**2.d0)*gamma*(gamma-1.d0))&
      !             +Q(0,0,:,:)*(in_G(1,:,:)**2.d0+in_G(2,:,:)**2.d0+in_G(3,:,:)**2.d0)*0.5d0!Et
    endsubroutine inflow
    subroutine outflow(UVWT,dUVWTx,Vx,dVx,dUVWTy,Vy,dVy)
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Vx,Vy,dUVWTx,dUVWTy,UVWT,dVx,dVy
      !無反射流出条件の時の条件を設定するsubroutine
      !境界のdVx,dVyに条件を設定するのでdVの計算ができた後に境界値のみ上書きをする
      !τ11=Vx(1),τ12=Vx(2),τ13=Vx(3)
      !τ21=Vy(1),τ22=Vy(2),τ23=Vy(3)
      !dτ12/dx,dτ13/dx,dτ21/dy,dτ23/dy,dq/dx,dq/dyの全てが0(境界のみ)
      !x方向右側の条件設定
      dVx(2,Nx,:,:) = 0.d0!dτ12/dx
      dVx(3,Nx,:,:) = 0.d0!dτ13/dx
      dVx(4,Nx,:,:) = dVx(1,Nx,:,:)*UVWT(1,Nx,:,:)+dUVWTx(1,Nx,:,:)*Vx(1,Nx,:,:)&
                    +dUVWTx(2,Nx,:,:)*Vx(2,Nx,:,:)+dUVWTx(3,Nx,:,:)*Vx(3,Nx,:,:)
      !dτ11/dx*u+τ11*du/dx+dτ12/dx(=0)*v+τ12*dv/dx+dτ13/dx(=0)*w+τ13*dw/dx+dq/dx(=0)
      !=dτ11/dx*u+τ11*du/dx+τ12*dv/dx+τ13*dw/dx
      !y方向左側の条件設定
      dVy(1,:,0,:) = 0.d0!dτ21/dy
      dVy(3,:,0,:) = 0.d0!dτ23/dy
      dVy(4,:,0,:) = dVy(2,:,0,:)*UVWT(2,:,0,:)+dUVWTy(1,:,0,:)*Vy(1,:,0,:)&
                    +dUVWTy(2,:,0,:)*Vy(2,:,0,:)+dUVWTy(3,:,0,:)*Vy(3,:,0,:)
      !dτ21/dy(=0)*u+τ21*du/dy+dτ22/dy*v+τ22*dv/dy+dτ23/dy(=0)*w+τ23*dw/dy+dq/dy(=0)
      !τ21*du/dy+dτ22/dy*v+τ22*dv/dy+τ23*dw/dy
      !y方向右側の条件設定
      dVy(1,:,Ny,:) = 0.d0
      dVy(3,:,Ny,:) = 0.d0
      dVy(4,:,Ny,:) = dVy(2,:,Ny,:)*UVWT(2,:,Ny,:)+dUVWTy(1,:,Ny,:)*Vy(1,:,Ny,:)&
                    +dUVWTy(2,:,Ny,:)*Vy(2,:,Ny,:)+dUVWTy(3,:,Ny,:)*Vy(3,:,Ny,:)
      !dq/dy
    endsubroutine outflow
    !buffer領域の設定subroutine
    !x方向
    !ここでまず計算に必要なU(x)とσ(x)を定義している
    !係数等はタサイさんの資料を参考にしている
    subroutine buffer_x(c_infty,Ux,sigma_x,zeta_fx)
      integer i
      double precision c_infty,Xmax,Xmin,x1
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Ux,sigma_x
      double precision,dimension(0:Nx) :: zeta_fx
      double precision,parameter ::alpha_u=1.5d0,alpha_sigma=1.125d0,beta_r=0.01d0,beta_l=0.01d0
      Xmax = Lx;Xmin = 0.d0
      do i = 0,Nx
        x1 = zeta_fx(i)
        !Uxは本来0:Nxの一次の配列で十分だが座標変換の際に他の配列と計算する際に配列が揃っていないと:で省略して計算できないので
        !無駄に同じ値を入れて4次の配列にしている
        Ux(:,i,:,:) = alpha_u*c_infty*(dtanh(dble(atanh(beta_r/alpha_u-1.d0))*(x1-Xmax)/(-Wrx))&
        -dtanh(dble(atanh(beta_l/alpha_u-1.d0))*(x1-Xmin)/Wlx))
        if(x1<(Xmax-Wrx)) then
          sigma_x(:,i,:,:) = 0.d0!流出部にのみBufferをつけるのでx左側もσは0となる
        else
          sigma_x(:,i,:,:) = alpha_sigma*c_infty*((x1-(Xmax-Wrx))/Wrx)**3.d0
        endif
      enddo
       !dx不要かも？
     Ux(:,0:NUx,:,:) = 0.d0!x左側のBufferを取るためにWlxの範囲のUxを確実に0に設定している
!      open(100,file="ux-check.csv")
!      do i=0,Nx
!      write(100,*) zeta_fx(i),",",Ux(0,i,0,0)
!      enddo
!      close(100)
      !格子伸長が入っているためほぼ計算領域の真ん中になる値を調べて代入した
    endsubroutine buffer_x
    !y方向
    subroutine buffer_y(c_infty,Uy,sigma_y,zeta_fy)
      integer i
      double precision Ymax,Ymin,c_infty,y1
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Uy,sigma_y
      double precision,dimension(0:Ny) :: zeta_fy
      double precision,parameter ::alpha_u=1.15d0,alpha_sigma=1.125d0,beta_r=0.01d0,beta_l=0.01d0
      Ymax = Ly;Ymin = -Ly
      !格子伸長を行うので新しい座標ζ_yを用いてUyとsigma_yを設定する
      do i = 0,Ny
        y1 = zeta_fy(i)
        Uy(:,:,i,:) = alpha_u*c_infty*(dtanh(dble(atanh(beta_r/alpha_u-1.d0))*(y1-Ymax)/(-Wry))&
        -dtanh(dble(atanh(beta_l/alpha_u-1.d0))*(y1-Ymin)/Wly))
        if(y1<(Wly+Ymin)) then
          sigma_y(:,:,i,:) = alpha_sigma*c_infty*((-y1+Ymin+Wly)/Wly)**3.d0
        elseif((y1>=(Wly+Ymin)).and.(y1<(Ymax-Wry))) then
          sigma_y(:,:,i,:) = 0.d0
        elseif(y1>=(Ymax-Wry)) then
          sigma_y(:,:,i,:) = alpha_sigma*c_infty*((y1-(Ymax-Wry))/Wry)**3.d0
        endif
      enddo
    endsubroutine buffer_y
    !ζ,dζ/dxの定義subroutine
    subroutine lattice_x(dx,zeta_fx,dzeta_inx)
      integer i
      double precision,dimension(0:Nx) :: zeta_fx
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: dzeta,dzeta_inx
      double precision dx,x,width,a1,a2,b1
      dzeta=0.d0;width=10.8d0;a1=1d0/14d0;a2=7.d0;b1=1.d0/1.4d0
      !widthは格子間隔を細かくする範囲。この式では-width<=x<=widthの範囲で適用される
      !a2は粗い所と細かい所の境界の傾きの大きさを設定している
      !a1はどの程度の格子数の差をつけるかを設定する係数。このパラメータは条件によって適宜調整する
      do i= 0,Nx
        x = dx*dble(i)
       zeta_fx(i) = b1 * ((1.7d0*x) - a1 * &
        (-dlog(dcosh(a2*(x - width))) + dlog(dcosh(a2*(x + width)))))
        dzeta(:,i,:,:) = b1 * (1.7d0 - (a1*a2) * &
        (-dtanh(a2*(x - width)) + dtanh(a2*(x + width))))
      enddo
      dzeta_inx = 1.d0/dzeta
    endsubroutine lattice_x
    !ζ,dζ/dyの定義subroutine
    subroutine lattice_y(dy,zeta_fy,dzeta_iny)
      integer i
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: dzeta,dzeta_iny
      double precision,dimension(0:Ny) :: zeta_fy
      double precision dy,y,width,a1,a2,b1
      dzeta=0.d0;width=3.d0;a1=1d0/14d0;a2=7d0;b1=1.d0/1.4d0
      !widthは格子間隔を細かくする範囲。この式では-width<=y<=widthの範囲で適用される
      !a2は粗い所と細かい所の境界の傾きの大きさを設定している
      !a1はどの程度の格子数の差をつけるかを設定する係数
      do i= 0,Ny
        y = -Ly + dy*dble(i)
        zeta_fy(i) = b1 * ((1.7d0*y) - a1 * &
        (-dlog(dcosh(a2*(y - width))) + dlog(dcosh(a2*(y + width)))))
        dzeta(:,:,i,:) = b1 * (1.7d0 - (a1*a2) * &
        (-dtanh(a2*(y - width)) + dtanh(a2*(y + width))))
      enddo
      dzeta_iny = 1.d0/dzeta
    endsubroutine lattice_y
    !作成したdx/dζをdF/dyなどに掛けて微分変換を行うsubroutine
    subroutine combine(dzeta_in,dF,dFzeta)
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: dzeta_in,dFzeta,dF
      dFzeta = dF * dzeta_in
    endsubroutine combine
    subroutine Q_boundary(Q)
      integer i
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Q
        Q(:,0,:,:) = Q(:,1,:,:)
        Q(:,:,0,:) = Q(:,:,1,:)
        Q(:,Nx,:,:) = Q(:,Nx-1,:,:)
        Q(:,:,Ny,:) = Q(:,:,Ny-1,:)
    endsubroutine Q_boundary
end module threedim

    program main
      use threedim
      implicit none
      character(len = 16) filename
      !時間更新毎に出力ファイルを変更するためのファイル名設定
      !NSCBCでdrho,du,dpを計算する必要があるので計算しやすいようにG行列を作成しG(rho,u,p)
      !となるように設定した。すなわちdGは密度、速度、圧力のそれぞれの微分項を含む行列となる
      !x方向
      !Runge-Kutta法のためQ1,Q2用の値を設定
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: G,Q,Q0,Q1,Q2,Qn,Fpx,Fmx,xp,xm,oldG
      !y方向
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Fpy,Fmy,yp,ym
      !z方向
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Fpz,Fmz,zp,zm
      double precision,dimension(0:Nx,0:Ny,0:Nz-1) :: myu
      !A,L,U,x,y,bはF+とF-用に分けるため別々で定義
      !x方向
      double precision,allocatable,dimension(:,:) :: LUmx,LUpx
      !y方向
      double precision,allocatable,dimension(:,:) :: LUmy,LUpy
      !z方向
      double precision,allocatable,dimension(:,:) :: LUmz,LUpz
      !粘性項の計算に使う行列(x方向)
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Vx,dVx,UVWT,dUVWTx
      double precision,allocatable,dimension(:,:) :: LUccsx
      !y方向
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Vy,dVy,dUVWTy
      double precision,allocatable,dimension(:,:) :: LUccsy
      !z方向
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Vz,dVz,dUVWTz
      double precision,allocatable,dimension(:,:) :: LUccsz
      ! NSCBC用
      double precision,dimension(0:3,0:Ny,0:Nz-1) :: in_G
      !x方向
      double  precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: dGx,dFx
      double precision  dx,x,pNx_infty
      !y方向
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: dGy,dFy
      double precision dy,y,p0y_infty,pNy_infty
      !z方向
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: dFz,dGz
      double precision dz,z
      integer i,j,k,M,ii,jj,kk
      double precision t0,t1,theta
      double precision c_infty
      double precision,dimension(0:Ny) :: ur,Tu
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: Ux,sigma_x,Uy,sigma_y,dQx,dQy
      double precision,dimension(0:4,0:Nx,0:Ny,0:Nz-1) :: dzeta_iny,dzeta_inx
      double precision,dimension(0:Nx) :: zeta_fx
      double precision,dimension(0:Ny) :: zeta_fy
      double precision,dimension(0:Nx,0:Ny,0:Nz-1) :: omega_1,omega_2,omega_3,dp!渦度と圧力変動差を入れる配列
      !計算にかかる時間をCPU時間で計測する
      call cpu_time(t0)

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
      dy = 2.d0*Ly /dble(Ny)
      dz = Lz / dble(Nz)
!一応ゼロクリア
      G=0.d0;Q=0.d0;Qn=0.d0;Q0=0.d0;Q1=0.d0;Q2=0.d0
      pNx_infty=0.d0;p0y_infty=0.d0;pNy_infty=0.d0;in_G=0.d0;ur=0.d0;Tu=0.d0
      Ux=0.d0;sigma_x=0.d0;Uy=0.d0;sigma_y=0.d0;zeta_fy=0.d0;dzeta_iny=0.d0
      zeta_fx=0.d0;dzeta_inx=0.d0;omega_1=0.d0;omega_2=0.d0;omega_3=0.d0;dp=0.d0;oldG=0.d0

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
      G(0,:,:,:) = 1.d0!ρ
      G(1,:,:,:) = 0.d0!u
      G(2,:,:,:) = 0.d0!v
      G(3,:,:,:) = 0.d0!w
      G(4,:,:,:) = G(0,:,:,:)*Temp/((Ma**2.d0)*gamma)!p
      !Bufferの計算のための初期値を用いて無限遠方での音速を定義
      c_infty = sqrt(Temp/Ma**2.d0)
      !Buffer領域の計算に使うUx,Uy,sigma_x,sigma_yの計算
      call buffer_x(c_infty,Ux,sigma_x,zeta_fx)
      call buffer_y(c_infty,Uy,sigma_y,zeta_fy)
    !流入条件
    !x=0の軸上にのみ流入条件を適用することでここからどんどん流入が起こる
    do i =0,Nz-1
      in_G(0,:,i) = 1.d0/Tu(:)!密度ρは理想気体状態方程式に従うから
      in_G(1,:,i) = ur(:)! NSCBC流入条件で使う流入値のみを保存する配列
      in_G(2,:,i) = 0.d0! NSCBC流入条件で使う流入値のみを保存する配列
      in_G(3,:,i) = 0.d0! NSCBC流入条件で使う流入値のみを保存する配列
    enddo
!      !Doループ内でin_G(2,:,:)をDirichlet条件で設定して撹乱を導入しているのでここでは設定しない
      !初期値の出力
      !まずt=0はループ外で個別に作成
      open(10, file = "result_3D/parameter000000.d")
      !もちろん出力もζ_y座標系とζ_x座標系で行う
      !      do k=0,Nz-1
      !        z = dz*dble(k)
              z = dz*dble(Nz/2)
              do i = 0,Ny
                do j = 0,Nx
                  write(10,'(6f24.16)') zeta_fx(j),zeta_fy(i),z,G(0,j,i,Nz/2),omega_3(j,i,Nz/2),dp(j,i,Nz/2)/dt
                enddo
                write(10,*)
              enddo
      !        write(10,*)
      !      enddo
!      open(20,file = "result_3D/1pressure.d")
!      write(20,'(1I1,1f24.16)') 0,G(3,162,Ny/2,Nz/2)!(23,7)を指定しているが実際は(22.89,6.97)にずれてしまう
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
        !v方向から流入させる撹乱のためのsin波の流入条件を設定。時間変動させている
!        theta = 2.d0*pi*dble(M)*dt
!        in_G(2,:,:)=0.d0
!        !計算高速化
!        do i = 0,Ny
!          !if(zeta_fy(i) < -b) then
!          !  in_G(2,i,:) = 0.d0
!          if((zeta_fy(i) >= -b).and.(zeta_fy(i) < b)) then
!          !elseif((zeta_fy(i) >= -b).and.(zeta_fy(i) < b)) then
!            in_G(2,i,:) = A2*sin(T2*theta)!1秒で1周期になる。1/4*θなら4秒で1周期
!          !elseif(zeta_fy(i) >= b) then
!          !  in_G(2,i,:) = 0.d0
!          endif
!        enddo
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
                dFx = xm + xp
                dFy = ym + yp
                dFz = zm + zp
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
        call NSCBC_x(G,dGx,dFx,pNx_infty)
        !y方向
        call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
        call NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
        !無反射流出条件の際の境界での粘性項の条件を設定
        call outflow(UVWT,dUVWTx,Vx,dVx,dUVWTy,Vy,dVy)
        !Buffer領域の計算
        !計算に必要なdQ/dx,dQ/dyをCCSで導出。今まではdQ/dtしか求めていなかった
        !dQ/dx,dQ/dy自体はdG/dx,dG/dyを組み合わせて作ることができるのでそうして作成すると微分せずに済み計算の短縮に繋がる
        !しかしEtの微分だけは足し算が多くなってしまい遅くなるのでEtはそれで新規に微分、Qのそれ以外の値は組み合わせて作るのが最善
        !だがこの場合Etは単体で定義していないので今回はそれができない。そのためQを直接微分する
        call dif_x(ccs_sigma,dx,Q,dQx,LUccsx,dzeta_inx)
        call dif_y(ccs_sigma,dy,Q,dQy,LUccsy,dzeta_iny)
       !計算して求めたdF,dVそしてBuffer領域の計算のための値などを組み合わせ、代入してdQ/dtを求める
          Q1(:,:,:,:) = Q(:,:,:,:) + c*dt*(dVx(:,:,:,:)+dVy(:,:,:,:)+dVz(:,:,:,:)-dFx(:,:,:,:)&
              &-dFy(:,:,:,:)-dFz(:,:,:,:)-sigma_x(:,:,:,:)*(Q(:,:,:,:)-Q0(:,:,:,:))&
              -sigma_y(:,:,:,:)*(Q(:,:,:,:)-Q0(:,:,:,:))-Ux(:,:,:,:)*dQx(:,:,:,:)-Uy(:,:,:,:)*dQy(:,:,:,:))
      !call Q_boundary(Q1)
      !i=0で流入条件させるのでその部分のQ1を上書きして流入させ続ける
      call inflow(Q1,in_G)!dirichlet条件で流入部を固定
      !==========
      !NSCBC_xを上書きしてNeumannにしてしまう
      !===========
      !call Neumann(Q1)
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
              dFx = xm + xp
              dFy = ym + yp
              dFz = zm + zp
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
        call NSCBC_x(G,dGx,dFx,pNx_infty)
        !y方向
        call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
        call NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
        call outflow(UVWT,dUVWTx,Vx,dVx,dUVWTy,Vy,dVy)
        !Buffer領域の計算
        call dif_x(ccs_sigma,dx,Q1,dQx,LUccsx,dzeta_inx)
        call dif_y(ccs_sigma,dy,Q1,dQy,LUccsy,dzeta_iny)
        Q2(:,:,:,:) = (0.75d0)*Q(:,:,:,:) +(0.25d0) * Q1(:,:,:,:) + c* (dt*0.25d0)&
            &*(dVx(:,:,:,:)+dVy(:,:,:,:)+dVz(:,:,:,:)-dFx(:,:,:,:)-dFy(:,:,:,:)&
            -dFz(:,:,:,:)-sigma_x(:,:,:,:)*(Q1(:,:,:,:)-Q0(:,:,:,:))&
            -sigma_y(:,:,:,:)*(Q1(:,:,:,:)-Q0(:,:,:,:))-Ux(:,:,:,:)*dQx(:,:,:,:)-Uy(:,:,:,:)*dQy(:,:,:,:))
!        call Q_boundary(Q2)
        call inflow(Q2,in_G)
        !==========
        !NSCBC_xを上書きしてNeumannにしてしまう
        !===========
        !call Neumann(Q2)
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
              dFx = xm + xp
              dFy = ym + yp
              dFz = zm + zp
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
        !x方向のNSCBCの計算
        call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
        call NSCBC_x(G,dGx,dFx,pNx_infty)
        !y方向
        call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
        call NSCBC_y(G,dGy,dFy,pNy_infty,p0y_infty)
        call outflow(UVWT,dUVWTx,Vx,dVx,dUVWTy,Vy,dVy)
        !Buffer領域の計算
        call dif_x(ccs_sigma,dx,Q2,dQx,LUccsx,dzeta_inx)
        call dif_y(ccs_sigma,dy,Q2,dQy,LUccsy,dzeta_iny)
            Qn(:,:,:,:)=Q(:,:,:,:)/3.d0+(2.d0/3.d0)*Q2(:,:,:,:)+c*&
            &((2.d0*dt)/3.d0)*(dVx(:,:,:,:)+dVy(:,:,:,:)+dVz(:,:,:,:)-dFx(:,:,:,:)&
            -dFy(:,:,:,:)-dFz(:,:,:,:)-sigma_x(:,:,:,:)*(Q2(:,:,:,:)-Q0(:,:,:,:))&
            -sigma_y(:,:,:,:)*(Q2(:,:,:,:)-Q0(:,:,:,:))-Ux(:,:,:,:)*dQx(:,:,:,:)-Uy(:,:,:,:)*dQy(:,:,:,:))
!        call Q_boundary(Qn)
        call inflow(Qn,in_G)
        !==========
        !NSCBC_xを上書きしてNeumannにしてしまう
        !===========
        !call Neumann(Qn)
        call rho_u_p(G,Qn)
!        if (mod(M,p_output) == 0) then
!          write(20,'(1I7,1f24.16)') M,G(3,162,Ny/2,Nz/2)
!        endif
         if(mod(M,output_count) == 0) then!dt=1.d-4で0.01秒刻みで出力するためにMの条件を設定
         !渦度用のdGの計算
         dGx=0.d0;dGy=0.d0;dGz=0.d0
         call dif_x(ccs_sigma,dx,G,dGx,LUccsx,dzeta_inx)
         call dif_y(ccs_sigma,dy,G,dGy,LUccsy,dzeta_iny)
         call dif_z(ccs_sigma,dz,G,dGz,LUccsz)
           omega_1(:,:,:) = dGy(3,:,:,:) - dGz(2,:,:,:)
           omega_2(:,:,:) = dGz(1,:,:,:) - dGx(3,:,:,:)
           omega_3(:,:,:) = dGx(2,:,:,:) - dGy(1,:,:,:)
           call rho_u_p(oldG,Q)
           dp(:,:,:) = (G(4,:,:,:) - oldG(4,:,:,:))/dt
           write(filename, '(i6.6)') M
           !Mの計算毎に出力ファイル名を変更して出力する
           !i5.5で5桁分の数字を表示できるのでdt=1.d-5以下で計算するならここも変更が必要
           open(10, file = "result_3D/parameter"//trim(filename)//".d")
           !do kk= 0,Nz-1
             !z=dz*dble(kk)
             z=dz*dble(Nz/2)
             do ii = 0,Ny
               do jj = 0,Nx
                 write(10,'(6f24.16)') zeta_fx(jj),zeta_fy(ii),z,G(0,jj,ii,Nz/2),omega_3(jj,ii,Nz/2),dp(jj,ii,Nz/2)
               enddo
               write(10,*)
             enddo
             !write(10,*)
             !一度に全てを出力する際にはデータの切れ目として空白を一行挿入しなくてはいけない
           !enddo
           write(10,'(2A1,1I7)') "#","M",M
           write(10,'(7A10)')"#","x","y","z","rho","vorticity","dp/dt"
           close(10)
         endif
        !計算が破綻している場合に計算を終了させるプログラム
        do k = 0,Nz-1
            do j = 0,Ny
              do i = 0,Nx
                  if(isnan(Qn(0,i,j,k))) then
!                  if(isnan(Q1(0,i,j,k))) then
                    !渦度用のdGの計算
                    oldG=0;dGx=0.d0;dGy=0.d0;dGz=0.d0
                    call rho_u_p(oldG,Q)
                    call dif_x(ccs_sigma,dx,oldG,dGx,LUccsx,dzeta_inx)
                    call dif_y(ccs_sigma,dy,oldG,dGy,LUccsy,dzeta_iny)
                    call dif_z(ccs_sigma,dz,oldG,dGz,LUccsz)
                    omega_1(:,:,:) = dGy(3,:,:,:) - dGz(2,:,:,:)
                    omega_2(:,:,:) = dGz(1,:,:,:) - dGx(3,:,:,:)
                    omega_3(:,:,:) = dGx(2,:,:,:) - dGy(1,:,:,:)
                    !dp(:,:,:) = (G(4,:,:,:) - oldG(4,:,:,:))/dt Gは計算破綻時間には常にNaNなのでdpは計算不能
                    write(filename, '(i6.6)') M
                    !Mの計算毎に出力ファイル名を変更して出力する
                    !i5.5で5桁分の数字を表示できるのでdt=1.d-5以下で計算するならここも変更が必要
!                    open(10, file = "result_3D/parameter"//trim(filename)//".d")
                    do kk= 0,Nz-1
!                      z=dz*dble(Nz/2)
                      z=dz*dble(kk)
                      write(filename, '(i6.6)') kk
                      open(10, file = "result_3D/parameter"//trim(filename)//".d")
                      do ii = 0,Ny
                        do jj = 0,Nx
!                          write(10,'(6f24.16)') zeta_fx(jj),zeta_fy(ii),z,oldG(0,jj,ii,Nz/2),omega_3(jj,ii,Nz/2)
                          write(10,'(6f24.16)') zeta_fx(jj),zeta_fy(ii),z,oldG(0,jj,ii,kk),omega_3(jj,ii,kk)
                        enddo
                        write(10,*)
                        !一度に全てを出力する際にはデータの切れ目として空白を一行挿入しなくてはいけない
                      enddo
                      write(10,'(7A10)')"#","x","y","z","rho","vorticity"
                      close(10)
                    enddo
!                    write(10,'(2A1,1I7)') "#","M",M
!                    write(10,'(7A10)')"#","x","y","z","rho","vorticity"
!                    close(10)
                    write(*,*) "x=",i,"y=",j,"z=",k,"M=",M
                    call cpu_time(t1)
                    write(*,'("Time required = ",i3,"min",f4.1,"sec")') &
                    &int((((t1-t0) - mod(t1-t0,60.d0)) /60.d0)), mod(t1-t0,60.d0)
                    stop "rho becomes NAN"
                  endif
                enddo
              enddo
          enddo
!              Q = Q1!オイラー法の時間の更新
        !RK法の時間の更新
              Q = Qn
        if(mod(M,1) == 0) then
          write(*,*) "M=",M!計算に時間がかかるので進行状況の確認用に出力
        endif
      enddo
      call cpu_time(t1)
      write(*,'("Time required = ",i3,"min",f4.1,"sec")') &
      &int((((t1-t0) - mod(t1-t0,60.d0)) /60.d0)), mod(t1-t0,60.d0)
!      write(20,'(2A1,1f24.16,1A1,1f24.16)') "#","x",zeta_fx(Nx/2),"y",zeta_fy(Ny/2),"z",dz*dble(Nz/2)
!      write(20,'(3A15)')"#","M","Local Pressure"
!      close(20)
    end program main
