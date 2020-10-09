!
program iryuu6
implicit none
integer i,j,k
double precision dz,pi,dzinv,Usum,Lsum,sigma
double precision,parameter :: alpha3 = 0.25d0
double precision,parameter :: alpha5 = 1.d0/3.d0
double precision,parameter :: Dalpha = 1.d0  !5次精度のDCSとなるための係数設定
double precision,parameter :: psigma = -0.25d0
double precision,parameter :: msigma = 0.25d0
double precision,parameter :: ccs_sigma = 0.d0
double precision,parameter :: b = 1.d0!Jet半径は1で固定してしまう
double precision,parameter :: Cz = 2.d0*b !z軸の幅の設定
double precision,parameter :: Wrz = 0.5d0*b!Buffer領域y方向右側の幅
double precision,parameter :: Wlz = Wrz!Buffer領域y方向左側の幅
double precision,parameter :: Lz = 2.d0*Cz+Wrz+Wlz!z方向長さ.計算領域の中心を0にする
integer,parameter::Nz=100

double  precision,parameter :: ra = 14.d0/9.d0, rb = 1.d0/9.d0&
&,da = 4.d0 / 9.d0,db = 2.d0 / 9.d0 !5次精度のDCSとなるための係数設定

double precision,allocatable::L(:,:),U(:,:),A(:,:),LU(:,:),Fz(:)
double precision,allocatable::y(:),Q(:),z(:),RHS_z(:),x(:)
double precision,allocatable::D2(:),D4(:),D6(:),D8(:)

allocate(Fz(0:Nz),y(0:Nz),Q(0:Nz),z(0:Nz),RHS_z(0:Nz))
allocate(A(0:Nz,0:Nz),L(0:Nz,0:Nz),U(0:Nz,0:Nz))
allocate(LU(-2:2,0:Nz))
allocate(D2(2:Nz-2),D4(2:Nz-2),D6(2:Nz-2),D8(2:Nz-2))
allocate(x(0:Nz))

  sigma = psigma
  pi=dacos(-1d0)
  dz=Lz/dble(Nz)

  A=0d0;L=0d0;U=0d0;LU=0d0
  Q=0d0;Fz=0d0;y=0d0;x=0d0
  D2=0.d0;D4=0.d0;D6=0.d0;D8=0.d0
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
  do i = 2,Nz-2
    A(i,i-1) = (1.d0 - Dalpha * sigma) * alpha5
    A(i,i) = 1.d0
    A(i,i+1) = (1.d0 + Dalpha * sigma) * alpha5
  enddo
  !$omp end parallel do

  !Nz-1行目 3次精度DCS
  A(Nz-1,Nz-2) = alpha3 * (1.d0 - Dalpha * sigma)
  A(Nz-1,Nz-1) = 1.d0
  A(Nz-1,Nz) = alpha3 * (1.d0 + Dalpha * sigma)
  !Nz行目 片側DCS
  A(Nz,Nz-1) = 3.d0
  A(Nz,Nz) = 1.d0
  !LU分解
  !まずL,Uの初期値を設定
  !===========並列化しない=====================================================
  !L(i,0)を求める際のU(0,0)だが、doループ内で求める値なので、並列化するとNaNになってしまう
  !U(0,i)とL(i,0),L(i,i)でDoループを分けることも可能だが、並列化の利点よりもループを増やす欠点の方が
  !大きいと判断し、ここでは並列化は行わない。
  do i = 0,Nz
    U(0,i) = A(0,i)
    L(i,0) = A(i,0) / U(0,0)
    L(i,i) = 1.d0
  enddo
  !===========並列化しない=====================================================

  !========並列化不可能======================================================
  !Uは行ごとに、Lは列ごとに求めていく。
  !ただしUの２行目、 Lの２列目、Uの3行目、Lの３列目といった順番
  do i = 1,Nz !初期条件の結果を利用してUのi列、Lのi行の順に求めていく
    do j = i,Nz !Uのi行の列要素を求めていく。i<=jを考慮する
      Usum = 0.d0
      do k = 0, i-1
        Usum = Usum + L(i,k) * U(k,j)
      enddo
      U(i,j) = A(i,j) - Usum
    enddo
    !U(i,j)はi=0,Nzでj=i,Nzで回しながら行要素を求めるが
    !Lはj=0,Nzでi=j+1,Nzで回しながら列要素を求めなければならない
    !しかしそれではdoループをi=0,Nzとj=0,Nzで分けなければならない
    !LはL(j,i)として考えるとi=0,Nzの一つのdo文でU,L両方を定義できる
    !その結果Uの２行目、 Lの２列目、Uの3行目、Lの３列目といった順に求めるコードが書ける
    !そうでなければNzの数だけコードを書かなければならず非現実的になってしまう
    !それがCroutのアルゴリズム
    do j = i+1,Nz !Lのi列の行要素を求めていく。j<iを考慮
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
  LU(0,Nz) = U(Nz,Nz)
  LU(1,Nz) = 0.d0
  !$omp parallel do
  do i = 1,Nz
    LU(-1,i) = L(i,i-1)
  enddo
  !$omp end parallel do

  !$omp parallel do
  do j = 0,1
    do i = 0,Nz-1
      LU(j,i) = U(i,i+j)
    enddo
  enddo
  !$omp end parallel do


  do k=0,Nz
    z(k)=-(Lz/2.d0) + dz*dble(k)
    Fz(k)=sin(2.d0*pi*z(k))
  end do

  dzinv = 1.d0/dz

  !片側DCS,3次精度DCSも入れた非周期条件の際のbの設定
  !片側DCSの右辺設定
 RHS_z(0) = ((-17.d0/6.d0)*Fz(0)+1.5d0*(Fz(1)+&
                  Fz(2))-Fz(3)/6.d0)*dzinv
 RHS_z(Nz)=((1.d0/6.d0)*Fz(Nz-3)-1.5d0*(Fz(Nz-2)+&
                  Fz(Nz-1))+(17.d0/6.d0)*Fz(Nz))*dzinv

 !3次精度DCSの右辺設定
 RHS_z(1)=((1.5d0)*(-Fz(0)+Fz(2))*(0.5d0*dzinv))+sigma*&
              ((Fz(0)-2.d0*Fz(1)+Fz(2))*(0.5d0*dzinv))
 RHS_z(Nz-1)=((1.5d0)*(-Fz(Nz-2)+Fz(Nz))*(0.5d0*dzinv))&
    +sigma*((Fz(Nz-2)-2.d0*Fz(Nz-1)+Fz(Nz))*(0.5d0*dzinv))

    !$omp parallel do
     do k=2,Nz-2
       D2(k) = (-Fz(k-1)+Fz(k+1)) * (0.5d0*dzinv)
       D4(k) = (-Fz(k-2)+Fz(k+2)) * (0.25d0*dzinv)

       D6(k) = (Fz(k-1)+Fz(k+1)- 2.d0* Fz(k)) * dzinv
       D8(k) = (Fz(k-2)+Fz(k+2)- 2.d0* Fz(k)) * (0.25d0*dzinv)
     enddo
   !$omp end parallel do
   !==========RHSにはD2~D8が必要で、Doループが大きいので、分割して並列化する
   !$omp parallel do
   do k =2,Nz-2
      RHS_z(k)=ra*D2(k)+rb*D4(k)+sigma*(da*D6(k)+db*D8(k))
   enddo
   !$omp end parallel do
  !前進代入法、後退代入法の計算サブルーチン(z方向)
    !前進代入
           y(0) = RHS_z(0)!例外の境界値
    !=============前後のyを使うので、並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
        do k = 1,Nz!y(0)でk=0は上で定義したので残りの1〜Nzを定義する
          !これはfill-inのない通常部分
          y(k) = RHS_z(k) - LU(-1,k)*y(k-1)
        enddo
    !=============並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝====＝＝＝======
      !後退代入
      !fill-inの計算の範囲外なので別で計算
             x(Nz) = y(Nz) / LU(0,Nz)!例外の境界値
    !=============前後のxを使うので、並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
      !x(k)でk=Nzは定義したので残りのk=0~Nz-1を定義する
      do k = Nz-1, 0, -1!後退するので-1ずつ進む
        x(k) = (y(k) - LU(1,k)*x(k+1)) / LU(0,k)
      enddo
    !=============並列化不可能＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝====＝＝＝======
  open(10,file ='testz_NPB.csv')
  do k=0,Nz
  write(10,*) z(k),',',Fz(k),',',x(k)/(2.d0*pi)
  end do
  close(10)
end program iryuu6
