! ***************  JET  *****************************************
! *****  TOSHI program with viscidterms for 2D eigen value  *****
! ***************************************************************

      PROGRAM EIGEN
      IMPLICIT REAL*8 (A-H,M,O-Z)
      INTEGER LDA,NK,JJ,ISB,ISB2
      PARAMETER(N=128,LDA=5*N,NK=60)
      COMMON /WORKSP/RWKSP
      REAL*8 bpara,Dm
      REAL*8 M1,U1,U2,T1,T2,A,Q, &
             ALFA0,ALFA,DALFA,BETA, &
             X,V,VX,PAI,RE,PR,GAMMA &
             OMEGAI,OMEGA2,CR,C
      REAL*8 TEST(0:N-1),TES(0:LDA-1,0:LDA-1),RWKSP(412820), &
             TE(0:LDA-1,0:LDA-1),MAR(0:LDA-1,0:LDA-1), &
             MAI(0:LDA-1,0:LDA-1),IVAL(LDA),RVAL(LDA),CRE(LDA)

      COMPLEX*16 AMAIN(0:LDA-1,0:LDA-1),CIJ(0:N-1,0:N-1), &
                 DIJ(0:N-1,0:N-1), &
                 U(0:N-1),DUDX(0:N-1),T(0:N-1),DTDX(0:N-1), &
                 ROU(0:N-1),DROUDX(0:N-1),CI,EVAL(LDA), &
                 MYU(0:N-1),DMYUDT(0:N-1),DDMDT2(0:N-1), &
                 DDUDX2(0:N-1),DDTDX2(0:N-1)

!********* DATA FILE **********************************************
      OPEN(10,FILE='DM2.4Re1000y128a.dat',STATUS='NEW')
      OPEN(11,FILE='DM2.4Re1000y128b.dat',STATUS='NEW')
      OPEN(30,FILE='DM2.4Re1000y128_all.dat',STATUS='NEW')
!********* INPUT INITIAL CONDITION ********************************
      CALL IWKIN(825640)

      A=1.0d0

      GAMMA=1.4d0
      PR=1.0d0
!     --- alfa no hani  -----
      ALFA0=0.1d0
      ALFA=ALFA0
      DALFA=0.1d0
      BETA=0.0d0
!     ----------------
      T1=1.12d0
      T2=1.0d0
      U1=1.0d0
      U2=0.0d0
!     ----------------
      bpara=0.08d0
      Dm=1.0d0
!     ----------------
      CI=(0.0d0,1.0d0)
!     ----------------
      M1=2.4d0
       C=U1/M1
!     -----------

      RE=1000.0d0

      PAI=4.0d0*dATAN(1.0d0)

! -------------------------------
       write (*,997) '======================================='
       write (*,997) ' Mach:',M1
       write (*,997) '   Re:',Re
       write (*,997) 'Alpha:',ALFA
       write (*,997) ' Beta:',BETA
       write (*,997) '   T1:',T1
       write (*,997) '   T2:',T2
       write (*,997) '   U1:',U1
       write (*,997) '   U2:',U2
       write (*,997) '======================================='
  997 format (a8,F8.2)
! -------------------------------

      X=1.0d-16
      V=1.0d-16
      VX=1.0E10
!******************************************************************
!     ---- K=0-----
               U( 0 )= U2
            DUDX( 0 )= 0.0d0
          DDUDX2( 0 )= 0.0d0
               T( 0 )= T2
            DTDX( 0 )= 0.0d0
          DDTDX2( 0 )= 0.0d0
             ROU( 0 )= 1.0d0/T(0)
          DROUDX( 0 )=-T(0)**(-2.0)*DTDX(0)
             MYU( 0 )=             T(0)**( 2.0/3.0)
          DMYUDT( 0 )= 2.0d0/3.0d0*T(0)**(-1.0/3.0)
          DDMDT2( 0 )=-2.0d0/9.0d0*T(0)**(-4.0/3.0)
!     ---- K=N/2-----
               U(N/2)= U1
            DUDX(N/2)= 0.0d0
          DDUDX2(N/2)= 0.0d0
               T(N/2)= T1
            DTDX(N/2)= 0.0d0
          DDTDX2(N/2)= 0.0d0
             ROU(N/2)= 1.0d0/T(N/2)
          DROUDX(N/2)=-T(N/2)**(-2.0)*DTDX(N/2)
             MYU(N/2)=             T(N/2)**( 2.0/3.0)
          DMYUDT(N/2)= 2.0d0/3.0d0*T(N/2)**(-1.0/3.0)
          DDMDT2(N/2)=-2.0d0/9.0d0*T(N/2)**(-4.0/3.0)

      DO 10 K=N/2+1,N-1
                    y=-A/dTAN(PAI*DBLE(K)/DBLE(N))
               U(K)= U1*(1.0d0 - dTanh( (y/Dm - Dm/y)/(4.0d0*bpara) ) ) &
                                                                 /2.0d0
            DUDX(K)=-(1/Dm + Dm/y**2)*U1 &
                 *( 1.0d0/dCosh(  (y/Dm - Dm/y)/(4.0d0*bpara) ) )**2 &
                                                     /(8.0d0*bpara)
          DDUDX2(K)= Dm*U1 &
                 *( 1.0d0/dCosh( ( y/Dm - Dm/y)/(4.0d0*bpara) ) )**2 &
                                                     /(4.0d0*bpara*y**3) &
                 +( 1.0d0/Dm  + Dm/y**2)**2*U1 &
                 *( 1.0d0/dCosh( ( y/Dm - Dm/y)/(4.0d0*bpara) ) )**2 &
                 *Tanh((  y/Dm  - Dm/y)/(4.0d0*bpara))/(16.0d0*bpara**2)
               T(K)= M1*M1*(GAMMA-1.0d0)/2.0d0*(U(K)-(U(K)**2.0d0)) &
                  + (T1-T2)*U(K) &
                  + T2
            DTDX(K)= M1*M1*(GAMMA-1.0d0)/2.0d0 &
                 *(  DUDX(K)-2.0d0*   U(K)*DUDX(K)) &
                  + (T1-T2)*DUDX(K)
          DDTDX2(K)= M1*M1*(GAMMA-1.0d0)/2.0d0 &
                 *(DDUDX2(K)-2.0d0*DUDX(K)*DUDX(K)-2.0d0*U(K)*DDUDX2(K))
             ROU(K)= 1.0d0/T(K)
          DROUDX(K)=-T(K)**(-2.0)*DTDX(K)
             MYU(K)=             T(K)**( 2.0/3.0)
          DMYUDT(K)= 2.0d0/3.0d0*T(K)**(-1.0/3.0)
          DDMDT2(K)=-2.0d0/9.0d0*T(K)**(-4.0/3.0)
   10 CONTINUE

      DO 11 K=N/2+1,N-1
             KK=N-K
               U(KK)=       U(K)
            DUDX(KK)= -  DUDX(K)
          DDUDX2(KK)=  DDUDX2(K)
               T(KK)=       T(K)
            DTDX(KK)= -  DTDX(K)
          DDTDX2(KK)=  DDTDX2(K)
             ROU(KK)=     ROU(K)
          DROUDX(KK)= -DROUDX(K)
             MYU(KK)=     MYU(K)
          DMYUDT(KK)=  DMYUDT(K)
          DDMDT2(KK)=  DDMDT2(K)
   11 CONTINUE
!     -----------------------------------
      DO 12 K=0,N-1
            IF (DABS(DBLE(     U(K))).LT.1.0d-20)       U(K)=0.0d0
            IF (DABS(DBLE(  DUDX(K))).LT.1.0d-20)    DUDX(K)=0.0d0
            IF (DABS(DBLE(DDUDX2(K))).LT.1.0d-20)  DDUDX2(K)=0.0d0
            IF (DABS(DBLE(     T(K))).LT.1.0d-20)       T(K)=0.0d0
            IF (DABS(DBLE(  DTDX(K))).LT.1.0d-20)    DTDX(K)=0.0d0
            IF (DABS(DBLE(DDTDX2(K))).LT.1.0d-20)  DDTDX2(K)=0.0d0
            IF (DABS(DBLE(   ROU(K))).LT.1.0d-20)     ROU(K)=0.0d0
            IF (DABS(DBLE(DROUDX(K))).LT.1.0d-20)  DROUDX(K)=0.0d0
            IF (DABS(DBLE(   MYU(K))).LT.1.0d-20)     MYU(K)=0.0d0
            IF (DABS(DBLE(DMYUDT(K))).LT.1.0d-20)  DMYUDT(K)=0.0d0
            IF (DABS(DBLE(DDMDT2(K))).LT.1.0d-20)  DDMDT2(K)=0.0d0
   12 CONTINUE
!     -----------------------------------

      open (20,file='mean_M2.4.dat')
      DO 13 K=1,N-1
            y=-A/dTAN(PAI*DBLE(K)/DBLE(N))
            write (20,999)  K,y &
                 ,DBLE(  U(K)),DBLE(  DUDX(K)),DBLE(DDUDX2(K)) &
                 ,DBLE(  T(K)),DBLE(  DTDX(K)),DBLE(DDTDX2(K)) &
                 ,DBLE(ROU(K)),DBLE(DROUDX(K)) &
                 ,DBLE(MYU(K)),DBLE(DMYUDT(K)),DBLE(DDMDT2(K))
   13 CONTINUE
      close(20)

  999 format (I5,12E17.8)

!************************D NO KEISAN****************************
      DO 20 I=0,N-1
      DO 30 J=0,N-1

      IF (I.EQ.J) THEN

      CIJ(I,J)=DCMPLX((DBLE(1-N)*dSIN((2.0d0*PAI*I/N))/N/A))

      ELSE

      IF (IABS(I-J).EQ.N/2) THEN

      CIJ(I,J)=DCMPLX((-dSIN((2.0d0*PAI*J/N)) &
                       *dSIN((DBLE(N-1)*PAI*(I-J)/N)) &
                       /N/dSIN((PAI*(I-J)/N))/A))

      ELSE

      IF (MOD(IABS(I-J),2).EQ.0) THEN

      CIJ(I,J)=DCMPLX((((1.0d0-dCOS((2.0d0*PAI*J/N))) &
                       /2.0d0/dTAN((PAI*(I-J)/N)) &
                       -dSIN((2.0d0*PAI*J/N)) &
                       *dSIN((DBLE(N-1)*PAI*(I-J)/N)) &
                       /N/dSIN((PAI*(I-J)/N)))/A))

      ELSE

      CIJ(I,J)=DCMPLX((((dCOS((2.0d0*PAI*J/N))-1.0d0) &
                       /2.0d0/dTAN((PAI*(I-J)/N)) &
                       -dSIN((2.0d0*PAI*J/N)) &
                       *dSIN((DBLE(N-1)*PAI*(I-J)/N)) &
                       /N/dSIN((PAI*(I-J)/N)))/A))

      ENDIF
      ENDIF
      ENDIF

   30 CONTINUE
   20 CONTINUE

!************************DD NO KEISAN***************************
      DO 25 I=0,N-1
      DO 35 J=0,N-1

      IF (I.EQ.J) THEN

      DIJ(I,J)=DCMPLX(-1.0d0/(2.0d0*N*A*A) &
                *(dCOS(4.0d0*PAI*J/N)-4.0d0*dCOS(2.0d0*PAI*J/N)+3.0d0) &
                *1.0d0/3.0d0*(N/2-1)*(N/2)*(N-1) &
                -DBLE(N-1)/(N*A*A) &
                *(dCOS(4.0d0*PAI*I/N)-dCOS(2.0d0*PAI*I/N)))

      ELSE

      IF (IABS(I-J).EQ.N/2) THEN

      DIJ(I,J)=DCMPLX(-1.0d0/(4.0d0*A*A) &
                *(dCOS(4.0d0*PAI*J/N)-4.0d0*dCOS(2.0d0*PAI*J/N)+3.0d0) &
                *DBLE(1-N/2) &
                -1.0d0/(N*A*A) &
                *(dCOS(4.0d0*PAI*J/N)-dCOS(2.0d0*PAI*J/N)) &
                *dSIN((N-1)*PAI*(I-J)/N) &
                /dSIN(PAI*(I-J)/N))

      ELSE

      IF (MOD(IABS(I-J),2).EQ.0) THEN

      DIJ(I,J)=DCMPLX(-1.0d0/(4.0d0*A*A) &
                *(dCOS(4.0d0*PAI*J/N)-4.0d0*dCOS(2.0d0*PAI*J/N)+3.0d0) &
                *(1.0d0/(dSIN(PAI*(I-J)/N)*dSIN(PAI*(I-J)/N))-DBLE(N/2)) &
                -3.0d0/(4.0d0*A*A)/dTAN(PAI*(I-J)/N) &
                *(2.0d0*dSIN(2.0d0*PAI*J/N)-dSIN(4.0d0*PAI*J/N)) &
                -1.0d0/(N*A*A)*(dCOS(4.0d0*PAI*J/N)-dCOS(2.0d0*PAI*J/N)) &
                *dSIN(DBLE(N-1)*PAI*(I-J)/N) &
                /dSIN(PAI*(I-J)/N))

      ELSE

      DIJ(I,J)=DCMPLX(1.0d0/(4.0d0*A*A) &
                *(dCOS(4.0d0*PAI*J/N)-4.0d0*dCOS(2.0d0*PAI*J/N)+3.0d0) &
                *(1.0d0/(dSIN(PAI*(I-J)/N)*dSIN(PAI*(I-J)/N))-DBLE(N/2)) &
                +3.0d0/(4.0d0*A*A)/dTAN(PAI*(I-J)/N) &
                *(2.0d0*dSIN(2.0d0*PAI*J/N)-dSIN(4.0d0*PAI*J/N)) &
                -1.0d0/(DBLE(N)*A*A) &
                *(dCOS(4.0d0*PAI*J/N)-dCOS(2.0d0*PAI*J/N)) &
                *dSIN(DBLE(N-1)*PAI*(I-J)/N) &
                /dSIN(PAI*(I-J)/N))

      ENDIF
      ENDIF
      ENDIF

   35 CONTINUE
   25 CONTINUE

!**************** GYOLETU NO KEISAN ************************************
      DO 600 J=1,NK
      ALFA=DALFA*(J-1)+ALFA0
              write (*,*) J,ALFA

      DO 100 K=0,N-1
      DO 200 L=0,N-1

      IF (K.EQ.L) THEN

      AMAIN(5*K,5*K)=U(K)*ALFA
      AMAIN(5*K,5*K+1)=ROU(K)*ALFA
      AMAIN(5*K,5*K+2)=-CI*(DROUDX(K)+ROU(K)*CIJ(K,K))
      AMAIN(5*K,5*K+3)=ROU(K)*BETA
      AMAIN(5*K,5*K+4)=DCMPLX(0.0d0)

      AMAIN(5*K+1,5*K)=T(K)*ALFA/(ROU(K)*GAMMA*M1*M1)
      AMAIN(5*K+1,5*K+1)=U(K)*ALFA &
                        -CI/(ROU(K)*RE) &
                        *(MYU(K)*(4.0d0/3.0d0*ALFA*ALFA &
                        -DIJ(K,K)+BETA*BETA) &
                        -DMYUDT(K)*DTDX(K)*CIJ(K,K))
      AMAIN(5*K+1,5*K+2)=-CI*DUDX(K) &
                        -ALFA/(ROU(K)*RE) &
                        *(1.0d0/3.0d0*MYU(K)*CIJ(K,K) &
                        +DMYUDT(K)*DTDX(K))
      AMAIN(5*K+1,5*K+3)=DCMPLX(0.0d0) &
                        -CI*MYU(K)*ALFA*BETA &
                        /(3.0d0*ROU(K)*RE)
      AMAIN(5*K+1,5*K+4)=ALFA/(GAMMA*M1*M1) &
                        +CI/(ROU(K)*RE) &
                        *(DUDX(K)*DDMDT2(K)*DTDX(K)/T(K) &
                        +DMYUDT(K)*DDUDX2(K)/T(K) &
                        +DUDX(K)/T(K)*DMYUDT(K)*CIJ(K,K) &
                        -DUDX(K)*DMYUDT(K)*DTDX(K)/(T(K)*T(K)))

      AMAIN(5*K+2,5*K)=-CI*(DTDX(K)+T(K)*CIJ(K,K)) &
                        /(GAMMA*ROU(K)*M1*M1)
      AMAIN(5*K+2,5*K+1)=DCMPLX(0.0d0) &
                        -ALFA/(ROU(K)*RE) &
                        *(1.0d0/3.0d0*MYU(K)*CIJ(K,K) &
                        -2.0d0/3.0d0*DMYUDT(K)*DTDX(K))
      AMAIN(5*K+2,5*K+2)=ALFA*U(K) &
                        -CI/(ROU(K)*RE) &
                        *(MYU(K)*(ALFA*ALFA &
                        -4.0d0/3.0d0*DIJ(K,K) &
                        +BETA*BETA) &
                        -4.0d0/3.0d0*DMYUDT(K)*DTDX(K) &
                        *CIJ(K,K))
      AMAIN(5*K+2,5*K+3)=DCMPLX(0.0d0) &
                        -BETA/(ROU(K)*RE) &
                        *(1.0d0/3.0d0*MYU(K)*CIJ(K,K) &
                        -2.0d0/3.0d0*DMYUDT(K)*DTDX(K))
      AMAIN(5*K+2,5*K+4)=-CI*(DROUDX(K)+ROU(K)*CIJ(K,K)) &
                        /(GAMMA*ROU(K)*M1*M1) &
                        -ALFA/(ROU(K)*RE)*DUDX(K)/T(K) &
                        *DMYUDT(K)

      AMAIN(5*K+3,5*K)=BETA*T(K)/(GAMMA*ROU(K)*M1*M1)
      AMAIN(5*K+3,5*K+1)=DCMPLX(0.0d0) &
                        -CI*ALFA*BETA*MYU(K) &
                        /(3.0d0*ROU(K)*RE)
      AMAIN(5*K+3,5*K+2)=DCMPLX(0.0d0) &
                        -BETA/(ROU(K)*RE) &
                        *(1.0d0/3.0d0*MYU(K)*CIJ(K,K) &
                        +DMYUDT(K)*DTDX(K))
      AMAIN(5*K+3,5*K+3)=ALFA*U(K) &
                        -CI/(ROU(K)*RE) &
                        *(MYU(K)*(ALFA*ALFA &
                        -DIJ(K,K) &
                        +4.0d0/3.0d0*BETA*BETA) &
                        -DMYUDT(K)*DTDX(K)*CIJ(K,K))
      AMAIN(5*K+3,5*K+4)=BETA/(GAMMA*M1*M1)

      AMAIN(5*K+4,5*K)=DCMPLX(0.0d0)
      AMAIN(5*K+4,5*K+1)=T(K)*(GAMMA-1.0d0)*ALFA &
                        +2.0d0*CI*GAMMA*(GAMMA-1.0d0)*M1*M1*MYU(K) &
                        /(ROU(K)*RE)*DUDX(K)*CIJ(K,K)
      AMAIN(5*K+4,5*K+2)=-CI*(DTDX(K)+(GAMMA-1.0d0)*T(K)*CIJ(K,K)) &
                        -2.0d0*GAMMA*(GAMMA-1.0d0)*M1*M1*MYU(K) &
                        /(ROU(K)*RE)*ALFA*DUDX(K)
      AMAIN(5*K+4,5*K+3)=T(K)*(GAMMA-1.0d0)*BETA
      AMAIN(5*K+4,5*K+4)=ALFA*U(K) &
                        +CI*GAMMA*(GAMMA-1.0d0)*M1*M1 &
                        /(ROU(K)*RE)*DMYUDT(K) &
                        *DUDX(K)*DUDX(K)/T(K) &
                        -CI*GAMMA/(ROU(K)*RE*PR) &
                        *(MYU(K)*(ALFA*ALFA &
                        -DIJ(K,K)+BETA*BETA) &
                        -DMYUDT(K)*(1.0d0+1.0d0/T(K))*DTDX(K)*CIJ(K,K) &
                        -DDMDT2(K)*DTDX(K)*DTDX(K)/T(K) &
                        -DMYUDT(K)*DDTDX2(K)/T(K) &
                        +DMYUDT(K)*DTDX(K)*DTDX(K)/(T(K)*T(K)))

      ELSE

      AMAIN(5*K,5*L+2)=-CI*ROU(K)*CIJ(K,L)
      AMAIN(5*K+1,5*L+1)=CI/(ROU(K)*RE)*(MYU(K)*DIJ(K,L) &
                        +DMYUDT(K)*DTDX(K)*CIJ(K,L))
      AMAIN(5*K+1,5*L+2)=-1.0d0/(3.0d0*ROU(K)*RE) &
                        *MYU(K)*ALFA*CIJ(K,L)
      AMAIN(5*K+1,5*L+4)=CI/(ROU(K)*RE)*DUDX(K)/T(K) &
                        *DMYUDT(K)*CIJ(K,L)
      AMAIN(5*K+2,5*L)=-CI*T(K)*CIJ(K,L)/(ROU(K)*GAMMA*M1*M1)
      AMAIN(5*K+2,5*L+1)=-1.0d0/(3.0d0*ROU(K)*RE)*ALFA &
                        *MYU(K)*CIJ(K,L)
      AMAIN(5*K+2,5*L+2)=4.0*CI/(3.0d0*ROU(K)*RE) &
                        *(MYU(K)*DIJ(K,L) &
                        +DMYUDT(K)*DTDX(K)*CIJ(K,L))
      AMAIN(5*K+2,5*L+3)=-1.0d0/(3.0d0*ROU(K)*RE) &
                        *BETA*MYU(K)*CIJ(K,L)
      AMAIN(5*K+2,5*L+4)=-CI*CIJ(K,L)/(GAMMA*M1*M1)
      AMAIN(5*K+3,5*L+2)=-1.0d0/(3.0d0*ROU(K)*RE) &
                        *BETA*MYU(K)*CIJ(K,L)
      AMAIN(5*K+3,5*L+3)=CI/(ROU(K)*RE)*(MYU(K)*DIJ(K,L) &
                        +DMYUDT(K)*DTDX(K)*CIJ(K,L))
      AMAIN(5*K+4,5*L+1)=2.0d0*CI/(ROU(K)*RE)*GAMMA*(GAMMA-1.0d0) &
                        *M1*M1*MYU(K)*DUDX(K)*CIJ(K,L)
      AMAIN(5*K+4,5*L+2)=-CI*(GAMMA-1.0d0)*T(K)*CIJ(K,L)
      AMAIN(5*K+4,5*L+4)=CI*GAMMA/(ROU(K)*RE*PR) &
                        *(MYU(K)*DIJ(K,L) &
                        +DMYUDT(K)*(1.0d0+1.0d0/T(K))*DTDX(K)*CIJ(K,L))

      ENDIF
  200 CONTINUE
  100 CONTINUE
!***************GYOLETU NO HOSEI*****************************
      DO 50 K=0,LDA-1
      DO 60 L=0,LDA-1

      TES(K,L)=DBLE(AMAIN(K,L))
      IF (ABS(TES(K,L)).LT.V) THEN
      MAR(K,L)=0.0d0
      ELSE
      MAR(K,L)=TES(K,L)
      ENDIF

      TE(K,L)=DIMAG(AMAIN(K,L))
      IF (ABS(TE(K,L)).LT.V) THEN
      MAI(K,L)=0.0d0
      ELSE
      MAI(K,L)=TE(K,L)
      ENDIF

      AMAIN(K,L)=DCMPLX((MAR(K,L)),(MAI(K,L)))
   60 CONTINUE
   50 CONTINUE
!**********************KOYUCHI NO KEISAN********************************
      CALL DEVLCG (LDA,AMAIN,LDA,EVAL)
!****************W(I) NO CHUSHUTU**********************************
      DO 500 I=1,LDA
      IVAL(I)=DIMAG(EVAL(I))
      RVAL(I)=DBLE(EVAL(I))
      IF(RVAL(I).NE.0.0d0) THEN
      CRE(I)=(RVAL(I))/ALFA
      ELSE
      CRE(I)=1.E3
      ENDIF
  500 CONTINUE

!************subsonic mode******************************************

      OMEGAI=0.0d0
      CR=1.E3

      DO 510 I=1,LDA
         WRITE(30,210) I,ALFA,CRE(I),IVAL(I),RVAL(I)
         IF (CRE(I).GT.U2.AND.CRE(I).LT.U1) THEN
            IF (IVAL(I).GT.OMEGAI) THEN
               OMEGAI=IVAL(I)
               ISB=I
            ENDIF
         ENDIF
  510 CONTINUE

      OMEGA2=0.0d0
      DO 511 I=1,LDA
         IF (CRE(I).GT.U2.AND.CRE(I).LT.U1) THEN
            IF ((IVAL(I).GT.OMEGA2).and.(I.NE.ISB)) THEN
               OMEGA2=IVAL(I)
               ISB2=I
            ENDIF
         ENDIF
  511 CONTINUE

      WRITE(10,210)ISB ,ALFA,CRE(ISB ),IVAL(ISB ),RVAL(ISB )
      WRITE(11,210)ISB2,ALFA,CRE(ISB2),IVAL(ISB2),RVAL(ISB2)

  210 FORMAT(I5,4E25.16)

  600 CONTINUE

      CLOSE(10,STATUS='KEEP')
      CLOSE(11,STATUS='KEEP')
      CLOSE(30,STATUS='KEEP')

      STOP
      END
