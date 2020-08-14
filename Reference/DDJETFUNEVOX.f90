! ***************  JET  ****************************************
! *********** TOSHI program for eigenfunction  *****************
! **************************************************************

      PROGRAM EIGEN
      IMPLICIT REAL*8 (A-H,M,O-Z)
      INTEGER LDA,ISB,NN,ISB2,IW,ID
      PARAMETER(N=128,LDA=5*N,NN=3,LDEVEC=5*N)
      CHARACTER NAME(6)*6
      COMMON /WORKSP/RWKSP

      REAL*8 M1,U1,U2,T1,T2,A,Q, &
             ALFA,BETA, &
             X,V,VX,PAI,RE,PR,GAMMA &
             OMEGAI,OMEGAR,CR,C,OMEGA2
      REAL*8 TEST(0:N-1),TES(0:LDA-1,0:LDA-1),RWKSP(412820), &
             TE(0:LDA-1,0:LDA-1),MAR(0:LDA-1,0:LDA-1), &
             MAI(0:LDA-1,0:LDA-1),IVAL(LDA),MIVAL(NN),MCRE(NN), &
             CRE(LDA),RVAL(LDA)
      REAL*8 EX(0:N-1),IVECA(N),RVECA(N),AMP(N), &
             ISP(3*N-3),IAMP(N-1),ISO(N-1)
      COMPLEX*16 AMAIN(0:5*N-1,0:5*N-1),CIJ(0:N-1,0:N-1), &
                 DIJ(0:N-1,0:N-1), &
                 U(0:N-1),DUDX(0:N-1),T(0:N-1),DTDX(0:N-1), &
                 ROU(0:N-1),DROUDX(0:N-1),CI,EVAL(LDA), &
                 P(0:N-1), &
                 MYU(0:N-1),DMYUDT(0:N-1),DDMDT2(0:N-1), &
                 DDUDX2(0:N-1),DDTDX2(0:N-1), &
                 EVEC(LDA,LDA)

      DATA NAME(1),NAME(2),NAME(3),NAME(4),NAME(5),NAME(6) &
          /'den','velo1','velo2','velo3' &
          ,'temp','press'/

!********* DATA FILE **********************************************
      OPEN(10,FILE='mM2.4Re1000a0.85th51b0am1funx.dat',STATUS='NEW')
      OPEN(20,FILE='mM2.4Re1000a0.85th51RHOmx.dat',STATUS='NEW')
      OPEN(30,FILE='mM2.4Re1000a0.85th51Umx.dat',STATUS='NEW')
      OPEN(40,FILE='mM2.4Re1000a0.85th51Pmx.dat',STATUS='NEW')

      OPEN(51,FILE='M2.4Re1000a0.85th51rx.dat',STATUS='NEW')
      OPEN(52,FILE='M2.4Re1000a0.85th51ux.dat',STATUS='NEW')
      OPEN(53,FILE='M2.4Re1000a0.85th51vx.dat',STATUS='NEW')
      OPEN(54,FILE='M2.4Re1000a0.85th51wx.dat',STATUS='NEW')
      OPEN(55,FILE='M2.4Re1000a0.85th51Tx.dat',STATUS='NEW')
      OPEN(56,FILE='M2.4Re1000a0.85th51px.dat',STATUS='NEW')

      OPEN(61,FILE='M2.4Re1000a0.85th51rxb.dat',STATUS='NEW')
      OPEN(62,FILE='M2.4Re1000a0.85th51uxb.dat',STATUS='NEW')
      OPEN(63,FILE='M2.4Re1000a0.85th51vxb.dat',STATUS='NEW')
      OPEN(64,FILE='M2.4Re1000a0.85th51wxb.dat',STATUS='NEW')
      OPEN(65,FILE='M2.4Re1000a0.85th51Txb.dat',STATUS='NEW')
      OPEN(66,FILE='M2.4Re1000a0.85th51pxb.dat',STATUS='NEW')
!********* INPUT INITIAL CONDITION ********************************
      CALL IWKIN(825640)
      PAI=4.0d0*dATAN(1.0d0)
      A=1.0d0

      GAMMA=1.4d0
      PR=1.0d0
!     --- alfa no hani  -----
!      ALFA=2.00d0
!      BETA=2.00d0
       alfa =    0.85d0
       beta =    alfa*dtan(51.0d0/180.0d0*PAI)
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
!c     --onsoku---(ex MachN=2 then AC=0.5)
!c      C=1.0d0
!c      M1=U2/C
!     -----------
      M1=2.4d0
       C=U1/M1
!     -----------

      RE=1000.0d0

!      Q=2.0d0

!      CI=(0.0d0,1.0d0)

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
  997   format (a8,F8.2)
! -------------------------------
      X=1.0d-16
      V=1.0d-16
      VX=1.0E10
!******************************************************************
!     ---- K=0-----
               U(0)= U2
            DUDX(0)= 0.0d0
          DDUDX2(0)= 0.0d0
               T(0)= T2
            DTDX(0)= 0.0d0
          DDTDX2(0)= 0.0d0
             ROU(0)= 1.0d0/T(0)
          DROUDX(0)=-T(0)**(-2.0)*DTDX(0)
             MYU(0)=             T(0)**( 2.0/3.0)
          DMYUDT(0)= 2.0d0/3.0d0*T(0)**(-1.0/3.0)
          DDMDT2(0)=-2.0d0/9.0d0*T(0)**(-4.0/3.0)
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

      open (21,file='test_all.dat')
      DO 13 K=1,N-1
            y=-A/dTAN(PAI*DBLE(K)/DBLE(N))
            write (21,998)  K,y,DBLE(U(K)),DBLE(DUDX(K)),DBLE(DDUDX2(K)) &
                            ,DBLE(T(K)),DBLE(DTDX(K)),DBLE(DDTDX2(K)) &
                            ,DBLE(ROU(K)),DBLE(DROUDX(K)) &
                           ,DBLE(MYU(K)),DBLE(DMYUDT(K)),DBLE(DDMDT2(K))
   13 CONTINUE
      close(21)

  998 format (I5,12E17.8)

      DO 850 K=1,N-1
            y=-A/dTAN(PAI*DBLE(K)/DBLE(N))
      WRITE(20,999) 1,K-N/2,Y,DBLE(ROU(K))
      WRITE(30,999) 1,K-N/2,Y,DBLE(U(K))
      WRITE(40,999) 1,K-N/2,Y,DBLE(P(K))

  850 CONTINUE

      CLOSE(20,STATUS='KEEP')
      CLOSE(30,STATUS='KEEP')
      CLOSE(40,STATUS='KEEP')

!************************D NO KEISAN****************************
      DO 20 I=0,N-1
      DO 30 J=0,N-1

      IF (I.EQ.J) THEN

      CIJ(I,J)=CMPLX((DBLE(1-N)*dSIN((2.0d0*PAI*DBLE(I)/DBLE(N))) &
                      /DBLE(N)/A))

      ELSE

      IF (IABS(I-J).EQ.N/2) THEN

      CIJ(I,J)=CMPLX((-dSIN((2.0d0*PAI*DBLE(J)/DBLE(N))) &
                       *dSIN((DBLE(N-1)*PAI*DBLE(I-J)/DBLE(N))) &
                       /DBLE(N)/dSIN((PAI*DBLE(I-J)/DBLE(N)))/A))

      ELSE

      IF (MOD(IABS(I-J),2).EQ.0) THEN

      CIJ(I,J)=CMPLX((((1.0d0-dCOS((2.0d0*PAI*DBLE(J)/DBLE(N)))) &
                       /2.0d0/dTAN((PAI*DBLE(I-J)/DBLE(N))) &
                       -dSIN((2.0d0*PAI*DBLE(J)/DBLE(N))) &
                       *dSIN((DBLE(N-1)*PAI*DBLE(I-J)/DBLE(N))) &
                       /DBLE(N)/dSIN((PAI*DBLE(I-J)/DBLE(N))))/A))

      ELSE

      CIJ(I,J)=CMPLX((((dCOS((2.0d0*PAI*DBLE(J)/DBLE(N)))-1.0d0) &
                       /2.0d0/dTAN((PAI*DBLE(I-J)/DBLE(N))) &
                       -dSIN((2.0d0*PAI*DBLE(J)/DBLE(N))) &
                       *dSIN((DBLE(N-1)*PAI*DBLE(I-J)/DBLE(N))) &
                       /DBLE(N)/dSIN((PAI*DBLE(I-J)/DBLE(N))))/A))

      ENDIF
      ENDIF
      ENDIF

   30 CONTINUE
   20 CONTINUE

!************************DD NO KEISAN***************************
      DO 25 I=0,N-1
      DO 35 J=0,N-1

      IF (I.EQ.J) THEN

      DIJ(I,J)=CMPLX(-1.0d0/(2.0d0*DBLE(N)*A*A) &
                *(dCOS(4.0d0*PAI*DBLE(J)/DBLE(N)) &
                -4.0d0*dCOS(2.0d0*PAI*DBLE(J)/DBLE(N))+3.0d0) &
                *1.0d0/3.0d0*DBLE(N/2-1)*DBLE(N/2)*DBLE(N-1) &
                -DBLE(N-1)/(DBLE(N)*A*A) &
                *(dCOS(4.0d0*PAI*DBLE(I)/DBLE(N)) &
                -dCOS(2.0d0*PAI*DBLE(I)/DBLE(N))))

      ELSE

      IF (IABS(I-J).EQ.N/2) THEN

      DIJ(I,J)=CMPLX(-1.0d0/(4.0d0*A*A) &
                *(dCOS(4.0d0*PAI*DBLE(J)/DBLE(N)) &
                -4.0d0*dCOS(2.0d0*PAI*DBLE(J)/DBLE(N))+3.0d0) &
                *DBLE(1-N/2) &
                -1.0d0/(DBLE(N)*A*A) &
                *(dCOS(4.0d0*PAI*DBLE(J)/DBLE(N)) &
                -dCOS(2.0d0*PAI*DBLE(J)/DBLE(N))) &
                *dSIN(DBLE(N-1)*PAI*DBLE(I-J)/DBLE(N)) &
                /dSIN(PAI*DBLE(I-J)/DBLE(N)))

      ELSE

      IF (MOD(IABS(I-J),2).EQ.0) THEN

      DIJ(I,J)=CMPLX(-1.0d0/(4.0d0*A*A) &
                *(dCOS(4.0d0*PAI*DBLE(J)/DBLE(N)) &
                -4.0d0*dCOS(2.0d0*PAI*DBLE(J)/DBLE(N))+3.0d0) &
                *(1.0d0/(dSIN(PAI*DBLE(I-J)/DBLE(N)) &
                *dSIN(PAI*DBLE(I-J)/DBLE(N)))-DBLE(N/2)) &
                -3.0d0/(4.0d0*A*A)/dTAN(PAI*DBLE(I-J)/DBLE(N)) &
                *(2.0d0*dSIN(2.0d0*PAI*DBLE(J)/DBLE(N)) &
                -dSIN(4.0d0*PAI*DBLE(J)/DBLE(N))) &
                -1.0d0/(DBLE(N)*A*A)*(dCOS(4.0d0*PAI*DBLE(J)/DBLE(N)) &
                -dCOS(2.0d0*PAI*DBLE(J)/DBLE(N))) &
                *dSIN(DBLE(N-1)*PAI*DBLE(I-J)/DBLE(N)) &
                /dSIN(PAI*DBLE(I-J)/DBLE(N)))

      ELSE

      DIJ(I,J)=CMPLX(1.0d0/(4.0d0*A*A) &
                *(dCOS(4.0d0*PAI*DBLE(J)/DBLE(N)) &
                -4.0d0*dCOS(2.0d0*PAI*DBLE(J)/DBLE(N))+3.0d0) &
                *(1.0d0/(dSIN(PAI*DBLE(I-J)/DBLE(N)) &
                *dSIN(PAI*DBLE(I-J)/DBLE(N)))-DBLE(N/2)) &
                +3.0d0/(4.0d0*A*A)/dTAN(PAI*DBLE(I-J)/DBLE(N)) &
                *(2.0d0*dSIN(2.0d0*PAI*DBLE(J)/DBLE(N)) &
                -dSIN(4.0d0*PAI*DBLE(J)/DBLE(N))) &
                -1.0d0/(DBLE(N)*A*A)*(dCOS(4.0d0*PAI*DBLE(J)/DBLE(N)) &
                -dCOS(2.0d0*PAI*DBLE(J)/DBLE(N))) &
                *dSIN(DBLE(N-1)*PAI*DBLE(I-J)/DBLE(N)) &
                /dSIN(PAI*DBLE(I-J)/DBLE(N)))

      ENDIF
      ENDIF
      ENDIF

   35 CONTINUE
   25 CONTINUE

!**************** GYOLETU NO KEISAN ************************************

      DO 100 K=0,N-1
      DO 200 L=0,N-1

      IF (K.EQ.L) THEN

      AMAIN(5*K,5*K)=ALFA*U(K)
      AMAIN(5*K,5*K+1)=ALFA*ROU(K)
      AMAIN(5*K,5*K+2)=-CI*(DROUDX(K)+ROU(K)*CIJ(K,K))
      AMAIN(5*K,5*K+3)=BETA*ROU(K)
      AMAIN(5*K,5*K+4)=(0.0d0,0.0d0)


      AMAIN(5*K+1,5*K)=ALFA*T(K)/(GAMMA*ROU(K)*M1*M1)
      AMAIN(5*K+1,5*K+1)=ALFA*U(K) &
                        -CI/(ROU(K)*RE) &
                        *(MYU(K)*(4.0d0/3.0d0*ALFA*ALFA &
                        -DIJ(K,K)+BETA*BETA) &
                        -DMYUDT(K)*DTDX(K)*CIJ(K,K))
      AMAIN(5*K+1,5*K+2)=-CI*DUDX(K) &
                        -ALFA/(ROU(K)*RE) &
                        *(1.0d0/3.0d0*MYU(K)*CIJ(K,K) &
                        +DMYUDT(K)*DTDX(K))
      AMAIN(5*K+1,5*K+3)= &
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
      AMAIN(5*K+2,5*K+1)= &
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
      AMAIN(5*K+2,5*K+3)= &
                        -BETA/(ROU(K)*RE) &
                        *(1.0d0/3.0d0*MYU(K)*CIJ(K,K) &
                        -2.0d0/3.0d0*DMYUDT(K)*DTDX(K))
      AMAIN(5*K+2,5*K+4)=-CI*(DROUDX(K)+ROU(K)*CIJ(K,K)) &
                        /(GAMMA*ROU(K)*M1*M1) &
                        -ALFA/(ROU(K)*RE)*DUDX(K)/T(K) &
                        *DMYUDT(K)

      AMAIN(5*K+3,5*K)=BETA*T(K)/(GAMMA*ROU(K)*M1*M1)
      AMAIN(5*K+3,5*K+1)= &
                        -CI*ALFA*BETA*MYU(K) &
                        /(3.0*ROU(K)*RE)
      AMAIN(5*K+3,5*K+2)= &
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

      AMAIN(5*K+4,5*K)=(0.0d0,0.0d0)
      AMAIN(5*K+4,5*K+1)=T(K)*(GAMMA-1.0d0)*ALFA &
                        +2.0d0*CI*GAMMA*(GAMMA-1.0d0)*M1*M1*MYU(K) &
                        /(ROU(K)*RE)*DUDX(K)*CIJ(K,K)
      AMAIN(5*K+4,5*K+2)=-CI*(DTDX(K)+(GAMMA-1)*T(K)*CIJ(K,K)) &
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
      AMAIN(5*K+2,5*L+2)=4.0d0*CI/(3.0d0*ROU(K)*RE) &
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
!      IF (ABS(TES(K,L)).LT.V) THEN
!      MAR(K,L)=0.0d0
!      ELSE
      MAR(K,L)=TES(K,L)
!      ENDIF

      TE(K,L)=DIMAG(AMAIN(K,L))
!      IF (ABS(TE(K,L)).LT.V) THEN
!      MAI(K,L)=0.0d0
!      ELSE
      MAI(K,L)=TE(K,L)
!      ENDIF

      AMAIN(K,L)=CMPLX((MAR(K,L)),(MAI(K,L)))
   60 CONTINUE
   50 CONTINUE
!**********************KOYUCHI NO KEISAN***************************

      CALL DEVCCG (LDA,AMAIN,LDA,EVAL,EVEC,LDEVEC)

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


!***************subsonic mode**************************************
      OMEGAI=0.0d0
      ISB=1

      DO 510 I=1,LDA
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

        write (*,*)  'ISB ',ISB ,IVAL(ISB )
        write (*,*)  'ISB2',ISB2,IVAL(ISB2)

      MIVAL(1)=OMEGAI
      MCRE(1)=CRE(ISB)

      WRITE(10,210)M1,ALFA,BETA
  210 FORMAT(/,1X,'M1=',F4.1,2X,'ALFA=',F5.2,2X,'BETA=',F5.2/)

      WRITE(10,213)MIVAL(1),MCRE(1)
  213 FORMAT(/,1X,'OMEGAi=',G10.3,2X,'Cr=',G10.3/)


!***** Eigen function *******************************************

      DO 1000 IW=1,2
         IF (IW.EQ.1) THEN
            ISB=ISB
            ID=0
         ELSE
            ISB=ISB2
            ID=10
         ENDIF

      DO 777 J=1,5

      DO 535 I=1,N

!---[subsonic mode]---

              RVECA(I)=DBLE(EVEC((I-1)*5+J,ISB))
              IVECA(I)=DIMAG(EVEC((I-1)*5+J,ISB))
              AMP(I)=DSQRT((RVECA(I))**2.0d0+(IVECA(I))**2.0d0)

  535 CONTINUE


      DO 536 I=2,N
        IAMP(I-1)=AMP(I)
  536 CONTINUE

      DO 537 I=1,N-1

              ISP(I)=IAMP(I)
              ISP(N-1+I)=RVECA(I+1)
              ISP(2*N-2+I)=IVECA(I+1)

!      write (*,*) J,I

!      if (j.EQ.4) then
!           ISO(I)=0.0d0
!         else
!           if (I.EQ.N/2) then
!              ISO(I)=0.0d0
!           else
!              ISO(I)=datan2(ISP(2*N-2+I),ISP(N-1+I))
!           endif
!      endif
  537 CONTINUE

!***** Y ********************************************************
      DO 230 K=1,N-1
              IF (K.EQ.N/2) THEN
                EX(K)=0.0d0
              ELSE
                EX(K)=-A/(dTAN(PAI*DBLE(K)/DBLE(N)))
              ENDIF
  230 CONTINUE

!***** Out put ***************************************************
      if (j.EQ.1) then
      DO 901 I=1,N-1
        WRITE(51+ID,999) ISB,I-N/2,EX(I),ISP(I),ISP(N-1+I),ISP(2*N-2+I)
  901 CONTINUE

      elseif (j.EQ.2) then
      DO 902 I=1,N-1
        WRITE(52+ID,999) ISB,I-N/2,EX(I),ISP(I),ISP(N-1+I),ISP(2*N-2+I)
  902 CONTINUE

      elseif (j.EQ.3) then
      DO 903 I=1,N-1
        WRITE(53+ID,999) ISB,I-N/2,EX(I),ISP(I),ISP(N-1+I),ISP(2*N-2+I)
  903 CONTINUE

      elseif (j.EQ.4) then
      DO 904 I=1,N-1
        WRITE(54+ID,999) ISB,I-N/2,EX(I),ISP(I),ISP(N-1+I),ISP(2*N-2+I)
  904 CONTINUE

      elseif (j.EQ.5) then
      DO 905 I=1,N-1
        WRITE(55+ID,999) ISB,I-N/2,EX(I),ISP(I),ISP(N-1+I),ISP(2*N-2+I)
  905 CONTINUE

      endif

  999 format (2I5,5E25.16)

  777 CONTINUE


!%%%%%%%%%%%%%%%%%%%%       PRESSURE       %%%%%%%%%%%%%%%%%%%%%

      DO 321 I=1,N
!---[subsonic mode]---

      RVECA(I)=DBLE( &
                EVEC((I-1)*5+1,ISB)*T(I-1) &
               +EVEC((I-1)*5+5,ISB)*ROU(I-1))
      IVECA(I)=DIMAG( &
                EVEC((I-1)*5+1,ISB)*T(I-1) &
               +EVEC((I-1)*5+5,ISB)*ROU(I-1))
      AMP(I)=DSQRT( ((RVECA(I))**2.0d0) + ((IVECA(I))**2.0d0) )

  321 CONTINUE

      DO 800 I=2,N
      IAMP(I-1)=AMP(I)
  800 CONTINUE

      DO 810 I=1,N-1

!---[subsonic mode]---

      ISP(I)=IAMP(I)
      ISP(N-1+I)=RVECA(I+1)
      ISP(2*N-2+I)=IVECA(I+1)

!           if (I.eq.N/2) then
!              ISO(I)=0.0d0
!           else
!              ISO(I)=datan2(ISP(2*N-2+I),ISP(N-1+I))
!           endif

  810 CONTINUE

      DO 831 I=1,N-1
        WRITE(56+ID,999) 1,I-N/2,EX(I),ISP(I),ISP(N-1+I),ISP(2*N-2+I)
  831 CONTINUE

 1000 CONTINUE

!-------------------------------------------------------------

      CLOSE(51,STATUS='KEEP')
      CLOSE(52,STATUS='KEEP')
      CLOSE(53,STATUS='KEEP')
      CLOSE(54,STATUS='KEEP')
      CLOSE(55,STATUS='KEEP')
      CLOSE(56,STATUS='KEEP')

      CLOSE(61,STATUS='KEEP')
      CLOSE(62,STATUS='KEEP')
      CLOSE(63,STATUS='KEEP')
      CLOSE(64,STATUS='KEEP')
      CLOSE(65,STATUS='KEEP')
      CLOSE(66,STATUS='KEEP')

      STOP
      END

