! *************** WAKE *****************************************
! *********** TOSHI program for eigenfunction (WAKE) ***********
! **************************************************************

!!!-----(f77...COMMON) parameter definition--------------------------------------------------------------!
module Para_
      REAL(8):: RE,Re_d,s_para,Ly
      integer, PARAMETER:: N=450
      integer, PARAMETER:: By=10000
end module
subroutine parapara()
	use Para_
      Re_d=1000.0d0
      Ly=50.0d0
      s_para=4.0d0
end subroutine parapara
!!!------------------------------------------------------------------------------------------------------!
      PROGRAM EIGEN
use Para_
!      IMPLICIT REAL*8 (A-H,M,O-Z)
      IMPLICIT NONE
!      INTEGER :: N,LDA,ISB,NN,LDEVEC,ISB2,IW,ID
      INTEGER :: ISB,ISB2,IW,ID
!      integer, PARAMETER:: N=256,LDA=5*(N+1),NN=3,LDEVEC=5*(N+1)
!      integer, PARAMETER:: N=64
      integer, PARAMETER:: LDA=5*(N+1),NN=3,LDEVEC=5*(N+1)
      INTEGER :: I,J,K,L,M,KK
      CHARACTER :: NAME(6)*6
!      COMMON /WORKSP/RWKSP
      REAL(8):: y(0:N)
      REAL(8):: M1,U1,U2,T1,T2,A,Q, &
             & ALFA,BETA, &
             & X,V,VX,PI,PR,GAMMA, &
             & OMEGAI,OMEGAR,CR,C,OMEGA2
      REAL(8):: TEST(0:N),TES(0:LDA-1,0:LDA-1),RWKSP(412820), &
             & TE(0:LDA-1,0:LDA-1),MAR(0:LDA-1,0:LDA-1), &
             & MAI(0:LDA-1,0:LDA-1),IVAL(0:LDA-1),MIVAL(NN),MCRE(NN), &
             & CRE(0:LDA-1),RVAL(0:LDA-1)
      REAL(8):: EX(0:N),IVECA(0:N),RVECA(0:N),AMP(0:N), &
             & ISP(0:N,3),IAMP(0:N),ISO(0:N)
!!!外部サブルーチン
		integer :: INFO
		integer,parameter ::LDVL=N, LDVR=N
 		real(8) ::RWORK(8*N)
 		integer ::LWORK=100*N
		complex(kind(0d0))::WORK(3*N), VL(LDVL, N), VR(LDVR,N)
!!!!!


      complex(kind(0d0)):: EALPHA(0:LDA-1)
      complex(kind(0d0)):: EBETA(0:LDA-1)
      complex(kind(0d0)):: BMAIN(0:LDA-1,0:LDA-1)

      complex(kind(0d0)):: AMAIN(0:5*(N+1)-1,0:5*(N+1)-1), &
                & EIJ(0:N,0:N),CIJ(0:N,0:N),DIJ(0:N,0:N), &
                     & U(0:N),  DUDX(0:N),T(0:N),DTDX(0:N), &
                   & ROU(0:N),DROUDX(0:N),CI, &
                     & P(0:N), &
                   & MYU(0:N),DMYUDT(0:N),DDMDT2(0:N), &
                & DDUDX2(0:N),DDTDX2(0:N), &
                & EVAL(0:LDA-1), &
                & EVEC(0:LDA-1,0:LDA-1)
      REAL(8):: se,s,ss,coef(0:N)
!*************************
      REAL(8):: U_min,U_max
!      REAL(8):: RE,s_para,Ly
!      COMMON /mUmU/RE,s_para,Ly
!*************************
      INTEGER :: FR1,FR2,FR3,FR4,FR5,FR6
      INTEGER :: FM1,FM2,FM3
      INTEGER :: FA1,FA2,FA3
      INTEGER :: FB1,FB2,FB3
      CHARACTER(50):: FNAME1,FNAME2,FNAME3,FNAME4
      CHARACTER(50):: NAMEA1,NAMEA2,NAMEA3,NAMEA4,NAMEA5,NAMEA6
      CHARACTER(50):: NAMEB1,NAMEB2,NAMEB3,NAMEB4,NAMEB5,NAMEB6
      DATA NAME(1),NAME(2),NAME(3),NAME(4),NAME(5),NAME(6) &
          /'den','velo1','velo2','velo3' &
          ,'temp','press'/

call parapara()


      PI=4.0d0*dATAN(1.0d0)
      A=1.0d0

      GAMMA=1.4d0
      PR=1.0d0
 !     Ly=50.0d0
!*************************
!*************************
 !     s_para=4.0d0
       U_min=0.1d0
       U_max=0.7d0
!*************************
!*************************
!     --- alfa no hani  -----
!      ALFA=1.0d0
 !     BETA=0.0d0
!     ----------------
      T1=1.0d0
      T2=1.0d0
      U1=0.0d0
      U2=1.0d0
!     ----------------
      CI=(0.0d0,1.0d0)
!     ----------------
      M1=0.01d0
       C=U1/M1
!     -----------

!      RE=5772.0d0
!      RE=10000.0d0


!     =======================================
        open(20,file='para_fun.txt',status='old')
          read (20,*) RE
          read (20,*) M1
          read (20,*) ALFA
          read (20,*) BETA
        close(20)



! -------------------------------
       write (*,997) '======================================='
       write (*,997) '   Re:',Re
       write (*,997) ' Mach:',M1
       write (*,997) 'Alpha:',ALFA
       write (*,997) ' Beta:',BETA
       write (*,997) '   T1:',T1
       write (*,997) '   T2:',T2
       write (*,997) '   U1:',U1
       write (*,997) '   U2:',U2
       write (*,997) '======================================='
  997 format (a8,F9.3)
! -------------------------------

      X=1.0d-16
      V=1.0d-16
      VX=1.0E10

!********* DATA FILE **********************************************
         FA3=(INT(ALFA*100.0d0)               )/100
         FA2=(INT(ALFA*100.0d0)-FA3*100       )/10
         FA1=(INT(ALFA*100.0d0)-FA3*100-FA2*10)/1

         FB3=(INT(BETA*100.0d0)               )/100
         FB2=(INT(BETA*100.0d0)-FB3*100       )/10
         FB1=(INT(BETA*100.0d0)-FB3*100-FB2*10)/1

         FM3=(INT(M1*100.0d0)               )/100
         FM2=(INT(M1*100.0d0)-FM3*100       )/10
         FM1=(INT(M1*100.0d0)-FM3*100-FM2*10)/1

         FR6= INT(Re)                            /100000
         FR5=(INT(Re)-FR6*100000                )/10000
         FR4=(INT(Re)-FR6*100000-FR5*10000                 )/1000
         FR3=(INT(Re)-FR6*100000-FR5*10000-FR4*1000        )/100
         FR2=(INT(Re)-FR6*100000-FR5*10000-FR4*1000-FR3*100)/10
         FR1=(INT(Re)-FR6*100000-FR5*10000-FR4*1000-FR3*100-FR2*10)/1

         FNAME1= 'mM'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  '_b0am1funx.dat'
         FNAME2= 'mM'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  '_RHOmx.dat'
         FNAME3= 'mM'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  '_Umx.dat'
         FNAME4= 'mM'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  '_Pmx.dat'


         NAMEA1=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'rx.dat'
         NAMEA2=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'ux.dat'
         NAMEA3=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'vx.dat'
         NAMEA4=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'wx.dat'
         NAMEA5=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'Tx.dat'
         NAMEA6=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'px.dat'

         NAMEB1=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'rxb.dat'
         NAMEB2=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'uxb.dat'
         NAMEB3=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'vxb.dat'
         NAMEB4=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'wxb.dat'
         NAMEB5=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'Txb.dat'
         NAMEB6=  'M'//char(FM3+48)//'.'//char(FM2+48)//char(FM1+48)// &
                 'Re'//char(FR6+48)//char(FR5+48)//char(FR4+48) &
                     //char(FR3+48)//char(FR2+48)//char(FR1+48)// &
                  'a'//char(FA3+48)//'.'//char(FA2+48)//char(FA1+48)// &
                  'b'//char(FB3+48)//'.'//char(FB2+48)//char(FB1+48)// &
                  'pxb.dat'


      OPEN(10,FILE=FNAME1)
      OPEN(20,FILE=FNAME2)
      OPEN(30,FILE=FNAME3)
      OPEN(40,FILE=FNAME4)

      OPEN(51,FILE=NAMEA1)
      OPEN(52,FILE=NAMEA2)
      OPEN(53,FILE=NAMEA3)
      OPEN(54,FILE=NAMEA4)
      OPEN(55,FILE=NAMEA5)
      OPEN(56,FILE=NAMEA6)

      OPEN(61,FILE=NAMEB1)
      OPEN(62,FILE=NAMEB2)
      OPEN(63,FILE=NAMEB3)
      OPEN(64,FILE=NAMEB4)
      OPEN(65,FILE=NAMEB5)
      OPEN(66,FILE=NAMEB6)
!********* INPUT INITIAL CONDITION ********************************
      CALL CAL_mF(U,DUDX,DDUDX2,EIJ,CIJ,DIJ,y)
!******************************************************************
      DO K=0,N
!!ポアズイユ
!            y(K)=dcos(PI*DBLE(K)/DBLE(N))
 !              U(K)=  1.0d0 - y(K)*y(K)
!            DUDX(K)= -2.0d0*dcos(PI*DBLE(K)/DBLE(N))
!          DDUDX2(K)= -2.0d0
               T(K)= 1.0d0+0.5d0*(M1**2)*(GAMMA-1.0d0) &
                                *(0.0d0-(U(K)**2.0d0))
            DTDX(K)= -0.5d0*(M1**2)*(GAMMA-1.0d0) &
                                 *2.0d0*U(K)*DUDX(K)
          DDTDX2(K)= -0.5d0*(M1**2)*(GAMMA-1.0d0) &
                                 *2.0d0*(DUDX(K)**2) &
                     -0.5d0*(M1**2)*(GAMMA-1.0d0) &
                                 *2.0d0*U(K)*DDUDX2(K)
             ROU(K)= 1.0d0/T(K)
          DROUDX(K)=-T(K)**(-2.0)*DTDX(K)
             MYU(K)=             T(K)**( 2.0/3.0)
          DMYUDT(K)= 2.0d0/3.0d0*T(K)**(-1.0/3.0)
          DDMDT2(K)=-2.0d0/9.0d0*T(K)**(-4.0/3.0)
	end do
!     -----------------------------------
      DO K=0,N
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
	end do
!     -----------------------------------

      open (21,file='test_all.dat')
      DO K=0,N
!            y=dcos(PI*K/DBLE(N))
!            y(K)=dcos(PI*K/DBLE(N))
         write (21,998)  K,y(K),DBLE(U(K)),DBLE(DUDX(K)),DBLE(DDUDX2(K)) &
                            ,DBLE(T(K)),DBLE(DTDX(K)),DBLE(DDTDX2(K)) &
                            ,DBLE(ROU(K)),DBLE(DROUDX(K)) &
                           ,DBLE(MYU(K)),DBLE(DMYUDT(K)),DBLE(DDMDT2(K))
	end do
      close(21)

  998 format (I5,12E17.8)

      DO K=0,N
!            y(K)=dcos(PI*K/DBLE(N))
!            y(K)=dcos(PI*K/DBLE(N))
      WRITE(20,999) 1,K-N/2,Y(K),DBLE(ROU(K))
      WRITE(30,999) 1,K-N/2,Y(K),DBLE(  U(K))
      WRITE(40,999) 1,K-N/2,Y(K),DBLE(  P(K))

	end do

      CLOSE(20,STATUS='KEEP')
      CLOSE(30,STATUS='KEEP')
      CLOSE(40,STATUS='KEEP')



!**************** GYOLETU NO KEISAN ************************************


      DO K=0,N
      DO L=0,N

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
                        *(DUDX(K)*DDMDT2(K)*DTDX(K) &
                        +DMYUDT(K)*DDUDX2(K) &
                        +DUDX(K)*DMYUDT(K)*CIJ(K,K) )
!     &                  -DUDX(K)*DMYUDT(K)*DTDX(K)/(T(K)*T(K)))!"/(T(K)*T(K))"これいらなくない？

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
                        -ALFA/(ROU(K)*RE)*DUDX(K) &
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
                        *DUDX(K)*DUDX(K)/T(K) & !"/T(K)"がなぜあるか不明
                        -CI*GAMMA/(ROU(K)*RE*PR) &
                        *(MYU(K)*(ALFA*ALFA &
                        -DIJ(K,K)+BETA*BETA) &
                        -DMYUDT(K)*(2.0d0)*DTDX(K)*CIJ(K,K) &!論文と並びが違うだけ(2)
                        -DDMDT2(K)*DTDX(K)*DTDX(K) &!論文と並びが違うだけ(1)
                        -DMYUDT(K)*DDTDX2(K) )
!     &                  +DMYUDT(K)*DTDX(K)*DTDX(K)/(T(K)*T(K)))!これの存在意味がわからない>だからコメントアウトしてるのかも？

      ELSE

      AMAIN(5*K,5*L+2)=-CI*ROU(K)*CIJ(K,L)
      AMAIN(5*K+1,5*L+1)=CI/(ROU(K)*RE)*(MYU(K)*DIJ(K,L) &
                        +DMYUDT(K)*DTDX(K)*CIJ(K,L))
      AMAIN(5*K+1,5*L+2)=-1.0d0/(3.0d0*ROU(K)*RE) &
                        *MYU(K)*ALFA*CIJ(K,L)
      AMAIN(5*K+1,5*L+4)=CI/(ROU(K)*RE)*DUDX(K) &
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
                        +DMYUDT(K)*(2.0d0)*DTDX(K)*CIJ(K,L))

      ENDIF

	end do
	end do

!*******************BOUNDARY CONDITIONS*********************
      DO K=0,LDA-1
       DO L=0,LDA-1
        BMAIN(K,L)=0.0d0
	end do
	end do
      DO K=0,LDA-1
        BMAIN(K,K)=1.0d0
	end do

      K=0
      DO L=0,LDA-1
!        AMAIN(5*K+0,L)=0.0d0
        AMAIN(5*K+1,L)=0.0d0
        AMAIN(5*K+2,L)=0.0d0
        AMAIN(5*K+3,L)=0.0d0
        AMAIN(5*K+4,L)=0.0d0
	end do
!        BMAIN(5*K+0,5*K+0)=0.0d0
        BMAIN(5*K+1,5*K+1)=0.0d0
        BMAIN(5*K+2,5*K+2)=0.0d0
        BMAIN(5*K+3,5*K+3)=0.0d0
        BMAIN(5*K+4,5*K+4)=0.0d0

      DO L=0,N
!        AMAIN(5*L+0,5*K+0)=EIJ(K,L)
!        AMAIN(5*L+1,5*K+1)=EIJ(K,L)
!        AMAIN(5*L+2,5*K+2)=EIJ(K,L)
!        AMAIN(5*L+3,5*K+3)=EIJ(K,L)
!        AMAIN(5*L+4,5*K+4)=EIJ(K,L)

        AMAIN(5*K+1,5*L+1)=EIJ(K,L)
        AMAIN(5*K+2,5*L+2)=EIJ(K,L)
        AMAIN(5*K+3,5*L+3)=EIJ(K,L)
        AMAIN(5*K+4,5*L+4)=EIJ(K,L)
	end do

      IF (0.EQ.0) THEN
      K=1
      DO L=0,LDA-1
!        AMAIN(5*K+0,L)=0.0d0
        AMAIN(5*K+1,L)=0.0d0
        AMAIN(5*K+2,L)=0.0d0
        AMAIN(5*K+3,L)=0.0d0
        AMAIN(5*K+4,L)=0.0d0
	end do
!        BMAIN(5*K+0,5*K+0)=0.0d0
        BMAIN(5*K+1,5*K+1)=0.0d0
        BMAIN(5*K+2,5*K+2)=0.0d0
        BMAIN(5*K+3,5*K+3)=0.0d0
        BMAIN(5*K+4,5*K+4)=0.0d0

      DO L=0,N
!        AMAIN(5*L+0,5*K+0)=EIJ(K,L)
!        AMAIN(5*L+1,5*K+1)=EIJ(K,L)
!        AMAIN(5*L+2,5*K+2)=EIJ(K,L)
!        AMAIN(5*L+3,5*K+3)=EIJ(K,L)
!        AMAIN(5*L+4,5*K+4)=EIJ(K,L)

        AMAIN(5*K+1,5*L+1)=CIJ(K-1,L)
        AMAIN(5*K+2,5*L+2)=CIJ(K-1,L)
        AMAIN(5*K+3,5*L+3)=CIJ(K-1,L)
        AMAIN(5*K+4,5*L+4)=CIJ(K-1,L)
	end do
      ENDIF


      K=N
      DO L=0,LDA-1
!        AMAIN(5*K+0,L)=0.0d0
        AMAIN(5*K+1,L)=0.0d0
        AMAIN(5*K+2,L)=0.0d0
        AMAIN(5*K+3,L)=0.0d0
        AMAIN(5*K+4,L)=0.0d0
	end do

!        BMAIN(5*K+0,5*K+0)=0.0d0
        BMAIN(5*K+1,5*K+1)=0.0d0
        BMAIN(5*K+2,5*K+2)=0.0d0
        BMAIN(5*K+3,5*K+3)=0.0d0
        BMAIN(5*K+4,5*K+4)=0.0d0

      DO L=0,N
!        AMAIN(5*L+0,5*K+0)=EIJ(K,L)
!        AMAIN(5*L+1,5*K+1)=EIJ(K,L)
!        AMAIN(5*L+2,5*K+2)=EIJ(K,L)
!        AMAIN(5*L+3,5*K+3)=EIJ(K,L)
!        AMAIN(5*L+4,5*K+4)=EIJ(K,L)

        AMAIN(5*K+1,5*L+1)=EIJ(K,L)
        AMAIN(5*K+2,5*L+2)=EIJ(K,L)
        AMAIN(5*K+3,5*L+3)=EIJ(K,L)
        AMAIN(5*K+4,5*L+4)=EIJ(K,L)
	end do

      IF (0.EQ.0) THEN
      K=N-1
      DO L=0,LDA-1
!        AMAIN(5*K+0,L)=0.0d0
        AMAIN(5*K+1,L)=0.0d0
        AMAIN(5*K+2,L)=0.0d0
        AMAIN(5*K+3,L)=0.0d0
        AMAIN(5*K+4,L)=0.0d0
	end do

!        BMAIN(5*K+0,5*K+0)=0.0d0
        BMAIN(5*K+1,5*K+1)=0.0d0
        BMAIN(5*K+2,5*K+2)=0.0d0
        BMAIN(5*K+3,5*K+3)=0.0d0
        BMAIN(5*K+4,5*K+4)=0.0d0

      DO L=0,N
!        AMAIN(5*L+0,5*K+0)=EIJ(K,L)
!        AMAIN(5*L+1,5*K+1)=EIJ(K,L)
!        AMAIN(5*L+2,5*K+2)=EIJ(K,L)
!        AMAIN(5*L+3,5*K+3)=EIJ(K,L)
!        AMAIN(5*L+4,5*K+4)=EIJ(K,L)

        AMAIN(5*K+1,5*L+1)=CIJ(K+1,L)
        AMAIN(5*K+2,5*L+2)=CIJ(K+1,L)
        AMAIN(5*K+3,5*L+3)=CIJ(K+1,L)
        AMAIN(5*K+4,5*L+4)=CIJ(K+1,L)
	end do
      ENDIF
!***************GYOLETU NO HOSEI*****************************
      DO K=0,LDA-1
      DO L=0,LDA-1

      TES(K,L)=DBLE(AMAIN(K,L))
      IF (ABS(TES(K,L)).LT.V) THEN !LTは<を意味する
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

      AMAIN(K,L)=CMPLX((MAR(K,L)),(MAI(K,L)))
	end do
	end do
!**********************KOYUCHI NO KEISAN***************************

!      CALL DEVCCG (LDA,AMAIN,LDA,EVAL,EVEC,LDEVEC)
!      CALL DGVCCG (LDA,AMAIN,LDA,BMAIN,LDA,EALPHA,EBETA,EVEC,LDA)
	call ZGEGV( 'N', 'V', LDA, AMAIN, LDA, BMAIN, LDA, EALPHA, EBETA, &
       & 0, LDVL, EVEC, LDEVEC, WORK, LWORK, RWORK, INFO )

      DO I = 0,LDA-1
         IF (EBETA(I).EQ.0.0d0) THEN
            EVAL(I)=999.9d0
         ELSE
            EVAL(I)=(EALPHA(I)/EBETA(I))
         ENDIF
	end do

!****************W(I) NO CHUSHUTU**********************************
      DO I=0,LDA-1
         IVAL(I)=DIMAG(EVAL(I))
         RVAL(I)= DBLE(EVAL(I))
         IF(RVAL(I).NE.0.0d0) THEN
            CRE(I)=(RVAL(I))/ALFA
         ELSE
            CRE(I)=1.E3
         ENDIF
	end do


!***************subsonic mode**************************************
      OMEGAI=-0.5d0
      ISB=1

      DO I=0,LDA-1
         IF (CRE(I).GT.U_min.AND.CRE(I).LT.U_max) THEN
            IF ((IVAL(I).GT.OMEGAI).AND.(IVAL(I).LT.1.0d0)) THEN
               OMEGAI=IVAL(I)
               ISB=I
            ENDIF
         ENDIF
	end do

      OMEGA2=-0.5d0
      DO I=0,LDA-1
         IF (CRE(I).GT.U_min.AND.CRE(I).LT.U_max) THEN
            IF ((IVAL(I).GT.OMEGA2).AND.(IVAL(I).LT.1.0d0) &
                                              .and.(I.NE.ISB)) THEN
               OMEGA2=IVAL(I)
               ISB2=I
            ENDIF
         ENDIF
	end do

        write (*,*)  'ISB ',ISB ,IVAL(ISB )
        write (*,*)  'ISB2',ISB2,IVAL(ISB2)

      MIVAL(1)=OMEGAI
       MCRE(1)=CRE(ISB)

      WRITE(10,210)M1,ALFA,BETA
  210 FORMAT(/,1X,'M1=',F4.1,2X,'ALFA=',F5.2,2X,'BETA=',F5.2/)

      WRITE(10,213)MIVAL(1),MCRE(1)
  213 FORMAT(/,1X,'OMEGAi=',G10.3,2X,'Cr=',G10.3/)


!***** Eigen function *******************************************

      DO IW=1,2
         IF (IW.EQ.1) THEN
            ISB=ISB
            ID=0
         ELSE
            ISB=ISB2
            ID=10
         ENDIF

      DO J=0,4

      DO I=0,N

!---[subsonic mode]---

              RVECA(I)= DBLE(EVEC(I*5+J,ISB))
              IVECA(I)=DIMAG(EVEC(I*5+J,ISB))
                AMP(I)=DSQRT((RVECA(I))**2.0d0+(IVECA(I))**2.0d0)

	end do


      DO I=0,N
        IAMP(I)=AMP(I)
	end do

      DO I=0,N

              ISP(I,1)= IAMP(I)
              ISP(I,2)=RVECA(I)
              ISP(I,3)=IVECA(I)

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
	end do

!***** Y ********************************************************
      DO K=0,N
              IF (K.EQ.N/2) THEN
                EX(K)=0.0d0
              ELSE
                EX(K)=dcos(PI*DBLE(K)/DBLE(N))
              ENDIF
	end do

!***** Out put ***************************************************
      if (j.EQ.0) then
      DO I=0,N
        WRITE(51+ID,999) ISB,I-N/2,Y(I),ISP(I,1),ISP(I,2),ISP(I,3)
	end do

      elseif (j.EQ.1) then
      DO I=0,N
        WRITE(52+ID,999) ISB,I-N/2,Y(I),ISP(I,1),ISP(I,2),ISP(I,3)
	end do

      elseif (j.EQ.2) then
      DO I=0,N
        WRITE(53+ID,999) ISB,I-N/2,Y(I),ISP(I,1),ISP(I,2),ISP(I,3)
	end do

      elseif (j.EQ.3) then
      DO I=0,N
        WRITE(54+ID,999) ISB,I-N/2,Y(I),ISP(I,1),ISP(I,2),ISP(I,3)
	end do

      elseif (j.EQ.4) then
      DO I=0,N
        WRITE(55+ID,999) ISB,I-N/2,Y(I),ISP(I,1),ISP(I,2),ISP(I,3)
	end do

      endif

  999 format (2I5,5E25.16)

	end do


!%%%%%%%%%%%%%%%%%%%%       PRESSURE       %%%%%%%%%%%%%%%%%%%%%

      DO I=0,N
!---[subsonic mode]---

      RVECA(I)=DBLE( &
                EVEC(I*5+0,ISB)*  T(I) &
               +EVEC(I*5+4,ISB)*ROU(I))
      IVECA(I)=DIMAG( &
                EVEC(I*5+0,ISB)*  T(I) &
               +EVEC(I*5+4,ISB)*ROU(I))
      AMP(I)=DSQRT( ((RVECA(I))**2.0d0) + ((IVECA(I))**2.0d0) )

	end do

      DO I=0,N
      IAMP(I)=AMP(I)
	end do

      DO I=0,N

!---[subsonic mode]---

      ISP(I,1)= IAMP(I)
      ISP(I,2)=RVECA(I)
      ISP(I,3)=IVECA(I)

!           if (I.eq.N/2) then
!              ISO(I)=0.0d0
!           else
!              ISO(I)=datan2(ISP(2*N-2+I),ISP(N-1+I))
!           endif

	end do

      DO I=0,N
        WRITE(56+ID,999) 1,I-N/2,Y(I),ISP(I,1),ISP(I,2),ISP(I,3)
	end do

	end do

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

!      STOP
 !     END
contains
!///////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE CAL_mF(ZU,ZdUdy,ZddUdy,ZEIJ,ZCIJ,ZDIJ,y)
use Para_
      IMPLICIT NONE
      INTEGER :: I,J,K,L,M!,N
!      integer, PARAMETER:: N=128
!      INTEGER By
!      integer, PARAMETER:: By=10000
      REAL(8)::      EIJ(0:N,0:N), CIJ(0:N,0:N), DIJ(0:N,0:N)
      complex(kind(0d0)):: ZEIJ(0:N,0:N),ZCIJ(0:N,0:N),ZDIJ(0:N,0:N)
      REAL(8):: se,s,ss,coef(0:N)
      REAL(8):: PI
      REAL(8):: sum1,sum2,sum3
      REAL(8)::      y(0:N)
      REAL(8)::    Eta(0:N)
      REAL(8)::  dEtdy(0:N)
      REAL(8):: ddEtdy(0:N)
      REAL(8)::      U(0:N)
      REAL(8)::   dUdy(0:N)
      REAL(8)::  ddUdy(0:N)
      complex(kind(0d0)):: ZU(0:N)
      complex(kind(0d0)):: ZdUdy(0:N)
      complex(kind(0d0)):: ZddUdy(0:N)
      REAL(8)::     HU(0:By)
call parapara()!      REAL(8):: Re_d,s_para,Ly
!      COMMON /mUmU/Re_d,s_para,Ly
!***********************************************************************
      PI=4.0d0*DATAN(1.0d0)
!     s_para=4.0d0
!      Re_d=1000.0d0
!      Ly=50.0d0
!*************************
      CALL CAL_HU(HU)
!*************************
      DO I=0,N-1
            Eta(I)=DCOS(DBLE(I)*PI/DBLE(N))
              y(I)=-s_para*DLOG((Eta(I)+1.0d0)/2.0d0)
          dEtdy(I)=-1.0d0/s_para*(Eta(I)+1.0d0)
         ddEtdy(I)= 1.0d0/(s_para**2)*(Eta(I)+1.0d0)
!         write(*,*)I,Eta(I),y(I)
	end do
          I=N
            Eta(I)=-0.99999d0
              y(I)=-s_para*DLOG((Eta(I)+1.0d0)/2.0d0)
          dEtdy(I)=-1.0d0/s_para*(Eta(I)+1.0d0)
         ddEtdy(I)= 1.0d0/(s_para**2)*(Eta(I)+1.0d0)
!*************************
      CALL CAL_mU(y,HU,U)
!*************************
!************************E,D NO KEISAN**********************************
       pi=4.0d0*datan(1.0d0)
       coef(0)=2.0d0
       coef(N)=2.0d0
       do j=1,N-1
       coef(j)=1.0d0!吉野論文のC_lのこと
	end do

       do m=0,N
       do l=0,N
      se=0.0d0
       s=0.0d0
       do k=0,N
       se=se+dcos(PI*DBLE(m)/DBLE(N)*DBLE(k)) &
            *dcos(PI*DBLE(k)/DBLE(N)*DBLE(l))/coef(l)/coef(k)
       do j=k+1,N,2
        s=s+dble(j)*dcos(pi*DBLE(j)/DBLE(N)*DBLE(l)) &
                   *dcos(pi*DBLE(k)/DBLE(N)*DBLE(m)) &
                   /coef(j)/coef(k)/coef(l)
	end do
	end do
       Eij(m,l)=2.0d0*se/dble(N)!一体なんのか不明
       Cij(m,l)=4.0d0* s/dble(N)!吉野論文のc_ijのこと。Σの計算はsで済ませている。ユーキのc_ijとは計算式が異なる。
	end do
	end do

!************************DD NO KEISAN***********************************
       do m=0,N
       do l=0,N
       s=0.0d0
       do k=0,N
       do j=k+2,N,2
       s=s+dble(j)*dble(j*j-k*k)*dcos(pi*DBLE(j)/DBLE(N)*DBLE(l)) &
                                *dcos(pi*DBLE(k)/DBLE(N)*DBLE(m)) &
         /coef(j)/coef(k)/coef(l)
	end do
	end do
       Dij(m,l)=2.0d0*s/dble(N)!吉野論文のd_ijのこと。Σの計算はsで済ませている。ユーキのd_ijとは計算式が異なる。
	end do
	end do
!***********************************************************************
      DO I=0,N
         DO J=0,N
            Dij(I,J)=Dij(I,J)*(dEtdy(I)**2)+Cij(I,J)*ddEtdy(I)
            Cij(I,J)=Cij(I,J)* dEtdy(I)
            ZEij(I,J)=Eij(I,J)
            ZCij(I,J)=Cij(I,J)
            ZDij(I,J)=Dij(I,J)
	end do
	end do



      DO I=0,N
         sum1=0.0d0
         sum2=0.0d0
         sum3=0.0d0
         DO J=0,N
            sum1=sum1+Eij(I,J)*U(J)
            sum2=sum2+Cij(I,J)*U(J)
            sum3=sum3+Dij(I,J)*U(J)
	end do
               U(I)=sum1
            dUdy(I)=sum2
           ddUdy(I)=sum3
	end do
          ddUdy(0)=0.0d0


      DO I=0,N
             ZU(I)=    U(I)
          ZdUdy(I)= dUdy(I)
         ZddUdy(I)=ddUdy(I)
	end do

      open(50,file='mean_flow.dat')
      DO I=0,N
         write(50,999) I,y(I),U(I),dUdy(I),ddUdy(I)
	end do
      close(50)

  999 FORMAT(I5,5E25.16E3)

      RETURN
      END SUBROUTINE CAL_mF
!///////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE CAL_HU(U)
use Para_
      IMPLICIT NONE
!      INTEGER By
      integer, PARAMETER:: By=10000
      INTEGER:: Iy
      REAL(8):: dy,Leta,deta,a,Re_x,delta,x_x
      REAL(8):: Eta(0:By)
      REAL(8)::   y(0:By)
      REAL(8)::   U(0:By)
      REAL(8)::   V(0:By)
      REAL(8)::   F(0:By),FA,FB,FC
      REAL(8):: dF1(0:By),dF1A,dF1B,dF1C
      REAL(8):: dF2(0:By),dF2A,dF2B,dF2C
      REAL(8):: dF3(0:By),dF3A,dF3B,dF3C
call parapara()!      REAL(8):: Re_d,s_para,Ly
!      COMMON /mUmU/Re_d,s_para,Ly
!--------------------------
!      Re_d=1000.0d0
      Re_x=(Re_d/1.7208d0)**2.0d0
       x_x=Re_x/Re_d
!--------------------------
!      Ly=50.0d0
      dy=Ly/DBLE(By)

      Leta=Ly*DSQRT(Re_x)/x_x
      deta=Leta/DBLE(By-1)


         F(0)=0.0d0
       dF1(0)=0.0d0
       dF2(0)=0.33d0
             a=1.0d-3

         DO Iy=0,By
            Eta(Iy)=DBLE(Iy)*deta
	end do

  600 CONTINUE

         DO Iy=0,By-1
!--------- f'''=-f''f/2
!-------1--------
             dF3(Iy) = -dF2(Iy)*F(Iy)/2.0d0

               FA =   F(Iy)+dF1(Iy)*deta/2.0d0
             dF1A = dF1(Iy)+dF2(Iy)*deta/2.0d0
             dF2A = dF2(Iy)+dF3(Iy)*deta/2.0d0
             dF3A = -dF2A*FA/2.0d0
!-------2--------
               FB =   F(Iy)+dF1A*deta/2.0d0
             dF1B = dF1(Iy)+dF2A*deta/2.0d0
             dF2B = dF2(Iy)+dF3A*deta/2.0d0
             dF3B = -dF2B*FB/2.0d0
!-------3--------
               FC =   F(Iy)+dF1B*deta
             dF1C = dF1(Iy)+dF2B*deta
             dF2C = dF2(Iy)+dF3B*deta
             dF3C = -dF2C*FC/2.0d0

               F(Iy+1) =   F(Iy) &
                 +(dF1(Iy)+2.0d0*dF1A+2.0d0*dF1B+dF1C)/6.0d0*deta
             dF1(Iy+1) = dF1(Iy) &
                 +(dF2(Iy)+2.0d0*dF2A+2.0d0*dF2B+dF2C)/6.0d0*deta
             dF2(Iy+1) = dF2(Iy) &
                 +(dF3(Iy)+2.0d0*dF3A+2.0d0*dF3B+dF3C)/6.0d0*deta

!          write (*,*) Iy,DBLE(Iy)*dy,dF1(Iy)
	end do

      IF (dF1(By).GT.1.0d0) THEN
         dF2(0)= dF2(0)-a
         a=a*0.1d0
         dF2(0)= dF2(0)+a
         GOTO 600
      ELSEIF (DABS(dF1(By)-1.0d0).LT.1.0d-14) THEN
         GOTO 700
      ELSE
         dF2(0)= dF2(0)+a
         GOTO 600
      ENDIF


  700 CONTINUE
         open (50,file='check_HU.dat')
         DO Iy=0,By
           y(Iy)=Eta(Iy)/DSQRT(Re_x)*x_x
           U(Iy)=dF1(Iy)
           V(Iy)=0.5d0*DSQRT(1.0d0/(Re_x*x_x))*(Eta(Iy)*dF1(Iy)-F(Iy))
           write (50,*) y(Iy),U(Iy)
	end do
         close(50)

      RETURN
      END  SUBROUTINE CAL_HU

!///////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////
      SUBROUTINE CAL_mU(y,HU,U)
use Para_
      IMPLICIT NONE
      INTEGER :: Iy,Jy,ID,JD
!      integer, PARAMETER:: N=128
!      INTEGER By
!      integer, PARAMETER:: By=10000
      REAL(8)::  Y(0:N),YY(0:By)
      REAL(8)::  dy
      REAL(8)::    U(0:N)
      REAL(8)::   HU(0:By)
!      REAL(8):: Re_d,s_para,Ly
!      COMMON /mUmU/Re_d,s_para,Ly
call parapara()
! -----------------------------------------------------------
      dy=Ly/DBLE(By)
! -----------------------------------------------------------
      do Iy=0,By
         YY(Iy) = DBLE(Iy)*dy
	end do
!c -----------------------------------------------------------

       DO Iy=0,N
      DO Jy=0,By-1
          IF (Y(Iy).GE.YY(Jy).and.Y(Iy).LE.YY(Jy+1)) THEN
               U(Iy)= (Y(Iy)-YY(Jy))/(YY(Jy+1)-YY(Jy))*(HU(Jy+1)-HU(Jy)) &
                      +HU(Jy)
          ENDIF
          IF (Y(Iy).GE.Ly) THEN
               U(Iy)=1.0d0
          ENDIF
	end do
	end do


      RETURN
      END  SUBROUTINE CAL_mU
end PROGRAM EIGEN
