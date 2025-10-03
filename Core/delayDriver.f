C * * * * * * * * * * * * * * * 
C --- DRIVER FOR RADAR5 
C * * * * * * * * * * * * * * *
C       Use make command (relevant to Makefile)

        IMPLICIT REAL*8 (A-H,O-Z)

        INTEGER, PARAMETER :: DP=kind(1D0)
C --->  PARAMETERS FOR RADAR5 (FULL JACOBIAN) <---
        INTEGER, PARAMETER :: ND=8
        INTEGER, PARAMETER :: NRDENS=1
        INTEGER, PARAMETER :: NGRID=1
        INTEGER, PARAMETER :: NLAGS=1
        INTEGER, PARAMETER :: NJACL=1
        INTEGER, PARAMETER :: MXST=60000
        INTEGER, PARAMETER :: LWORK=30
        INTEGER, PARAMETER :: LIWORK=30
        INTEGER, PARAMETER :: NRPAR=60000
        REAL(kind=DP), dimension(ND) :: Y
        REAL(kind=DP), dimension(NGRID+1) :: GRID
        REAL(kind=DP), dimension(LWORK) :: WORK
        INTEGER, dimension(LIWORK) :: IWORK
        INTEGER, dimension(NRDENS+1) :: IPAST
        REAL(kind=DP), dimension(NRPAR) :: RPAR
        REAL(kind=DP), dimension(3) :: IPAR
        INTEGER, dimension(22) :: ISTAT
        REAL(kind=DP), dimension(:), allocatable :: PAR(:)
        REAL(kind=DP), dimension(:), allocatable :: INIT(:)
        REAL(kind=DP), dimension(:), allocatable :: TIME(:)
        REAL(kind=DP), dimension(:), allocatable :: PDATA(:)
        REAL(kind=DP), dimension(:), allocatable :: RDATA(:)
        REAL(kind=DP), dimension(:), allocatable :: PTHDATA(:)
        REAL(kind=DP), dimension(:), allocatable :: DPCDTDATA(:)
        REAL(kind=DP), dimension(:), allocatable :: DPTHDTDATA(:)
        REAL(kind=DP), dimension(:), allocatable :: PDSPL(:)
        REAL(kind=DP), dimension(:), allocatable :: PTHSPL(:)
        REAL(kind=DP), dimension(:), allocatable :: DPCDTSPL(:)
        REAL(kind=DP), dimension(:), allocatable :: DPTHDTSPL(:)
        REAL(kind=DP), dimension(:), allocatable :: RSPL(:)
        INTEGER :: NPAR,NINIT,NTIME,NPDATA,NPTHDATA,NDPCDTDATA
        INTEGER :: NDPTHDTDATA
        INTEGER :: NRDATA
        REAL(kind=DP) :: Delay,dt
        EXTERNAL  FCN,PHI,ARGLAG,JFCN,JACLAG,SOLOUT

        
C ------ FILE TO OPEN ----------
        OPEN(9,FILE='sol.out')
        OPEN(10,FILE='cont.out')
        REWIND 9
        REWIND 10
        
        OPEN(31,file="Pars.txt")
        READ(31,*) NPAR
        ALLOCATE(PAR(NPAR))
        READ(31,*) PAR
        CLOSE(31)

        OPEN(32,file="Init.txt")
        READ(32,*) NINIT
        ALLOCATE(INIT(NINIT))
        READ(32,*) INIT
        CLOSE(32)

        OPEN(33,file="Time.txt")
        READ(33,*) NTIME
        ALLOCATE(TIME(NTIME))
        READ(33,*) TIME
        CLOSE(33)

        OPEN(34,file="dt.txt")
        READ(34,*) dt
        CLOSE(34)

        OPEN(35,file="Pdata.txt")
        READ(35,*) NPDATA
        ALLOCATE(PDATA(NPDATA))
        READ(35,*) PDATA
        CLOSE(35)

        OPEN(36,file="Pthdata.txt")
        READ(36,*) NPTHDATA
        ALLOCATE(PTHDATA(NPTHDATA))
        READ(36,*) PTHDATA
        CLOSE(36)

        OPEN(37,file="Rdata.txt")
        READ(37,*) NRDATA
        ALLOCATE(RDATA(NRDATA))
        READ(37,*) RDATA
        CLOSE(37)
        
        OPEN(38,file="dPcdtdata.txt")
        READ(38,*) NDPCDTDATA
        ALLOCATE(DPCDTDATA(NDPCDTDATA))
        READ(38,*) DPCDTDATA
        CLOSE(38)
        
        OPEN(39,file="dPthdtdata.txt")
        READ(39,*) NDPTHDTDATA
        ALLOCATE(DPTHDTDATA(NDPTHDTDATA))
        READ(39,*) DPTHDTDATA
        CLOSE(39)

! --- ALLOCATE SPACE FOR INTERPOLANTS        
        ALLOCATE(PDSPL(NPDATA))
        ALLOCATE(PTHSPL(NPTHDATA))
        ALLOCATE(RSPL(NRDATA))
        ALLOCATE(DPCDTSPL(NDPCDTDATA))
        ALLOCATE(DPTHDTSPL(NDPTHDTDATA))
        CALL spline(TIME,PDATA,NPDATA,0,0,PDSPL)
        CALL spline(TIME,PTHDATA,NPTHDATA,0,0,PTHSPL)
        CALL spline(TIME,DPCDTDATA,NDPCDTDATA,0,0,DPCDTSPL)
        CALL spline(TIME,DPTHDTDATA,NDPTHDTDATA,0,0,DPTHDTSPL)
        CALL spline(TIME,RDATA,NRDATA,0,0,RSPL)
        
        
C --- Build Parameter Vector ----
        RPAR(1) = dt
        DO I=1,NPAR
           RPAR(I+1) = PAR(I)
        END DO
        DO I=1,NINIT
           RPAR(I+1+NPAR) = INIT(I)
        END DO
        DO I=1,NTIME
           RPAR(I+1+NPAR+NINIT) = TIME(I)
        END DO
        DO I=1,NPDATA
           RPAR(I+1+NPAR+NINIT+NTIME) = PDATA(I)
        END DO
        DO I=1,NPTHDATA
           RPAR(I+1+NPAR+NINIT+2*NTIME) = PTHDATA(I)
        END DO
        DO I=1,NPDATA
           RPAR(I+1+NPAR+NINIT+3*NTIME) = PDSPL(I)
        END DO
        DO I=1,NPTHDATA
           RPAR(I+1+NPAR+NINIT+4*NTIME) = PTHSPL(I)
        END DO
        DO I=1,NRDATA
           RPAR(I+1+NPAR+NINIT+5*NTIME) = RDATA(I)
        END DO
        DO I=1,NRDATA
           RPAR(I+1+NPAR+NINIT+6*NTIME) = RSPL(I)
        END DO
        DO I=1,NDPCDTDATA
           RPAR(I+1+NPAR+NINIT+7*NTIME) = DPCDTDATA(I)
        END DO
        DO I=1,NDPCDTDATA
           RPAR(I+1+NPAR+NINIT+8*NTIME) = DPCDTSPL(I)
        END DO
        DO I=1,NDPTHDTDATA
           RPAR(I+1+NPAR+NINIT+9*NTIME) = DPTHDTDATA(I)
        END DO
        DO I=1,NDPTHDTDATA
           RPAR(I+1+NPAR+NINIT+10*NTIME) = DPTHDTSPL(I)
        END DO
        
C --- Vector of lengths of each input
        IPAR(1) = NPAR
        IPAR(2) = NINIT
        IPAR(3) = NTIME
        
C --- DIMENSION OF THE SYSTEM
        N=ND
C --- COMPUTE THE JACOBIAN USING FINITE DIFFERENCES
        IJAC=0
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXLPICIT FORM
        IMAS=0
        MLMAS=N
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- Delay
        Delay = PAR(28)
C --- Initial Conditions        
        Y = INIT        
C --- INITIAL VALUES 
        X=TIME(1)
C --- ENDPOINT OF INTEGRATION
        XEND=TIME(NTIME)
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.0D-6
C       ATOL=RTOL*1.0D0
        ATOL = RTOL*1.0D0
C --- INITIAL STEP SIZE
        H=1.0D-4
C --- DEFAULT VALUES FOR PARAMETERS
        DO I=1,20
           IWORK(I)=0
           WORK(I)=0.0D0
        END DO
        
C       IWORK(8) = 1
C       IWORK(11) = 2
        
C --- WORKSPACE FOR PAST 
        IWORK(12)=MXST
C --- THE SEVENTH COMPONENT USES RETARDED ARGUMENTS
        IWORK(15)=NRDENS
        IPAST(1)=7
C ---  SET THE PRESCRIBED GRID-POINTS
        DO I=1,NGRID
          GRID(I)=Delay*I
        END DO
C --- WORKSPACE FOR GRID
        IWORK(13)=NGRID

C --- CALL OF THE SUBROUTINE RADAR5   
        CALL RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JFCN,IJAC,MLJAC,MUJAC,
     &                  JACLAG,NLAGS,NJACL,
     &                  IMAS,SOLOUT,IOUT,
     &                  WORK,IWORK,RPAR,IPAR,IDID,
     &                  GRID,IPAST,DUMMY,MLMAS,MUMAS)

!C --- PRINT FINAL SOLUTION SOLUTION
        WRITE (6,*) X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8)
        WRITE(6,*)' ***** TOL=',RTOL,' ****'
        WRITE(6,*) 'SOLUTION IS TABULATED IN FILES: sol.out & cont.out'

        END PROGRAM

C----------------------------------------------------------------------
        SUBROUTINE SOLOUT (NR,XOLD,X,HSOL,Y,CONT,LRC,N,
     &                     RPAR,IPAR,IRTRN)
C ----- PRINTS THE DISCRETE OUTPUT AND THE CONTINUOUS OUTPUT
C       AT EQUIDISTANT OUTPUT-POINTS
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP) :: XSTEP,XEND
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(LRC) :: CONT
        REAL(kind=DP), dimension(60000) :: RPAR
        REAL(kind=DP), dimension(3) :: IPAR
        REAL(kind=DP), dimension(:), ALLOCATABLE :: TIME(:)
        EXTERNAL PHI
C       XOUT IS USED FOR THE DENSE OUTPUT
        COMMON /INTERN/XOUT

        XSTEP = RPAR(1)

        NPAR  = IPAR(1)
        NINIT = IPAR(2)
        NTIME = IPAR(3)
        
        ALLOCATE(TIME(NTIME))
        DO I=1,NTIME
           TIME(I) = RPAR(1+NPAR+NINIT+I)
        END DO

        XEND = TIME(NTIME)
        
        WRITE (9,*) X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8)

        IF (NR.EQ.1) THEN
           WRITE (10,*) X,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8)
           XOUT=XSTEP
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (10,*) XOUT,CONTR5(1,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(2,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(3,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(4,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(5,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(6,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(7,N,XOUT,CONT,X,HSOL),
     &                           CONTR5(8,N,XOUT,CONT,X,HSOL)

              XOUT=XOUT+XSTEP
              GOTO 10
           END IF
        END IF

        IF (X .EQ. XEND) THEN
                     WRITE (10,*) XEND,CONTR5(1,N,XEND,CONT,X,HSOL),
     &                           CONTR5(2,N,XEND,CONT,X,HSOL),
     &                           CONTR5(3,N,XEND,CONT,X,HSOL),
     &                           CONTR5(4,N,XEND,CONT,X,HSOL),
     &                           CONTR5(5,N,XEND,CONT,X,HSOL),
     &                           CONTR5(6,N,XEND,CONT,X,HSOL),
     &                           CONTR5(7,N,XEND,CONT,X,HSOL),
     &                           CONTR5(8,N,XEND,CONT,X,HSOL)

        END IF
             
        
        RETURN
        END
C
        FUNCTION ARGLAG(IL,X,Y,RPAR,IPAR,PHI,PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(8) :: Y
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        REAL(kind=DP), dimension(60000) :: RPAR
        INTEGER, dimension(3) :: IPAR
        
        ARGLAG=X-RPAR(1+28) !Value of the deviating argument (delay); has to do with step-size

        RETURN
        END
C----------------------------------------------------------------------
        SUBROUTINE FCN(N,X,Y,F,ARGLAG,PHI,RPAR,IPAR,
     &                  PAST,IPAST,NRDS)
        IMPLICIT REAL*8 (A-H,K,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), PARAMETER :: PI=3.1415926536
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(N) :: F
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        REAL(kind=DP), dimension(60000) :: RPAR
        REAL(kind=DP), dimension(:), allocatable :: PAR(:)
        REAL(kind=DP), dimension(:), allocatable :: TIME(:)
        REAL(kind=DP), dimension(:), allocatable :: PD(:)
        REAL(kind=DP), dimension(:), allocatable :: PTH(:)
        REAL(kind=DP), dimension(:), allocatable :: R(:)
        REAL(kind=DP), dimension(:), allocatable :: PDSPL(:)
        REAL(kind=DP), dimension(:), allocatable :: PTHSPL(:)
        REAL(kind=DP), dimension(:), allocatable :: RSPL(:)
        REAL(kind=DP), dimension(:), allocatable :: DPCDT(:)
        REAL(kind=DP), dimension(:), allocatable :: DPCDTSPL(:)
        REAL(kind=DP), dimension(:), allocatable :: DPTHDT(:)
        REAL(kind=DP), dimension(:), allocatable :: DPTHDTSPL(:)
        REAL(kind=DP), dimension(3) :: IPAR
        REAL(kind=DP) :: Pc, Pthor, Resp, dPcdtV, dPthdtV, dPadt
        REAL(kind=DP) :: ewc, Pa, ewa, nm, Gp, Gs
        REAL(kind=DP) :: Gr, Htilde, Hs_new, taupb_new, osc
        EXTERNAL PHI,ARGLAG
    
C       It's dPcdtV and dPthdtV so the variable names aren't the same.

        dt = RPAR(1)
        
        NPAR  = IPAR(1)
        NINIT = IPAR(2)
        NTIME = IPAR(3)
        
C       Get values from stored information in RPAR
        
        ALLOCATE(PAR(NPAR))
        DO I=1,NPAR
           PAR(I) = RPAR(1+I)
        END DO
        ALLOCATE(TIME(NTIME))
        DO I=1,NTIME
           TIME(I) = RPAR(1+NPAR+NINIT+I)
        END DO
        ALLOCATE(PD(NTIME))
        DO I=1,NTIME
           PD(I) = RPAR(1+NPAR+NINIT+NTIME+I)
        END DO
        ALLOCATE(PTH(NTIME))
        DO I=1,NTIME
           PTH(I) = RPAR(1+NPAR+NINIT+2*NTIME+I)
        END DO
        ALLOCATE(PDSPL(NTIME))
        DO I=1,NTIME
           PDSPL(I) = RPAR(1+NPAR+NINIT+3*NTIME+I)
         END DO
        ALLOCATE(PTHSPL(NTIME))
        DO I=1,NTIME
           PTHSPL(I) = RPAR(1+NPAR+NINIT+4*NTIME+I)
        END DO
        ALLOCATE(R(NTIME))
        DO I=1,NTIME
           R(I) = RPAR(1+NPAR+NINIT+5*NTIME+I)
        END DO
        ALLOCATE(RSPL(NTIME))
        DO I=1,NTIME
           RSPL(I) = RPAR(1+NPAR+NINIT+6*NTIME+I)
        END DO
        ALLOCATE(DPCDT(NTIME))
        DO I=1,NTIME
           DPCDT(I) = RPAR(1+NPAR+NINIT+7*NTIME+I)
        END DO
        ALLOCATE(DPCDTSPL(NTIME))
        DO I=1,NTIME
           DPCDTSPL(I) = RPAR(1+NPAR+NINIT+8*NTIME+I)
        END DO
        ALLOCATE(DPTHDT(NTIME))
        DO I=1,NTIME
           DPTHDT(I) = RPAR(1+NPAR+NINIT+9*NTIME+I)
        END DO
        ALLOCATE(DPTHDTSPL(NTIME))
        DO I=1,NTIME
           DPTHDTSPL(I) = RPAR(1+NPAR+NINIT+10*NTIME+I)
        END DO
        
C       FIRST DELAY
        CALL LAGR5(1,X,Y,ARGLAG,PAST,ALPHA1,IPOS1,RPAR,IPAR,
     &       PHI,IPAST,NRDS)


        Y7L1=YLAGR5(7,ALPHA1,IPOS1,PHI,RPAR,IPAR,
     &       PAST,IPAST,NRDS)
        
!     --- INTERPOLATE DATA
C       Interpolate because X may not be at dt interval
        CALL splint(TIME,PD,PDSPL,NTIME,X,Pc,dt)
        CALL splint(TIME,PTH,PTHSPL,NTIME,X,Pthor,dt)
        CALL splint(TIME,R,RSPL,NTIME,X,Resp,dt)
        CALL splint(TIME,DPCDT,DPCDTSPL,NTIME,X,dPcdtV,dt)
        CALL splint(TIME,DPTHDT,DPTHDTSPL,NTIME,X,dPthdtV,dt)
        
        D     = PAR(1)
        B     = PAR(2)
        A     = PAR(3)
        KP    = PAR(4)
        Kb    = PAR(5)
        Kp    = PAR(6)
        Kr    = PAR(7)
        Ks    = PAR(8)
        tau   = PAR(9)
        taub  = PAR(10)
        taup  = PAR(11)
        taur  = PAR(12)
        taus  = PAR(13)
        tauH  = PAR(14)
        qw    = PAR(15)
        qp    = PAR(16)
        qs    = PAR(17)
        xiw   = PAR(18)
        xip   = PAR(19)
        xis   = PAR(20)
        HI    = PAR(21)
        Hp    = PAR(22)
        Hr    = PAR(23)
        Hs    = PAR(24)
        HIa   = PAR(25)
        tchar = PAR(26)
        respAmp = PAR(27)
C       Ds (Delay) is 28th parameter
        
        IF (X > (tchar + 1)) THEN
            HI = HIa
        END IF
        
C       Auxilirary Parameters        
        ewc = 1-SQRT((1+EXP(-qw*(Y(1)-xiw)))/(A+EXP(-qw*(Y(1)-xiw))))

        Pa  = Pc - Pthor
        ewa = 1-SQRT((1+EXP(-qw*(Y(2)-xiw)))/(A+EXP(-qw*(Y(2)-xiw))))
        
        dPadt = dPcdtV - dPthdtV

        nm = B*(ewc - Y(3)) + (1 - B)*(ewa - Y(4))
        
        Gp = 1/(1 + EXP(-qp*(nm - xip)))
        Gs = 1/(1 + EXP(qs*(nm - xis)))
        
        !osc = respAmp * COS(0 + respGamma*X)
        
        !Gr = (-1*Resp + 2*respMean)*respAmp/(1 + EXP(qr*(Pthor - xir)))
        !Gr = 1/(1 + EXP(qr*(Pthor - xir)))
        Gr = respAmp * Resp

        Htilde = HI*(1 - Hp*Y(5) + Hr*Y(6) + Hs*Y(7))
        
C       State Equations        
        F(1) = (-Y(1) + Pc - KP*dPcdtV)/tau
        F(2) = (-Y(2) + Pa - KP*dPadt)/tau
        F(3) = (-Y(3) + Kb*ewc)/taub
        F(4) = (-Y(4) + Kb*ewa)/taub
        F(5) = (-Y(5) + Kp*Gp)/taup
        F(6) = (-Y(6) + Kr*Gr)/taur
        !F(6) = 0.D0
        F(7) = (-Y7L1 + Ks*Gs)/taus !Only equation with delay
        F(8) = (-Y(8) + Htilde)/tauH
        
        RETURN
        END
C------------------------------------------------------------------------
        SUBROUTINE JFCN(N,X,Y,DFY,LDFY,ARGLAG,PHI,RPAR,IPAR,
     &                  PAST,IPAST,NRDS)
C ----- STANDARD JACOBIAN OF THE EQUATION
        IMPLICIT REAL*8 (A-H,K,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(LDFY,N) :: DFY
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        REAL(kind=DP), dimension(60000) :: RPAR
        REAL(kind=DP), dimension(3) :: IPAR
        EXTERNAL PHI,ARGLAGM
        
C       Dummy function because it uses internal finite-difference

        RETURN
        END
C----------------------------------------------------------------------
        SUBROUTINE JACLAG(N,X,Y,DFYL,ARGLAG,PHI,IVE,IVC,IVL,
     &                    RPAR,IPAR,PAST,IPAST,NRDS)
C ----- JACOBIAN OF DELAY TERMS IN THE EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(N) :: Y
        REAL(kind=DP), dimension(1) :: DFYL
        REAL(kind=DP), dimension(1) :: PAST
        INTEGER, dimension(1) :: IPAST
        REAL(kind=DP), dimension(60000) :: RPAR
        REAL(kind=DP), dimension(3) :: IPAR
        INTEGER, dimension(1) :: IVE,IVC,IVL
        EXTERNAL PHI
        
        IVL(1)=1 !Number of relevant delay
        IVE(1)=7 !Number of relevant equation
        IVC(1)=7 !Number of relevant component
        DFYL(1)= -1/RPAR(1+13) !Derivative of delay component 

        RETURN
        END
C-----------------------------------------------------------------------
        FUNCTION PHI(I,X,RPAR,IPAR)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER, PARAMETER :: DP=kind(1D0)
        REAL(kind=DP), dimension(60000) :: RPAR
        REAL(kind=DP), dimension(3) :: IPAR
        REAL(kind=DP), dimension(8) :: INIT
        
C       Initial value function, only called if there is a delay
C       Contains initial values for each state; up to delay for 7th eq.

        NPAR = IPAR(1)
        NINIT = IPAR(2)
        INIT = RPAR(1+NPAR+1:1+NPAR+NINIT)

        SELECT CASE (I)
        CASE (1)
            PHI=INIT(1)
        CASE (2) 
            PHI=INIT(2)
        CASE (3) 
            PHI=INIT(3)
        CASE (4)
           PHI=INIT(4)
        CASE (5)
           PHI=INIT(5)
        CASE (6)
           PHI=INIT(6)
        CASE (7)
           PHI=INIT(7)
        CASE (8)
           PHI=INIT(8)
        END SELECT
        
        RETURN
        END

