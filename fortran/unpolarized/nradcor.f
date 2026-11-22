C
C     PROGRAMME WHITNEY ADAPTE AU SUN
C
C     IN CHECKING THE PROGRAM OUT LOOK AT FUPPER,FLOWER,CFACT,
C     BBB,CCC,ETAS,ETAP,ATT,PA,WEXC,CRS,TRAD,EI,EF AND PRT CAREFULLY
C
C     Valeur trouvee dans le programme: VCOUL=-28.7 MeV
C     Only needed for Coulumb Sum Rule. Otherwise set to zero.
C     IF you use VCOUL check that Ep includes correction.
C
C     NF IS THE NUMBER OF PREVIOUSLY ANALYZED SPECTRA OR THE NUMBER OF
C     MODE SPECTRA INTRODUCED.  USING NF = NE - 1 RESTORES TO PREVIOUS
C     MODE FOR ALMOST THE SAME THING
C
C     SIMPSON INTEGRATION:
C       - NUM: MIN NUMBER OF STEPS
C       - DEL: MAX ENERGY STEP
C
C     CUTE IS THE CUT OFF IN THE RADIATION UNFOLDING
C     THE ANSWER SHOULD BE INDEPENDENT OF THE VALUE
C     OF CUTE. CHECK THIS BY VARING CUTE FROM ABOVE TO
C     BELOW DE(NE). NOTE FOR VERY SMALL VALUES OF CUTE
C     NUM MAY DEED BE INCREASED (CUTE =5,10).
C
C-----K. Slifer Temple University 1/29/03
C     Cleaned up original code to remove some obsolete options and
C     confusing GOTO labels.  Changed input format.  Modified definition
C     of CFACT, FES, and FEP to use the formalism of S. Stein et al. 
C     Phys. Rev. D 12 7 (1975).  Equation A82.  
C     The original code is apparently based on Guthrie Miller's Thesis, 
C     "Inelastic Electron Scattering at Large Angles", 
C     SLAC-R-129 (January 1971). Equation A.V.5
C
C     The two approaches agree quite well in the resonance region or
C     until about nu=Es/3.  After that there is a 10-20% disagreement
C     at JLAB typical incident energies.
C     
C     To use radcor.f to radiate a born xs instead of to unfold, set 
C     NITT to 1 and enter only one incident spectra in input file.  The 
C     result which is stored in variable CRS will be output to screen.  
C     The original code provided this option by setting IFF to 1, but it 
C     has been broken.  Now IFF should remain 0  

C------------------------------------------------------------------------
      PROGRAM RADCOR
      PARAMETER(NPTS=1500)     ! Max number of points per spectra
      PARAMETER( NE0=13  )     ! Max number of incident energies.
      CHARACTER TITLE*7

      REAL XNU2(NPTS,NE0), XFACTOR(NPTS,NE0)
      INTEGER NFINITE(NE0)


      REAL EP(NE0,  NPTS)      ! Scattered energy including collision loss
      REAL EP0(NE0, NPTS)      ! Scattered energy
      REAL ECR(NE0, NPTS)      ! Initial spectra 
      REAL EER(NE0, NPTS)      ! Statistical error on initial spectra
      REAL ESY(NE0, NPTS)      ! Systematic  error on initial spectra

      REAL PX(NE0   ,NPTS)     ! E-E'
      REAL XSI(NE0  ,NPTS)     ! initial spectra
      REAL XSF(NE0  ,NPTS)     ! final spectra
      REAL STATF(NE0,NPTS)     ! final statistical error
      REAL SYSTF(NE0,NPTS)     ! final systematic error
      REAL STATI(NE0,NPTS)    
      REAL SYSTI(NE0,NPTS)
      REAL ERRTI(NE0,NPTS)
      REAL RSYSTRC(NE0,NPTS)

      REAL YYA(NPTS)           ! for spline fit
      REAL YYB(NPTS)           ! for spline fit
      REAL YYC(NPTS)           ! for spline fit
      REAL YYD(NPTS)           ! for spline fit
      REAL PX_TMP(NPTS)        ! for spline fit
      REAL XSF_TMP(NPTS)       ! for spline fit
      REAL STATF_TMP(NPTS)     ! for spline fit
      REAL SYSTF_TMP(NPTS)     ! for spline fit
      REAL RSYSTRC_TMP(NPTS)   ! 

      REAL XNUMIN(NE0)         ! minimum energy transfer
      REAL XNUMAX(NE0)         ! maximum energy transfer

      DIMENSION  TX(NE0)       ! Target thickness in radiation lengths.
      DIMENSION TXC(NE0)       ! TXC(ITG) = 0.5*TX(ITG)
      DIMENSION TIW(NE0)       ! Thickness before scattering in radiation lengths.
      DIMENSION TFW(NE0)       ! Thickness after  scattering in radiation lengths.
      DIMENSION RC_SYS(NE0)       ! Estimate of relative error on Peaking approx.


      DIMENSION  JSTRT(NE0)
      DIMENSION RNORME(NE0)    ! Renormalization factor for experimental xs
      DIMENSION     ES(NE0)    ! Incident energies.
      DIMENSION    ES0(NE0)
      DIMENSION     DE(NE0)    ! Binsize for each energy
      DIMENSION     NP(NE0)    ! Number of points for each spectra
      DIMENSION  WMOTT(NE0)
      DIMENSION    ETA(NE0)
   
      DIMENSION DELTAS(NE0)    ! Collision loss before scattering
      DIMENSION DELTAP(NE0)    ! Collision loss after scattering

      DIMENSION   AR(NE0,NPTS)
      DIMENSION   AT(NE0,NPTS)
      DIMENSION   PA(NE0,NPTS) ! (ES(I) - EP(I,J))/ES(I)   Y
      DIMENSION ECRA(NE0,NPTS) ! Used to test convergence of iterative procedure.

      DIMENSION   EV(NPTS)     ! Used in plotting procedure
      DIMENSION  SIG(NPTS)     ! Used in plotting procedure
      DIMENSION ESIG(NPTS)     ! Used in plotting procedure
      DIMENSION WORK(NPTS)     ! Used in plotting procedure
      CHARACTER*24 TIMESTR

      COMMON ALPH,AP,AT,BSWZW,BSZ,C1,C2,C3,C4,C5,C6,C7,D1,EF,EI
      COMMON ES,FME,FMESQ,FMOTT,FMT,GAMT,NE,NP,PA
      COMMON PI,SCZ,SCZW,SINSQ,TH,TRAD,TXCO,COT,Z,Tb,Ta,SPENCE
      COMMON INTERNAL,EXTERNAL
      COMMON/REVISE/QSQ
      LOGICAL  PLOTIT,THEORY,RADIATE
      LOGICAL INTERNAL,EXTERNAL,COLLISIONAL_LOSS
      LOGICAL FINITE_ACC
      EXTERNAL FES,FEP
 
      !These are the E94010 collisional losses before (S) and after (P) scattering.
      !Set to zero if unsure.
      DATA DELTAS/0.4,0.5,0.5,0.5,0.5,0.5,7*0.0/
      DATA DELTAP/2.4,3.3,2.2,3.3,2.5,2.2,7*0.0/
C----------------------------------------------------------------
C     Some Functions

      RL(XS,XZ,XA,XSC)=1./(24.092/XA*ALPH*2.818**2*XZ*(XZ+
     . XS)*XSC*1.E-3) ! Not used anymore

      SE(XZ,XSC)=ALOG(1194./(XZ**(2./3.)))/XSC

      B(XS,XZ,XSC)=4./3.+(XZ+1.)/((XZ+XS)*9.*XSC)

      SC(XZ)=ALOG(184.15/XZ**(1./3.))-1.202*(ALPH*XZ)**2
C----------------------------------------------------------------
      PI    = 3.1415926536 
      ALPH  = 1./137.035989
      FME   = 0.510999     
      AMASS = 931.494      
      HBARC = 197.32705    

      AP    = ALPH/PI
      FMESQ = FME*FME
C--------------------- 
      PLOTIT = .FALSE. ! Outputs graphical result.
      THEORY = .FALSE. ! PERMITS READING IN OF MODEL GENERATED LOWER ENERGY SPECTRA

      INTERNAL         =.TRUE.  ! Turn on/off internal R.C.
      EXTERNAL         =.TRUE.  ! Turn on/off external R.C.

                                ! These two only work for E94010
      COLLISIONAL_LOSS =.TRUE. ! Turn on/off collision loss. Values specific to E94010
      FINITE_ACC       =.TRUE. ! Added this option 06/25/04. Values spicific to E94010

      RADIATE          =.FALSE. ! Use the code in reverse, that is, to radiate a born xs 
                                ! instead of unfolding a born xs from measured data.



      VCOUL  = 0.0              ! Old CSR experiment stuff. Leave zero
      NF     = 0
      NUM    = 80
      DEL    = 1.               ! Step size in the simpson integration.
      CUTE   = 10.              ! "Delta" parameter from Stein.
C----------------------------------------------------------------

      OPEN (1,FILE='ele-input',status='old')
      OPEN (99,FILE='output/E94010_xs_unfolded.dat',STATUS='UNKNOWN')

      !CALL FDATE(TIMESTR)             ! This seems to work only on solaris
      !WRITE(99,'(A,A)') '# ',TIMESTR

      WRITE(99,'(A)') "# Energy in MeV.  cross sections in pb/MeV-sr.    "
      WRITE(99,'(A)') "# XSi : before RC"
      WRITE(99,'(A)') "# XSf : after  RC"
    
      !---------------------------------------------------------------------------------
      IF(FINITE_ACC) THEN
        WRITE(99,'(A)') "# XSfc: after  RC, with finite acceptance correction"
        WRITE(99,'(A)') "#"
        WRITE(99,'(A,A5,x,A6,5A12)') "#","E0","E-E'","XSi", "XSf","XSfc", "STAT", "SYST"
      ELSE
        WRITE(99,'(A)') "#"
        WRITE(99,'(A,A5,x,A6,4A12)') "#","E0","E-E'","XSi", "XSf", "STAT", "SYST"
      ENDIF

      WRITE(99,'(A,A72)')'#','------------------------------------------------------------------------'

      DO I = 1,6
        CALL GET_FINITE_COR(I,XNU2,XFACTOR,NFINITE)
      ENDDO

10    CONTINUE

      DO I=1,NE0
         NP(I)    = 0
         DE(I)    = 0.
         ES(I)    = 0.
         ES0(I)   = 0.
         ETA(I)   = 0.
         WMOTT(I) = 0.
      ENDDO
      CALL ZERO(AR,    NPTS*NE0)
      CALL ZERO(EER,   NPTS*NE0) 
      CALL ZERO(ESY,   NPTS*NE0)
      CALL ZERO(EP,    NPTS*NE0)
      CALL ZERO(AT,    NPTS*NE0)
      CALL ZERO(PA,    NPTS*NE0)
      CALL ZERO(ECR,   NPTS*NE0)
      CALL ZERO(ECRA,  NPTS*NE0)
     
      IFF = 0  ! Leave at zero
      EL  = 0.
      DDEE= 0.

C-----------------------------------------------------------------------
C
C     NITT         : NOMBRE D'ITERATIONS
C     NE           : Number of spectra to analyze.
C     Z,ZW         : Charge of target & windows.
C     A,AW         : Mass # of target & windows.
C     THETA        : Scattering angle in degrees.
C     TTHETA       : Target angle in degrees.
C
C                    All thicknesses should be in radiation lengths.
C                    -----------------------------------------------
C     TX           : Total target thickness.  Set this to zero.  
C                    TIF and TIW now take target thickness into account.
C     TIW          : Total thickness before scattering including 1/2 target thickness.
C     TFW          : Total thickness after  scattering including target contribution.
C
C     SYS          : Relative error on unfolding procedure.
C                   
C     RNORME       : Renormalization factor. Can leave at 1.0
C
C                    JUST GIVE ZW AND AW REALISTIC NUMBERS
C------------------------------------------------------------------------

      READ (1,*) NITT             
      READ (1,*) NE                
      READ (1,*) Z, A, THETA       
      READ (1,*) ZW,AW,TTHETA      
      DO K=1,NE
        READ(1,*) TX(K),TIW(K),TFW(K),RC_SYS(K),RNORME(K)
      ENDDO
      
      IF(IBIN .EQ. 0) IBIN = 1
     
      IBIN = 1     

      FMT  = A*AMASS
      TH   = THETA*PI/180.
      THT  = TTHETA*PI/180.
      SINSQ= SIN(TH/2.)**2
      COT  = COS(TH)

444   CONTINUE

      SCZ   = SC(Z)
      SCZW  = SC(ZW)
      S     = SE(Z,SCZ)
      SW    = SE(ZW,SCZW)
C     XO    = RL(S,Z,A,SCZ)      ! Not used anymore
C     XW    = RL(SW,ZW,AW,SCZW)  ! Not used anymore
      BSZ   = B(S,Z,SCZ)
      BSWZW = B(SW,ZW,SCZW)

      FMOTT = (HBARC*ALPH*COS(TH/2.)/2./SINSQ)**2
      FMOTT = FMOTT*1.E7                                ! MeV^2-nb

      DO ITG = 1, NE
        TXC(ITG) = 0.5*TX(ITG)
      ENDDO

      IF(IFF.EQ.1) GO TO 427

C---------------------------------------------------------------
      DO I=1,NE

        READ (1,'(A7)') TITLE

C
C
C ==>	Introduction du potentiel coulombien dans la cinematique
C
C	- EMX IS THE MAXIMUM FINAL ENERGY
C	- NP(I) ARE THE NUMBER OF SECONDARY POINTS FOR EACH ENERGY
C	- DE(I)IS THE ENERGY SPACING BETWEEN DATA POINTS FOR EACH ENERGY
C	- ES(I) ARE THE NE NUMBER OF INCIDENT ENERGIES IN MEV
C	  VARING FROM LOW TO HIGH
C

        READ (1,*)  ES0(I),EMX,DE(I),NP(I)
        IF (COLLISIONAL_LOSS) THEN
          ES(I)=ES0(I)-DELTAS(I)
        ELSE
          ES(I)=ES0(I)
        ENDIF

C       This is if you want individual output files for each spectra
        !-----------------------------------------------------------
C       OPEN(10+I,FILE=TITLE//'_radcor.out',STATUS='unknown')
C       OPEN(30+I,FILE=TITLE//'f.dat',STATUS='unknown')
C       WRITE(10+I,'(A,4F7.1)') "#",ES(I),EMX,DE(I),NP(I)
C       WRITE(30+I,'(A,4F7.1)') "#",ES(I),EMX,DE(I),NP(I)

        ES(I)=ES(I)-VCOUL
        EMX=EMX-VCOUL



        WMOTT(I) = FMOTT/(ES(I)+VCOUL)**2
        ETA(I)   = 1.+2.*ES(I)*SINSQ/FMT
        NPI      = NP(I)
        IF(NPI .GT. (NPTS-1)) THEN
          WRITE(6,'(A)') 'TOO MANY DATA POINTS. INCREASE PARAMETER NPTS'
          STOP
        ENDIF

C
C	- EP(I,J) ARE THE SECONDARY ENERGIES IN MEV
C	- ECR(I,J) ARE THE CROSS SECTIONS AT THESE ENERGIES
C	- EER(I,J) ARE THE UNCERTAINTIES IN THE VALUES -- AT THE MOMENT
C	  NOTHING IS DONE WITH THEM 
C       
CSL     - ESY(I,J) Systematic error on initial spectra. Can be zero if unknown.
C

        DO J=1,NPI
          READ (1,*) XNU,ECR(I,J),EER(I,J),ESY(I,J)
          IF(J.EQ.1)   XNUMIN(I)=XNU
          IF(J.EQ.NPI) XNUMAX(I)=XNU
          XNU=XNU+VCOUL
          EP0(I,J)=ES0(I)-XNU
          IF (COLLISIONAL_LOSS) THEN
              EP(I,J)=EP0(I,J)+DELTAP(I)
          ELSE
              EP(I,J)=EP0(I,J)
          ENDIF
        ENDDO

        DO J=1,NPI
           ECR(I,J) = ECR(I,J)*RNORME(I)
           EER(I,J) = EER(I,J)*RNORME(I)
C           EP(I,J) = EMX-(J-1)*DE(I)

           !This tries to fill in any regions where spectra are not continuous 
           !Probably not the safest way to solve this problem if you have it.
           IF( (J.GE.2).AND.(J.LE.(NPI-1)) ) THEN  ! Fill in the holes.
              IF(ECR(I,J) .LE. 0.) THEN
                 J1 = J - 1
                 J2 = J + 1
                 IF( (ECR(I,J1).NE.0.) .AND. (ECR(I,J2).NE.0.) ) THEN
                    ECR(I,J) =(ECR(I,J1) + ECR(I,J2)    )/2.
                    EER(I,J) = EER(I ,J1) + EER(I, J2)
                 ENDIF
               ENDIF
           ENDIF

C          PURPOSE OF THIS IS TO GIVE BETTER INTERPOLATION TRAJECTORIES USING
C          ALL AVAILABLE DATA AT ALL EXCITATION ENERGIES
C          PA(I,J) ARE THE EXCITATION ENERGIES

           AR(I,J) = ECR(I,J)/WMOTT(I)          ! Holds the initial spectra
           AT(I,J) = AR(I,J)                    ! Serves as model for R.C.
           PA(I,J) = (ES(I) - EP(I,J))/ES(I)
           IF(THEORY .AND. (I .LE.NF))   AT(I,J) = AR(I,J)
        ENDDO ! J

      ENDDO ! I=1,NE
C------------------------------------------------------------------------


      NPI=NP(1)

      DO J=1,NPI
        AT(1,J)=AR(1,J)
      ENDDO
CSL--------------------------------------------    
      THR   = THETA*PI/180.
      ARG   = COS(THR/2.)**2

      SPENCE= PI**2/6.-LOG(ARG)*LOG(1.-ARG)
      DO  NSP=1,50
        SPENCE = SPENCE-ARG**NSP/FLOAT(NSP)**2
      ENDDO
CSL-------------------------------------------


427   CONTINUE
      DO 98  NUMBER = 1, NITT  ! ---MAIN  LOOP---------------------------
        DO 96 I=1,NE           ! ---INNER LOOP---------------------------

        IF(IFF.EQ.1) GO TO 428

        TXCO = TXC(I)
        C5   = 1.-ALPH/PI*28./9.
        C6   = 13.*ALPH/PI/6.

        C3   = TIW(I)
        C4   = TFW(I)
        C1   = BSZ*TXCO + BSWZW*C3
        C2   = BSZ*TXCO + BSWZW*C4
        GAMT = GAMMA(1. + C1 + C2)

        IF(THEORY .AND. (I .LE.NF))  GO TO 96

        EI   = ES(I)
        ETAS = ETA(I)

        WRITE (6,*)EI,NUMBER

        NPI  = NP(I)
        JSZ  = 0                    !JSTRT(I)     
        EJSZ = EI/ETAS  - 3.

428     CONTINUE

CSL--------------------LOOP OVER DATA POINTS-------------------------------------
        DO J=1,NPI
          IF(IFF.EQ.1)GO TO 430
          EF = EP(I,J)
430       IF(IFF.NE.1) GO TO 33
          !I  = NE
          EF = EL-DDEE
429       EF = EF+DDEE
          !IF(EF.GT.EU) GO TO 400

33        CONTINUE
          ETAP = 1./(1.-2.*EF*SINSQ/FMT)
          C7   = SQRT(ETAS*ETAP)
          QSQ  = 4.*EI*EF*SINSQ
          TRAD = 2.*ALPH/PI*(ALOG(QSQ/FMESQ)-1.)
          D1   = SQRT(EI*EF)
CSL          
          XR=ETAS*ETAP

          BBB    = 0.
          CCC    = 0.
          FUPPER = EI-ETAS*ETAP*CUTE
          FLOWER = EF*ETAP


          IF(FUPPER .LE. 0.) THEN
            WRITE (6,'(A,G10.4)') ' FUPPER=',FUPPER
            STOP
          ENDIF

          IF(FLOWER .LE. 0.) THEN 
             WRITE (6,'(A,G10.4)') ' FLOWER = ',FLOWER
             STOP
          ENDIF
          IF(FLOWER-FUPPER.GE.0.0) THEN
            !WRITE(6,*) "DOES THIS EVER HAPPEN?",EI,EF
            FUPPER=EI/ETAS
          ELSE
            CALL SIM(FLOWER,FUPPER,DEL,NUM,FES,BBB)
            FUPPER=EI/ETAS
          ENDIF


C
C       INSERT TO KEEP PROGRAM FROM DOUBLE SUBRTACTING THE ELASTIC RAD TAIL
C

        !EJSZ = EI/ETAS - 3.
        !IF(FUPPER .GT. EJSZ) FUPPER = EJSZ   ! This is an old option that was 
                                              ! mistakenly left on for E94010 analysis.
                                              ! The effect is tiny but these lines must
                                              ! be uncommented to reproduce "official"
                                              ! E94010 results.

          FLOWER = EF+CUTE

          IF(FLOWER-FUPPER.GE.0.0) THEN
            CONTINUE
          ELSE
            CALL SIM(FLOWER,FUPPER,DEL,NUM,FEP,CCC)
          ENDIF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CSL ----The original CFACT          
CSL       CFACT=(C7*CUTE/D1)**(TRAD+C1+C2)*(C5+C6*ALOG(QSQ/FMESQ))/GAMT
CSL ----Exact from Miller thesis.
CSL       CFACT=(XR*CUTE/EI)**(TRAD/2.+C1) * (CUTE/EF)**(TRAD/2.+C2)
CSL       CFACT=CFACT* (C5+C6*ALOG(QSQ/FMESQ))/GAMT   ! no different
CSL ----Changed to this new formalism from Stein.
CSL     Equation numbers are from the Stein paper
C----------------------------------------------------------------------
          DEL   = CUTE
          xb    = 4./3.  ! A45

          IF (EXTERNAL) THEN
            Tb    = TXC(I) + TIW(I)
            Ta    = TXC(I) + TFW(I) 
          ELSE
            Tb    = 0.0
            Ta    = 0.0
          ENDIF

          XM      = 931.49 ! Mass of the nucleon
          XMT     = FMT    ! Mass of the target
          EMASS   = FME

          QMS     = 4.*EI*EF*SIN(THR/2.)**2
          

          IF (INTERNAL) THEN
            TR    = (2.*ALPH/PI)*(LOG(QMS/EMASS**2)-1.) 
            TR    = TR/2./xb                            ! = trad (A57)
           
            D2 = 13.*(LOG(QMS/EMASS**2)-1.)/12.-17./36. ! this term dominates D2
            D2 = D2 +0.5*(PI**2/6.-SPENCE)
            D2 = D2 -1./4.*( LOG( EI/(EF) ) )**2        ! Correct. to peak. appr.
            D2 = D2*(2.*ALPH/PI)
            D2 = D2+0.5772*xb*(Tb+Ta)                   ! Here D2= F-1
            xF = (1.+D2)
          ELSE
            TR = 0.0
            xF = 1.0+0.5772*xb*(Tb+Ta)                  ! (A44)
          ENDIF

          Tpb = tr + Tb
          Tpa = tr + Ta

          R    = ( XMT+EI*(1-COS(THR)) )/( XMT-(EF)*(1-COS(THR)) ) ! (A83)
          xeta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )    ! (A46)
          xi   = (PI*EMASS/2./ALPH)*(Ta+Tb)
          xi   = xi/( (Z+xeta)*LOG(183.*Z**(-1./3.)) )             ! (A52)

          CFACT = xF
          CFACT = CFACT*( (R*DEL/EI)**(xb*Tpb))*((DEL/EF)**(xb*Tpa) )
          CFACT = CFACT*(1. - xi/DEL/( 1.-xb*(Tpb+Tpa)) )  ! not much effect
CSL-----     


          WEXC = (EI - EF)/EI
          ATT  = TERP(I,WEXC)
          CRS  = CFACT*WMOTT(I)*ATT+BBB+CCC          

          IF(IFF.EQ.1) THEN
            PRT  = ATT*WMOTT(I)
            GO TO 429
          ENDIF
          IF (RADIATE.AND.I.LT.NE) THEN             ! If radiating, don't modify
            CONTINUE                                ! the lower energy cross section model.
          ELSE
            AT(I,J)= (AR(I,J)-(BBB+CCC)/WMOTT(I))/CFACT  ! Update the model
          ENDIF

          IF(J .LT. JSZ) AT(I, J) = AR(I, J)        ! HERE

          WR1 = WMOTT(I)*AR(I,J)                       ! Initial spectra
          WR2 = WMOTT(I)*AT(I,J)                       ! Newest model
C
C
C ==>	  POUR LES TESTS DE CALAGE EN ENERGIE:
C
C	  - WR2 = ECR(I,J)
C	  - EER(I,J) NON MODIFIES
C
C         MO AND TSAI ESTIMATE ABOUT 4%ERROR IN CORRECTION FROM NEGLECTED
C         DIAGRAMS.  1% ERROR IN NUMERICAL INTEGRAL.  TWO ADDED ALGEBRAICALLY
C	  SINCE THEY ARE NOT STATISTICALLY INDEPENDENT
C
C
          ECR(I,J) = WR2
          IF(ECR(I,J) .LT. 0.) THEN 
             ECR(I,J) = 0.
          ENDIF
          DC3      = .05*.05*(BBB*BBB + CCC*CCC)
C         EER(I,J) = SQRT(EER(I,J)**2 + DC3)
          REER     = EER(I,J) / WMOTT(I)
          PX(I,J)  = ES(I)-EP(I,J)

          IF(NUMBER.EQ.(NITT-1)) THEN
            ECRA(I,J) = WR2           ! to study convergence.
          ELSEIF(NUMBER.EQ.NITT) THEN

            IF( RADIATE .AND. I.EQ.NE ) THEN      ! Output a radiated cross section.
              XNUR = PX(I,J)
              XWR  = GET_W(EI,THETA,XNUR)
              WRITE (6,'(3F7.1,2E12.4)') EI,XNUR,XWR,WMOTT(I)*AR(I,J),CRS
            ENDIF

            XSI(I,J) = WR1                  ! Initial XS
            XSF(I,J) = WR2                  ! Final XS
            SCALE    = XSF(I,J)/XSI(I,J)    ! To adjust statistical error
                                            ! 
            STATI(I,J) = EER(I,J)           ! Initial statistical error 
            SYSTI(I,J) = ESY(I,J)           ! Initial systematic  error 
            ERRTI(I,J) = SQRT( (STATI(I,J))**2 + (SYSTI(I,J))**2 ) ! Total init. error

            STATF(I,J) = STATI(I,J) * SCALE                        ! Final stat. error 
            SYSTRC     = ABS( XSI(I,J)-XSF(I,J) )*RC_SYS(I)        ! R.C. systematic. 
                                                                   ! 
            SYSTF(I,J) = ( SYSTI(I,J) * SCALE )**2 + (SYSTRC)**2   ! 
            SYSTF(I,J) = SQRT ( SYSTF(I,J) )                       ! Final syst. error
            RSYSTRC(I,J) = SYSTRC/XSF(I,J)     *100.0              ! 

          ENDIF
      ENDDO ! --------J, Loop over data points----------------------------- 


      IS  = I
      NPI = NP(IS)

      IF(IFF.EQ.1) GO TO 98

      IF(.NOT.PLOTIT) GO TO 688

      CALL GETROW(EP  , EV  , IS  , NPI)
      CALL GETROW(ECR , SIG , IS  , NPI)
      CALL GETROW(EER , ESIG, IS  , NPI)
      CALL ADDVEC(ESIG, SIG , WORK, NPI)


      NYS = 10.1/DE(IS)


      CALL YSCALE(WORK, NPI, 100, DVY, FVY, NYS)
      CALL PRPLT(EV, SIG, ESIG, NPI, DVY, FVY)

688   CONTINUE
96    CONTINUE  !----------INNER LOOP------------------------------------
98    CONTINUE  !----------MAIN  LOOP------------------------------------



      IF (RADIATE) STOP  ! No need to go through with unfolding procedure




CSL----This option has been broken.  Leave IFF=0.
C
C	- IFF = 1: RADIATE THE DATA FOR THE GIVEN TARGET AND
C		   WINDOW THICKNESSES. DO THIS FOR THE MAX
C		   INCIDENT ENERGY
C	- IFF = 2: START PROGRAM OVER WITH NEW DATA (GOTO 10)
C	- IFF = 0: WHICH CAN BE BLANK CARD, ... STOP
C
C

      EL   = EP(IS, NPI-1)
      EU   = EP(IS, 3)
      DDEE = DE(IS)
      
      IF    (IFF.EQ.0) THEN
        GO TO 80
      ELSEIF(IFF.EQ.1) THEN
        GO TO 444
      ELSEIF(IFF.EQ.2) THEN
        GO TO 10
      ELSE
        GO TO 80
      ENDIF

72    IE=1
      GO TO 78
78    CONTINUE


80    CONTINUE
      DO I=1,NE
      
        YP1=0
        YPN=0
        DO J=1,NP(I)
          PX_TMP(J)      =   PX(I,J)
          XSF_TMP(J)     =  XSF(I,J)
          STATF_TMP(J)   = STATF(I,J)
          SYSTF_TMP(J)   = SYSTF(I,J)
          RSYSTRC_TMP(J) = RSYSTRC(I,J)
        ENDDO
        
        ! Subroutine "SPLINE" must be called before "SPLINT" can be used.
        ! This is all to shift the spectra back to the original x-bins
        ! from before the collisional loss.

        CALL SPLINE(PX_TMP,  XSF_TMP   ,NP(I),YP1,YPN,YYA )
        CALL SPLINE(PX_TMP,STATF_TMP   ,NP(I),YP1,YPN,YYB )
        CALL SPLINE(PX_TMP,SYSTF_TMP   ,NP(I),YP1,YPN,YYC )
        CALL SPLINE(PX_TMP,RSYSTRC_TMP ,NP(I),YP1,YPN,YYD )

        DO J=1,NP(I)-1
          IF (COLLISIONAL_LOSS) THEN
            X = PX(I,J) + ( DELTAS(I) + DELTAP(I) )  ! original nu value
          ELSE
            X = PX(I,J)
          ENDIF
          CALL SPLINT(PX_TMP,  XSF_TMP   ,YYA ,NP(I),X,Y   ) ! shift final xs
          CALL SPLINT(PX_TMP,STATF_TMP   ,YYB ,NP(I),X,YE  ) ! shift final stat error
          CALL SPLINT(PX_TMP,SYSTF_TMP   ,YYC ,NP(I),X,YSE ) ! shift final syst error
          CALL SPLINT(PX_TMP,RSYSTRC_TMP ,YYD ,NP(I),X,RRCSE)! shift final RCsyst error

          !WRITE (30+I,'(F6.1,4E12.4)') X,XSI(I,J),Y,YE,YSE

          ! Apply correction for finite Acceptance effect.
          ! These values are specific to E94010 and shouldn't be used otherwize.

          IF (FINITE_ACC) THEN 
            CALL INTERP(NFINITE(I),XNU2(1,I),XFACTOR(1,I),DUMMY,X,FIN_COR,DM,15,-1,IERR)
            WRITE (99,'(F6.1,x,F6.1,5E12.4)') ES0(I),X,XSI(I,J),Y,Y*FIN_COR,YE,YSE
          ELSE
            WRITE (99,'(F6.1,x,F6.1,4E12.4)') ES0(I),X,XSI(I,J),Y,YE,YSE
          ENDIF

        ENDDO
      ENDDO
C----------------------------------------------------------------

      STOP
      END


C----------------------------------------------------------------
C----------------------------------------------------------------
C
      SUBROUTINE GETROW(A, V, INDX, NJ)
      PARAMETER(NPTS=1500,NE0=13)
      DIMENSION A(NE0,NPTS), V(NPTS)

      DO 1  J = 1, NJ
         V(J) = A(INDX, J)
1     CONTINUE

2     FORMAT(/, 10G12.4)
3     FORMAT(///)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ADDVEC(X, Y, Z, N)
      PARAMETER(NPTS=1500,NE0=13) 

      DIMENSION X(NPTS), Y(NPTS), Z(NPTS)

         DO 100  I = 1, N
100   Z(I) = X(I) + Y(I)

2     FORMAT(/,10G12.4      )
3     FORMAT(/////)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE YSCALE(V, NV, NPLT, DVY, FVY, ISTART)
      PARAMETER(NPTS=1500,NE0=13) 

      DIMENSION V(NPTS)

      PLT = NPLT
      YMIN = 0.0
      YMAX =-100.
      FVY = 0.

      DO 100 I = ISTART, NV
      IF(V(I) .GE. YMAX) YMAX = V(I)
100   IF(V(I) .LE. YMIN) YMIN = V(I)

      IF(YMIN .GE. 0.) YMIN = 0.
      IF(YMIN .LT. 0.) FVY = YMIN
      
      DVY = (YMAX - YMIN)/PLT

      WRITE (2,911) DVY, YMIN, YMAX
911   FORMAT(/,5X,'DVY = ',G14.4,' YMIN = ',G14.4,' YMAX = ',G14.4)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRPLT(X, Y, D, N, DVY, FVY)

      DIMENSION X(1), Y(1), D(1)
      DIMENSION ALINE(100), EDGE(11)
      DATA STAR/'*'/, BL/' '/, BAR/'-'/

      N20 = N - 20

      DO 100  I = 1, 11
100   EDGE(I) = FVY + (I - 1)*10.*DVY

      !WRITE (2,400) EDGE
400   FORMAT(//, 30X, '  UNFOLDED CROSS-SECTS    '///, 2X,
     +  11E10.3, 1X// 11X, 20('+----'))
       N20 = N

      DO 200  I = 1, N
      	DO 600  J = 1, 100
600   	ALINE(J) = BL
		K = (Y(I) - FVY)/DVY
		K = K + 1
		KD = (D(I)/DVY) + 1
		KM = K - KD
		KP = K + KD

		IF(K .GT. 100) K = 100
		IF(KP .GT. 100) KP = 100
		IF(KM .GT. 100) KM = 100
		IF(KM.LE.0)KM=1

			DO 700  JJ = KM, KP
700			ALINE(JJ) = BAR
			ALINE(K) = STAR

			WRITE (2,500) X(I), ALINE
500			FORMAT(1X, F8.2, ' I', 100A1)

200   CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION TERP(I,WEXC)
      PARAMETER(NPTS=1500,NE0=13)
      DIMENSION ES(NE0),NP(NE0)
      DIMENSION AT(NE0, NPTS), PA(NE0, NPTS)
      COMMON ALPH,AP,AT,BSWZW,BSZ,C1,C2,C3,C4,C5,C6,C7,D1,EF,EI
      COMMON ES,FME,FMESQ,FMOTT,FMT,GAMT,NE,NP,PA
      COMMON PI,SCZ,SCZW,SINSQ,TH,TRAD,TXCO,COT,Z,Tb,Ta,SPENCE

      TERP = 0.0
      NAR=NP(I)

      IF(WEXC.GE.PA(I,1)) GO TO 10

      TERP=0.

      RETURN

10    IF(WEXC.LT.PA(I,NAR)) GO TO 20

      TERP=AT(I,NAR)

      RETURN

20    J=NAR/2
      MAX=NAR
      MIN=0

21    IF(WEXC.LT.PA(I,J)) GO TO 30

      IF(WEXC.LE.PA(I,J+1)) GO TO 33

      MIN=1
      J=J+(MAX-J+1)/2
      GO TO 21

30    MAX=J
      J=J-(J-MIN+1)/2
      GO TO 21

33    IF(J.LE. 2) GO TO 50

40    X1 = PA(I, J - 1)
      Y1 = AT(I, J - 1)
      X2 = PA(I,J)
      Y2 = AT(I,J)
      X3 = PA(I, J+1)
      Y3 = AT(I, J + 1)

C
C ==>	PARABOLIC FIT FOR  X1<X2<X3
C

      XX=WEXC-X2
      XP=X3-X2
      XM=X1-X2
      AA=Y2
      CC=(Y3-Y1)*(XP+XM)/(XP-XM)-(Y3+Y1-2.*Y2)
      CC=CC/XP/XM/2.
      BB=(Y3-Y1)/(XP-XM)-CC*(XP+XM)
      TERP=AA+BB*XX+CC*XX**2

      RETURN

50    TERP=(AT(I,J)*(PA(I,J+1)-WEXC)+AT(I,J+1)*(WEXC-PA(I,J)))/
     . (PA(I,J+1)-PA(I,J) + 0.00001)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ZERO(X, N)

      DIMENSION X(N)

      DO 1 I = 1, N
1     X(I) = 0.

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TTERP(X, YT, AX, AY, N)
      PARAMETER(NPTS=1500,NE0=13) 

      DIMENSION AX(NPTS), AY(NPTS)
C
C       - X = X VALUE
C       - YT = INTERPOLATED Y VALUE
C       - AX, AY ARE X,Y ARRAYS WITH SMALLEST RPT, SMALLEST X FIRST.
C       - N = NO OF CHANNELS
C
C


C
C ==>	WE TEST IF X IS OUT OF RANGE, AND IF SO MAKE LINEAR INTERP.
C

      IF(X .LE. AX(N-1)  .OR. X.EQ. AX(N)) GO TO 5

      YT = (AY(N)-AY(N-1))*(X- AX(N-1))/(AX(N) - AX(N-1)) + AY(N-1)

      RETURN

5     CONTINUE

      IF(X .GE. AX(1) ) GO TO 4

      YT = (AY(2) - AY(1))*(X - AX(1))/(AX(2) - AX(1)) + AY(1)

      RETURN

4     CONTINUE

      NN = N - 1

      DO 1  I = 2, NN
         IF(X .LT. AX(I) ) GO TO 2
         IF(X .EQ. AX(I)) GO TO 3
1     CONTINUE

3     YT = AY(I)

      RETURN

2     X1 = AX(I - 1)

C
C ==>	X LIES BETWEEN X2 AND X3
C

      X2 = AX(I)
      X3 = AX(I + 1)
      Y1 = AY(I - 1)
      Y2 = AY(I)
      Y3 = AY(I + 1)
      X0=((Y1-Y2)*(X3*X3-X2*X2)+(X1*X1-X2*X2)*(Y2-Y3))/
     .(2.*((X2-X1)*(Y3-Y2)+(X3-X2)*(Y1-Y2)) + .0001)
      AA=(Y1-Y2)/(X1*X1-X2*X2+2.*X0*(X2-X1) + .00001)
      Y0=Y1-AA*(X1-X0)**2
      YT=(X- X0)**2*AA + Y0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION FTCS(EO,EB)
      PARAMETER(NPTS=1500,NE0=13)

      DIMENSION ES(NE0),NP(NE0)
      DIMENSION AT(NE0, NPTS), PA(NE0, NPTS)
      COMMON ALPH,AP,AT,BSWZW,BSZ,C1,C2,C3,C4,C5,C6,C7,D1,EF,EI
      COMMON ES,FME,FMESQ,FMOTT,FMT,GAMT,NE,NP,PA
      COMMON  PI,SCZ,SCZW,SINSQ,TH,TRAD,TXCO,COT,Z,Tb,Ta,SPENCE

      EO=EO
      EB=EB

      IF(EB.LT.EO/(1.+2.*EO*SINSQ/FMT)) GO TO 10

      FTCS=0.

      RETURN

10    CONTINUE

      EO=EO
      EB=EB
      WEXC = (EO - EB)/ EO

      IF(EO.GT.ES(1)) GO TO 20

      FTCS=TERP(1,WEXC)*FMOTT/EO/EO

      RETURN

20    IF(ES(NE).GE.EO) GO TO 21

      WRITE (2,96)
96    FORMAT(1X,22HES(NE) IS LESS THAN EO)

      STOP

21    I=1

22    IF(ES(I+1)-EO) 23,24,24

23    I=I+1

      GO TO 22

24    CONTINUE

      F=TERP(I,WEXC)
      G=TERP(I+1,WEXC)
      H=(F*(ES(I+1)-EO)+G*(EO-ES(I)))/(ES(I+1)-ES(I))
      FTCS=H*FMOTT/EO/EO

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION V(EO,EV)
      PARAMETER(NPTS=1500,NE0=13)

      DIMENSION ES(NE0),NP(NE0)
      DIMENSION AT(NE0, NPTS), PA(NE0, NPTS)
      COMMON ALPH,AP,AT,BSWZW,BSZ,C1,C2,C3,C4,C5,C6,C7,D1,EF,EI
      COMMON ES,FME,FMESQ,FMOTT,FMT,GAMT,NE,NP,PA
      COMMON PI,SCZ,SCZW,SINSQ,TH,TRAD,TXCO,COT,Z,Tb,Ta,SPENCE
      COMMON/REVISE/QSQ
C
C       VARIOUS FORMS OF V MAY BE USED
C       SEE MO AND TSAI OR MILLER:
C
C       - V WITH .2 AND 0.7 IS FROM MILLER
C       - V WITH .75 FROM MO-TSAI
C


      X=EV/EO
      V = TRAD/2.*( X - 0.2*(1. - X) + 0.7*(1. - X)**2)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION FMIN(X,Y)

      IF(X-Y) 1,2,2

1     FMIN=X

      RETURN

2     FMIN=Y

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION WM(EO,E,B)

      WM=(1.+(E/EO)**2-(E/EO)*(2.-B))/(EO-E)

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SIM(FL,FU,DEL,N,F,ANS)

      NX = ((FU - FL)/DEL/2.)
      NX = 2*NX

      IF(NX.LT.N) NX=(N/2)*2

      DT=(FU-FL)/FLOAT(NX)
      SI = F(FL)
      SI = SI + F(FU)
      XX = FL + DT
      SI = SI + F(XX) *4.
      AA=FL
      II=0
5     II=II+2

      IF (II.GE.NX)  GO TO 7

      AA =FLOAT (II)*DT
      AA=AA+FL
      BB=AA+DT
      SI = SI + 2.*F(AA) + 4.*F(BB)
      GO TO 5

7      ANS = DT * SI/3.0
8      FORMAT (' FL,FU,DEL,DT,N,NX=', 6(G11.4,1X), ' ANS ', G11.4)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION GAMMA(X)

       Y = X-1.

      GAMMA=1.-.57486*Y+.95123*Y*Y-.69986*Y**3
     . +.42455*Y**4-.10107*Y**5

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION FES(VAR)
      PARAMETER(NPTS=1500,NE0=13)

      DIMENSION ES(NE0),NP(NE0)
      DIMENSION AT(NE0, NPTS), PA(NE0, NPTS)
      COMMON ALPH,AP,AT,BSWZW,BSZ,C1,C2,C3,C4,C5,C6,C7,D1,EF,EI
      COMMON ES,FME,FMESQ,FMOTT,FMT,GAMT,NE,NP,PA
      COMMON PI,SCZ,SCZW,SINSQ,TH,TRAD,TXCO,COT,Z,Tb,Ta,SPENCE
      COMMON INTERNAL,EXTERNAL
      LOGICAL INTERNAL,EXTERNAL
CSL      ORIGINAL CODE
CSL      X=D1/3.
CSL      Y=(EI-VAR)/C7
CSL      FK=FMIN(X,Y)
CSL      EO=EI
CSL      EB=EF
CSL      FES=FTCS(VAR,EB)*(V(EO,VAR)/(EO-VAR)+C3*WM(EO,VAR,
CSL     . BSWZW)+TXCO*WM(EO,VAR,BSZ))*(FK/D1)**(TRAD+C1+C2)*
CSL     . (C5+C6*ALOG(4.*EB*VAR*SINSQ/FMESQ))/GAMT

CSL       MILLER FORM
CSL       EO=EI
CSL       EB=EF
CSL       X1 = EI/3.
CSL       X2 = EI-VAR
CSL       K1 = FMIN(X1,X2)
CSL
CSLCC      EPMAX = VAR/(1.0 + 2.0*VAR/FMT*SINSQ)
CSL       EPMAX = EI/(1.0 + 2.0*EI/FMT*SINSQ)
CSL       R     = C7**2
CSL
CSL       Y1 = EF/3.
CSLCC       Y2 = (EPMAX-EF)/R
CSL       Y2 = (EI-VAR)/R  ! stein
CSL       K1p= FMIN(Y1,Y2)
CSL
CSL       FES = FTCS(VAR,EB) * (V(EO,VAR)/(EO-VAR) +
CSL     .       C3*WM(EO,VAR,BSWZW) + TXCO * WM(EO,VAR,BSZ)  ) *
CSL     .       (K1/EI)**(TRAD/2.+C1) * (K1p/EF)**(TRAD/2.+C2) *
CSL     .       (C5 + C6*ALOG(4.*EB*VAR*SINSQ/FMESQ))/GAMT

CSL   S. Stein et al. Phys. Rev. D 12 7
      xb    = 4./3.  ! A45
      XMT   = FMT    ! Mass of the target
      EMASS = FME
      THR   = TH
      ARG   = COS(THR/2.)**2

      xeta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )
      xi  = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi  = xi/( (Z+xeta)*LOG(183.*Z**(-1./3.)) )
      R   = ( XMT+EI*(1-COS(THR)) )/( XMT-(EF)*(1-COS(THR)) )
      QMS1  = 4.*(VAR)*(EF)*SIN(THR/2.)**2    
      
      IF (INTERNAL) THEN
        tr1   = 1./xb*(ALPH/PI)*(LOG(QMS1/EMASS**2)-1.)

        D2    = 13.*(LOG(QMS1/EMASS**2)-1.)/12.-17./36.
        D2    = D2 - 1./4.*( LOG( EI/(EF) ) )**2       ! Corr. to peak. approx.
        D2    = D2 + 0.5*(PI**2/6.-SPENCE)
        D2    = D2 * (2.*ALPH/PI)
        D2    = D2 + 0.5772*xb*(Tb+Ta)
      ELSE
C        write(6,*) "DEBUG: FES, INTERNAL FALSE"
        tr1   = 0.0
        D2    = 0.5772*xb*(Tb+Ta)
      ENDIF
      Tpb   = tr1 + Tb
      Tpa   = tr1 + Ta
      
      FES   = ( xb*Tpb/(EI-VAR) ) *phi((EI-VAR)/(EI))   
      FES   = FES + xi/(2.*(EI-VAR)**2)             
      FES   = FES * FTCS(VAR,EF)*(1.+D2)            
      FES   = FES * ( (EI-VAR)/(EF*R) )**(xb*Tpa)    
      FES   = FES * ( (EI-VAR)/(EI)   )**(xb*Tpb)   


      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CSL
      function phi(x)
      IMPLICIT REAL (A-H,O-Z)
      phi=1.0-x+3./4.*x**2
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION FEP(VAR)
      PARAMETER(NPTS=1500,NE0=13)

      DIMENSION  ES(NE0),NP(NE0)
      DIMENSION AT(NE0, NPTS), PA(NE0, NPTS)
      COMMON ALPH,AP,AT,BSWZW,BSZ,C1,C2,C3,C4,C5,C6,C7,D1,EF,EI
      COMMON ES,FME,FMESQ,FMOTT,FMT,GAMT,NE,NP,PA
      COMMON PI,SCZ,SCZW,SINSQ,TH,TRAD,TXCO,COT,Z,Tb,Ta,SPENCE
      COMMON INTERNAL,EXTERNAL
      LOGICAL INTERNAL,EXTERNAL
CSL ORIGINAL CODE
CSL     X=D1/3.
CSL      Y=(VAR-EF)*C7
CSL     FK=FMIN(X,Y)
CSL      EO=EI
CSL      EB=EF
CSL      FEP=FTCS(EO,VAR)*(V(VAR,EB)/(VAR-EB)+C4*WM(VAR,EB,
CSL     . BSWZW)+TXCO*WM(VAR,EB,BSZ))*(FK/D1)**(TRAD+C1+C2)*
CSL     . (C5+C6*ALOG(4.*EO*VAR*SINSQ/FMESQ))/GAMT

CSL Miller
CSL      EO=EI
CSL      EB=EF
CSL      R     = C7**2
CSL      EPMAX = EI/(1.0 + 2.0*EI/FMT*SINSQ)
CSL
CSL      X1    = EI/3.
CSLCC     X2    = R*(EPMAX-VAR)
CSLCC     X2    = R*(EPMAX-EF)
CSL      X2    = R*(VAR-EF)   ! stein
CSL      K1    = FMIN(X1,X2)
CSL
CSL
CSL      Y1    = EF/3.
CSLCC     Y2    = (EPMAX-VAR)/R
CSLCC     Y2    = (EPMAX-EF)/R
CSL      Y2    = (VAR-EF)/R  ! stein
CSL      K1p   = FMIN(Y1,Y2)
CSL
CSL      FEP   = FTCS(EO,VAR) * (V(VAR,EB)/(VAR-EB) +
CSL    .        C4*WM(VAR,EB,BSWZW) + TXCO*WM(VAR,EB,BSZ) ) *
CSLCC     .        (K1/EI)**(TRAD/2.+C1)*(K1p/EF)**(TRAD/2.+C2)*
CSL     .        (K1/EI)**(TRAD/2.+C1)*(K1p/VAR)**(TRAD/2.+C2)*  ! stein
CSL     .        (C5+C6*ALOG(4.*EO*VAR*SINSQ/FMESQ))/GAMT
CSL

CSL   S. Stein et al. Phys. Rev. D 12 7

      xb    = 4./3.  ! A45
      XMT   = FMT    ! Mass of the target
      EMASS = FME
      THR   = TH
      ARG   = COS(THR/2.)**2

      xeta  = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )
      xi    = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi    = xi/( (Z+xeta)*LOG(183.*Z**(-1./3.)) )
      R     = ( XMT+EI*(1-COS(THR)) )/( XMT-(EF)*(1-COS(THR)) )

      QMS2  = 4.*EI*(VAR)*SIN(THR/2.)**2  

      IF (INTERNAL) THEN
        tr2   = 1./xb*(ALPH/PI)*(LOG(QMS2/EMASS**2)-1.)
        D2    = 13.*(LOG(QMS2/EMASS**2)-1.)/12.-17./36.
        D2    = D2 - 1./4.*( LOG( EI/(EF) ) )**2 !KS. Correction to peak. approx.
        D2    = D2 + 0.5*(PI**2/6.-SPENCE)
        D2    = D2 * (2.*ALPH/PI)
        D2    = D2 + 0.5772*xb*(Tb+Ta)
      ELSE
C        write(6,*) "DEBUG: FEP, INTERNAL FALSE"
        tr2   = 0.0
        D2    = 0.5772*xb*(Tb+Ta)
      ENDIF

      Tpb   = tr2 + Tb
      Tpa   = tr2 + Ta

      FEP    = ( xb*Tpa/(VAR-EF) ) *phi((VAR-EF)/(VAR))
      FEP    = FEP + xi/(2.*(VAR-EF)**2)
      FEP    = FEP * FTCS(EI,VAR)*(1.+D2)
      FEP    = FEP * ( (VAR-EF)/(VAR) )**(xb*Tpa)
      FEP    = FEP * ( (VAR-EF)*R/(EI))**(xb*Tpb)
      RETURN
      END

CSL-----------------------------------------------------------------------------
CSL----------SUBROUTINES THAT ARE NO LONGER USED--------------------------------
CSL-----------------------------------------------------------------------------

C
C ==>   SUBROUTINE AVGPRP COLLAPSES AN ARRAY OF LENGTH NOLD INTO ONE
C       OF LENGTH NNEW: E=ENERGY, S=CROSS-SEC, DS=ERROR, DE=DELTA-E,
C       IBIN=AVERAGING INTERVAL
C

      SUBROUTINE AVGPRP(E, S, DS, DE, IBIN, NOLD, NNEW)

      DIMENSION  E(1), S(1), DS(1)

      L = NOLD/IBIN
      L = L*IBIN
      INEW = 1

C
C
C ==>   WE CHOOSE E(INEW) TO BE LEFT EDGE TO PRESERVE TOTAL COMPATA
C       WITH ALL OTHER PROGRAMS WHICH WILL DEPEND UPON THIS PROGRAM.
C

      DO 100  I = 1,L, IBIN
         SX = 0.
         DSX = 0.
         E(INEW) = E(I)
         M = IBIN
         K = I + IBIN - 1

      DO 102  J = I, K
         IF(S(J) .LE. 0.) M = M - 1
         IF(S(J) .LT. 0.) S(J) = 0.
         SX = SX + S(J)
         IF(S(J) .LE. 0.) DS(J) = 0.
         DSX = DSX + DS(J)**2
102   CONTINUE

C
C
C ==>   TO PREVENTS DIVIDE CHECK IF ALL BINS EMPTY
C       THE INDEX M GIVES THE CORRECT AVERAGE IN A NEW BIN WITH COM
C       BINS CONTAINING ZERO'S, BUT THE ENERGY POINT IS NOT MOVED. SEE AB
C
C

         IF(M .LE. 0) M = 1

         AM =  M
         S(INEW) = SX/AM
         DS(INEW) = SQRT(DSX)/AM
         INEW = INEW + 1

100   CONTINUE

      NNEW = INEW - 1
      DE = IBIN*DE
C
C
C ==>   WE ERASE CONTENTS OF THE ORIGINAL BINS BEYOND INEW WHICH IS APPROX NO
C
      DO 200  I = INEW, NOLD
              E(I) = 0.
              DS(I) = 0.
              S(I) = 0.
200   CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PUTROW(A, V, INDX, NJ)
      PARAMETER(NPTS=1500,NE0=13)

      DIMENSION A(NE0,NPTS), V(NPTS)

      DO 1  J = 1, NJ
1     A(INDX, J) = V(J)

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C
C ==>	SUBROUTINE TAILPL REMOVE RADIATIVE TAIL OF THE ELASTIC PEAK
C
C	RADIATIVE TAIL SUBTRACTION BEGINS ONLY FOR EXCITED ENERGY
C	GREATER THAN 3 MEV:
C
C	- IK = NUMBER OF THE PARTICULAR SPECTRUM
C	- NPI = NUMBER OF POINTS IN THAT DATA SPECTRUM
C	- ES = BEAM ENERGY
C	- ECR = THE CROSS-SECT MATRIX
C	- EER = ERROR IN CROSS SECTION
C
C	AN ESTIMATED ERROR FOR THE TAIL IS CALCULATED BY TAKING THE
C	RELATIVE STRENGTH OF THE BORIE TERM RELATIVE TO THE REST OF THE
C	PEAKING APPROX AND SAYING THAT THAT FRACTION IS THE UNCERTAINTY
C	IN THE TAIL. NOT MUCH BUT SOMETHING
C

      SUBROUTINE TAILPL(ECR, EP, ES, SINSQ, IK, NPI, JSTRT,DE,EER,Z,A)
      PARAMETER(NPTS=1500,NE0=13)

      DIMENSION ECR(NE0, NPTS), EP(NE0, NPTS), ETAIL(200), SGTAIL(200)
      DIMENSION EER(NE0, NPTS)

      CALL ZERO(ETAIL, 200)
      CALL ZERO(SGTAIL,200)
      WRITE(6,*) "DEBUG IN TAILPL"
C
C ==>	NPTAIL = NUMBER OF POINTS IN TAIL CALC (MUST BE < OR = TO 80)
C

      READ (1,508) NPTAIL
508   FORMAT(I5, F10.5)

      WRITE (2,513) NPTAIL
513   FORMAT(/,' NPTAIL = ', I5)

      IF(NPTAIL.GT.200) WRITE (2,401)

401   FORMAT (1X,'  NPTAIL IS GREATER THAN 200')

      BORIE = 3.1415926*Z/137.
      SINE = SQRT(SINSQ)
      COS2 = 1. - SINSQ
      CERR = SINE*(1. - SINE)/COS2
      BORIE = BORIE * CERR

      DO 509  ITL = 1, NPTAIL
         NTL = NPTAIL - ITL + 1
         READ (1,510) ETAIL(NTL), SGTAIL(NTL)
510      FORMAT(6G10.4)
509   CONTINUE

      WRITE (2,400) (ETAIL(I), SGTAIL(I), I = 1, NPTAIL)
400   FORMAT (1X,6G10.4)

      ESHIFT = ES/(1. + 2.*ES*SINSQ/931./ A)
      ESHIFT = ES - ESHIFT
      JSTRT = 1

      DO 511  JS = 1, NPI
         EX = EP(IK, JS)-DE/2.

C
C ==>	SETS TO LEFT EDGE OF BIN. AGREES WITH CONVENTIONS
C

         IF(EX .LE. ES - ESHIFT - 3.) GO TO 512
         JSTRT = JS

C
C
C	TELLS MAIN PROGRAM WHERE TO START UNFOLDING AFTER ELASTIC PK
C	REPLACES, EFFECTIVELY, NEXT CARD, NOW INOPERATIVE
C	ECR(IK, JS) = 0.
C	SETS ALL CROSS-SECTION ABOVE END OF ELASTIC PK = 0. SO WE DON'T
C	DOUBLE-COUNT WHEN PULLING TAIL
C

         GO TO 511

512      CALL TTERP(EX, RTL, ETAIL, SGTAIL, NPTAIL)

         IF(ECR(IK, JS) .GT. 0.)
     .ECR(IK, JS) = ECR(IK, JS) - RTL

C
C
C ==>	RTL IS THE RETURNED INTERPOLATED VALUE OF RADIATIVE TAIL
C

         IF(ECR(IK, JS) .LT. 0.) ECR(IK, JS) = 0.
         EER(IK, JS) = EER(IK, JS) + RTL * BORIE

511   CONTINUE

      RETURN
      END

C---------------------------------------------------------------------------
C     Some new subroutines that I added
C---------------------------------------------------------------------------

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)  ! from Numerical recipes
      IMPLICIT NONE
      INTEGER N,NMAX
      REAL YP1,YPN,X(N),Y(N),Y2(N)
      PARAMETER (NMAX=1500)
      INTEGER I,K
      REAL P,QN,SIG,UN,U(NMAX)
      IF (YP1.GT..99E30) THEN
         Y2(1) = 0.
         U(1)  = 0.
      ELSE
         Y2(1) = -0.5
         U(1)  = ( 3./(X(2)-X(1)) ) *
     &         ( ( Y(2)-Y(1) )/( X(2)-X(1) ) - YP1 )
      ENDIF
      DO I=2,N-1
         SIG  = ( X(I)-X(I-1) ) / ( X(I+1)-X(I-1) )
         P    = SIG*Y2(I-2) + 2
         Y2(I)= (SIG-1)/P
         U(I) = (6.*( (Y(I+1)-Y(I)   )/( X(I+1) - X(I) )
     &              - (Y(I)  -Y(I-1) )/( X(I)   - X(I-1)))
     &              / (X(I+1)-X(I-1) )- SIG*U(I-1))/P
      ENDDO
     
      IF (YPN.GT..99E30) THEN
         QN = 0.
         UN = 0.
      ELSE
         QN = 0.5
         UN = ( 3./(X(N)-X(N-1)) ) *
     &        ( YPN - ( Y(N)-Y(N-1) ) / ( X(N)-X(N-1) ) )
      ENDIF
     
      Y2(N) = ( UN - QN*U(N-1) ) / ( QN*Y2(N-1)+1. )
      DO K=N-1,1,-1
         Y2(K) = Y2(K)*Y2(K+1) + U(K)
      ENDDO 
      RETURN
      END


      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      IMPLICIT NONE
      INTEGER N
      REAL X,Y,XA(N),Y2A(N),YA(N)
      INTEGER K,KHI,KLO
      REAL A,B,H
C      write(6,*) "debug:",N,XA(N-100)
      KLO = 1
      KHI = N
 1    IF (KHI-KLO.GT.1) THEN
        K = ( KHI+KLO ) / 2
        IF ( XA(K).GT.X ) THEN
          KHI = K
        ELSE
          KLO = K
        ENDIF
      GOTO 1
      ENDIF
      H = XA(KHI) - XA(KLO)
      IF ( H.EQ.0. ) THEN
C         write(6,*) KHI,XA(KHI),XA(KLO)
         PAUSE 'BAD XA INPUT IN SPLINT'
      ENDIF
      A = ( XA(KHI)-X ) / H
      B = ( X-XA(KLO) ) / H
      Y = A*YA(KLO) + B*YA(KHI) +
     &    ( (A**3-A)*Y2A(KLO) + (B**3-B)*Y2A(KHI) ) * (H**2)/6.
      RETURN
      END


C-----
      SUBROUTINE GET_FINITE_COR(M,XNU2,XFACTOR,NFINITE)
      PARAMETER(NL=1500,NE=13)
      INTEGER NFINITE(NE)
      REAL XNU2(NL,NE), XFACTOR(NL,NE)

      CHARACTER LINE1*80, CELL(6)*4
      DATA CELL/'8_SY','1_AR','2_JN','3_AR','4_2T','5_JN'/

      OPEN(UNIT=12,
     &FILE='input/'//CELL(M)//'_qfs_finite_acceptance.out',
     &STATUS='old')

      N=0
 10   CONTINUE
      READ(12,'(A)',ERR=20,END=20) LINE1
      N=N+1
      IF (LINE1(1:1).EQ.'#') THEN
        N=N-1
        GOTO 10
      ENDIF

      READ(LINE1,*) XNU2(N,M),XJUNK,XJUNK,XFACTOR(N,M)
      GOTO 10
 20   CONTINUE
      NFINITE(M) = N
      RETURN
      END


C     Seonho Choi's "Interpol_1dim"
C     -----------------------------
CSL   SUBROUTINE INTERP(N ,XP,YP,EYP,X,Y,EY,IFLAG,ITYPE_ERR,IERR)
      SUBROUTINE INTERP(NI,XI,YI,EYI,X,Y,EY,IFLAG,ITYPE_ERR,IERR)
      PARAMETER (DELTA=1.0E-3)
CSL   REAL XP,YP,EYP,X,Y,EY
      REAL XP,YP,EYP,X,Y,EY,XI,YI,EYI
      INTEGER N,ITYPE_ERR
CSL   DIMENSION XP(*),YP(*),EYP(*),NPT_MIN(20)
      DIMENSION XI(*),YI(*),EYI(*),NPT_MIN(20)
      DIMENSION XP(2000),YP(2000),EYP(2000)
      LOGICAL L_LEFT,L_RIGHT
      DATA NPT_MIN/2,3,3,4,4,5,6,7,8,11*0/

      IF (ITYPE_ERR.GE.0) THEN
        CALL STRIP(NI,XI,YI,EYI,N,XP,YP,EYP)
      ELSE
        DO KK=1,NI
           XP(KK)=XI(KK)
           YP(KK)=YI(KK)
          EYP(KK)=EYI(KK)
        ENDDO
        N=NI
      ENDIF
*
* IFLAG = 1 : linear interpolation
* IFLAG = 2 : quadratic interpolation (a x b c)
* IFLAG = 3 : quadratic interpolation (a b x c)
* IFLAG = 4 : quadratic interpolation (average if IFLAG=2 and IFLAG=3)
* IFLAG = 5 : cubic interpolation
* IFLAG = 6 : 4th order
* IFLAG = 7 : 5th order
* IFLAG = 8 : 6th order
* IFLAG = 9 : 7th order
*
* IFLAG > 10 : does the extrapolation, too
*
* ITYPE_ERR = -1 : no error calculation
* ITYPE_ERR =  0 : all the errors are independent (statistical error)
* ITYPE_ERR =  1 : all the errors are fully correlated (systematic error)
*
* IERR = 0  : normal exit (interpolation)
* IERR = 1  : normal exit (left extrapolation)
* IERR = 2  : normal exit (right extrapolation)
* IERR = -1 : left out of bound
* IERR = -2 : right out of bound
* IERR = -3 : only 1 data point
*
*
* Check if we have sufficient number of points for a given method
*
      JFLAG = MOD(IFLAG,10)
      IF (N.LT.NPT_MIN(JFLAG)) THEN
         IF (N.LT.2) THEN ! nothing is possible, return
            IERR = -3
            RETURN
         ELSE IF (N.EQ.2) THEN ! only linear interpolation is possible
            JFLAG = 1
         ELSE IF (N.EQ.3.AND.JFLAG.GT.2) THEN  ! only quadratic interpolation is possible
            JFLAG = 2
         ELSE IF (N.EQ.4.AND.JFLAG.GT.5) THEN
            JFLAG = 5
         ELSE IF (N.GE.5.AND.JFLAG.GT.6) THEN
            JFLAG = N + 1
         ENDIF
      ENDIF

      L_LEFT = .FALSE.
      L_RIGHT = .FALSE.
      IF ((X-XP(1))/X.LT.-DELTA) L_LEFT = .TRUE.
      IF ((X-XP(N))/X.GT.DELTA) L_RIGHT = .TRUE.
*
* No extrapolation when JFLAG < 10
*
      IF (IFLAG.LT.10) THEN
         IF (L_LEFT) THEN
            IERR = -1
            Y = YP(1)
            EY = EYP(1)
            RETURN
         ELSE IF (L_RIGHT) THEN
            IERR = -2
            Y = YP(N)
            EY = EYP(N)
            RETURN
         ENDIF
      ENDIF
*
* When one of XP(I) is close enough to X, no need to go further
*
      DO I = 1,N
         IF (ABS((X-XP(I))/X).LT.DELTA) THEN
            Y = YP(I)
            EY = EYP(I)
            IERR = 0
            RETURN
         ENDIF
      ENDDO

*
* Find the corresponding interval
*      
      J = 0
      IF (L_LEFT) THEN
         J = 0
      ELSE IF (L_RIGHT) THEN
         J = N
      ELSE
         J = -1
         DO I = 1,N-1
            IF (X.GT.XP(I).AND.X.LT.XP(I+1)) THEN
               J = I
            ENDIF
         ENDDO
      ENDIF
*
* If everything is OK, you should not have J = -1, but just in case
*
      IF (J.EQ.-1) THEN
         IERR = -4
         RETURN
      ENDIF
*
* At boundaries, force other methods for quadratic interpolation
*
      IF (J.EQ.0.OR.J.EQ.1) THEN
         IF (JFLAG.EQ.3.OR.JFLAG.EQ.4) JFLAG = 2
      ELSE IF (J.EQ.N-1.OR.J.EQ.N) THEN
         IF (JFLAG.EQ.2.OR.JFLAG.EQ.4) JFLAG = 3
      ENDIF

      IF (JFLAG.EQ.1) THEN
         K = J
         L = K + 1
      ELSE IF (JFLAG.EQ.2) THEN
         K = J
         L = K + 2
      ELSE IF (JFLAG.EQ.3) THEN
         K = J - 1
         L = J + 1
      ELSE IF (JFLAG.EQ.4.OR.JFLAG.EQ.5) THEN
         K = J - 1
         L = J + 2
      ELSE IF (JFLAG.EQ.6) THEN
         K = J - 1
         L = J + 3
      ELSE IF (JFLAG.EQ.7) THEN
         K = J - 2
         L = J + 3
      ELSE IF (JFLAG.EQ.8) THEN
         K = J - 2
         L = J + 4
      ELSE IF (JFLAG.EQ.9) THEN
         K = J - 3
         L = J + 4
      ENDIF
*
* Consider boundary conditions
*
      IF (K.LT.1) THEN
         L = L + (1 - K)
         K = 1
      ENDIF
      IF (L.GT.N) THEN
         K = K - (L - N)
         L = N
      ENDIF

      IF (JFLAG.EQ.1) THEN ! linear interpolation
         A = (XP(K+1)-X)/(XP(K+1)-XP(K))
         B = (X-XP(K))/(XP(K+1)-XP(K))
         Y = A*YP(K)+B*YP(K+1)
         IF (ITYPE_ERR.EQ.-1) THEN
            EY = -1.0
         ELSE IF (ITYPE_ERR.EQ.0) THEN
            EY = (A*EYP(K))**2+(B*EYP(K+1))**2
            EY = SQRT(EY)
         ELSE IF (ITYPE_ERR.EQ.1) THEN
            EY = A*EYP(K)+B*EYP(K+1)
         ENDIF
      ELSE IF (JFLAG.EQ.2.OR.JFLAG.EQ.3) THEN
         CALL INTERPOL_PARA(XP(K),YP(K),EYP(K),X,Y,EY,ITYPE_ERR)
      ELSE IF (JFLAG.EQ.4) THEN
         CALL INTERPOL_PARA(XP(K),YP(K),EYP(K),X,Y1,EY1,ITYPE_ERR)
         CALL INTERPOL_PARA(XP(K+1),YP(K+1),EYP(K+1),X,Y2,EY2,ITYPE_ERR)
         IF (EY1.GT.0.0) THEN
            W1 = 1.0/EY1**2
         ELSE
            W1 = 1.0
         ENDIF
         IF (EY2.GT.0.0) THEN
            W2 = 1.0/EY2**2
         ELSE
            W2 = 0.0
         ENDIF
         SW = 1.0/(W1 + W2)
         Y = (W1*Y1+W2*Y2)*SW
         EY = SQRT(SW)
      ELSE IF (JFLAG.EQ.5) THEN
         CALL INTERPOL_CUBIC(XP(K),YP(K),EYP(K),X,Y,EY,ITYPE_ERR)
      ELSE
         NDEG = JFLAG - 2
         CALL INTERPOL_NDEG(NDEG,XP(K),YP(K),EYP(K),X,Y,EY,ITYPE_ERR)
      ENDIF
      IERR = 0
      RETURN
      END

      SUBROUTINE INTERPOL_PARA(XP,YP,EYP,X,Y,EY,ITYPE_ERR)
      REAL XP,YP,EYP,X,Y
      DIMENSION XP(*),YP(*),EYP(*)
      INTEGER ITYPE_ERR

      REAL A,B,C,D1,D2,D3,D12,D13,D23

      D1 = X - XP(1)
      D2 = X - XP(2)
      D3 = X - XP(3)

      D12 = XP(1) - XP(2)
      D13 = XP(1) - XP(3)
      D23 = XP(2) - XP(3)

      A = D1*D2/(D13*D23)
      B = D2*D3/(D12*D13)
      C = D3*D1/(-D23*D12)

      Y = A*YP(3)+B*YP(1)+C*YP(2)
      IF (ITYPE_ERR.EQ.-1) THEN
         EY = -1.0
      ELSE IF (ITYPE_ERR.EQ.0) THEN
         EY = (A*EYP(3))**2 + (B*EYP(1))**2 + (C*EYP(2))**2
         EY = SQRT(EY)
      ELSE IF (ITYPE_ERR.EQ.1) THEN
         EY = A*EYP(3) + B*EYP(1) + C*EYP(2)
      ENDIF
      RETURN
      END

      SUBROUTINE INTERPOL_CUBIC(XP,YP,EYP,X,Y,EY,ITYPE_ERR)
      REAL XP,YP,EYP,X,Y,EY
      DIMENSION XP(*),YP(*),EYP(*),C(4)
      INTEGER ITYPE_ERR

      D1 = X - XP(1)
      D2 = X - XP(2)
      D3 = X - XP(3)
      D4 = X - XP(4)

      D12 = XP(1) - XP(2)
      D13 = XP(1) - XP(3)
      D14 = XP(1) - XP(4)
      D23 = XP(2) - XP(3)
      D24 = XP(2) - XP(4)
      D34 = XP(3) - XP(4)

      C(1) = (D2*D3*D4)/(D12*D13*D14)
      C(2) = (D1*D3*D4)/(-D12*D23*D24)
      C(3) = (D1*D2*D4)/(D13*D23*D34)
      C(4) = (D1*D2*D3)/(-D14*D24*D34)

      Y = 0.0
      EY = 0.0
      DO I = 1,4
         Y = Y + C(I)*YP(I)
         IF (ITYPE_ERR.EQ.-1) THEN
            EY = -1.0
         ELSE IF (ITYPE_ERR.EQ.0) THEN
            EY = EY + (C(I)*EYP(I))**2
         ELSE IF (ITYPE_ERR.EQ.1) THEN
            EY = EY + C(I)*EYP(I)
         ENDIF
      ENDDO
      IF (ITYPE_ERR.EQ.0) THEN
         EY = SQRT(EY)
      ENDIF
      RETURN
      END

      SUBROUTINE INTERPOL_NDEG(NDEG,XP,YP,EYP,X,Y,EY,ITYPE_ERR)
      REAL XP,YP,EYP,X,Y,EY
      INTEGER NDEG,ITYPE_ERR
      DIMENSION XP(*),YP(*),EYP(*),D(100),DIJ(100,100),C(100)

      NP = NDEG + 1

      DO I = 1,NP
         D(I) = X - XP(I)
         DO J = I+1,NP
            DIJ(I,J) = XP(I) - XP(J)
            DIJ(J,I) = -DIJ(I,J)
         ENDDO
      ENDDO

      DO I = 1,NP
         C(I) = 1.0
         DO J = 1,NP
            IF (J.NE.I) C(I) = C(I)*D(J)/DIJ(I,J)
         ENDDO
      ENDDO

      Y = 0.0
      EY = 0.0
      DO I = 1,NP
         Y = Y + C(I)*YP(I)
         IF (ITYPE_ERR.EQ.-1) THEN
            EY = -1.0
         ELSE IF (ITYPE_ERR.EQ.0) THEN
            EY = EY + (C(I)*EYP(I))**2
         ELSE IF (ITYPE_ERR.EQ.1) THEN
            EY = EY + C(I)*EYP(I)
         ENDIF
      ENDDO
      IF (ITYPE_ERR.EQ.0) THEN
         EY = SQRT(EY)
      ENDIF
      RETURN
      END


      SUBROUTINE STRIP(NIN,XIN,YIN,EIN,NOUT,XOUT,YOUT,EOUT)
      PARAMETER(NPTS      = 2000  )
      REAL XIN(*),YIN(*),EIN(*),XOUT(*),YOUT(*),EOUT(*)
      K = 0
      DO J = 1,NIN
        IF ( (EIN(J).GT.0) ) THEN
          K       = K + 1
          XOUT(K) = XIN(J)
          YOUT(K) = YIN(J)
          EOUT(K) = EIN(J)
        ENDIF
      ENDDO
      NOUT = K

      ISTART=NOUT+1
      DO J = ISTART,NIN
        XOUT(J) = 0.0
        YOUT(J) = 0.0
        EOUT(J) = -1.E6
      ENDDO

      RETURN
      END
      
      REAL FUNCTION GET_W(E0,TH,XNU)
      XM  = 938.272
      XMT = XM
      PI  = ACOS(-1.)
      THR = TH*PI/180.0
      S2  = ( SIN(THR/2.0) )**2

      EP = E0-XNU
      Q2 = 4.0*E0*EP*S2
      W2 = XM**2 - Q2 + 2.0*XM*XNU
      GET_W  = SQRT(W2)
      RETURN
      END

