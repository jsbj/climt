      SUBROUTINE GFLUXV(DTDEL,TDEL,TAUCUM,WDEL,CDEL,UBAR0,F0PI,RSF,
     *                  BTOP,BSURF,FMIDP,FMIDM,DIFFV,FLUXUP,FLUXDN)

C  GCM2.0  Feb 2003
C
C ?? COMBINE INTO ONE GLFLUX ROUTINE, NOT VIS AND IR
C 
C  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITIONS
C  FOR THE VISIBLE  FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
C  THE LEVELS. THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
C  MEASURED FROM THE TOP OF EACH LAYER.  (DTAU) TOP OF EACH LAYER HAS  
C  OPTICAL DEPTH TAU(N).IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
C  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
C  FROM TOP TO BOTTOM. SEE C.P. MCKAY, TGM NOTES.
C THIS SUBROUTINE DIFFERS FROM ITS IR COUNTERPART IN THAT HERE WE SOLVE FOR 
C THE FLUXES DIRECTLY USING THE GENERALIZED NOTATION OF MEADOR AND WEAVOR
C J.A.S., 37, 630-642, 1980.
C THE TRI-DIAGONAL MATRIX SOLVER IS DSOLVER AND IS DOUBLE PRECISION SO MANY 
C VARIABLES ARE PASSED AS SINGLE THEN BECOME DOUBLE IN DSOLVER
C
C NLL           = NUMBER OF LEVELS (NAYER + 1) THAT WILL BE SOLVED
C NAYER         = NUMBER OF LAYERS (NOTE DIFFERENT SPELLING HERE)
C WAVEN         = WAVELENGTH FOR THE COMPUTATION
C DTDEL(NLAYER) = ARRAY OPTICAL DEPTH OF THE LAYERS
C TDEL(NLL)     = ARRAY COLUMN OPTICAL DEPTH AT THE LEVELS
C WDEL(NLEVEL)  = SINGLE SCATTERING ALBEDO
C CDEL(NLL)     = ASYMMETRY FACTORS, 0=ISOTROPIC
C UBARV         = AVERAGE ANGLE, 
C UBAR0         = SOLAR ZENITH ANGLE
C F0PI          = INCIDENT SOLAR DIRECT BEAM FLUX
C RSF           = SURFACE REFLECTANCE
C BTOP          = UPPER BOUNDARY CONDITION ON DIFFUSE FLUX
C BSURF         = REFLECTED DIRECT BEAM = (1-RSFI)*F0PI*EDP-TAU/U
C FP(NLEVEL)    = UPWARD FLUX AT LEVELS
C FM(NLEVEL)    = DOWNWARD FLUX AT LEVELS
C FMIDP(NLAYER) = UPWARD FLUX AT LAYER MIDPOINTS
C FMIDM(NLAYER) = DOWNWARD FLUX AT LAYER MIDPOINTS
C added Dec 2002
C DIFFV         = downward diffuse solar flux at the surface
C 
      implicit none
      integer NL
      PARAMETER (NL=101)  ! MUST BE LARGER THAN NLL

      include "grid.h"
      include "radinc.h"

      REAL*8 EM, EP
      REAL*8 W0(L_NLAYRAD), COSBAR(L_NLAYRAD), DTAU(L_NLAYRAD)
      REAL*8 TAU(L_NLEVRAD), WDEL(L_NLAYRAD), CDEL(L_NLAYRAD)
      REAL*8 DTDEL(L_NLAYRAD), TDEL(L_NLEVRAD)
      REAL*8 FMIDP(L_NLAYRAD), FMIDM(L_NLAYRAD)
      REAL*8 LAMDA(NL), ALPHA(NL),XK1(NL),XK2(NL)
      REAL*8 G1(NL), G2(NL), G3(NL), GAMA(NL), CP(NL), CM(NL), CPM1(NL)
      REAL*8 CMM1(NL), E1(NL), E2(NL), E3(NL), E4(NL), EXPTRM(NL)
      REAL*8 FLUXUP, FLUXDN
      REAL*8 TAUCUM(L_LEVELS), taucump(L_LEVELS+2)

      integer NAYER, L, K
      real*8  ubar0, f0pi, rsf, btop, bsurf, g4, denom, am, ap
      real*8  taumax, taumid, cpmid, cmmid
      real*8  diffv

C======================================================================C
 
      NAYER  = L_NLAYRAD
      TAUMAX = L_TAUMAX    !Default is 35.0

C  Delta function updates.  Use scaled values, calculated below.      
c     DO L=1,L_NLAYRAD
c       W0(L)     = WDEL(L)
c       COSBAR(L) = CDEL(L)
c       DTAU(L)   = DTDEL(L)
c       TAU(L)    = TDEL(L)
c     END DO
c     TAU(L_NLEVRAD) = TDEL(L_NLEVRAD)

C  Scaled (Delta function) values

      TAU(1)=TDEL(1)*(1.-WDEL(1)*CDEL(1)**2)
      TAUCUMP(1) = 0.0
      TAUCUMP(2) = TAUCUM(2)*(1.-WDEL(1)*CDEL(1)**2)
      TAUCUMP(3) = TAUCUM(2)+ (TAUCUM(3)-TAUCUM(2))*
     *             (1.-WDEL(1)*CDEL(1)**2)
      DO L=1,L_NLAYRAD-1
        W0(L)        = WDEL(L)*(1.-CDEL(L)**2)/
     *                 (1.-WDEL(L)*CDEL(L)**2)
        COSBAR(L)    = CDEL(L)/(1.+CDEL(L))
        DTAU(L)      = DTDEL(L)*(1.-WDEL(L)*CDEL(L)**2)
        TAU(L+1)     = TAU(L)+DTAU(L)
        K            = 2*(L+1)
        TAUCUMP(K)   = TAU(L+1)
        TAUCUMP(K+1) = TAUCUMP(k)+(TAUCUM(K+1) - TAUCUM(k))*
     *                 (1.-WDEL(L)*CDEL(L)**2)
      END DO

C  Bottom layer

      L=L_NLAYRAD
      W0(L)      = WDEL(L)*(1.-CDEL(L)**2)/(1.-WDEL(L)*CDEL(L)**2)
      COSBAR(L)  = CDEL(L)/(1.+CDEL(L))
      DTAU(L)    = DTDEL(L)*(1.-WDEL(L)*CDEL(L)**2)
      TAU(L+1)   = TAU(L)+DTAU(L)
      K          = 2*(L+1)
      TAUCUMP(K) = TAU(L+1)

C     WE GO WITH THE QUADRATURE APPROACH HERE.  THE "SQRT(3)" factors
C     ARE THE UBARV TERM.

      DO L=1,L_NLAYRAD
        
        ALPHA(L)=SQRT( (1.0-W0(L))/(1.0-W0(L)*COSBAR(L)) )

C       SET OF CONSTANTS DETERMINED BY DOM 

        G1(L)    = (SQRT(3.0)*0.5)*(2.0- W0(L)*(1.0+COSBAR(L)))
        G2(L)    = (SQRT(3.0)*W0(L)*0.5)*(1.0-COSBAR(L))
        G3(L)    = 0.5*(1.0-SQRT(3.0)*COSBAR(L)*UBAR0)
        LAMDA(L) = SQRT(G1(L)**2 - G2(L)**2)
        GAMA(L)  = (G1(L)-LAMDA(L))/G2(L)
      END DO

      DO L=1,L_NLAYRAD
        G4    = 1.0-G3(L)
        DENOM = LAMDA(L)**2 - 1./UBAR0**2
 
C       THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C       THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C       THE SCATTERING GOES TO ZERO
C       PREVENT THIS WITH AN IF STATEMENT

        IF ( DENOM .EQ. 0.) THEN
          DENOM=1.E-10
        END IF

        AM = F0PI*W0(L)*(G4   *(G1(L)+1./UBAR0) +G2(L)*G3(L) )/DENOM
        AP = F0PI*W0(L)*(G3(L)*(G1(L)-1./UBAR0) +G2(L)*G4    )/DENOM

C       CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
C       AT THE TOP OF THE LAYER, THAT IS LOWER   OPTICAL DEPTH TAU(L)
 
        CPM1(L) = AP*EXP(-TAU(L)/UBAR0)
        CMM1(L) = AM*EXP(-TAU(L)/UBAR0)

C       CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
C       BOTTOM OF THE LAYER.  THAT IS AT HIGHER OPTICAL DEPTH TAU(L+1)

        CP(L) = AP*EXP(-TAU(L+1)/UBAR0)
        CM(L) = AM*EXP(-TAU(L+1)/UBAR0)

      END DO
 
C     NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
C     FOR THE TRIDIAGONAL ROTATED LAYERED METHOD

      DO L=1,L_NLAYRAD
        EXPTRM(L) = MIN(TAUMAX,LAMDA(L)*DTAU(L))  ! CLIPPED EXPONENTIAL
        EP = EXP(EXPTRM(L))

        EM        = 1.0/EP
        E1(L)     = EP+GAMA(L)*EM
        E2(L)     = EP-GAMA(L)*EM
        E3(L)     = GAMA(L)*EP+EM
        E4(L)     = GAMA(L)*EP-EM
      END DO
 
      CALL DSOLVER(NAYER,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,
     *             BSURF,RSF,XK1,XK2)

C     NOW WE CALCULATE THE FLUXES AT THE MIDPOINTS OF THE LAYERS.

C  For original, unscaled, version, use TAUCUM, not TAUCUMP 
      DO L=1,L_NLAYRAD-1
c     EXPTRM(L) = MIN(TAUMAX,LAMDA(L)*(TAUCUM(2*L+1)-TAUCUM(2*L)))
      EXPTRM(L) = MIN(TAUMAX,LAMDA(L)*(TAUCUMP(2*L+1)-TAUCUMP(2*L)))
 
        EP = EXP(EXPTRM(L))

        EM    = 1.0/EP
        G4    = 1.0-G3(L)
        DENOM = LAMDA(L)**2 - 1./UBAR0**2

C       THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C       THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C       THE SCATTERING GOES TO ZERO
C       PREVENT THIS WITH A IF STATEMENT

        IF ( DENOM .EQ. 0.) THEN
          DENOM=1.E-10
        END IF

        AM = F0PI*W0(L)*(G4   *(G1(L)+1./UBAR0) +G2(L)*G3(L) )/DENOM
        AP = F0PI*W0(L)*(G3(L)*(G1(L)-1./UBAR0) +G2(L)*G4    )/DENOM

C       CPMID AND CMMID  ARE THE CPLUS AND CMINUS TERMS EVALUATED
C       AT THE MIDDLE OF THE LAYER.

        TAUMID   = TAUCUMP(2*L+1)
c       TAUMID   = TAUCUM(2*L+1)

        CPMID = AP*EXP(-TAUMID/UBAR0)
        CMMID = AM*EXP(-TAUMID/UBAR0)

        FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
        FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID
 
C       ADD THE DIRECT FLUX TO THE DOWNWELLING TERM

        FMIDM(L)= FMIDM(L)+UBAR0*F0PI*EXP(-MIN(TAUMID,TAUMAX)/UBAR0)
   
      END DO
 
C     FLUX AT THE Ptop layer

      EP    = 1.0
      EM    = 1.0
      G4    = 1.0-G3(1)
      DENOM = LAMDA(1)**2 - 1./UBAR0**2

C     THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C     THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C     THE SCATTERING GOES TO ZERO
C     PREVENT THIS WITH A IF STATEMENT

      IF ( DENOM .EQ. 0.) THEN
        DENOM=1.E-10
      END IF

      AM = F0PI*W0(1)*(G4   *(G1(1)+1./UBAR0) +G2(1)*G3(1) )/DENOM
      AP = F0PI*W0(1)*(G3(1)*(G1(1)-1./UBAR0) +G2(1)*G4    )/DENOM

C     CPMID AND CMMID  ARE THE CPLUS AND CMINUS TERMS EVALUATED
C     AT THE MIDDLE OF THE LAYER.

      CPMID  = AP
      CMMID  = AM

      FLUXUP = XK1(1)*EP + GAMA(1)*XK2(1)*EM + CPMID
      FLUXDN = XK1(1)*EP*GAMA(1) + XK2(1)*EM + CMMID

C     ADD THE DIRECT FLUX TO THE DOWNWELLING TERM

      fluxdn = fluxdn+UBAR0*F0PI*EXP(-MIN(TAUCUMP(2),TAUMAX)/UBAR0)
c     fluxdn = fluxdn+UBAR0*F0PI*EXP(-MIN(TAUCUM(1),TAUMAX)/UBAR0)
 
C     This is for the "special" bottom layer, where we take
C     DTAU instead of DTAU/2.

      L     = L_NLAYRAD 
      EXPTRM(L) = MIN(TAUMAX,LAMDA(L)*(TAUCUMP(L_LEVELS)-
     *                                 TAUCUMP(L_LEVELS-1)))
c     EXPTRM(L) = MIN(TAUMAX,LAMDA(L)*(TAUCUM(L_LEVELS)-
c    *                                 TAUCUM(L_LEVELS-1)))

      EP    = EXP(EXPTRM(L))
      EM    = 1.0/EP
      G4    = 1.0-G3(L)
      DENOM = LAMDA(L)**2 - 1./UBAR0**2

C     THERE IS A POTENTIAL PROBLEM HERE IF W0=0 AND UBARV=UBAR0
C     THEN DENOM WILL VANISH. THIS ONLY HAPPENS PHYSICALLY WHEN 
C     THE SCATTERING GOES TO ZERO
C     PREVENT THIS WITH A IF STATEMENT

      IF ( DENOM .EQ. 0.) THEN
        DENOM=1.E-10
      END IF

      AM = F0PI*W0(L)*(G4   *(G1(L)+1./UBAR0) +G2(L)*G3(L) )/DENOM
      AP = F0PI*W0(L)*(G3(L)*(G1(L)-1./UBAR0) +G2(L)*G4    )/DENOM

C     CPMID AND CMMID  ARE THE CPLUS AND CMINUS TERMS EVALUATED
C     AT THE MIDDLE OF THE LAYER.

      TAUMID   = MIN(TAUCUMP(L_LEVELS),TAUMAX)
c     TAUMID   = MIN(TAUCUM(L_LEVELS),TAUMAX)
      CPMID    = AP*EXP(-MIN(TAUMID,TAUMAX)/UBAR0)
      CMMID    = AM*EXP(-MIN(TAUMID,TAUMAX)/UBAR0)

      FMIDP(L) = XK1(L)*EP + GAMA(L)*XK2(L)*EM + CPMID
      FMIDM(L) = XK1(L)*EP*GAMA(L) + XK2(L)*EM + CMMID

C  Save the diffuse downward flux for TEMPGR calculations

      DIFFV = FMIDM(L)
 
C     ADD THE DIRECT FLUX TO THE DOWNWELLING TERM

      FMIDM(L)= FMIDM(L)+UBAR0*F0PI*EXP(-MIN(TAUMID,TAUMAX)/UBAR0)

      RETURN
      END
