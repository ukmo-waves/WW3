!/ ------------------------------------------------------------------- /
      MODULE PPLOOKMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |        Chris Bunney, UKMO         |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Jan-2016 |
!/                  +-----------------------------------+
!/
!/    13-Jan-2016 : Initial version. (C. Bunney, UKMO)  ( version 4.18 )
!/
!  1. Purpose :
!
!     Module to map PP header entries to fields for use when outputting
!     data to Met Office PP or FieldFile format.
!
!  2. Method :
!
!     Data is read from a text file and mapped to each field using a
!     key derived from the field group (IFI), section (IFJ), component
!     number (ICOMP) and partition number (IPART).
!     A simple linear lookup table is employed - it is envisaged that 
!     the cost of this is magitudes smaller than the acual field
!     I/O so it is not worth writing anything "clever"!
!
!     The ppheaders.txt file should have one line for each entry encoded:
!
!        IFI IFJ ICOM  IPART   Stash FieldCode FieldType Packing Scale
! 
!  4. Subroutines used :
!
!     PP_LOOKUP_INIT    Loads data from the ppheaders.txt file into
!                       the PP_HDR structure array.
!     PP_LLOOKUP        Finds entry in PP_HDR array using key and
!                       sets the relavent fields in the IHEAD and 
!                       RHEAD PP header arrays.
!
!
!/ ------------------------------------------------------------------- /

      USE W3ODATMD, ONLY: NDSO, NDSE

      PUBLIC

      ! Lookup table for PP header information.
      TYPE PP_HDR_TYPE
        INTEGER  :: key ! Combination of IFI, IFJ, ICMP, IPART
        INTEGER  :: stash, fldcode, fldtype
        INTEGER  :: packing
        REAL     :: scalefac
        ! TYPE(PP_HDR), POINTER :: next
      END TYPE PP_HDR_TYPE

      INTEGER, PARAMETER :: MAX_PP_HDR = 200
      INTEGER            :: NHDR

      TYPE(PP_HDR_TYPE), ALLOCATABLE :: PP_HDR(:)

      CONTAINS 

      SUBROUTINE PP_LOOKUP_INIT()
      ! Populates pp header info from file called ppheaders.txt

      IMPLICIT NONE

      INTEGER :: ierr
      INTEGER :: ifi, ifj, icmp, ipart, st, fc, ft, pck
      REAL    :: fac
      CHARACTER :: check

      ! Allocate space:
      IF(.NOT. ALLOCATED(PP_HDR)) ALLOCATE(PP_HDR(MAX_PP_HDR))
      NHDR = 0

      ! Open file and read headers
      OPEN(unit=70, file='ppheaders.txt', status='old', iostat=ierr)
      IF(ierr .NE. 0) THEN
         WRITE(NDSE,*) "Error opening ppheaders.txt"
         WRITE(NDSE,*) "Stashcodes will not be correctly set"
         RETURN
      ENDIF
      
      ! Read line by line, skipping comment ($) lines:
      DO
        READ(70,'(A1)',end=100) check
        IF(check .EQ. '$') CYCLE
        BACKSPACE(70)


        READ(70, *, end=100) ifi, ifj, icmp, ipart, st, fc, ft, pck, fac
        NHDR = NHDR + 1
        PP_HDR(NHDR)%key = ipart + icmp * 10 + ifj * 1E3 + ifi * 1E5
        PP_HDR(NHDR)%stash = st
        PP_HDR(NHDR)%fldcode = fc
        PP_HDR(NHDR)%fldtype = ft
        PP_HDR(NHDR)%packing = pck
        PP_HDR(NHDR)%scalefac = fac
      ENDDO

      WRITE(NDSO,*) "Read ",NHDR, " entries from ppheaders.txt"
 100  CONTINUE

      RETURN
 200  PRINT*,'Failed to read from ppheaders.txt'
      CALL EXIT(10)

      END SUBROUTINE PP_LOOKUP_INIT


      SUBROUTINE PP_LOOKUP(PPIHEAD, PPRHEAD, IFI, IFJ, ICMP, IPART)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: PPIHEAD(64)
      REAL,    INTENT(INOUT) :: PPRHEAD(64)
      INTEGER, INTENT(IN) :: IFI, IFJ, ICMP, IPART
      INTEGER             :: key, i

      ! Key is encoded as: IIJJCP
      !  wjere I=IFI, J=IFJ, C=ICMP, P=IPART
      key = ipart + icmp * 10 + ifj * 1E3 + ifi * 1E5

      DO I=1,NHDR
         IF(key .EQ. PP_HDR(I)%key) THEN
            PPIHEAD(42) = PP_HDR(I)%stash
            PPIHEAD(23) = PP_HDR(I)%fldcode
            PPIHEAD(32) = PP_HDR(I)%fldtype
            PPRHEAD(64) = PP_HDR(I)%scalefac

            IF(PP_HDR(I)%packing .NE. 0) THEN
               PPIHEAD(21) = 1
               PPRHEAD(51) = REAL(PP_HDR(I)%packing)
            ELSE 
               PPIHEAD(21) = 0
               PPRHEAD(51) = 0.0
            ENDIF

            RETURN
         ENDIF
      ENDDO

      WRITE(NDSE,*) 'WARNING STASH NOT FOUND FOR KEY', key
      WRITE(NDSE,*) 'WILL SET STASH CODE TO BE EQUAL TO KEY'
      PPIHEAD(42) = key
      PPIHEAD(23) = 0
      PPIHEAD(32) = 0
      PPRHEAD(64) = 1.0
      PPIHEAD(21) = 0
      PPRHEAD(51) = 0.0

      RETURN

      END SUBROUTINE PP_LOOKUP


      END MODULE PPLOOKMD


!  --------------------------------------------------
!        1                          Forcing Fields
!   -------------------------------------------------
!  T  T  1     1   DW         DPT   Water depth.
!  T  T  1     2   C[X,Y]     CUR   Current velocity.
!  T  T  1     3   UA         WND   Wind speed.
!  T  T  1     4   AS         AST   Air-sea temperature difference.
!  T  T  1     5   WLV        WLV   Water levels.
!  T  T  1     6   ICE        ICE   Ice concentration.
!  T  T  1     7   IBG        IBG   Iceberg-induced damping
!  T  T  1     8   D50        D50   Median sediment grain size
!   -------------------------------------------------
!        2                          Standard mean wave Parameters
!   -------------------------------------------------
!  T  T  2     1   HS         HS    Wave height.
!  T  T  2     2   WLM        LM    Mean wave length.
!  T  T  2     3   T02        T02   Mean wave period (Tm02).
!  T  T  2     4   T0M1       T0M1  Mean wave period (Tm0,-1).
!  T  T  2     5   T01        T01   Mean wave period (Tm01).
!  T  T  2     6   FP0        FP    Peak frequency.
!  T  T  2     7   THM        DIR   Mean wave direction.
!  T  T  2     8   THS        SPR   Mean directional spread.
!  T  T  2     9   THP0       DP    Peak direction.
!  T  T  2    10   HIG        HIG   Infragravity height
!   -------------------------------------------------
!        3                          Spectral Parameters (first 5)
!   -------------------------------------------------
!  F  F  3     1   Ef         EF    Wave frequency spectrum
!  F  F  3     2   th1m       TH1M  Mean wave direction from a1,b2
!  F  F  3     3   sth1m      STH1M Directional spreading from a1,b2
!  F  F  3     4   th2m             Mean wave direction from a2,b2
!  F  F  3     5   sth2m            Directional spreading from a2,b2
!  F  F  3     6   WN         WN    Wavenumber array
!   -------------------------------------------------
!        4                          Spectral Partition Parameters 
!   -------------------------------------------------
!  T  T  4     1   PHS        PHS   Partitioned wave heights.
!  T  T  4     2   PTP        PTP   Partitioned peak period.
!  T  T  4     3   PLP        PLP   Partitioned peak wave length.
!  T  T  4     4   PDIR       PDIR  Partitioned mean direction.
!  T  T  4     5   PSI        PSPR  Partitioned mean directional spread.
!  T  T  4     6   PT01       PT01  Partitioned mean period (Tm01).
!  T  T  4     7   PT02       PT02  Partitioned mean period (Tm02).
!  T  T  4     8   PTM1       PTM1  Partitioned mean period (Tm-10).
!  T  T  4     9   PWS        PWS   Partitionned wind sea fraction.
!  T  T  4    10   PWST       TWS   Total wind sea fraction.
!  T  T  4    11   PNR        PNR   Number of partitions.
!   -------------------------------------------------
!        5                          Atmosphere-waves layer
!   -------------------------------------------------
!  T  T  5     1   UST        UST   Friction velocity.
!  F  T  5     2   CHARN      CHA   Charnock parameter
!  F  T  5     3   CGE        CGE   Energy flux
!  F  T  5     4   PHIAW      FAW   Air-sea energy flux
!  F  T  5     5   TAUWI[X,Y] TAW   Net wave-supported stress
!  F  T  5     6   TAUWN[X,Y] TWA   Negative part of the wave-supported stress
!  F  F  5     7   WHITECAP   WCC   Whitecap coverage
!  F  F  5     8   WHITECAP   WCF   Whitecap thickness
!  F  F  5     9   WHITECAP   WCH   Mean breaking height
!  F  F  5    10   WHITECAP   WCM   Whitecap moment
!   -------------------------------------------------
!        6                          Wave-ocean layer 
!   -------------------------------------------------
!  F  F  6     1   S[XX,YY,XY] SXY  Radiation stresses.
!  F  F  6     2   TAUO[X,Y]  TWO   Wave to ocean momentum flux
!  F  F  6     3   BHD        BHD   Bernoulli head (J term) 
!  F  F  6     4   PHIOC      FOC   Wave to ocean energy flux
!  F  F  6     5   TUS[X,Y]   TUS   Stokes transport
!  F  F  6     6   USS[X,Y]   USS   Surface Stokes drift
!  F  F  6     7   [PR,TP]MS  P2S   Second-order sum pressure 
!  F  F  6     8   US3D       USF   Spectrum of surface Stokes drift
!  F  F  6     9   P2SMS      P2L   Micro seism  source term
!   -------------------------------------------------
!        7                          Wave-bottom layer 
!   -------------------------------------------------
!  F  F  7     1   ABA        ABR   Near bottom rms amplitides.
!  F  F  7     2   UBA        UBR   Near bottom rms velocities.
!  F  F  7     3   BEDFORMS   BED   Bedforms
!  F  F  7     4   PHIBBL     FBB   Energy flux due to bottom friction 
!  F  F  7     5   TAUBBL     TBB   Momentum flux due to bottom friction
!   -------------------------------------------------
!        8                          Spectrum parameters
!   -------------------------------------------------
!  F  F  8     1   MSS[X,Y]   MSS   Mean square slopes
!  F  F  8     2   MSC[X,Y]   MSC   Spectral level at high frequency tail
!   -------------------------------------------------
!        9                          Numerical diagnostics  
!   -------------------------------------------------
!  T  T  9     1   DTDYN      DTD   Average time step in integration.
!  T  T  9     2   FCUT       FC    Cut-off frequency.
!  T  T  9     3   CFLXYMAX   CFX   Max. CFL number for spatial advection. 
!  T  T  9     4   CFLTHMAX   CFD   Max. CFL number for theta-advection. 
!  F  F  9     5   CFLKMAX    CFK   Max. CFL number for k-advection. 
!   -------------------------------------------------
!        10                         User defined          
!   -------------------------------------------------
!  F  F  10    1              U1    User defined #1. (requires coding ...)
!  F  F  10    2              U2    User defined #1. (requires coding ...)
!   -------------------------------------------------
