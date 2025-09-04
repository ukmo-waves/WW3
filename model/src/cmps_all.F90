! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   SUBROUTINE COEX,COEX2,CMPS,XPND,INSTIN,EXTRIN -----------------
!LL
!LL   PURPOSE:   PACK TO AND UNPACK FROM WGDOS FORMAT
!LL
!LL   (Note that for optimal performance the routines INSTIN and
!LL    EXTRIN are inline expanded (enable option 8 on fpp)
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2  19/04/93  Code for ne real missing data indicator (TCJ).
!LL   4.2  Nov. 96   T3E migration: WHENNE,WHENEQ replaced by
!LL                   portable fortran code.   S.J.Swarbrick
!LL   4.3  Apr. 97 T3E migration: Calling the functions CRI2IBM/IBM2CRI
!LL                directly replaces the calls to the (now unsupported)
!LL                routines USICTI, USICTC, USSCTC, USSCTI. Ian Edmond
!LL   5.0 05/03/99 Remove DEF,C98_1A. D. Robinson.
!LL   5.1  23/03/00 Add declaration for CRI2IBM. S.D.Mullerworth
!LL   5.2  15/11/00 Allow use in reconfiguration. P.Selwood.
!LL
!LL   5.1 May. 00  The max number of bits per value used in the wgdos
!LL                packed data is 31. This implies that the maximum
!LL                difference from the base value which can be stored
!LL                is 2.147483E+09. A check is now added to prevent
!LL                values greater than this being stored.  Ian Edmond
!LL   5.3  07/06/01 Declare cri2ibm.  A van der Wal
!     5.5  14/02/03 Add Implicit None + appropriate declarations where
!                   necessary. Tidy up some declarations. T.White
!LL   5.5  28/02/03 Insert code for portable data conversion routines
!LL                 to replace Cray-specific CRI2IBM etc.     P.Dando
!     6.0  10/09/03 Conversion of portable data conversion routines
!                   (IEEE2IBM etc) into functions with error return
!                   codes matching those of CRAY routines.   P.Dando
!     6.0  25/09/03 Optimisation for NEC.
!                   Bob Carruthers, Gerrit v.d. Velde, Paul Selwood.
!     6.1  23/08/04 Add FLDCALC def. D Robinson
!     6.1  13/09/04 Fix bug with insertion of zero and missing data
!                   bit maps (filling to end of next 32-bit word)
!                   in CMPS routine when row_length is a multiple of
!                   64.                                    P.Dando
!     6.2  26/10/05 Added function ishift which behaves in a standard
!                   way (as expected by coex) to fix WGDOS errors in
!                   ifc compiled runs.                     T. Edwards
!     6.2  24/11/05 Optimisation for NEC. Two new subroutines added -
!                   CMPS_ALL and XPND_ALL. J-C.Rioual(NEC)/D.Robinson.
!LL
!LL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
!LL   STANDARD B, VERSION 2, DATED 18/01/90
!LL
!LL  Logical component number: S72
!LL
!LL   SYSTEM TASK: P7
!LL   OWNED BY P J SMITH
!LL
!LL   DOCUMENTATION:  ??????
!LLEND-------------------------------------------------------------





!     This is a portable replacement for ishft, which behaves as the UM
!     expects it too when shifting more than 64 bits.





      SUBROUTINE CMPS_ALL(FIELD,ICOMP64,N,IX,IY,NUM,ISC,RMDI,           &
     &                    ICODE,CMESSAGE)

!     Better vectorized version of CMPS compressing the whole field
!     at once.
      IMPLICIT NONE

!     Subroutine arguments

      INTEGER(KIND=4), INTENT(IN)  :: n, ix, iy, isc
      INTEGER(KIND=4), INTENT(OUT) :: num
      REAL(KIND=4),    INTENT(IN)  :: field(ix,iy)
      REAL(KIND=4),    INTENT(IN)  :: rmdi
      INTEGER(KIND=8), INTENT(OUT) :: icomp64(n)
      INTEGER(KIND=4)              :: icode
      CHARACTER            :: cmessage*80

!     Local variables

      INTEGER(KIND=8) :: i, j, i2, iword, istart, npoint, nshft, ival
      INTEGER :: isig, iexp, iman
      INTEGER :: is1, is2, is3
      INTEGER :: nbits_bmap, nwords_bmap, nwords_data
      INTEGER :: nbits_pack, nvals_pack

      REAL(KIND=8)    :: aprec, bprec

      LOGICAL :: obtmis, obtzer

      INTEGER(KIND=8), DIMENSION (iy*(2*ix+2)+4) :: icomp  ! should be plenty
      INTEGER(KIND=8), DIMENSION (iy)      :: ibase, imax, ibit
      INTEGER(KIND=8), DIMENSION (iy)      :: nmiss, nzero, ibm
      INTEGER(KIND=8), DIMENSION (2*ix+32) :: itmp

      REAL(KIND=8),    DIMENSION(ix)       :: atmp
      REAL(KIND=8),    DIMENSION(iy)       :: base, fmax

      INTEGER(KIND=8), PARAMETER :: MASK32 = Z'FFFFFFFF'
      INTEGER(KIND=8), PARAMETER :: MASK_EXPT_IBM = Z'7F000000'
      INTEGER(KIND=8), PARAMETER :: MASK_SIGN_IBM = Z'80000000'

! GENERAL REMARK:
! All lengths and alignments in WGDOS packing are for 32-bit words,
! so life gets much easier when we treat the packed data as 32-bit
! words.
! So we gather all compressed data in the low 32 bits of icomp
! and compress this array at the end to 64 bits

! Scale factor

      aprec = 2.**isc
      bprec = 1./aprec

!     Find minimum and maximum value for every row,
!     count number of missing and zero numbers

      base(:)  = 999999.
      fmax(:)  = -HUGE(0.0)
      nmiss(:) = 0
      nzero(:) = 0

      DO j=1,iy
        DO i=1,ix
          IF(field(i,j)/=rmdi) THEN
            base(j) = MIN(base(j),field(i,j))
            fmax(j) = MAX(fmax(j),field(i,j))
            IF(field(i,j)==0) nzero(j) = nzero(j)+1
          ELSE
            nmiss(j) = nmiss(j)+1
          ENDIF
        ENDDO
      ENDDO

!     ROUND BASE TO PRECISION REQUIRED

      DO j=1,iy
        ibase(j) = nint(base(j)*bprec)
        base(j)  = ibase(j)*aprec
      ENDDO

!     IBM floating point representation of base:

      DO j=1,iy
        IF(base(j)==0) THEN
          ibm(j) = 0
        ELSE
          isig = sign(1.0_8,base(j))
          base(j) = abs(base(j))
          iexp = exponent(base(j)) + 256
          iman = fraction(base(j)) * 16777216.
          i = MOD(iexp,4)
          IF(i==1) iman = ISHFT(iman,-3)
          IF(i==2) iman = ISHFT(iman,-2)
          IF(i==3) iman = ISHFT(iman,-1)
          iexp = (iexp+3)/4
          ibm(j) = IOR(IAND(ISHFT(iexp,24),MASK_EXPT_IBM),iman)
          IF(isig<0) ibm(j)=IOR(ibm(j),MASK_SIGN_IBM)
        ENDIF
      ENDDO

!     Find maximum scaled value

      IMAX(:) = 0
      DO j=1,iy
        IF(nmiss(j)<ix) THEN
          imax(j) = MAX(NINT(fmax(j)*bprec)-ibase(j), 0)
        ENDIF
      ENDDO

!     Check if this may be represented with 31 bits

      IF (ANY(imax(:)>2147483647)) THEN
        ICODE = 2
        CMESSAGE='COEX: Unable to WGDOS pack to this accuracy'
        GOTO 999
      ENDIF

!     FIND NUMBER OF BITS REQUIRED TO STORE MAX DIFFERENCE

      ibit(:) = 0
      i2 = 1
      DO i=1,32
        DO j=1,iy
          IF(imax(j)>=i2) ibit(j) = i
        ENDDO
        i2 = i2*2
      ENDDO

!     Fill data into output array

      iword = 4 ! first 3 words filled at end

      DO j=1,iy

        istart = iword
        iword = iword+2 ! the two header words are inserted at end

        ! Check which bitmaps we need

        obtmis = (nmiss(j)>0)
        ! Check if it is worthwile to use zero-bitmap:
        obtzer = (ibit(j)*nzero(j) > ix)

        ! Set itmp with the bitmap pattern

        nbits_bmap = 0
        IF (obtmis) THEN
          itmp(1:ix) = MERGE(1,0,field(:,j)==rmdi)
          nbits_bmap = ix
        ENDIF

        IF (obtzer) THEN
          itmp(nbits_bmap+1:nbits_bmap+ix) = MERGE(1,0,field(:,j)/=0)
          nbits_bmap = nbits_bmap+ix
        ENDIF

        ! Insert bitmap

        IF(nbits_bmap>0) THEN

          ! add 1's to the end since bitmaps should be padded with 1's

          itmp(nbits_bmap+1:nbits_bmap+31) = 1

          ! The number of words to be used for bitmap

          nwords_bmap = (nbits_bmap+31)/32

          ! Compress itmp

          ! Combine 4 contiguous 1-bit-itmes to one 4-bit-item
          DO i=1,nwords_bmap*8
            itmp(i) = IOR(IOR(ISHFT(itmp(4*i-3),3),                     &
     &                        ISHFT(itmp(4*i-2),2)),                    &
     &                    IOR(ISHFT(itmp(4*i-1),1),itmp(4*i)))
          ENDDO
          ! Combine 4 contiguous 4-bit-itmes to one 16-bit-item
          DO i=1,nwords_bmap*2
            itmp(i) = IOR(IOR(ISHFT(itmp(4*i-3),12),                    &
     &                        ISHFT(itmp(4*i-2), 8)),                   &
     &                    IOR(ISHFT(itmp(4*i-1), 4),itmp(4*i)))
          ENDDO
          ! Combine 2 contiguous 16-bit-itmes to the final destination
          DO i=1,nwords_bmap
            icomp(iword+i-1) = IOR(ISHFT(itmp(2*i-1),16),itmp(2*i))
          ENDDO

          iword = iword + nwords_bmap

        ENDIF

        ! Insert data

        IF(ibit(j)>0) THEN

          ! Get rid of missing values

          IF(obtmis) THEN
            npoint = 0
            DO i=1,ix
              IF(field(i,j)/=rmdi) THEN
                npoint = npoint+1
                atmp(npoint) = field(i,j)
              ENDIF
            ENDDO
          ELSE
            npoint = ix
            atmp(:) = field(:,j)
          ENDIF

          ! Get rid of zero values

          IF(obtzer) THEN
            ival = npoint
            npoint = 0
!CDIR NODEP
            DO i=1,ival
              IF(atmp(i)/=0) THEN
                npoint = npoint+1
                atmp(npoint) = atmp(i)
              ENDIF
            ENDDO
          ENDIF

          ! Number of words used for the compressed data

          nwords_data = (npoint*ibit(j)+31)/32

          ! Scale and find difference from base

          DO i=1,npoint
            itmp(i) = MAX(NINT(atmp(i)*bprec)-ibase(j), 0)
          ENDDO

          ! As long as ibit(j) is <=16 we can combine two contiguous
          ! items to one with the double number of bits, halfing the
          ! number of words to be compressed

          nbits_pack = ibit(j)
          nvals_pack = npoint

          DO WHILE (nbits_pack <= 16)
            itmp(nvals_pack+1) = 0 ! for odd numbers
            DO i=1,(nvals_pack+1)/2
              itmp(i) = IOR(ISHFT(itmp(2*i-1),nbits_pack),itmp(2*i))
            ENDDO
            nbits_pack = 2*nbits_pack
            nvals_pack = (nvals_pack+1)/2
          ENDDO

          IF(nbits_pack == 32) THEN

            ! This is the case if ibit(j) is 1, 2, 4, 8 or 16
            ! We have not much to do, just copy itmp

            DO i=1,nwords_data
              icomp(iword+i-1) = itmp(i)
            ENDDO

          ELSE

            ! Shift every value in itmp to the left and append
            ! the bits of the following 2 words

            is1 = 64-nbits_pack   ! amount to shift itmp(i)
            is2 = 64-2*nbits_pack ! amount to shift itmp(i+1)
            is3 = 64-3*nbits_pack ! amount to shift itmp(i+2)

            itmp(nvals_pack+1) = 0
            itmp(nvals_pack+2) = 0

            DO i=1,nvals_pack
              itmp(i) = IOR(IOR(ISHFT(itmp(i  ),is1),                   &
     &                          ISHFT(itmp(i+1),is2)),                  &
     &                          ISHFT(itmp(i+2),is3))
            ENDDO

            ! Now itmp contains enough data so that we can cut out
            ! the compressed data words

            DO i=1,nwords_data

              ! Word which contains compressed data word:
              ival = itmp(((i-1)*32)/nbits_pack + 1)

              ! Number of bits we have to shift to the left
              ! so that we have the compressed data word left packed:
              nshft = MOD((i-1)*32,nbits_pack)

              ! Normally we could now code:
              !   ival = ISHFT(ival,nshft)
              ! but since vector-shift-by-vector is not
              ! implemented in H/W we have to do:

              if(IAND(nshft,16)/=0) ival = ISHFT(ival,16)
              if(IAND(nshft, 8)/=0) ival = ISHFT(ival, 8)
              if(IAND(nshft, 4)/=0) ival = ISHFT(ival, 4)
              if(IAND(nshft, 2)/=0) ival = ISHFT(ival, 2)
              if(IAND(nshft, 1)/=0) ival = ISHFT(ival, 1)

              ! Shift to the right half and mask out upper bits
              ! (for the case that ISHFT does an arithmetic shift)

              icomp(iword+i-1) = IAND(ISHFT(ival,-32),MASK32)

            ENDDO
          ENDIF

          iword = iword + nwords_data

        ENDIF

        ! Now insert the header for this row:
        ! First word of compressed data: IBM representation of base
        ! Second word of compressed data:
        ! 16 bits: ibit(j) + flags
        ! 16 bits: number of words of data following

        icomp(istart) = ibm(j)
        IF(obtzer) ibit(j) = ibit(j) + 128
        IF(obtmis) ibit(j) = ibit(j) + 32
        icomp(istart+1) = IOR(ISHFT(ibit(j),16),iword-istart-2)

      ENDDO

!     Fill first 3 words of compressed data

      num = iword-1
      icomp(1) = num
      icomp(2) = IAND(isc,MASK32)
      icomp(3) = IOR(ISHFT(ix,16),iy)

!     Compress to 64 bit words

      IF ((num+1)/2 > n) THEN
        ICODE = 2
        CMESSAGE='COEX: Dimension of ICOMP too small'
        GOTO 999
      ENDIF

      ! ChrisB: I get the feeling that we could not do this bit and instead
      ! return the icomp array, which is aleady in 32 bit format....
      icomp(num+1) = 0
      DO i=1,(num+1)/2
        icomp64(i) = IOR(ISHFT(icomp(2*i-1),32),IAND(icomp(2*i),MASK32))
        !icomp64(i) = icomp(i) ! ChrisB: assume icomp64 is 32bit
      ENDDO


 999  CONTINUE
      RETURN

      END SUBROUTINE CMPS_ALL

