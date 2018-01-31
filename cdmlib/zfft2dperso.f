C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2011, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@cs.tsukuba.ac.jp
C
C
C     2-D COMPLEX FFT ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     CALL ZFFT2D(A,NX,NY,IOPT)
C
C     A(NX,NY) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
C       ------------------------------------
C         NX = (2**IP) * (3**IQ) * (5**IR)
C         NY = (2**JP) * (3**JQ) * (5**JR)
C       ------------------------------------
C     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
C          = -1 FOR FORWARD TRANSFORM
C          = +1 FOR INVERSE TRANSFORM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE ZFFT2D(A,NX,NY,IOPT)
      IMPLICIT REAL*8 (A-H,O-Z)
C The maximum supported number of processors is 65536.
      PARAMETER (MAXNPU=65536)
C The maximum supported 2-D transform length is 65536.
      PARAMETER (NDA2=65536)
C The maximum supported 3-D transform length is 4096.
      PARAMETER (NDA3=4096)
C The parameter NBLK is a blocking parameter.
      PARAMETER (NBLK=16)
C The parameter NB is a blocking parameter for NVIDIA GPUs.
      parameter (NB=128)
C The parameter NP is a padding parameter to avoid cache conflicts in
C the FFT routines.
      PARAMETER (NP=8)
C Size of L2 cache
      PARAMETER (L2SIZE=2097152)
      COMPLEX*16 A(*)
      COMPLEX*16 B((NDA2+NP)*NBLK),C(NDA2)
      COMPLEX*16 WX(NDA2),WY(NDA2)
      DIMENSION LNX(3),LNY(3)
      DATA NX0,NY0/0,0/
      SAVE NX0,NY0,WX,WY
C
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
C
      IF (NX .NE. NX0) THEN
         CALL SETTBL(WX,NX)
         NX0=NX
      END IF
      IF (NY .NE. NY0) THEN
         CALL SETTBL(WY,NY)
         NY0=NY
      END IF
C     
      IF (IOPT .EQ. 1) THEN
!$OMP PARALLEL DO
!DIR$ VECTOR ALIGNED
        DO 10 I=1,NX*NY
          A(I)=DCONJG(A(I))
   10   CONTINUE
      END IF
C
!$OMP PARALLEL PRIVATE(B,C)
      CALL ZFFT2D0(A,B,C,WX,WY,NX,NY,LNX,LNY)
!$OMP END PARALLEL
C
      IF (IOPT .EQ. 1) THEN
!$OMP PARALLEL DO
!DIR$ VECTOR ALIGNED
        DO 20 I=1,NX*NY
          A(I)=DCONJG(A(I))
   20   CONTINUE
      ELSE
         DN=1.0D0/DBLE(NX*NY)
         DO 30 I=1,NX*NY
            A(I)=A(I)*DN
 30      CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE ZFFT2D0(A,B,C,WX,WY,NX,NY,LNX,LNY)
      IMPLICIT REAL*8 (A-H,O-Z)
C The maximum supported number of processors is 65536.
      PARAMETER (MAXNPU=65536)
C The maximum supported 2-D transform length is 65536.
      PARAMETER (NDA2=65536)
C The maximum supported 3-D transform length is 4096.
      PARAMETER (NDA3=4096)
C The parameter NBLK is a blocking parameter.
      PARAMETER (NBLK=16)
C The parameter NB is a blocking parameter for NVIDIA GPUs.
      parameter (NB=128)
C The parameter NP is a padding parameter to avoid cache conflicts in
C the FFT routines.
      PARAMETER (NP=8)
C Size of L2 cache
      PARAMETER (L2SIZE=2097152)
      COMPLEX*16 A(NX,*)
      COMPLEX*16 B(NY+NP,*),C(*)
      COMPLEX*16 WX(*),WY(*)
      DIMENSION LNX(*),LNY(*)
C
!$OMP DO
      DO 70 II=1,NX,NBLK
        DO 30 JJ=1,NY,NBLK
          DO 20 I=II,MIN0(II+NBLK-1,NX)
!DIR$ VECTOR ALIGNED
            DO 10 J=JJ,MIN0(JJ+NBLK-1,NY)
              B(J,I-II+1)=A(I,J)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
        DO 40 I=II,MIN0(II+NBLK-1,NX)
          CALL FFT235(B(1,I-II+1),C,WY,NY,LNY)
   40   CONTINUE
        DO 60 J=1,NY
!DIR$ VECTOR ALIGNED
          DO 50 I=II,MIN0(II+NBLK-1,NX)
            A(I,J)=B(J,I-II+1)
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE
!$OMP DO
      DO 80 J=1,NY
        CALL FFT235(A(1,J),C,WX,NX,LNX)
   80 CONTINUE
      RETURN
      END
