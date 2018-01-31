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
C     3-D COMPLEX FFT ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     CALL ZFFT3D(A,NX,NY,NZ,IOPT)
C
C     A(NX,NY,NZ) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
C     NZ IS THE LENGTH OF THE TRANSFORMS IN THE Z-DIRECTION (INTEGER*4)
C       ------------------------------------
C         NX = (2**IP) * (3**IQ) * (5**IR)
C         NY = (2**JP) * (3**JQ) * (5**JR)
C         NZ = (2**KP) * (3**KQ) * (5**KR)
C       ------------------------------------
C     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
C          = -1 FOR FORWARD TRANSFORM
C          = +1 FOR INVERSE TRANSFORM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE ZFFT3D(A,NX,NY,NZ,IOPT)
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
      COMPLEX*16 B((NDA3+NP)*NBLK),C(NDA3)
      COMPLEX*16 WX(NDA3),WY(NDA3),WZ(NDA3)
      DIMENSION LNX(3),LNY(3),LNZ(3)
      DATA NX0,NY0,NZ0/0,0,0/
      SAVE NX0,NY0,NZ0,WX,WY,WZ
C
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
      CALL FACTOR(NZ,LNZ)
C
      IF (NX .NE. NX0) THEN
        CALL SETTBL(WX,NX)
        NX0=NX
      END IF
      IF (NY .NE. NY0) THEN
        CALL SETTBL(WY,NY)
        NY0=NY
      END IF
      IF (NZ .NE. NZ0) THEN
        CALL SETTBL(WZ,NZ)
        NZ0=NZ
      END IF
C
      IF (IOPT .EQ. 1) THEN
!$OMP PARALLEL DO
!DIR$ VECTOR ALIGNED
        DO 10 I=1,NX*NY*NZ
          A(I)=DCONJG(A(I))
   10   CONTINUE
      END IF
C
!$OMP PARALLEL PRIVATE(B,C)
      CALL ZFFT3D0(A,B,B,C,WX,WY,WZ,NX,NY,NZ,LNX,LNY,LNZ)
!$OMP END PARALLEL
C
      IF (IOPT .EQ. 1) THEN     
!$OMP PARALLEL DO
!DIR$ VECTOR ALIGNED
        DO 20 I=1,NX*NY*NZ
          A(I)=DCONJG(A(I))
   20   CONTINUE
      ELSE
         DN=1.0D0/DBLE(NX*NY*NZ)
         DO 30 I=1,NX*NY*NZ
            A(I)=A(I)*DN
 30      CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE ZFFT3D0(A,BY,BZ,C,WX,WY,WZ,NX,NY,NZ,LNX,LNY,LNZ)
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
      COMPLEX*16 A(NX,NY,*)
      COMPLEX*16 BY(NY+NP,*),BZ(NZ+NP,*),C(*)
      COMPLEX*16 WX(*),WY(*),WZ(*)
      DIMENSION LNX(*),LNY(*),LNZ(*)
C
!$OMP DO
      DO 80 J=1,NY
        DO 70 II=1,NX,NBLK
          DO 30 KK=1,NZ,NBLK
            DO 20 I=II,MIN0(II+NBLK-1,NX)
!DIR$ VECTOR ALIGNED
              DO 10 K=KK,MIN0(KK+NBLK-1,NZ)
                BZ(K,I-II+1)=A(I,J,K)
   10         CONTINUE
   20       CONTINUE
   30     CONTINUE
          DO 40 I=II,MIN0(II+NBLK-1,NX)
            CALL FFT235(BZ(1,I-II+1),C,WZ,NZ,LNZ)
   40     CONTINUE
          DO 60 K=1,NZ
!DIR$ VECTOR ALIGNED
            DO 50 I=II,MIN0(II+NBLK-1,NX)
              A(I,J,K)=BZ(K,I-II+1)
   50       CONTINUE
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE
!$OMP DO
      DO 170 K=1,NZ
        DO 150 II=1,NX,NBLK
          DO 110 JJ=1,NY,NBLK
            DO 100 I=II,MIN0(II+NBLK-1,NX)
!DIR$ VECTOR ALIGNED
              DO 90 J=JJ,MIN0(JJ+NBLK-1,NY)
                BY(J,I-II+1)=A(I,J,K)
   90         CONTINUE
  100       CONTINUE
  110     CONTINUE
          DO 120 I=II,MIN0(II+NBLK-1,NX)
            CALL FFT235(BY(1,I-II+1),C,WY,NY,LNY)
  120     CONTINUE
          DO 140 J=1,NY
!DIR$ VECTOR ALIGNED
            DO 130 I=II,MIN0(II+NBLK-1,NX)
              A(I,J,K)=BY(J,I-II+1)
  130       CONTINUE
  140     CONTINUE
  150   CONTINUE
        DO 160 J=1,NY
          CALL FFT235(A(1,J,K),C,WX,NX,LNX)
  160   CONTINUE
  170 CONTINUE
      RETURN
      END
