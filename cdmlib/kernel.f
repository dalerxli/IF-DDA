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
C     RADIX-2, 3, 4, 5 AND 8 FFT KERNEL ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE FFT2(A,B,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(M,*),B(M,*)
      COMPLEX*16 C0,C1
C
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        C0=A(I,1)
        C1=A(I,2)
        B(I,1)=C0+C1
        B(I,2)=C0-C1
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT3A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(L,*),B(3,*),W(2,*)
      COMPLEX*16 C0,C1,C2,D0,D1,D2,W1,W2
      DATA C31/0.86602540378443865D0/C32/0.5D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        W1=W(1,J)
        W2=W(2,J)
        C0=A(J,1)
        C1=A(J,2)
        C2=A(J,3)
        D0=C1+C2
        D1=C0-C32*D0
        D2=(0.0D0,-1.0D0)*C31*(C1-C2)
        B(1,J)=C0+D0
        B(2,J)=W1*(D1+D2)
        B(3,J)=W2*(D1-D2)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT3B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(M,L,*),B(M,3,*),W(2,*)
      COMPLEX*16 C0,C1,C2,D0,D1,D2,W1,W2
      DATA C31/0.86602540378443865D0/C32/0.5D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        C0=A(I,1,1)
        C1=A(I,1,2)
        C2=A(I,1,3)
        D0=C1+C2
        D1=C0-C32*D0
        D2=(0.0D0,-1.0D0)*C31*(C1-C2)
        B(I,1,1)=C0+D0
        B(I,2,1)=D1+D2
        B(I,3,1)=D1-D2
   10 CONTINUE
      DO 30 J=2,L
        W1=W(1,J)
        W2=W(2,J)
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          C0=A(I,J,1)
          C1=A(I,J,2)
          C2=A(I,J,3)
          D0=C1+C2
          D1=C0-C32*D0
          D2=(0.0D0,-1.0D0)*C31*(C1-C2)
          B(I,1,J)=C0+D0
          B(I,2,J)=W1*(D1+D2)
          B(I,3,J)=W2*(D1-D2)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT4A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(L,*),B(4,*),W(3,*)
      COMPLEX*16 C0,C1,C2,C3,D0,D1,D2,D3,W1,W2,W3
C
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
        C0=A(J,1)
        C1=A(J,2)
        C2=A(J,3)
        C3=A(J,4)
        D0=C0+C2
        D1=C0-C2
        D2=C1+C3
        D3=(0.0D0,-1.0D0)*(C1-C3)
        B(1,J)=D0+D2
        B(2,J)=W1*(D1+D3)
        B(3,J)=W2*(D0-D2)
        B(4,J)=W3*(D1-D3)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT4B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(M,L,*),B(M,4,*),W(3,*)
      COMPLEX*16 C0,C1,C2,C3,D0,D1,D2,D3,W1,W2,W3
C
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        C0=A(I,1,1)
        C1=A(I,1,2)
        C2=A(I,1,3)
        C3=A(I,1,4)
        D0=C0+C2
        D1=C0-C2
        D2=C1+C3
        D3=(0.0D0,-1.0D0)*(C1-C3)
        B(I,1,1)=D0+D2
        B(I,2,1)=D1+D3
        B(I,3,1)=D0-D2
        B(I,4,1)=D1-D3
   10 CONTINUE
      DO 30 J=2,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          C0=A(I,J,1)
          C1=A(I,J,2)
          C2=A(I,J,3)
          C3=A(I,J,4)
          D0=C0+C2
          D1=C0-C2
          D2=C1+C3
          D3=(0.0D0,-1.0D0)*(C1-C3)
          B(I,1,J)=D0+D2
          B(I,2,J)=W1*(D1+D3)
          B(I,3,J)=W2*(D0-D2)
          B(I,4,J)=W3*(D1-D3)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT5A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(L,*),B(5,*),W(4,*)
      COMPLEX*16 C0,C1,C2,C3,C4,D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10
      COMPLEX*16 W1,W2,W3,W4
      DATA C51/0.95105651629515357D0/C52/0.61803398874989485D0/
     1     C53/0.55901699437494742D0/C54/0.25D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
        W4=W(4,J)
        C0=A(J,1)
        C1=A(J,2)
        C2=A(J,3)
        C3=A(J,4)
        C4=A(J,5)
        D0=C1+C4
        D1=C2+C3
        D2=C51*(C1-C4)
        D3=C51*(C2-C3)
        D4=D0+D1
        D5=C53*(D0-D1)
        D6=C0-C54*D4
        D7=D6+D5
        D8=D6-D5
        D9=(0.0D0,-1.0D0)*(D2+C52*D3)
        D10=(0.0D0,-1.0D0)*(C52*D2-D3)
        B(1,J)=C0+D4
        B(2,J)=W1*(D7+D9)
        B(3,J)=W2*(D8+D10)
        B(4,J)=W3*(D8-D10)
        B(5,J)=W4*(D7-D9)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT5B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(M,L,*),B(M,5,*),W(4,*)
      COMPLEX*16 C0,C1,C2,C3,C4,D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10
      COMPLEX*16 W1,W2,W3,W4
      DATA C51/0.95105651629515357D0/C52/0.61803398874989485D0/
     1     C53/0.55901699437494742D0/C54/0.25D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        C0=A(I,1,1)
        C1=A(I,1,2)
        C2=A(I,1,3)
        C3=A(I,1,4)
        C4=A(I,1,5)
        D0=C1+C4
        D1=C2+C3
        D2=C51*(C1-C4)
        D3=C51*(C2-C3)
        D4=D0+D1
        D5=C53*(D0-D1)
        D6=C0-C54*D4
        D7=D6+D5
        D8=D6-D5
        D9=(0.0D0,-1.0D0)*(D2+C52*D3)
        D10=(0.0D0,-1.0D0)*(C52*D2-D3)
        B(I,1,1)=C0+D4
        B(I,2,1)=D7+D9
        B(I,3,1)=D8+D10
        B(I,4,1)=D8-D10
        B(I,5,1)=D7-D9
   10 CONTINUE
      DO 30 J=2,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
        W4=W(4,J)
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          C0=A(I,J,1)
          C1=A(I,J,2)
          C2=A(I,J,3)
          C3=A(I,J,4)
          C4=A(I,J,5)
          D0=C1+C4
          D1=C2+C3
          D2=C51*(C1-C4)
          D3=C51*(C2-C3)
          D4=D0+D1
          D5=C53*(D0-D1)
          D6=C0-C54*D4
          D7=D6+D5
          D8=D6-D5
          D9=(0.0D0,-1.0D0)*(D2+C52*D3)
          D10=(0.0D0,-1.0D0)*(C52*D2-D3)
          B(I,1,J)=C0+D4
          B(I,2,J)=W1*(D7+D9)
          B(I,3,J)=W2*(D8+D10)
          B(I,4,J)=W3*(D8-D10)
          B(I,5,J)=W4*(D7-D9)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT8A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(L,*),B(8,*),W(7,*)
      COMPLEX*16 C0,C1,C2,C3,C4,C5,C6,C7,D0,D1,D2,D3,D4,D5,D6,D7
      COMPLEX*16 E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,W1,W2,W3,W4,W5,W6,W7
      DATA C81/0.70710678118654752D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
        W4=W(4,J)
        W5=W(5,J)
        W6=W(6,J)
        W7=W(7,J)
        C0=A(J,1)
        C1=A(J,2)
        C2=A(J,3)
        C3=A(J,4)
        C4=A(J,5)
        C5=A(J,6)
        C6=A(J,7)
        C7=A(J,8)
        D0=C0+C4
        D1=C0-C4
        D2=C2+C6
        D3=(0.0D0,-1.0D0)*(C2-C6)
        D4=C1+C5
        D5=C1-C5
        D6=C3+C7
        D7=C3-C7
        E0=D0+D2
        E1=D0-D2
        E2=D4+D6
        E3=(0.0D0,-1.0D0)*(D4-D6)
        E4=C81*(D5-D7)
        E5=(0.0D0,-1.0D0)*C81*(D5+D7)
        E6=D1+E4
        E7=D1-E4
        E8=D3+E5
        E9=D3-E5
        B(1,J)=E0+E2
        B(2,J)=W1*(E6+E8)
        B(3,J)=W2*(E1+E3)
        B(4,J)=W3*(E7-E9)
        B(5,J)=W4*(E0-E2)
        B(6,J)=W5*(E7+E9)
        B(7,J)=W6*(E1-E3)
        B(8,J)=W7*(E6-E8)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT8B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(M,L,*),B(M,8,*),W(7,*)
      COMPLEX*16 C0,C1,C2,C3,C4,C5,C6,C7,D0,D1,D2,D3,D4,D5,D6,D7
      COMPLEX*16 E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,W1,W2,W3,W4,W5,W6,W7
      DATA C81/0.70710678118654752D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        C0=A(I,1,1)
        C1=A(I,1,2)
        C2=A(I,1,3)
        C3=A(I,1,4)
        C4=A(I,1,5)
        C5=A(I,1,6)
        C6=A(I,1,7)
        C7=A(I,1,8)
        D0=C0+C4
        D1=C0-C4
        D2=C2+C6
        D3=(0.0D0,-1.0D0)*(C2-C6)
        D4=C1+C5
        D5=C1-C5
        D6=C3+C7
        D7=C3-C7
        E0=D0+D2
        E1=D0-D2
        E2=D4+D6
        E3=(0.0D0,-1.0D0)*(D4-D6)
        E4=C81*(D5-D7)
        E5=(0.0D0,-1.0D0)*C81*(D5+D7)
        E6=D1+E4
        E7=D1-E4
        E8=D3+E5
        E9=D3-E5
        B(I,1,1)=E0+E2
        B(I,2,1)=E6+E8
        B(I,3,1)=E1+E3
        B(I,4,1)=E7-E9
        B(I,5,1)=E0-E2
        B(I,6,1)=E7+E9
        B(I,7,1)=E1-E3
        B(I,8,1)=E6-E8
   10 CONTINUE
      DO 30 J=2,L
        W1=W(1,J)
        W2=W(2,J)
        W3=W(3,J)
        W4=W(4,J)
        W5=W(5,J)
        W6=W(6,J)
        W7=W(7,J)
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          C0=A(I,J,1)
          C1=A(I,J,2)
          C2=A(I,J,3)
          C3=A(I,J,4)
          C4=A(I,J,5)
          C5=A(I,J,6)
          C6=A(I,J,7)
          C7=A(I,J,8)
          D0=C0+C4
          D1=C0-C4
          D2=C2+C6
          D3=(0.0D0,-1.0D0)*(C2-C6)
          D4=C1+C5
          D5=C1-C5
          D6=C3+C7
          D7=C3-C7
          E0=D0+D2
          E1=D0-D2
          E2=D4+D6
          E3=(0.0D0,-1.0D0)*(D4-D6)
          E4=C81*(D5-D7)
          E5=(0.0D0,-1.0D0)*C81*(D5+D7)
          E6=D1+E4
          E7=D1-E4
          E8=D3+E5
          E9=D3-E5
          B(I,1,J)=E0+E2
          B(I,2,J)=W1*(E6+E8)
          B(I,3,J)=W2*(E1+E3)
          B(I,4,J)=W3*(E7-E9)
          B(I,5,J)=W4*(E0-E2)
          B(I,6,J)=W5*(E7+E9)
          B(I,7,J)=W6*(E1-E3)
          B(I,8,J)=W7*(E6-E8)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@cs.tsukuba.ac.jp
C
C
C     RADIX-2, 3, 4, 5 AND 8 FFT ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE FFT235(A,B,W,N,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
      DIMENSION IP(*)
C
      IF (IP(1) .NE. 1) THEN
        KP4=2-MOD(IP(1)+2,3)
        KP8=(IP(1)-KP4)/3
      ELSE
        KP4=0
        KP8=0
      END IF
C
      KEY=1
      J=1
      L=N
      M=1
      DO 10 K=1,KP8
        L=L/8
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT8(A,B,W(J),M,L)
          ELSE
            CALL FFT8(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT8(A,A,W(J),M,L)
          ELSE
            CALL FFT8(B,A,W(J),M,L)
          END IF
        END IF
        M=M*8
        J=J+L*7
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT5(A,B,W(J),M,L)
          ELSE
            CALL FFT5(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT5(A,A,W(J),M,L)
          ELSE
            CALL FFT5(B,A,W(J),M,L)
          END IF
        END IF
        M=M*5
        J=J+L*4
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT4(A,B,W(J),M,L)
          ELSE
            CALL FFT4(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT4(A,A,W(J),M,L)
          ELSE
            CALL FFT4(B,A,W(J),M,L)
          END IF
        END IF
        M=M*4
        J=J+L*3
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT3(A,B,W(J),M,L)
          ELSE
            CALL FFT3(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT3(A,A,W(J),M,L)
          ELSE
            CALL FFT3(B,A,W(J),M,L)
          END IF
        END IF
        M=M*3
        J=J+L*2
   40 CONTINUE
      IF (IP(1) .EQ. 1) THEN
        IF (KEY .GE. 0) THEN
          CALL FFT2(A,A,M)
        ELSE
          CALL FFT2(B,A,M)
        END IF
      END IF
      RETURN
      END
      SUBROUTINE FFT3(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT3A(A,B,W,L)
      ELSE
        CALL FFT3B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT4(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT4A(A,B,W,L)
      ELSE
        CALL FFT4B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT5(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT5A(A,B,W,L)
      ELSE
        CALL FFT5B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT8(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT8A(A,B,W,L)
      ELSE
        CALL FFT8B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE SETTBL(W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 W(*)
      DIMENSION IP(3)
C
      CALL FACTOR(N,IP)
C
      IF (IP(1) .NE. 1) THEN
        KP4=2-MOD(IP(1)+2,3)
        KP8=(IP(1)-KP4)/3
      ELSE
        KP4=0
        KP8=0
      END IF
C
      J=1
      L=N
      DO 10 K=1,KP8
        L=L/8
        CALL SETTBL0(W(J),8,L)
        J=J+L*7
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        CALL SETTBL0(W(J),5,L)
        J=J+L*4
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        CALL SETTBL0(W(J),4,L)
        J=J+L*3
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        CALL SETTBL0(W(J),3,L)
        J=J+L*2
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL0(W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 W(M-1,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(M)*DBLE(L))
      DO 20 J=1,L
!DIR$ VECTOR ALIGNED
        DO 10 I=1,M-1
          TEMP=PX*DBLE(I)*DBLE(J-1)
          W(I,J)=DCMPLX(DCOS(TEMP),DSIN(TEMP))
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL2(W,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 W(NX,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(NX)*DBLE(NY))
!$OMP PARALLEL DO PRIVATE(TEMP)
      DO 20 J=1,NY
!DIR$ VECTOR ALIGNED
        DO 10 I=1,NX
          TEMP=PX*DBLE(I-1)*DBLE(J-1)
          W(I,J)=DCMPLX(DCOS(TEMP),DSIN(TEMP))
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@cs.tsukuba.ac.jp
C
C
C     FACTORIZATION ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE FACTOR(N,IP)
      DIMENSION IP(*)
C
      IP(1)=0
      IP(2)=0
      IP(3)=0
      N2=N
      IF (MOD(N,2) .NE. 0 .AND. MOD(N,3) .NE. 0 .AND.
     1    MOD(N,5) .NE. 0) RETURN
   10 IF (N2 .LE. 1) RETURN
      IF (MOD(N2,2) .EQ. 0) THEN
        IP(1)=IP(1)+1
        N2=N2/2
        GO TO 10
      ELSE IF (MOD(N2,3) .EQ. 0) THEN
        IP(2)=IP(2)+1
        N2=N2/3
        GO TO 10
      ELSE IF (MOD(N2,5) .EQ. 0) THEN
        IP(3)=IP(3)+1
        N2=N2/5
        GO TO 10
      END IF
      RETURN
      END
      SUBROUTINE GETNXNY(N,NX,NY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IP(3),LNX(3),LNY(3)
C
      ISQRTN=IDINT(DSQRT(DBLE(N)))
      CALL FACTOR(N,IP)
      DO 10 I=1,3
        LNX(I)=0
   10 CONTINUE
      IRES=ISQRTN
      DO 40 K=0,(IP(3)+1)/2
        DO 30 J=0,(IP(2)+1)/2
          DO 20 I=0,(IP(1)+1)/2
            NX=(2**I)*(3**J)*(5**K)
            IF (NX .LE. ISQRTN) THEN
              IRES2=ISQRTN-NX
              IF (IRES2 .LT. IRES) THEN
                LNX(1)=I
                LNX(2)=J
                LNX(3)=K
                IRES=IRES2
              END IF
            END IF
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE
      DO 50 I=1,3
        LNY(I)=IP(I)-LNX(I)
   50 CONTINUE
      NX=(2**LNX(1))*(3**LNX(2))*(5**LNX(3))
      NY=(2**LNY(1))*(3**LNY(2))*(5**LNY(3))
      RETURN
      END
