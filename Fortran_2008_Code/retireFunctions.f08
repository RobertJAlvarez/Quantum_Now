!!!
!Extra functions for reference but never used
!!!
MODULE retireFunctions
  USE fortranFunctions, ONLY: PI, DBL, FMOD, ABSO, DIV
  USE arrayFunctions, ONLY: J2x2

  PRIVATE
  PUBLIC :: Class_SQR, Class_DIV, Class_DIAGNxN !Only TestingFunctions.f08 use it

  CONTAINS

  REAL(DBL) FUNCTION Class_DIV(AA,BB) !Returns A/B
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: AA, BB

    !Stop if denominator is 0
    IF (ABSO(BB) < 1.D-15) THEN
      STOP "Can't divide by 0"
    END IF

    !Calculate division where BB is always positive
    IF (BB < 0.D0) THEN
      Class_DIV = -AA*RECIPROCAL(-BB)
    ELSE
      Class_DIV = AA*RECIPROCAL(BB)
    END IF
  END FUNCTION Class_DIV

  REAL(DBL) FUNCTION Class_SQR(num) RESULT(xn)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num
    REAL(DBL) :: xn_1, xn_2 !xn_1 = x_(n-1), xn_2 = x_(n-2)
    INTEGER :: i

    xn = 1.D0
    xn_1 = -0.1D0
    xn_2 = -0.2D0

!     X ~=(num)^0.5
!  X + E = (num)^0.5
! (X + E)^2 = num
! X^2 + 2XE + E^2 = num, But E is very small, so E^2=0:
!      E = (num - X^2)/2X, this way:
!     X ~= num^0.5 = X + E -> X =~ X + E

    DO i=1, 100
      xn = xn + DIV(num - xn*xn, 2.D0*xn)
      !Exit loop if values are cycling
      IF(xn == xn_1) THEN
        EXIT
      ELSE IF (xn == xn_2) THEN
        xn = DIV(xn_2 + xn_1, 2.D0)
        EXIT
      END IF
      !Update x_ns
      xn_2 = xn_1 !x_(n-2) = x_(n-1)
      xn_1 = xn   !x_(n-1) = x_n
    END DO
  END FUNCTION Class_SQR

  REAL(DBL) FUNCTION RECIPROCAL(Z)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: Z
    REAL(DBL) :: RECIP

    RECIP = 1.D0 !Only stays as 1.0 if z = 1

    IF (Z > 1.D0)THEN                        !IF Z > 1
      RECIP = GT1DIVIDE(Z - 1.D0)
    ELSE IF(Z > 0.1D0 .AND. Z < 1.D0) THEN   !IF 0.1 < Z < 1
      RECIP = LT1DIVIDE(1.D0 - Z)
    ELSE IF (Z > 0.D0 .AND. Z <= 0.1D0) THEN !IF 0 < Z <= 0.1
      RECIP = TNYDIVIDE(Z)
    END IF

    RECIPROCAL = RECIP
  END FUNCTION RECIPROCAL

  REAL(DBL) FUNCTION DIVIDER(X) !DIVIDER (1/(1+x))   0<x<1
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: X
    REAL(DBL) :: TWOTHIRDS
    REAL(DBL) :: Y

    IF (X < 0.D0) STOP 'X LT 0'
    IF (X > 1.D0) STOP 'X GT 1'

    TWOTHIRDS = 0.666666666666667D0

    IF (X > 0.5D0) THEN
      Y = TWOTHIRDS*(X - 0.5D0)
      DIVIDER = DIVIDE(Y)*TWOTHIRDS
    ELSE
      DIVIDER = DIVIDE(X)
    END IF
  END FUNCTION DIVIDER

  REAL(DBL) FUNCTION DIVIDE(X)  !DIVIDE (1/(1+x))   0<x<0.5
    IMPLICIT NONE

    REAL(DBL), INTENT(IN) :: X
    INTEGER :: i, N = 50
    REAL(DBL) :: P, RECID

    IF (X < 0.D0) STOP 'X LT 0   DIVIDE'
    IF (X > 0.5D0) STOP 'X GT 0.5 DIVIDE'

    RECID = 0.D0
    P = 1.D0

    DO i=1, N
      RECID = RECID + P
      P = -P*X
    END DO

    DIVIDE = RECID
  END FUNCTION DIVIDE

  !new function for 1/Z, where 0<Z<1
  !1/(1-x)=(1+x+x**2+....x**N)
  REAL(DBL) FUNCTION LT1DIVIDE(X) RESULT(RECIP)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: X
    REAL(DBL) :: P
    INTEGER :: i, N

    IF (X < 0.D0) STOP 'X LT 0 LT1DIVIDE'
    IF (X > 1.D0) STOP 'X GT 1 LT1DIVIDE'

    N = 1000
    RECIP = 0.D0
    P = 1.D0

    DO i=1, N
      RECIP = RECIP + P
      P=P*X
    END DO
  END FUNCTION LT1DIVIDE

  !1/Z = 1/(0.5**N) + (Z-0.5**N)
  !    = (1./0.5**N) 1/(1+ (Z-0.5**N)/(0.5**N))
  !    = 2**N [1/(1 + (2^N(Z-0.5^N)))]
  REAL(DBL) FUNCTION TNYDIVIDE(Z) RESULT(RECIP)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: Z
    REAL(DBL) :: PP, PR, X
    INTEGER :: i

    RECIP = 0.D0
    PP = 1.D0
    PR = 1.D0
    DO i=0, 64
      IF(Z > PR .and. Z < PR*2.D0)THEN
        X = (Z-PR)*PP
        RECIP = PP*DIVIDER(X)
      END IF
      PR = PR*0.5D0
      PP = PP*2.D0
    END DO
  END FUNCTION TNYDIVIDE

  REAL(DBL) FUNCTION GT1DIVIDE(X) RESULT(RECIP)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: X
    REAL(DBL) :: ONED2N(0:64), ARE(0:64), TWON(0:64)
    REAL(DBL) :: Y
    INTEGER :: i

    ONED2N(0) = 1.D0
    DO i=1, 64
      ONED2N(i) = ONED2N(i-1)*0.5D0
    END DO

    ARE(0) = 0.5D0
    DO i=1, 64
      ARE(i) = ONED2N(i)*DIVIDE(ONED2N(i))
    END DO

    TWON(0) = 1.D0
    DO i=1, 64
      TWON(i) = TWON(i-1)*2.D0
    END DO

    RECIP = 0.D0
    IF (X .LT. 1.D0) THEN
      RECIP = DIVIDER(X)
    ELSE
      DO i=0, 63
        IF (X .GE. TWON(i) .AND. X .LT. TWON(i+1)) THEN
          Y = (X - TWON(i))*ARE(i)
          RECIP = DIVIDER(Y)
          RECIP = RECIP*ARE(i)
        END IF
      END DO
    END IF
  END FUNCTION GT1DIVIDE

  REAL(DBL) FUNCTION SINB(num)  !Bhaskara approximation
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num
    REAL(DBL) :: x
    REAL(DBL) :: temp

    x = FMOD(ABSO(num),PI)

    temp = (PI-x)*x
    SINB = (16.D0*temp) / (5.D0*PI*PI - 4.D0*temp)

    !Adjust for negative angles
    IF (num > 0.D0) THEN
      IF (FMOD(num,2.D0*PI) > PI) SINB = -SINB
    ELSE
      IF (FMOD(num,2.D0*PI) > -PI) SINB = -SINB
    END IF
  END FUNCTION SINB

!
!Only difference is when sorting idxs
!
  SUBROUTINE Class_DIAGNxN(NDH, NBS, HAM, UMT, PRD, SPC)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDH, NBS
    REAL(DBL), INTENT(INOUT) :: HAM(:,:)
    REAL(DBL), INTENT(OUT) :: UMT(:,:), PRD(:,:), SPC(:,:)

    REAL(DBL) :: H(2,2), O(2,2), E(2)!, V(2,2), T(2,2), D(2,2)  !J2X2
    REAL(DBL) :: OFF_MAX, ERRPREV, ERRNW, FOM(1000)    !fom = Figure of marite

    INTEGER :: i, j, k, l, m, n, iTry
    INTEGER :: MXIT, nPair, idx(2,1000)
    LOGICAL :: update(1000), new
    INTEGER :: mPair, no

    IF (NDH > 1000) STOP 'NDH must be smaller then 1000'

    IF (NBS > NDH) THEN
      WRITE(*,*) 'NDH must be larger than ', NBS
      STOP
    END IF

    UMT = 0.D0
    DO i=1, NBS
      UMT(i,i) = 1.D0
    END DO

    PRD = HAM
    MXIT = NBS*NBS*2

    DO iTry = 1, MXIT
      update = .TRUE.
      OFF_MAX = 0.D0
      ERRPREV = 0.D0
      nPair = 0

      DO i=1, NBS
        IF (update(i)) THEN
          new = .FALSE.
          DO j=1, NBS
!            IF (update(j)) THEN
            IF (j /= i) THEN
              ERRPREV = ERRPREV + PRD(i,j)*PRD(j,i)
              IF (ABSO(PRD(i,j)) >= OFF_MAX) THEN
                k = i
                l = j
                OFF_MAX = ABSO(PRD(i,j))
                new = .TRUE.
              END IF
            END IF
          END DO

          IF (new) THEN
            mPair = 0
            DO j=1, nPair
              no = 0
              IF (idx(1,j) == k) no=1
              IF (idx(2,j) == k) no=1
              IF (idx(1,j) == l) no=1
              IF (idx(2,j) == l) no=1
              IF (no == 0) THEN
                mPair = mPair + 1
                idx(1,mPair) = idx(1,j)
                idx(2,mPair) = idx(2,j)
                FOM(mPair) = FOM(j)
              ELSE
                update(idx(1,j)) = .TRUE.
                update(idx(2,j)) = .TRUE.
              END IF
            END DO

            nPair = mPair + 1
            DO j=nPair, 2, -1
              idx(1,j) = idx(1,j-1)
              idx(2,j) = idx(2,j-1)
              FOM(j) = FOM(j-1)
            END DO

            update(k) = .FALSE.
            update(l) = .FALSE.

            idx(1,1) = MIN(k,l)
            idx(2,1) = MAX(k,l)

            FOM(1) = OFF_MAX
          END IF
        END IF
      END DO

      DO i=1, nPair
        WRITE(*,*) idx(1,i), idx(2,i), PRD(idx(1,i),idx(2,i))
      END DO

      k = idx(1,1)
      l = idx(2,1)

      H(1,1) = PRD(k,k)
      H(1,2) = PRD(k,l)
      H(2,1) = PRD(l,k)
      H(2,2) = PRD(l,l)
      CALL J2X2(H, E, O)

!      WRITE(*,*) 'E and O values:'
!      DO i=1, 2
!        WRITE(*,*) E(i), (O(i,j),j=1,2)
!      END DO

      SPC = 0.D0
      DO i=1, NBS
        SPC(i,i) = 1.D0
      END DO

      SPC(k,k) = O(1,1)
      SPC(l,k) = O(2,1)
      SPC(l,l) = O(2,2)
      SPC(k,l) = O(1,2)

      !Get new unitary matrix
      PRD = 0.D0
      DO m=1, NBS
        IF (m /= k .AND. m /= l) THEN
          DO n=1, NBS
            PRD(n,m) = UMT(n,m)
          END DO
        END IF
      END DO

      DO n=1, NBS
        PRD(n,k) = PRD(n,k) + UMT(n,k)*O(1,1)
        PRD(n,k) = PRD(n,k) + UMT(n,l)*O(2,1)
        PRD(n,l) = PRD(n,l) + UMT(n,k)*O(1,2)
        PRD(n,l) = PRD(n,l) + UMT(n,l)*O(2,2)
      END DO

      UMT = PRD
      DO i=1, NBS
        DO k=1, NBS
          SPC(k,i) = 0.D0
          DO l=1, NBS
            SPC(k,i) = SPC(k,i) + UMT(l,i)*HAM(l,k) !l,k better than k,l base on memory allocation
          END DO
        END DO
      END DO

      PRD = 0.D0
      !Make new HAM
      DO i=1, NBS
        DO j=1, NBS
          DO k=1, NBS
            PRD(j,i) = PRD(j,i) + UMT(k,j)*SPC(k,i)
          END DO
        END DO
      END DO

      ERRNW = 0.D0
      DO i=1,NBS
        DO j=i+1, NBS
          ERRNW = ERRNW + PRD(i,j)*PRD(j,i)
        END DO
      END DO

      WRITE(*,30) iTry, ERRNW, ERRPREV, ERRPREV - OFF_MAX*OFF_MAX
      30 FORMAT(I3, 3G15.6)

      IF (ERRNW < 1.D-12) EXIT
    END DO

    IF (iTry == MXIT) THEN
      WRITE(*,*) 'Warning: No Convergence'
    END IF

    WRITE(*,*) iTry, NBS, float(iTry)/float(NBS*NBS), 'Diag Eff'

    HAM = PRD
  END SUBROUTINE Class_DIAGNxN
END MODULE retireFunctions

