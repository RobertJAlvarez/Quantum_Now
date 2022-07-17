!!!
!Extra functions for reference but never used
!!!
MODULE retireFunctions
  USE FortranFunctions, ONLY: PI, DBL, FMOD, ABSO, DIV

  PRIVATE
  PUBLIC :: Class_DIV !Only TestingFunctions.f08 use it

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
END MODULE retireFunctions

