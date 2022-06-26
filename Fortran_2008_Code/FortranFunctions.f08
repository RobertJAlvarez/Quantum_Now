MODULE FortranFunctions
  PRIVATE
  PUBLIC :: DBL, PI, ABSO, COSINE, SINE, DIV, SQR, FMOD
  PUBLIC :: Class_DIV !Only use at TestingFunctions.f08

  INTEGER, PARAMETER :: INT_QD = SELECTED_INT_KIND(18)  ! Use 64 bits, quadruple integer precision
  INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)  ! Use 64 bits, double real precision
  REAL(DBL), PARAMETER :: PI = 3.14159265358979D0

  CONTAINS

  !Return |N|%|D|
  REAL(DBL) FUNCTION FMOD(NN,DD)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: NN, DD !NN = numerator, DD = denominator
    REAL(DBL) :: N, D

    !Make N and D positive
    IF (NN < 0.0D0 .OR. DD < 0.0D0) THEN
      IF (NN < 0.0D0 .AND. DD < 0.0D0) THEN !If both are negative
        N = -NN
        D = -DD
        FMOD = 1.0D0
      ELSE IF (NN < 0.0D0) THEN !If N is negative and D positive
        N = -NN
        D = DD
        FMOD = -1.0D0
      ELSE                      !If N is positive and D negative
        N = NN
        D = -DD
        FMOD = -1.0D0
      END IF
    ELSE  !If both are positive
      N = NN
      D = DD
      FMOD = 1.0D0
    END IF

    DO WHILE (N > D)
      N = N - D
    END DO
    FMOD = FMOD * N
  END FUNCTION FMOD

  !Return absolute value of num
  REAL(DBL) FUNCTION ABSO(num)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num

    IF (num < 0.0D0) THEN
      ABSO = -num
    ELSE
      ABSO = num
    END IF
  END FUNCTION ABSO

  !Goldschmidt division, return N/D
  REAL(DBL) FUNCTION DIV(NN,DD)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: NN, DD !NN = numerator, DD = denominator
    REAL(DBL) :: F,N,D

    !Stop if denominator is 0
    IF (ABSO(DD) < 1.0D-15) THEN
      STOP "Can't divide by 0"
    END IF

    IF (DD < 0.0D0) THEN  !If DD < 0 multiply D and N by -1 so D > 0
      N = -NN
      D = -DD
    ELSE
      N = NN
      D = DD
    END IF

    DO WHILE(D > 1.0D0)   !Scale N and D so 0<D<1
      N = N*0.1D0
      D = D*0.1D0
    END DO

    DO WHILE(D + 1.0D-15 <= 1.0D0) !Make D = 1
      IF (D <= 0.1D0) THEN
        F = 10.0D0
      ELSE
        F = 2.0D0 - D
      END IF
      N = N*F
      D = D*F
    END DO
    DIV = N
  END FUNCTION DIV

  REAL(DBL) FUNCTION COSINE(num)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num

    COSINE = SINE(DIV(PI,2.0D0) - num)
  END FUNCTION COSINE

!Bhaskara approx:   temp = (PI-x)*x
!                   sin(x) = (16*temp) / (5*PI*PI - 4*temp)
!Second approx:     temp = (x/PI)*(x/PI - 1)
!                   sin(x) = (temp/10)*(36*temp - 31)
!Weight average: Bhaskara -> 0.385  Second -> 0.615
!
!sin(x) approx with weight average: temp = (x/PI)*(x/PI - 1)
! sin(x) = temp(2.21652(temp - 31/36) - 1.5372/(1.25 + temp))
  REAL(DBL) FUNCTION SINE(num)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num
    REAL(DBL) :: x
    REAL(DBL) :: temp

    x = FMOD(ABSO(num),PI)

    temp = DIV(x,PI)
    temp = temp*(temp - 1.0D0)
    SINE = temp*(2.21652D0*(temp - DIV(31.0D0,36.0D0)) - DIV(1.5372D0,1.25D0 + temp))

    !Adjust for negative angles and shift the graph down by 2.6D-5
    SINE = SINE + 2.6D-5
    IF (num > 0.0D0) THEN
      IF (FMOD(num,2.0D0*PI) > PI) SINE = -SINE
    ELSE
      IF (FMOD(num,2.0D0*PI) > -PI) SINE = -SINE
    END IF
  END FUNCTION SINE

  !Return square root of num
  REAL(DBL) FUNCTION SQR(num) !SQR = xn
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num
    REAL(DBL) :: xn_1, xn_2, EPS  !xn_1 = x_(n-1), xn_2 = x_(n-2)
    INTEGER :: i

    SQR = 1.0D0   !SQR needs a first value to start the iterations
    xn_1 = -0.1D0 !SQA, SQB & SQC are dummy variables to keep track
    xn_2 = -0.2D0 !of how SQR evolves through the iterations.

!      X =~(num)^0.5
!  X + E = (num)^0.5
! (X + E)^2 = num
!  X^2 + 2XE + E^2 = num, But E is very small, so E^2=0:
!        E = DIV(num - X^2,2X), this way:
!  X =~ num^0.5 = X + E -> X =~ X + E, this updates X to a new value.
! So by repeating this multiple times we will get num^0.5

    DO i=1, 100
      EPS = DIV(num - SQR*SQR, 2.0D0*SQR) !Get the value of E according to the current S
      SQR = SQR + EPS         !Update x_n to x_(n+1)
      IF(SQR == xn_1) THEN        !If the last two values are the same then the code stops.
        EXIT
      ELSE IF (SQR == xn_2) THEN  !If the ith and ith-2 value are the same, then the
        SQR = DIV(xn_2 + xn_1, 2.0D0) !code is cycling, so we stop it and average between this values
        EXIT
      END IF
      !Update dummy variables so they reflect the two previous values of SQR
      xn_2 = xn_1 !x_(n-2) = x_(n-1)
      xn_1 = SQR  !x_(n-1) = x_n
    END DO
    write(*,*) i
  END FUNCTION SQR

!!!
!Extra functions for reference but never used
!!!

  !Returns A/B
  REAL(DBL) FUNCTION Class_DIV(AA,BB)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: AA, BB

    !Stop if denominator is 0
    IF (ABSO(BB) < 1.0D-15) THEN
      STOP "Can't divide by 0"
    END IF

    !Calculate division where BB is always positive
    IF (BB < 0.0D0) THEN
      Class_DIV = -AA*RECIPROCAL(-BB)
    ELSE
      Class_DIV = AA*RECIPROCAL(BB)
    END IF
  END FUNCTION Class_DIV

  REAL(DBL) FUNCTION RECIPROCAL(Z)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: Z
    REAL(DBL) :: RECIP

    RECIP = 1.0D0 !Only stays as 1.0 if z = 1

    IF (Z > 1.0D0)THEN                        !IF Z > 1
      RECIP = GT1DIVIDE(Z - 1.0D0)
    ELSE IF(Z > 0.1D0 .AND. Z < 1.0D0) THEN   !IF 0.1 < Z < 1
      RECIP = LT1DIVIDE(1.0D0 - Z)
    ELSE IF (Z > 0.0D0 .AND. Z <= 0.1D0) THEN !IF 0 < Z <= 0.1
      RECIP = TNYDIVIDE(Z)
    END IF

    RECIPROCAL = RECIP
  END FUNCTION RECIPROCAL

  !DIVIDER (1/(1+x))   0<x<1
  REAL(DBL) FUNCTION DIVIDER(X)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: X
    REAL(DBL) :: TWOTHIRDS
    REAL(DBL) :: Y

    IF (X < 0.0D0) STOP 'X LT 0'
    IF (X > 1.0D0) STOP 'X GT 1'

    TWOTHIRDS = 0.666666666666667D0

    IF (X > 0.5D0) THEN
      Y = TWOTHIRDS*(X - 0.5D0)
      DIVIDER = DIVIDE(Y)*TWOTHIRDS
    ELSE
      DIVIDER = DIVIDE(X)
    END IF
  END FUNCTION DIVIDER

  !DIVIDE (1/(1+x))   0<x<0.5
  REAL(DBL) FUNCTION DIVIDE(X)
    IMPLICIT NONE

    REAL(DBL), INTENT(IN) :: X
    INTEGER :: i, N = 50
    REAL(DBL) :: P, RECID

    IF (X < 0.0D0) STOP 'X LT 0   DIVIDE'
    IF (X > 0.5D0) STOP 'X GT 0.5 DIVIDE'

    RECID = 0.0D0
    P = 1.0D0

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

    IF (X < 0.0D0) STOP 'X LT 0 LT1DIVIDE'
    IF (X > 1.0D0) STOP 'X GT 1 LT1DIVIDE'

    N = 1000
    RECIP = 0.0D0
    P = 1.0D0

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

    RECIP = 0.0D0
    PP = 1.0D0
    PR = 1.0D0
    DO i=0, 64
      IF(Z > PR .and. Z < PR*2.0D0)THEN
        X = (Z-PR)*PP
        RECIP = PP*DIVIDER(X)
      END IF
      PR = PR*0.5D0
      PP = PP*2.0D0
    END DO
  END FUNCTION TNYDIVIDE

  REAL(DBL) FUNCTION GT1DIVIDE(X) RESULT(RECIP)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: X
    REAL(DBL) :: ONED2N(0:64), ARE(0:64), TWON(0:64)
    REAL(DBL) :: Y
    INTEGER :: i

    ONED2N(0) = 1.0D0
    DO i=1, 64
      ONED2N(i) = ONED2N(i-1)*0.5D0
    END DO

    ARE(0) = 0.5D0
    DO i=1, 64
      ARE(i) = ONED2N(i)*DIVIDE(ONED2N(i))
    END DO

    TWON(0) = 1.0D0
    DO i=1, 64
      TWON(i) = TWON(i-1)*2.0D0
    END DO

    RECIP = 0.0D0
    IF (X .LT. 1.0D0) THEN
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

!Bhaskara approx: temp = (PI-x)*x
!                 sin(x) = (16*temp) / (5*PI*PI - 4*temp)
  REAL(DBL) FUNCTION SINB(num)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num
    REAL(DBL) :: x
    REAL(DBL) :: temp

    x = FMOD(ABSO(num),PI)

    temp = (PI-x)*x
    SINB = (16.0D0*temp) / (5.0D0*PI*PI - 4.0D0*temp)

    !Adjust for negative angles
    IF (num > 0.0D0) THEN
      IF (FMOD(num,2.0D0*PI) > PI) SINB = -SINB
    ELSE
      IF (FMOD(num,2.0D0*PI) > -PI) SINB = -SINB
    END IF
  END FUNCTION SINB
END MODULE FortranFunctions
