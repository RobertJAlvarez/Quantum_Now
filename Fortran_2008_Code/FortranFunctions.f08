MODULE FortranFunctions
  PRIVATE
  PUBLIC :: DBL, PI, ABSO, COSINE, SINE, DIV, SQR, FMOD
  PUBLIC :: Class_DIV, SQR2 !Only use at TestingFunctions.f08

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
    REAL*8, INTENT(IN) :: AA, BB

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

  REAL(DBL) FUNCTION RECIPROCAL(ZZ)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: ZZ
    REAL(DBL) :: RECID, Z

    Z = ZZ

    IF (Z > 1.0D0)THEN                        !IF Z > 1
      CALL GT1DIVIDE(Z - 1.0D0,RECID)
    ELSE IF(Z > 0.1D0 .AND. Z < 1.0D0) THEN   !IF 0.1 < Z < 1
      CALL LT1DIVIDE(1.0D0 - Z,RECID)
    ELSE IF (Z > 0.0D0 .AND. Z <= 0.1D0) THEN !IF 0 < Z <= 0.1
      CALL TNYDIVIDE(Z,RECID)
    ELSE                                      !IF Z = 1
      RECID = 1.0D0
    END IF

    RECIPROCAL = RECID
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

  !new subroutine for 1/Z, where 0<Z<1
  !1/(1-x)=(1+x+x**2+....x**N)
  SUBROUTINE LT1DIVIDE(X, RECID)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: X
    REAL(DBL), INTENT(OUT) :: RECID
    REAL(DBL) :: P
    INTEGER :: i, N

    IF (X < 0.0D0) STOP 'X LT 0 LT1DIVIDE'
    IF (X > 1.0D0) STOP 'X GT 1 LT1DIVIDE'

    N = 1000
    RECID = 0.0D0
    P = 1.0D0

    DO i=1, N
      RECID = RECID + P
      P=P*X
    END DO
  END SUBROUTINE LT1DIVIDE

  !1/Z = 1/(0.5**N) + (Z-0.5**N)
  !    = (1./0.5**N) 1/(1+ (Z-0.5**N)/(0.5**N))
  !    = 2**N [1/(1 + (2^N(Z-0.5^N)))]
  SUBROUTINE TNYDIVIDE(Z,RECIDUAL)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: Z
    REAL(DBL), INTENT(INOUT) :: RECIDUAL
    REAL(DBL) :: PP, PR, X
    INTEGER :: i

    PP = 1.0D0
    PR = 1.0D0
    DO i=0, 64
      IF(Z > PR .and. Z < PR*2.0D0)THEN
        X = (Z-PR)*PP
        RECIDUAL = PP*DIVIDER(X)
      END IF
      PR = PR*0.5D0
      PP = PP*2.0D0
    END DO
  END SUBROUTINE TNYDIVIDE

  SUBROUTINE GT1DIVIDE(X,RECID)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: X
    REAL(DBL), INTENT(OUT) :: RECID
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

    IF (X .LT. 1.0D0) THEN
      RECID = DIVIDER(X)
    ELSE
      DO i=0, 63
        IF (X .GE. TWON(i) .AND. X .LT. TWON(i+1)) THEN
          Y = (X - TWON(i))*ARE(i)
          RECID = DIVIDER(Y)
          RECID = ARE(i)*RECID
        END IF
      END DO
    END IF
  END SUBROUTINE GT1DIVIDE

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

  REAL(DBL) FUNCTION SQR2(pass_num)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: pass_num
    REAL(DBL) :: real_part, quot  ! Decimal part of pass_num -- quotient
    INTEGER(INT_QD) :: divi, rem  ! Divisor -- remainder
    INTEGER :: natural_part       ! Whole numbers in pass_num
    INTEGER :: n_dig_pairs        ! Save number of pair for integers
    INTEGER :: new_quot_dig       ! Save new quotient digit
    INTEGER :: i, two_dig         ! Dummy variable -- next two digits to add to rem, 
    INTEGER :: nat_r              ! Save natural_part in reverse pair, Ex: 1234 -> 3412

    natural_part = int(pass_num)
    real_part = pass_num - natural_part

    divi = 0
    quot = 0.0D0
    rem = 0
    n_dig_pairs = 0
    ! If pass_num is greater than 1
    IF (natural_part /= 0) THEN
      ! Make pair
      nat_r = 0
      DO WHILE (natural_part > 0)
        nat_r = nat_r*100 + MOD(natural_part,100)
        natural_part = natural_part/100
        n_dig_pairs = n_dig_pairs + 1
      END DO
      ! Calculate sqrt root
      DO i=1, n_dig_pairs
        two_dig = MOD(nat_r,100)                        ! Get first two digits of natural_part
        nat_r = nat_r/100                               ! Delete those first two digits
        rem = rem*100 + two_dig                         ! Add next 2 digits in the dividend
        new_quot_dig = new_divi(divi,rem)               ! Find new digit for quotient
        quot = quot*10.0D0 + new_quot_dig               ! Update quotient
        rem = rem - ((divi+new_quot_dig)*new_quot_dig)  ! Perform division
        divi = (divi + 2*new_quot_dig)*10               ! Update dividend
      END DO
    END IF
    ! Real part is always calculated for maximum precision
    DO i=n_dig_pairs+1, 14
      real_part = real_part*100
      two_dig = int(real_part)
      real_part = real_part - two_dig
      rem = rem*100 + two_dig                         ! Add next 2 digits in the dividend
      new_quot_dig = new_divi(divi,rem)               ! Find new digit for quotient
      quot = quot*10.0D0 + new_quot_dig               ! Update quotient
      rem = rem - ((divi+new_quot_dig)*new_quot_dig)  ! Perform division
      divi = (divi + 2*new_quot_dig)*10               ! Update dividend
    END DO

    DO i=n_dig_pairs+1, 14
      quot = quot * 0.1D0
    END DO

    SQR2 = quot
  END FUNCTION SQR2

  INTEGER FUNCTION new_divi(divisor,dividend)
    IMPLICIT NONE
    INTEGER(INT_QD), INTENT(IN) :: divisor, dividend
    INTEGER :: low, high, mid, best_dig
    INTEGER(INT_QD) :: temp

    low = 0
    high = 9
    best_dig = 0
    DO WHILE (low <= high)
      mid = (high + low)/2
      temp = (divisor + mid)*mid
      IF (temp == dividend) THEN
        best_dig = mid
        EXIT
      ELSE IF (temp < dividend) THEN
        best_dig = mid
        low = mid + 1
      ELSE
        high = mid - 1
      END IF
    END DO
    new_divi = best_dig
  END FUNCTION new_divi
END MODULE FortranFunctions
