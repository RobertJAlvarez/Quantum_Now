! File:   FortranFunctions.f08
! Brief:  Contain elementary functions
! 
! See line 11 for all public functions
!
! Author: Robert Alvarez
! Bugs:   No known bugs
!
MODULE FortranFunctions
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: DBL, PI, ABSO, COSINE, SINE, DIV, SQR, FMOD

  INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=15)  ! Use 64 bits, double real precision
  REAL(DBL), PARAMETER :: PI = 3.14159265358979D0

  CONTAINS

  !Last modification: Jun 7th, 2022
  !  Simplification on copying NN and DD
  REAL(DBL) FUNCTION FMOD(NN,DD)  !Return |N|%|D|
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: NN, DD !NN = numerator, DD = denominator
    REAL(DBL) :: N, D

    N = NN
    D = DD
    FMOD = 1.D0  !Save a factor of 1 or -1

    !Make N and D positive
    IF (N < 0.D0) THEN
      N = -N
      FMOD = -FMOD
    END IF
    IF (D < 0.D0) THEN
      D = -D
      FMOD = -FMOD
    END IF

    DO WHILE (N > D)
      N = N - D
    END DO
    FMOD = FMOD * N
  END FUNCTION FMOD

  !Last modification: July 14th, 2022
  !  Simplification on IF statement
  REAL(DBL) FUNCTION ABSO(num)  !Return absolute value of num
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num

    ABSO = num  !Assume num is positive
    IF (num < 0.D0) ABSO = -ABSO
  END FUNCTION ABSO

  !Last modification: Jun 20th, 2022
  !  Simplification on coping NN and DD
  REAL(DBL) FUNCTION DIV(NN,DD) RESULT(N) !Goldschmidt division, return N/D
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: NN, DD !NN = numerator, DD = denominator
    REAL(DBL) :: F, D

    !Stop if denominator is 0
    IF (ABSO(DD) < 1.D-15) STOP "Can't divide by 0"

    !If DD < 0 multiply D and N by -1 so D > 0
    IF (DD < 0.D0) THEN
      N = -NN
      D = -DD
    ELSE
      N = NN
      D = DD
    END IF

    DO WHILE(D > 1.D0)   !Scale N and D so 0<D<1
      N = N*0.1D0
      D = D*0.1D0
    END DO

    DO WHILE(D + 1.D-15 < 1.D0) !Make D = 1
      F = 2.D0 - D
      N = N*F
      D = D*F
    END DO
  END FUNCTION DIV

  !Last modification: July 15th, 2022
  !  Algorithm simplification
  REAL(DBL) FUNCTION SQR(num) RESULT(xn)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num
    INTEGER :: i

    xn = 1.D0

    DO i=1, 10
      xn = xn + DIV(num - xn*xn, 2.D0*xn)
    END DO
  END FUNCTION SQR

!Bhaskara approx: temp = (PI-x)*x
!                 sin(x) = (16*temp) / (5*PI*PI - 4*temp)
!Second approx:   temp = (x/PI)*(x/PI - 1)
!                 sin(x) = (temp/10)*(36*temp - 31)
!Weight average: Bhaskara -> 0.385  Second -> 0.615
!
!sin(x) approx with weight average: temp = (x/PI)*(x/PI - 1)
! sin(x) = temp(2.21652(temp - 31/36) - 1.5372/(1.25 + temp))

  !Last modification: January 12th, 2022
  !  Add another temp to only call DIV() once
  REAL(DBL) FUNCTION SINE(num)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num
    REAL(DBL) :: x, temp

    x = FMOD(ABSO(num),PI)

    temp = DIV(x,PI)
    temp = temp*(temp - 1.D0)
    SINE = temp*(2.21652D0*(temp - DIV(31.D0,36.D0)) - DIV(1.5372D0,1.25D0 + temp)) + 2.563D-5

    !Adjust for negative angles
    IF (num > 0.D0) THEN
      IF (FMOD(num,2.D0*PI) > PI) SINE = -SINE
    ELSE
      IF (FMOD(num,2.D0*PI) > -PI) SINE = -SINE
    END IF
  END FUNCTION SINE

  !Last modification: January 12th, 2022
  !  Creation date
  REAL(DBL) FUNCTION COSINE(num)  !Return cosine(x) by using sine(x) function
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: num

    COSINE = SINE(0.5D0*PI - num)
  END FUNCTION COSINE
END MODULE FortranFunctions
