! @file TestingFunctions.f08
! @brief A tester driver
! 
! This contains the basics to test every function/subroutine used by the Quantum Now software
!
! @author Robert Alvarez
! @bug No known bugs
!
PROGRAM testing_func
  USE FortranFunctions, ONLY: PI, DBL, FMOD, DIV, Class_DIV, SINE, COSINE
  USE ArrayFunctions
  USE Applications
  IMPLICIT NONE

  INTEGER :: input, NBS

DO input=-5,5 
NBS = NBS + 1
END DO
WRITE(*,*) NBS

  DO
    WRITE(*,*)
    WRITE(*,*) 'You can test for:'
    WRITE(*,*) 'FortranFunctions: 1. modulus, 2. division, 3. sine'
    WRITE(*,*) 'ArrayFUnctions: 4. Inverse '
    WRITE(*,*) 'Anything else to exit'
    READ(*,*) input
    IF (input == 1) THEN
      CALL modulus_test()
    ELSE IF (input == 2) THEN
      CALL div_comp()
    ELSE IF (input == 3) THEN
      CALL trigTest()
    ELSE IF (input == 4) THEN
      CALL inverseTest()
    ELSE
      WRITE(*,*) 'Have a nice day:)'
      EXIT
    END IF
  END DO

  CONTAINS

  ! Read arbitrary size matrix values from keyboard and send matrix to get inverted
  SUBROUTINE inverseTest()
    IMPLICIT NONE
    
    REAL(DBL), ALLOCATABLE :: A(:,:), B(:,:), nums(:)
    INTEGER :: i, n

    WRITE(*,*) 'What is the size of the matrix?'
    READ(*,*) n

    ALLOCATE(A(n,n), B(n,n), nums(n))

    DO i=1,n
      10 FORMAT (A,I2,A,I2)
      WRITE(*,10) 'Enter ',n,' numbers for row ',i
      READ(*,*) nums(:)
      A(i,:) = nums
    END DO

    B = INVERSE(A)
  END SUBROUTINE

  !Generate 10000 different numerator values from -720 to 720 and compare
  !modulus operation with a denominator of PI and test if value is correct
  SUBROUTINE modulus_test
    IMPLICIT NONE

    INTEGER :: i
    REAL(DBL) :: N, D, ans, ans1

    WRITE(*,*) 'modulus_test'

    D = PI
    DO i = 1,10000
      N = -10.0D0
      CALL newAngle(N)

      ans = MOD(N,D)
      ans1 = FMOD(N,D)
      IF (ABS(ans-ans1) > 1.0D-5) THEN
        WRITE(*,*)  'N = ', N, 'D = ', D
        WRITE(*,*) MOD(N,D)
        WRITE(*,*) FMOD(N,D)
      END IF
    END DO
    WRITE(*,*) 'None of the 10000 cases have an error greater than 1E-5'
  END SUBROUTINE modulus_test

  !Print error generated and time take between Class_DIV and DIV functions for the
  !denominator intervals of 0.0 < D < 0.1, 0.1 < D < 1.0, and 1.0 < D < 10000
  SUBROUTINE div_comp
    IMPLICIT NONE

    PROCEDURE(DIV), POINTER :: p
    REAL(DBL) :: factor(10)

    WRITE(*,*) 'Division functions comparition:'

    !Class_DIV with 0.0 < D < 0.1
    WRITE(*,*)
    WRITE(*,*) 'Function: Class_DIV - range 0.0 < D < 0.1'
    p => Class_DIV
    factor = [1.0D-1,1.0D-2,1.0D-3,1.0D-4,1.0D-5,1.0D-6,1.0D-7,1.0D-8,1.0D-9,1.0D-10]
    CALL auto_div(p, factor, 0.1, 1.0)

    !DIV with 0.0 < D < 0.1
    WRITE(*,*)
    WRITE(*,*) 'Function: DIV - range: 0.0 < D < 0.1'
    p => DIV
    CALL auto_div(p, factor, 0.1, 1.0)

    !Class_DIV with 0.1 < D < 1.0
    WRITE(*,*)
    WRITE(*,*) 'Function: Class_DIV - range 0.1 < D < 1.0'
    p => Class_DIV
    factor = [0.15D0,0.2D0,0.3D0,0.4D0,0.5D0,0.55D0,0.6D0,0.7D0,0.8D0,0.85D0]
    CALL auto_div(p, factor, 0.0, 0.1)

    !DIV with 0.1 < D < 1.0
    WRITE(*,*)
    WRITE(*,*) 'Function: DIV - range 0.1 < D < 1.0'
    p => DIV
    CALL auto_div(p, factor, 0.0, 0.1)

    !Class_DIV with 1.0 < D < 1,000
    WRITE(*,*)
    WRITE(*,*) 'Function: Class_DIV - range 1.0 < D < 1,600'
    p => Class_DIV
    factor = [1.0D1,2.0D1,4.0D1,6.0D1,8.0D1,1.0D2,2.0D2,4.0D2,6.0D2,8.0D2]
    CALL auto_div(p, factor, 0.0, 1.0)

    !DIV with 1.0 < D < 10,000
    WRITE(*,*)
    WRITE(*,*) 'Function: DIV - range 1.0 < D < 1,600'
    p => DIV
    CALL auto_div(p, factor, 0.0, 1.0)
  END SUBROUTINE div_comp

  !Perform a division between a random numerator and denominator number in the intervals
  !between low_b and up+b using the P_DIV function which could be Class_DIV or DIV
  SUBROUTINE auto_div(P_DIV, factor, low_b, up_b)
    IMPLICIT NONE

    PROCEDURE(DIV), POINTER, INTENT(IN) :: P_DIV
    REAL(DBL), INTENT(IN) :: factor(:)
    REAL, INTENT(IN) :: low_b, up_b

    INTEGER :: i, j, nTimes
    REAL(DBL) :: N, D, t, error

    j = 1
    error = 0.0D0
    nTimes = 100000
    t = GetTime()

    DO i=1, nTimes
      CALL RANDOM_NUMBER(N)
      CALL RANDOM_NUMBER(D)

      IF (DBLE(low_b) < D .AND. D < DBLE(up_b)) THEN
        IF (low_b == 0.1 .AND. up_b == 1.0) D = D*factor(j)
        IF (low_b == 0.0 .AND. up_b == 0.1) D = D+factor(j)
        IF (low_b == 0.0 .AND. up_b == 1.0) D = (D+1.0D0)*factor(j)
        j = j + 1
        IF (j >= 11) j=1
      END IF
      error = error + ABS(N/D - P_DIV(N,D))
    END DO

    WRITE(*,*) 'Time:', GetTime() - t
    WRITE(*,*) 'error = ', error
  END SUBROUTINE auto_div

  !Print angle, Fortran sine function for angle, our own sine function for angle, and the error.
  !Angle have equally step sizes of PI/36 from -PI to PI
  SUBROUTINE trigTest
    IMPLICIT NONE

    REAL(DBL) :: angle, add
    INTEGER :: i

    WRITE(*,*) 'trigTest:'

    add = PI/36.0D0
    angle = -PI

    WRITE(*,'(A,4(A9,A10))') ' ', 'angle', '', 'SINE', '', 'SIN', '', 'SIN error'
    DO i=0, 72
        WRITE(*,'(F15.10,4F20.16)') angle, SINE(angle), SIN(angle), SIN(angle)-SINE(angle)
        angle = angle + add
    END DO
  END SUBROUTINE trigTest

  !Generate and return a new number from -720 and 720
  SUBROUTINE newAngle(getNew)
    IMPLICIT NONE

    REAL(DBL), INTENT(INOUT) :: getNew

    IF (getNew >= 0.0D0) THEN
      CALL RANDOM_NUMBER(getNew)
      getNew = -getNew*720D0
    ELSE
      CALL RANDOM_NUMBER(getNew)
      getNew = getNew*720D0
    END IF
  END SUBROUTINE newAngle

  !Get current time
  REAL(DBL) FUNCTION GetTime()
    IMPLICIT NONE

    INTEGER i, timer_count_rate, timer_count_max
    CALL SYSTEM_CLOCK(i, timer_count_rate, timer_count_max)
    GetTime = DBLE(i) / DBLE(timer_count_rate)
  END FUNCTION GetTime
END PROGRAM testing_func

