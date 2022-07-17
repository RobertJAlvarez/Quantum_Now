! @file TestingFunctions.f08
! @brief A tester driver
! 
! This contains the basics to test every function/subroutine used by the Quantum Now software
!
! @author Robert Alvarez
! @bug No known bugs
!
PROGRAM testing_func
  USE FortranFunctions, ONLY: PI, DBL, FMOD, DIV, SINE, COSINE
  USE retireFunctions, ONLY: Class_DIV
  USE ArrayFunctions
  USE Applications
  IMPLICIT NONE

  INTEGER :: input

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
      N = -10.D0
      CALL newAngle(N)

      ans = MOD(N,D)
      ans1 = FMOD(N,D)
      IF (ABS(ans-ans1) > 1.D-5) THEN
        WRITE(*,*)  'N = ', N, 'D = ', D
        WRITE(*,*) MOD(N,D)
        WRITE(*,*) FMOD(N,D)
      END IF
    END DO
    WRITE(*,*) 'None of the 10000 cases have an error greater than 1E-5'
  END SUBROUTINE modulus_test
  
  !Print error generated and time take between Class_DIV and DIV functions for the
  !denominator with intervals of 0.0 < D < 0.1, 0.1 < D < 1.0, and 1.0 < D < 100000.0
  SUBROUTINE div_comp
    IMPLICIT NONE

    PROCEDURE(DIV), POINTER :: p
    REAL(DBL), DIMENSION(10) :: factor1, factor2, factor3
    INTEGER :: i, j

    factor1 = [1.D-1,1.D-2,1.D-3,1.D-4,1.D-5,1.D-6,1.D-7,1.D-8,1.D-9,1.D-10]
    factor2 = [0.15D0,0.2D0,0.3D0,0.4D0,0.5D0,0.55D0,0.6D0,0.7D0,0.8D0,0.85D0]
    factor3 = [1.D1,5.D1,1.D2,5.D2,1.D3,5.D3,1.D4,5.D4,1.D5,1.D5]

    DO i=1,3
      DO j=1,2
        WRITE(*,*)
        SELECT CASE (j)
        CASE (1)
          p => Class_DIV
          SELECT CASE (i)
          CASE (1)
            WRITE(*,*) 'Function: Class_DIV - range 0.0 < D < 0.1'
            CALL auto_div(p, factor1, 0.1, 1.0)
          CASE (2)
            WRITE(*,*) 'Function: Class_DIV - range 0.1 < D < 1.0'
            CALL auto_div(p, factor2, 0.0, 0.1)
          CASE (3)
            WRITE(*,*) 'Function: Class_DIV - range 1.0 < D < 100000.0'
            CALL auto_div(p, factor3, 0.0, 1.0)
          END SELECT
        CASE (2)
          p => DIV
          SELECT CASE (i)
          CASE (1)
            WRITE(*,*) 'Function: DIV - range 0.0 < D < 0.1'
            CALL auto_div(p, factor1, 0.1, 1.0)
          CASE (2)
            WRITE(*,*) 'Function: DIV - range 0.1 < D < 1.0'
            CALL auto_div(p, factor2, 0.0, 0.1)
          CASE (3)
            WRITE(*,*) 'Function: DIV - range 1.0 < D < 100000.0'
            CALL auto_div(p, factor3, 0.0, 1.0)
          END SELECT
        END SELECT
      END DO
    END DO
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
    error = 0.D0
    nTimes = 100000
    t = GetTime()

    DO i=1, nTimes
      CALL RANDOM_NUMBER(N)
      CALL RANDOM_NUMBER(D)

      IF (DBLE(low_b) < D .AND. D < DBLE(up_b)) THEN
        IF (low_b == 0.1 .AND. up_b == 1.0) D = D*factor(j)
        IF (low_b == 0.0 .AND. up_b == 0.1) D = D+factor(j)
        IF (low_b == 0.0 .AND. up_b == 1.0) D = (D+1.D0)*factor(j)
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

    add = PI/36.D0
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

    IF (getNew >= 0.D0) THEN
      CALL RANDOM_NUMBER(getNew)
      getNew = -getNew*7.D0
    ELSE
      CALL RANDOM_NUMBER(getNew)
      getNew = getNew*7.D0
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

