! file:   TestingFunctions.f08
! brief:  A tester driver
! 
! This contains the basics to test every function/subroutine used by the Quantum Now software
!
! author: Robert Alvarez
! bug:    No known bugs
!
PROGRAM testing_func
  USE FortranFunctions, ONLY: PI, DBL, FMOD, DIV, SINE
  USE retireFunctions, ONLY: Class_DIV, Class_DIAGNxN
  USE ArrayFunctions, ONLY: print_mtx, INVERSE, J2x2, JAC2BY2GEN, DIAGNxN, LEASTSQUARE
  USE Applications, ONLY: STARKDVR, RINGDVR, HMODVR, BOXDVR
  IMPLICIT NONE

  INTEGER :: input

  DO
    WRITE(*,*)
    WRITE(*,*) 'You can test for:'
    WRITE(*,*) 'FortranFunctions: 1. modulus, 2. division, 3. sine'
    WRITE(*,*) 'ArrayFunctions: 4. Inverse, 5. J2x2, 6. JAC2BY2GEN, 7. DIAGNxN, 8. LEASTSQUARE '
    WRITE(*,*) 'Applications: 9. STARKDVR, 10. RINGDVR, 11. BOXDVR, 12. HMODVR'
    WRITE(*,*) 'Anything else to exit'
    READ(*,*) input

    SELECT CASE (input)
    CASE (1)
      CALL modulus_test()
    CASE (2)
      CALL div_comp()
    CASE (3)
      CALL trigTest()
    CASE (4)
      CALL inverseTest()
    CASE (5)
      CALL J2x2_test()
    CASE (6)
      CALL gen_J2x2_test()
    CASE (7)
      CALL DIAGDVR()
    CASE (8)
      CALL LSA_test()
    CASE (9)
      CALL STARKDVR()
    CASE (10)
      CALL RINGDVR()
    CASE (11)
      CALL BOXDVR()
    CASE (12)
      CALL HMODVR()
    CASE DEFAULT
      WRITE(*,*) 'Have a nice day:)'
      EXIT
    END SELECT
  END DO

  CONTAINS
!
! Extra functions
!
  !Get current time
  REAL(DBL) FUNCTION GetTime()
    IMPLICIT NONE

    INTEGER i, timer_count_rate, timer_count_max
    CALL SYSTEM_CLOCK(i, timer_count_rate, timer_count_max)
    GetTime = DBLE(i) / DBLE(timer_count_rate)
  END FUNCTION GetTime

!
! FortranFunctions.f08
!
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
            CALL auto_div(p, factor1, 0.D0, 0.1D0,i)
          CASE (2)
            WRITE(*,*) 'Function: Class_DIV - range 0.1 < D < 1.0'
            CALL auto_div(p, factor2, 0.1D0, 1.D0,i)
          CASE (3)
            WRITE(*,*) 'Function: Class_DIV - range 1.0 < D < 100000.0'
            CALL auto_div(p, factor3, 1.D0, 1.D5,i)
          END SELECT
        CASE (2)
          p => DIV
          SELECT CASE (i)
          CASE (1)
            WRITE(*,*) 'Function: DIV - range 0.0 < D < 0.1'
            CALL auto_div(p, factor1, 0.D0, 0.1D0,i)
          CASE (2)
            WRITE(*,*) 'Function: DIV - range 0.1 < D < 1.0'
            CALL auto_div(p, factor2, 0.1D0, 1.D0,i)
          CASE (3)
            WRITE(*,*) 'Function: DIV - range 1.0 < D < 100000.0'
            CALL auto_div(p, factor3, 1.D0, 1.D5,i)
          END SELECT
        END SELECT
      END DO
    END DO
  END SUBROUTINE div_comp

  !Perform a division between a random numerator and denominator number in the intervals
  !between low_b and up+b using the P_DIV function which could be Class_DIV or DIV
  SUBROUTINE auto_div(P_DIV, factor, low_b, up_b, n_case)
    IMPLICIT NONE

    PROCEDURE(DIV), POINTER, INTENT(IN) :: P_DIV
    REAL(DBL), INTENT(IN) :: factor(:)
    REAL(DBL), INTENT(IN) :: low_b, up_b
    INTEGER, INTENT(IN) :: n_case

    INTEGER :: i, j, nTimes
    REAL(DBL) :: N, D, t, error

    j = 1
    error = 0.D0
    nTimes = 1000000
    t = GetTime()

    DO i=1, nTimes
      CALL RANDOM_NUMBER(N)
      CALL RANDOM_NUMBER(D)

      IF (DBLE(low_b) < D .OR. D > DBLE(up_b)) THEN
        SELECT CASE(n_case)
        CASE (1)
          D = D*factor(j) + 1.D-15
        CASE (2)
          D = D + factor(j)
        CASE (3)
          D = (D+1.D0)*factor(j)
        END SELECT
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

!
! ArrayFunctions.f08
!
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

  SUBROUTINE J2x2_test()
    IMPLICIT NONE
    !
  END SUBROUTINE J2x2_test
 
  SUBROUTINE gen_J2x2_test()
    IMPLICIT NONE
    !
  END SUBROUTINE gen_J2x2_test

  SUBROUTINE DIAGDVR()  !Diag driver
    IMPLICIT NONE
    REAL(DBL), ALLOCATABLE :: HAM(:,:), UMT(:,:), PRD(:,:)  !Hamiltonian, Unitary, Product
    REAL(DBL) :: G, X, P, TXR
    INTEGER :: i, j, NBS

    NBS = 7
    ALLOCATE(HAM(7,7), UMT(7,7), PRD(7,7))

    G = -1.D0   !Ground
    X = -0.5D0  !Exited
    P =  0.1D0  !Perturbation
    TXR = 0.2D0 !Transfer
    HAM = 0.D0

    HAM(1,1) = G
    HAM(2,2) = X
    HAM(3,3) = X
    HAM(4,4) = G
    HAM(5,5) = X
    HAM(6,6) = X
    HAM(7,7) = G

    HAM(1,2) = P
    HAM(2,3) = TXR
    HAM(3,4) = P
    HAM(4,5) = P
    HAM(5,6) = TXR
    HAM(6,7) = P

    DO i=1, NBS
      DO j=i, NBS
        HAM(j,i) = HAM(i,j)
      END DO
    END DO
    WRITE(*,*) 'Original Hamiltonian:'
    CALL print_mtx(HAM)

    CALL DIAGNxN(NBS, HAM, UMT, PRD)

    WRITE(*,*) 'Updated Hamiltonian:'
    CALL print_mtx(HAM)
  END SUBROUTINE DIAGDVR

  SUBROUTINE LSA_test()
    IMPLICIT NONE
    !
  END SUBROUTINE
END PROGRAM testing_func

