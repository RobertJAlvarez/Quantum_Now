! file:   TestingFunctions.f08
! brief:  A tester driver
! 
! This contains the basics to test every function/subroutine used by the Quantum Now software
!
! author: Robert Alvarez
! bug:    No known bugs
!
PROGRAM testingFunctions
  USE fortranFunctions, ONLY: PI, DBL, FMOD, DIV, SQR, SINE
  USE retireFunctions, ONLY: Class_DIV, Class_DIAGNxN
  USE arrayFunctions, ONLY: print_mtx, INVERSE, J2x2, JAC2BY2GEN, DIAGNxN, LEASTSQUARE
  USE applications, ONLY: STARKDVR, RINGDVR, BOXDVR, HMODVR, write_plot_instructions, open_plot
  IMPLICIT NONE

  INTEGER :: input

  DO
    WRITE(*,*)
    WRITE(*,*) 'You can test for:'
    WRITE(*,*) 'FortranFunctions: 1. Modulus, 2. Division, 3. Square root, 4. Sin'
    WRITE(*,*) 'ArrayFunctions: 5. Inverse, 6. J2x2, 7. JAC2BY2GEN, 8. DIAGNxN, 9. LEASTSQUARE'
    WRITE(*,*) 'Applications: 10. STARKDVR, 11. RINGDVR, 12. BOXDVR, 13. HMODVR'
    WRITE(*,*) 'Anything else to exit'
    READ(*,*) input

    SELECT CASE (input)
    CASE (1)
      CALL modulus_test()
    CASE (2)
      CALL div_comp()
    CASE (3)
      CALL sqr_test()
    CASE (4)
      CALL trig_test()
    CASE (5)
      CALL inverse_test()
    CASE (6)
      CALL solve_2by2()
    CASE (7)
      CALL solve_gen()
    CASE (8)
      CALL DIAGDVR()
    CASE (9)
      CALL LSA_test()
    CASE (10)
      CALL STARKDVR()
    CASE (11)
      CALL RINGDVR()
    CASE (12)
      CALL BOXDVR()
    CASE (13)
      CALL HMODVR()
    CASE DEFAULT
      WRITE(*,*) 'Have a nice day:)'
      EXIT
    END SELECT

    IF (input >= 10 .AND. input <= 13) THEN
      CALL write_plot_instructions()
      CALL open_plot()
    END IF
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
      IF (ABS(ans-ans1) > 1.D-11) THEN
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

  ! Calculate cumulative error between 1,000,000 square numbers from 0 to 1
  SUBROUTINE sqr_test()
    IMPLICIT NONE
    INTEGER :: i, j, nTimes
    REAL(DBL) :: N, t, error

    j = 1
    error = 0.D0
    nTimes = 1000000
    t = GetTime()

    DO i=1, nTimes
      CALL RANDOM_NUMBER(N)
      error = error + ABS(SQRT(N) - SQR(N))
    END DO

    WRITE(*,*) 'Time:', GetTime() - t
    WRITE(*,*) 'error = ', error
  END SUBROUTINE sqr_test

  !Print angle, Fortran sine function for angle, our own sine function for angle, and the error.
  !Angle have equally step sizes of PI/36 from -PI to PI
  SUBROUTINE trig_test
    IMPLICIT NONE

    REAL(DBL) :: angle, add
    INTEGER :: i

    WRITE(*,*) 'trig_test:'

    add = PI/36.D0
    angle = -PI

    WRITE(*,'(A,4(A9,A10))') ' ', 'angle', '', 'SINE', '', 'SIN', '', 'SIN error'
    DO i=0, 72
        WRITE(*,'(F15.10,4F20.16)') angle, SINE(angle), SIN(angle), SIN(angle)-SINE(angle)
        angle = angle + add
    END DO
  END SUBROUTINE trig_test

  SUBROUTINE newAngle(getNew)
    IMPLICIT NONE
    REAL(DBL), INTENT(INOUT) :: getNew

    IF (getNew >= 0.D0) THEN
      CALL RANDOM_NUMBER(GETNEW)
      getNew = -getNew*720.D0
    ELSE
      CALL RANDOM_NUMBER(GETNEW)
      getNew = getNew*720.D0
    END IF
  END SUBROUTINE newAngle

!
! ArrayFunctions.f08
!
  ! Read arbitrary size matrix values from keyboard and send matrix to get inverted
  SUBROUTINE inverse_test()
    IMPLICIT NONE
    
    REAL(DBL), ALLOCATABLE :: A(:,:), B(:,:)
    INTEGER :: i, n

    WRITE(*,*) 'What is the size of the matrix?'
    READ(*,*) n

    ALLOCATE(A(n,n), B(n,n))

    10 FORMAT (A,I2,A,I2)
    DO i=1,n
      WRITE(*,10) 'Enter ',n,' numbers for row ',i
      READ(*,*) A(i,:)
    END DO

    B = INVERSE(A)

    DEALLOCATE(A, B)
  END SUBROUTINE inverse_test

  SUBROUTINE solve_2by2()
    IMPLICIT NONE
    !

    !
  END SUBROUTINE solve_2by2

  SUBROUTINE solve_gen()
    IMPLICIT NONE
    REAL(DBL), DIMENSION(2,2) :: H, O, V
    REAL(DBL) :: E(2)

    INTEGER, PARAMETER :: seed = 71530
    INTEGER :: ND

    CALL SRAND(seed)
    WRITE(*,*) 'Dimension = ?'
    READ(*,*) ND

    BLOCK
      REAL(DBL), DIMENSION(ND,ND) :: HAM, OVR, EVC, TST, RST, UPD
      REAL(DBL) :: ERR
      INTEGER :: i, j, k, l, m, n, ITST, ITRY
 
      DO i=1, ND
        DO j=i, ND
          HAM(i,j) = DBLE(rand()) - 0.5D0
          OVR(i,j) = 0.3D0*(DBLE(rand()) - 0.D0)
          HAM(j,i) = HAM(i,j)
          OVR(j,i) = OVR(i,j)
        END DO
        HAM(i,i) = 5.D0*HAM(i,i)
        OVR(i,i) = 1.D0 + DBLE(rand())
        WRITE(*,'(12F10.3)') (OVR(i,j),j=1,ND)
      END DO

      EVC = 0.D0
      DO i=1, ND
        EVC(i,i) = 1.D0
      END DO

      WRITE(*,*) 'Initial Hamiltonian:'
      CALL print_mtx(HAM)

      DO ITRY=1, 100
        DO i=1, ND
          DO j=i+1, ND
            H = 0.D0
            O = 0.D0
            DO k=1, ND
              DO l=1, ND
                H(1,1) = H(1,1) + EVC(k,i)*EVC(l,i)*HAM(k,l)
                H(1,2) = H(1,2) + EVC(k,i)*EVC(l,j)*HAM(k,l)
                H(2,1) = H(2,1) + EVC(k,j)*EVC(l,i)*HAM(k,l)
                H(2,2) = H(2,2) + EVC(k,j)*EVC(l,j)*HAM(k,l)
                O(1,1) = O(1,1) + EVC(k,i)*EVC(l,i)*OVR(k,l)
                O(1,2) = O(1,2) + EVC(k,i)*EVC(l,j)*OVR(k,l)
                O(2,1) = O(2,1) + EVC(k,j)*EVC(l,i)*OVR(k,l)
                O(2,2) = O(2,2) + EVC(k,j)*EVC(l,j)*OVR(k,l)
              END DO
            END DO

            CALL JAC2BY2GEN(H,O,V,E)  ! H is never use again, O and E only after the last loop

            UPD = 0.D0
            DO k=1, ND
              UPD(k,k) = 1.D0
            END DO

            UPD(i,i) = V(1,1)
            UPD(i,j) = V(1,2)
            UPD(j,i) = V(2,1)
            UPD(j,j) = V(2,2)

            ! RST1(l,k) = V1_ml HAM(m,n) V1_nk
            ! RST2(q,p) = V2_lq RST!(l,k) V2_kp = V1_ml V2_lq [Hmn V1_nk V2_kp]_np
            ! V2_ji V1i
            TST = 0.D0
            DO n=1,ND
              DO l=1,ND
                DO k=1, ND
                  TST(n,l) = TST(n,l) + EVC(n,k)*UPD(k,l)
                END DO
              END DO
            END DO

            EVC = TST
            DO ITST=1,2
              IF (ITST == 1) TST = OVR
              IF (ITST == 2) TST = HAM
              DO k=1, ND
                DO l=1, ND
                  RST(l,k) = 0.D0
                  DO m=1, ND
                    DO n=1, ND
                      RST(l,k) = RST(l,k) + EVC(m,k)*EVC(n,l)*TST(m,n)
                    END DO
                  END DO
                END DO
!                WRITE(*,'(12F10.3)') (RST(l,k),l=1,ND)
              END DO
              ERR = 0.D0
              DO m=1, ND
                DO n=m+1,ND
                  ERR = ERR + RST(m,n)*RST(n,m)
                END DO
              END DO
!              WRITE(*,*) 'ERR:', SQRT(ERR)
            END DO
          END DO
        END DO
      END DO

      WRITE(*,'(/,2F12.6)') E(1), E(2)
      DO i=1, 2
        WRITE(*,*) E(i), (O(j,i), j=1,2)
      END DO
      WRITE(*,*) O(1,1)*O(1,2) + O(2,1)*O(2,2)
      WRITE(*,*) O(1,1)*O(1,1) + O(2,1)*O(2,1)
      WRITE(*,*) O(2,2)*O(2,2) + O(1,2)*O(1,2)

      WRITE(*,*) 'Resulting Hamiltonian:'
      CALL print_mtx(HAM)
    END BLOCK
  END SUBROUTINE solve_gen

  SUBROUTINE DIAGDVR()  !Diag driver
    IMPLICIT NONE
    REAL(DBL) :: HAM(7,7), UMT(7,7)  !Hamiltonian, Unitary
    REAL(DBL) :: G, X, P, TXR
    INTEGER :: i, j, NBS

    NBS = 7

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

    CALL DIAGNxN(HAM, UMT)

    WRITE(*,*) 'Diagonalized Hamiltonian:'
    CALL print_mtx(HAM)
  END SUBROUTINE DIAGDVR

  SUBROUTINE LSA_test()
    IMPLICIT NONE
    !
  END SUBROUTINE
END PROGRAM testingFunctions

