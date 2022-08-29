MODULE applications
  USE fortranFunctions, ONLY: DBL, PI, ABSO, COSINE, SINE, DIV, SQR
  USE arrayFunctions, ONLY: print_mtx, print_EV, DIAGNxN

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: STARKDVR, RINGDVR, BOXDVR, HMODVR, write_plot_instructions, open_plot

  CONTAINS

  SUBROUTINE STARKDVR()
    IMPLICIT NONE
    REAL(DBL), ALLOCATABLE, DIMENSION(:,:) :: HAM, UMT, DIP !Dipole
    REAL(DBL) :: EFIELD, t, TAU
    INTEGER :: i, j, NBS

    NBS = 3
    ALLOCATE(HAM(NBS,NBS), UMT(NBS,NBS), DIP(NBS,NBS))

    WRITE(*,*) 'WELCOME TO STARK DRIVER, EFIELD=?'
    READ(*,*) EFIELD

    DIP = 0.D0
    DIP(1,3) = 0.01D0
    DIP(2,3) = 0.03D0
    DO i=1,NBS
      DO j=i+1,NBS
        DIP(j,i) = DIP(i,j)
      END DO
    END DO

    HAM = 0.D0
    HAM(1,1) = -0.5D0
    HAM(2,2) = -0.125D0
    HAM(3,3) = -0.125D0

    HAM(1,3) = DIP(1,3)*EFIELD
    HAM(2,3) = DIP(2,3)*EFIELD
    DO i=1,NBS
      DO j=i+1,NBS
        HAM(j,i) = HAM(i,j)
      END DO
    END DO

    CALL DIAGNxN(HAM, UMT)
    CALL print_diag_mtx_info(HAM, UMT)

! AT TIME=0, OCCUPY THE 2S FUNCTION:
! |PHI_K (t)> = SUM_J exp(ie(J) t)*OVR(k,j)|PSI_J> = SUM_JL exp(iEjt)OVR(j,k)*UMT(k,l)|PHI_l> !E_j=LAMDA_J
    TAU = DIV(8.D0*PI,ABSO(HAM(1,1)))
    OPEN(12,FILE='PLOT')
    DO j=0,1000
      t = DBLE(j)*DIV(TAU,1.D3)
      CALL calc_something(HAM, UMT, DIP, t, 12)
    END DO
    CLOSE(12)
  
    DEALLOCATE(HAM, UMT, DIP)
  END SUBROUTINE STARKDVR

  SUBROUTINE RINGDVR()
    IMPLICIT NONE
    REAL(DBL), ALLOCATABLE, DIMENSION(:,:) :: HAM, UMT, DIP
    REAL(DBL) :: RINGSZ, E, aM, t, TAU
    INTEGER :: i, j, NBS, iAP

    NBS = 11
    ALLOCATE(HAM(NBS,NBS), UMT(NBS,NBS), DIP(NBS,NBS))

    WRITE(*,*) 'HOW LARGE IS YOUR RING?'
    READ(*,*)  RINGSZ
    WRITE(*,*) 'ELECTRIC(1) OR MAGNETIC(2)?'
    READ(*,*) iAP
    WRITE(*,*) 'HOW BIG THE FIELD?'
    READ(*,*) E

    aM = -5.D0
    HAM = 0.D0
    DO i=1,NBS 
      HAM(i,i) = DIV(aM*aM,RINGSZ*RINGSZ)
      aM = aM + 1.D0
    END DO

    DIP = 0.D0
    IF (iAP /= 2) THEN
      DO i=2,NBS
        DIP(i-1,i) = DIV(E*RINGSZ,2.D0)
        DIP(i,i-1) = DIP(i-1,i)
      END DO
    ELSE
      i = 0
!     E=-1*DIV(2.D0,RINGSZ*RINGSZ)
      DO j=-5,5
        i = i + 1
        DIP(i,i) = DBLE(j)*E 
      END DO
    END IF

    WRITE(*,*) "INITIAL HAM:"
    CALL print_mtx(HAM)

    DIP = DIP + HAM

    CALL DIAGNxN(HAM,UMT)
    CALL print_diag_mtx_info(HAM, UMT)

! AT TIME=0, OCCUPY THE 2S FUNCTION:
! |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(k,j)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
    TAU = DIV(8.D0*PI,ABSO(HAM(1,1)))
    OPEN(12,FILE='PLOT')
    DO j=-100,100000
      IF (j == 0)THEN
        HAM = DIP
        CALL DIAGNxN(HAM,UMT)
        WRITE(*,'(10F12.8)') (HAM(i,i),i=1,NBS)
      END IF
      t = DBLE(j)*DIV(TAU,5.D1)
      CALL calc_something(HAM, UMT, DIP, t, 12)
    END DO
    CLOSE(12)

    DEALLOCATE(HAM, UMT, DIP)
  END SUBROUTINE RINGDVR

  SUBROUTINE BOXDVR()
    IMPLICIT NONE
    REAL(DBL), ALLOCATABLE, DIMENSION(:,:) :: HAM, UMT, DIP
    REAL(DBL) :: twom, TAU, boxsz, t
    INTEGER :: i, j, NBS

    NBS = 9
    ALLOCATE(HAM(NBS,NBS), UMT(NBS,NBS), DIP(NBS,NBS))

    WRITE(*,*) 'WELCOME TO BOX DRIVER, HOW LARGE IS YOUR BOX?'
    READ(*,*) BOXSZ

    i = 0
    HAM = 0.D0
    DO j=-4,4
      i = i + 1
      TWOM = 2.D0*DBLE(j)*PI
      HAM(i,i) = DIV(TWOM,BOXSZ)*DIV(TWOM,BOXSZ)
    END DO

    WRITE(*,*) "INITIAL HAM:"
    CALL print_mtx(HAM)

    CALL DIAGNxN(HAM,UMT)
    CALL print_diag_mtx_info(HAM, UMT)

! AT TIME=0, OCCUPY THE 2S FUNCTION:
! |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(k,j)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
    !TAU = DIV(8.D0*PI,ABSO(HAM(1,1)))
    TAU = DIV(1.D0,ABSO(DIV(HAM(1,1),4.D0)))
    OPEN(12,FILE='PLOT')
    DO j=0,1000
      t = DBLE(j)*DIV(TAU,1.D3)
      CALL calc_something(HAM, UMT, DIP, t, 12)
    END DO
    CLOSE(12)

    DEALLOCATE(HAM, UMT, DIP)
  END SUBROUTINE BOXDVR

  SUBROUTINE HMODVR()
    IMPLICIT NONE
    REAL(DBL), ALLOCATABLE, DIMENSION(:,:) :: HAM, UMT, DIP
    REAL(DBL) :: rkmx, rkmn, emax, emin, energy, hmass
    REAL(DBL) :: efield, guess, scale, t, TAU
    INTEGER :: i, j, k, NBS

    NBS = 2
    ALLOCATE(HAM(NBS,NBS), UMT(NBS,NBS), DIP(NBS,NBS))

    WRITE(*,*) 'NRLMOL:',DIV(0.36005D0*1.6D-19,6.626D-34)

    WRITE(*,*) 'WELCOME TO HARMONIC OSCILLATOR DRIVER'
    WRITE(*,*) 'HCl Diatomic in and Electric Field' 
    WRITE(*,*) 'Clorine atom has infinite mass'

    energy = 8.88D13
    WRITE(*,'(A,ES9.2E2,A)') 'IR frequency of HCl molecule is', energy, ' Hz.'

    hmass = 1.6735575D-27
    WRITE(*,'(A,ES14.7E2,A)') 'The mass of the Hydrogen is', hmass, ' kg.'

    WRITE(*,*) 'Strength of Electric field?'
    READ(*,*) EFIELD  !!! EFIELD never used !!!

    DIP(1,2) = 1.D-2
    DO i=1,NBS
      DO j=i+1,NBS
        DIP(j,i) = DIP(i,j)
      END DO
    END DO

    rkmx = 1.D30
    rkmn = 0.D0
    emax = 1.D30
    emin = 0.D0

    DO k=1,30
      WRITE(*,*) 'Guess the spring constant in (Newtons/Meter)'
      WRITE(*,*) 'Guess should be between:',rkmn,' and ',rkmx
      READ(*,*) guess
      scale = SQR(DIV(guess,hmass))
      WRITE(*,'(4G15.6)') guess, hmass, DIV(guess,hmass), scale

      HAM = 0.D0
      HAM(1,1) = DIV(0.5D0*scale,2.D0*PI)
      HAM(2,2) = DIV(1.5D0*scale,2.D0*PI)
      HAM(1,2) = 0.D0 !DIP(1,2)
      HAM(2,1) = 0.D0 !DIP(2,1)
      DO i=1,NBS
        DO j=1,NBS
          HAM(i,j) = DIV(HAM(i,j),scale)
        END DO
      END DO

      WRITE(*,*) "INITIAL HAM:"
      CALL print_mtx(HAM)

      CALL DIAGNxN(HAM,UMT)

      DO i=1,NBS
        HAM(i,i) = HAM(i,i)*scale
      END DO

      HAM(3,3) = MAX(HAM(2,2),HAM(1,1)) - MIN(HAM(2,2),HAM(1,1))
      WRITE(*,"('frequency:',2G12.4)") HAM(3,3),energy
      IF (HAM(3,3) > energy) THEN
        WRITE(*,*) 'Too Big!'
        IF (HAM(3,3) < EMAX) THEN
          EMAX = HAM(3,3)
          RKMX = guess
        END IF
      END IF
      IF (HAM(3,3) < energy) THEN
        WRITE(*,*) 'Too SMALL!'
        IF (HAM(3,3) > EMIN) THEN
          EMIN = HAM(3,3)
          RKMN = guess
        END IF
      END IF
    END DO 

    CALL print_diag_mtx_info(HAM, UMT)

! AT TIME=0, OCCUPY THE 2S FUNCTION:
! |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(k,j)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
    TAU = DIV(8.D0*PI,ABSO(HAM(1,1)))
    OPEN(12,FILE='PLOT')
    DO j=0,1000
      t = DBLE(j)*DIV(TAU,1.D3)
      CALL calc_something(HAM, UMT, DIP, t, 12)
    END DO
    CLOSE(12)

    DEALLOCATE(HAM, UMT, DIP)
  END SUBROUTINE HMODVR

  SUBROUTINE print_diag_mtx_info(HAM, UMT)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: HAM(:,:), UMT(:,:)
    REAL(DBL) :: OVR(SIZE(HAM,1),SIZE(HAM,1))
    INTEGER :: NBS, i, j

    NBS = SIZE(HAM,1)

    WRITE(*,*) "UPDATED HAM"
    CALL print_mtx(HAM)

    !PROVE THAT THE INVERSE OF UMT IS THE TRANSPOSE OF UMT
    CALL print_EV([(HAM(i,i),i=1,NBS)], UMT)
    DO i=1,NBS
      DO j=1,NBS
        OVR(j,i) = UMT(i,j)
      END DO
    END DO

    WRITE(*,*) '2s Wavefunction in terms of new eigenstates'
    WRITE(*,*) (OVR(2,i),i=1,NBS)
  END SUBROUTINE print_diag_mtx_info

  SUBROUTINE calc_something(HAM, UMT, DIP, t, file_num)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: HAM(:,:), UMT(:,:), DIP(:,:), t
    INTEGER, INTENT(IN) :: file_num

    REAL(DBL) :: PRD(SIZE(HAM,1),SIZE(HAM,1))
    REAL(DBL) :: AR, AI, DR, DI, phs
    INTEGER :: i, j, k, NBS

    NBS = SIZE(HAM,1)

    DO k=1,NBS
      DO j=1,NBS
        AR = 0.D0
        AI = 0.D0
        DO i=1,NBS
          AR = AR + COSINE(HAM(i,i)*t)*UMT(k,i)*UMT(j,i)
          AI = AI + SINE(HAM(i,i)*t)*UMT(k,i)*UMT(j,i)
        END DO
        PRD(k,j) = AR*AR + AI*AI
      END DO
    END DO

! |phi_2t> = Sum_i u(2,i)exp(i eps_i t) |psi_i>
! <phi_2t| phi_2t> = sum_ij exp(i (eps_i-eps_j)t <psi_j|d|phi_i>*u(2,i)*u(2,j)
    DR = 0.D0
    DI = 0.D0
    DO i=1,NBS
      DO j=1,NBS
        phs = (HAM(i,i) - HAM(j,j))*t
        DR = DR + COSINE(phs)*UMT(2,i)*UMT(2,j)*DIP(i,j)
        DI = DI + SINE(phs)*UMT(2,i)*UMT(2,j)*DIP(i,j)
      END DO
    END DO
    WRITE(*,'(12F12.4)') t, (PRD(k,2),k=1,NBS), DR, DI
    WRITE(file_num,'(12F12.4)') t, (PRD(k,2),k=1,NBS), SQR(DR*DR + DI*DI)*1.D2 !Dipole
  END SUBROUTINE calc_something

  SUBROUTINE write_plot_instructions()
    IMPLICIT NONE
    OPEN(12,file='plot_directions')
    REWIND(12)
    WRITE(12,*)'gnuplot -p << END'
    WRITE(12,*)'set terminal png'
    WRITE(12,*)'set output "output.png"'
    WRITE(12,*)'p for [col=2:12] "PLOT" u 1:col w l'
    WRITE(12,*)'quit'
    WRITE(12,*)'END'
!         WRITE(12,12)
!12   FORMAT('p for [col=2:5] "PLOT" u 1:col w lp')    
    CLOSE(12)
  END SUBROUTINE write_plot_instructions

  SUBROUTINE open_plot()
    IMPLICIT NONE
    CALL SYSTEM('chmod +x plot_directions')
    CALL SYSTEM('./plot_directions')
    CALL SYSTEM('open output.png')
  END SUBROUTINE open_plot
END MODULE applications

