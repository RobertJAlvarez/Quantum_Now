MODULE Applications
  USE FortranFunctions, ONLY: DBL, PI, ABSO, COSINE, SINE, DIV, SQR
  USE ArrayFunctions, ONLY: DIAGNxN

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE STARKDVR()
    IMPLICIT NONE
    INTEGER, PARAMETER :: NDH = 10
    REAL(DBL), DIMENSION(NDH,NDH) :: HAM, OVR, UMT, PRD, DIP
    REAL(DBL) :: EFIELD, G, X, P, TXR, t, DIPOLE
    REAL(DBL) :: AI, AR, DI, DR, phs, tau
    INTEGER :: i, j, k, m, NBS

    WRITE(*,*) 'WELCOME TO STARK DRIVER, EFIELD=?'
    READ(*,*) EFIELD

    G = -1.D0
    X = -0.5D0
    P = 0.1D0
    TXR = 0.2D0
    HAM = 0.D0
    OVR = 0.D0
    DIP = 0.D0
    DIP(1,3) = 0.01D0
    DIP(2,3) = 0.03D0

    DO i=1,3
      DO j=i+1,3
        DIP(j,i) = DIP(i,j)
      END DO
    END DO

    DO i=1,NDH
      OVR(i,i) = 1.D0
    END DO

    HAM(1,1) = -0.5D0
    HAM(2,2) = -0.125D0
    HAM(3,3) = -0.125D0
    HAM(1,3) = DIP(1,3)*EFIELD
    HAM(3,1) = HAM(1,3)
    HAM(2,3) = DIP(2,3)*EFIELD
    HAM(3,2) = HAM(2,3)
    NBS = 3

    CALL DIAGNxN(NDH,NBS,HAM,UMT,PRD)

    WRITE(*,*) "UPDATED HAM"
    DO i=1,NBS
      WRITE(*,"(10F12.4)") (HAM(j,i),j=1,NBS)  !LAMDA_J
    END DO

! PROVE THAT THE INVERSE OF UMT IS THE TRANSPOSE OF UMT
    WRITE(*,*) 'EIGENVALUES AND EIGENVECTORS:'
    DO i=1,NBS
      WRITE(*,"(10F12.4)") HAM(i,i),(UMT(j,i),j=1,NBS)
    END DO
    DO j=1,NBS
      DO k=1,NBS
        OVR(k,j) = UMT(j,k)
      END DO
    END DO

    WRITE(*,*) '2s Wavefunction in terms of new eigenstates'
    WRITE(*,*) (OVR(2,k),k=1,3)
! AT TIME=0, OCCUPY THE 2S FUNCTION:
! |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(K,J)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
    TAU = DIV(8.D0*PI,ABSO(HAM(1,1)))
    DO M=0,1000
      t = DBLE(M)*DIV(TAU,1.D3)
      OPEN(12,FILE='PLOT')
      DO k=1,3
        DO j=1,3  
          AR = 0.D0
          AI = 0.D0
          DO i=1,3
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
      DO i=1,3
        DO j=1,3
          phs = (HAM(i,i) - HAM(j,j))*t
          DR = DR + COSINE(phs)*UMT(2,i)*UMT(2,j)*DIP(i,j)
          DI = DI + SINE(phs)*UMT(2,i)*UMT(2,j)*DIP(i,j) 
        END DO
      END DO
      DIPOLE = SQR(DR*DR + DI*DI)*1.D3
      WRITE(*,'(10F12.4)') t,(PRD(k,2),k=1,NBS),DR,DI             
      WRITE(12,'(10F12.4)') t,(PRD(k,2),k=1,NBS),DIPOLE               
    END DO
    CLOSE(12)
  
    OPEN(12,file='plot_directions')
    REWIND(12)
    WRITE(12,*)'gnuplot -p << END'
    WRITE(12,*)'set terminal png'
    WRITE(12,*)'set output "output.png"'
    WRITE(12,*)'p for [col=2:5] "PLOT" u 1:col w lp'
    WRITE(12,*)'quit'
    WRITE(12,*)'END'
!   WRITE(12,12)
!  12 FORMAT('p for [col=2:5] "PLOT" u 1:col w lp')    
    CLOSE(12)
    CALL SYSTEM('chmod +x plot_directions')
    CALL SYSTEM('./plot_directions')
    CALL SYSTEM('eog output.png')
    STOP
  END SUBROUTINE STARKDVR

  SUBROUTINE RINGDVR()
    IMPLICIT NONE
    INTEGER, PARAMETER :: NDH = 12
    REAL(DBL), DIMENSION(NDH,NDH) :: HAM, OVR, UMT, PRD, DIP
    REAL(DBL) :: RINGSZ, E, aM, phs, DIPOLE
    REAL(DBL) :: AI, AR, DI, DR, t, TAU
    INTEGER :: i, j, k, m, NBS, iAP

    WRITE(*,*) 'HOW LARGE IS YOUR RING?'
    READ(*,*)  RINGSZ
    WRITE(*,*) 'ELECTRIC(1) OR MAGNETIC(2)?'
    READ(*,*) iAP
    WRITE(*,*) 'HOW BIG THE FIELD?'
    READ(*,*) E

    DIP = 0.D0
    DO i=1,NDH
      OVR(i,i) = 1.D0
    END DO

    NBS = 0
    HAM = 0.D0
    aM = -5.D0
    DO m=-5,5 
      NBS = NBS + 1

      IF (NBS > 1) THEN
        DIP(NBS-1,NBS) = DIV(E*RINGSZ,2.D0)
        DIP(NBS,NBS-1) = DIP(NBS-1,NBS)
!       HAM(NBS-1,NBS) = DIV(E*RINGSZ,2.D0)
!       HAM(NBS,NBS-1) = HAM(NBS-1,NBS)
      END IF
      HAM(NBS,NBS) = DIV(aM*aM,RINGSZ*RINGSZ)
      aM = aM + 1.D0
    END DO

    IF (iAP == 2)THEN
      DIP = 0.D0
      NBS = 0
!     E=-1*DIV(2.D0,RINGSZ*RINGSZ)
      DO m=-5,5
        NBS = NBS + 1
        DIP(NBS,NBS) = DBLE(M)*E 
      END DO
    END IF

    DO i=1,NBS
      WRITE(*,'(11F7.2)') (HAM(i,j),j=1,NBS)
    END DO

    DIP = DIP + HAM
    CALL DIAGNxN(NDH,NBS,HAM,UMT,PRD)
    WRITE(*,*) "UPDATED HAM"
    DO I=1,NBS
      WRITE(*,'(10F12.4)') (HAM(J,I),J=1,NBS)  !LAMDA_J
    END DO
    WRITE(*,*) 'EIGENVALUES AND EIGENVECTORS:'
    DO I=1,NBS
      WRITE(*,'(10F12.4)') HAM(I,I),(UMT(J,I),J=1,NBS)
! PROVE THAT THE INVERSE OF UMT IS THE TRANPOSE OF UMT
    END DO
    DO J=1,NBS
      DO K=1,NBS
        OVR(K,J) = UMT(J,K)
      END DO
    END DO
    WRITE(*,*) '2s Wavefunction in terms of new eigenstates'
    WRITE(*,*) (OVR(2,K),K=1,NBS)
! AT TIME=0, OCCUPY THE 2S FUNCTION:
! |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(K,J)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
    TAU = DIV(8.D0*PI,ABSO(HAM(1,1)))
    OPEN(12,FILE='PLOT')
    DO M=-100,100000
      IF (M == 0)THEN
        HAM = DIP
        OVR = 0.D0
        DO I=1,NBS
          OVR(I,I) = 1.D0
        END DO
        CALL DIAGNxN(NDH,NBS,HAM,UMT,PRD)
        DO I=1,NBS
          WRITE(*,*) HAM(I,I)
        END DO
      END IF
      t = DBLE(M)*DIV(TAU,5.D1)
      DO K=1,NBS
        DO j=1,NBS
          AR = 0.D0
          AI = 0.D0
          DO I=1,NBS
            AR = AR + COSINE(HAM(i,i)*t)*UMT(k,i)*UMT(j,i) 
            AI = AI + SINE(HAM(i,i)*t)*UMT(k,i)*UMT(j,i) 
          END DO
          PRD(K,j) = AR*AR + AI*AI
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
      DIPOLE = SQR(DR*DR + DI*DI)*1.D2
      WRITE(*,'(12F12.4)') t,(PRD(K,2),K=1,NBS),DR,DI             
!     WRITE(12,"(3F12.4)")t,PRD(2,2),PRD(6,2)!,K=1,NBS)!,DIPOLE               
      WRITE(12,'(12F12.4)') t,(PRD(K,2),K=1,NBS)!,DIPOLE               
    END DO
    CLOSE(12)

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
    CALL SYSTEM('chmod +x plot_directions')
    CALL SYSTEM('./plot_directions')
    CALL SYSTEM('eog output.png')
    STOP
  END SUBROUTINE RINGDVR

  SUBROUTINE BOXDVR()
    IMPLICIT NONE
    INTEGER, PARAMETER :: NDH = 10
    REAL(DBL), DIMENSION(NDH,NDH) :: HAM, OVR, UMT, PRD, DIP
    REAL(DBL) :: alpha, a, twom, tau, dr, di, ai, ar, boxsz, dipole, phs, t
    INTEGER :: i, j, k, m, NBS

    WRITE(*,*) 'WELCOME TO BOX DRIVER, HOW LARGE IS YOUR BOX?'
    READ(*,*) BOXSZ
    Alpha = DIV(2*PI,BOXSZ)
    A = DBLE(M*M)*PI
    DO I=1,NDH
      OVR(I,I) = 1.D0
    END DO

    NBS = 0
    HAM = 0.D0
    DO M=-4,4
      NBS = NBS + 1
      TWOM = 2.D0*DBLE(M)*PI
      WRITE(*,*) M,TWOM
      HAM(NBS,NBS) = DIV(TWOM,BOXSZ)*DIV(TWOM,BOXSZ)
    END DO

    WRITE(*,*) "INITIAL HAM"
    DO I=1,NBS
      WRITE(*,'(10F12.4)') (HAM(J,I),J=1,NBS)  !LAMDA_J
    END DO
    CALL DIAGNxN(NDH,NBS,HAM,UMT,PRD)

    WRITE(*,*) "UPDATED HAM"
    DO I=1,NBS
      WRITE(*,'(10F12.4)') (HAM(J,I),J=1,NBS)  !LAMDA_J
    END DO

    WRITE(*,*) 'EIGENVALUES AND EIGENVECTORS:'
    DO I=1,NBS
      WRITE(*,'(10F12.4)') HAM(I,I),(UMT(J,I),J=1,NBS)
! PROVE THAT THE INVERSE OF UMT IS THE TRANPOSE OF UMT
    END DO
    DO J=1,NBS
      DO K=1,NBS
        OVR(K,J)=UMT(J,K)
      END DO
    END DO

    WRITE(*,*) '2s Wavefunction in terms of new eigenstates'
    WRITE(*,*) (OVR(2,K),K=1,3)
! AT TIME=0, OCCUPY THE 2S FUNCTION:
! |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(K,J)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
    TAU = DIV(8.D0*PI,ABSO(HAM(1,1)))
    TAU = DIV(1.D0,ABSO(DIV(HAM(1,1),4.D0)))
    DO M=0,1000
      t = DBLE(M)*DIV(TAU,1.D3)
      OPEN(12,FILE='PLOT')
      DO K=1,3
        DO J=1,3
          AR = 0.D0
          AI = 0.D0
          DO I=1,3
            AR = AR + COSINE(HAM(i,i)*t)*UMT(k,i)*UMT(j,i)
            AI = AI + SINE(HAM(i,i)*t)*UMT(k,i)*UMT(j,i)
          END DO
          PRD(K,J) = AR*AR + AI*AI
        END DO
      END DO
! |phi_2t> = Sum_i u(2,i)exp(i eps_i t) |psi_i>
! <phi_2t| phi_2t> = sum_ij exp(i (eps_i-eps_j)t <psi_j|d|phi_i>*u(2,i)*u(2,j)
      DR = 0.D0
      DI = 0.D0
      DO i=1,3
        DO j=1,3
          phs = (HAM(i,i) - HAM(j,j))*t
          DR = DR + COSINE(phs)*UMT(2,i)*UMT(2,j)*DIP(i,j)
          DI = DI + SINE(phs)*UMT(2,i)*UMT(2,j)*DIP(i,j)
        END DO
      END DO
      DIPOLE = SQR(DR*DR + DI*DI)*1.D3
      WRITE(*,'(10F12.4)') t,(PRD(K,2),K=1,NBS),DR,DI
      WRITE(12,'(10F12.4)') t,(PRD(K,2),K=1,NBS),DIPOLE
    END DO
    CLOSE(12)

    OPEN(12,file='plot_directions')
    REWIND(12)
    WRITE(12,*)'gnuplot -p << END'
    WRITE(12,*)'set terminal png'
    WRITE(12,*)'set output "output.png"'
    WRITE(12,*)'p for [col=2:5] "PLOT" u 1:col w lp'
    WRITE(12,*)'quit'
    WRITE(12,*)'END'
!         write(12,12)
!12   FORMAT('p for [col=2:5] "PLOT" u 1:col w lp')
    CLOSE(12)
    CALL SYSTEM('chmod +x plot_directions')
    CALL SYSTEM('./plot_directions')
    CALL SYSTEM('xdg-open output.png')
    STOP
  END SUBROUTINE BOXDVR

  SUBROUTINE HMODVR()
    IMPLICIT NONE
    INTEGER, PARAMETER :: NDH = 10
    REAL(DBL), DIMENSION(NDH,NDH) :: HAM, OVR, UMT, PRD, DIP
    REAL(DBL) :: rkmx, rkmn, emax, emin, energy, hmass, ai, ar, di, dr
    REAL(DBL) :: dipole, efield, guess, phs, scale, t, tau
    INTEGER :: i, j, k, m, NBS

    WRITE(*,*) 'NRLMOL:',DIV(0.36005D0*1.6D-19,6.626D-34)
    rkmx = 1.D30
    rkmn = 0.D0
    emax = 1.D30
    emin = 0.D0
    WRITE(*,*) 'WELCOME TO HARMONIC OSCILLATOR DRIVER, EFIELD=?'
    WRITE(*,*) 'HCl Diatomic in and Electric Field' 
    WRITE(*,*) 'Clorine atom has infinite mass'
    WRITE(*,*) 'IR frequency of HCl molecule is: 8.88*10^13 Hz' 
    energy = 8.88D13
    WRITE(*,*) 'energy:', energy
    WRITE(*,*) 'Please googlge this to determine if this is correct'
    WRITE(*,*) 'The mass of the Hydrogen is 1.67*10^{-27} kg'
    WRITE(*,*) 'Strength of Electric field?'
    hmass = 1.67D-27!*(35.)/(36.)
    WRITE(*,*) EField
    DO k=1,30
      WRITE(*,*) 'Guess the spring constant in (Newtons/Meter)'
      WRITE(*,*) 'Guess should be between:',rkmn,' and ',rkmx
      READ(*,*) guess
      DIP(1,2) = 1.D-2
      DIP(2,1) = DIP(1,2)
      OVR = 0.D0
      DO I=1,2
        DO J=I+1,2
          DIP(J,I) = DIP(I,J)
        END DO
      END DO
      DO I=1,NDH
        OVR(I,I) = 1.D13
      END DO
      HAM = 0.D0
      WRITE(*,'(4G15.6)') guess, hmass, DIV(guess,hmass),SQR(DIV(guess,hmass))
      scale = SQR(DIV(guess,hmass))
      HAM(1,1) = DIV(0.5D0*SQR(DIV(guess,hmass)),(2.D0*PI))
      HAM(2,2) = DIV(1.5D0*SQR(DIV(guess,hmass)),(2.D0*PI))
      HAM(1,2) = 0.D0 !DIP(1,2)
      HAM(2,1) = 0.D0 !DIP(2,1)
      NBS = 2
      DO I=1,NBS
        WRITE(*,*) (HAM(I,J),J=1,NBS)
      END DO
      DO I=1,NBS
        DO J=1,NBS
          HAM(I,J) = DIV(HAM(I,J),scale)
        END DO
      END DO
      CALL DIAGNxN(NDH,NBS,HAM,UMT,PRD)
      DO I=1,NBS
        HAM(I,I) = HAM(I,I)*scale
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
    WRITE(*,*) 'UPDATED HAM'
    DO I=1,NBS
      WRITE(*,'(10F12.4)') (HAM(J,I),J=1,NBS)  !LAMDA_J
    END DO
    WRITE(*,*) 'EIGENVALUES AND EIGENVECTORS:'
    DO I=1,NBS
      WRITE(*,'(10F12.4)') HAM(I,I),(UMT(J,I),J=1,NBS)
! PROVE THAT THE INVERSE OF UMT IS THE TRANPOSE OF UMT
    END DO
    DO J=1,NBS
      DO K=1,NBS
        OVR(K,J) = UMT(J,K)
      END DO
    END DO 
    WRITE(*,*) '2s Wavefunction in terms of new eigenstates'
    WRITE(*,*) (OVR(2,K),K=1,3)
! AT TIME=0, OCCUPY THE 2S FUNCTION:
! |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(K,J)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
    TAU = DIV(8.D0*PI,ABSO(HAM(1,1)))
    DO M=0,1000
      t = DBLE(M)*DIV(TAU,1.D3)
      OPEN(12,FILE='PLOT')
      DO k=1,3
        DO j=1,3  
          AR = 0.D0
          AI = 0.D0
          DO i=1,3
            AR = AR + COSINE(HAM(i,i)*t)*UMT(k,i)*UMT(j,i) 
            AI = AI + SINE(HAM(i,i)*t)*UMT(k,i)*UMT(j,i) 
          END DO
          PRD(k,j) = AR*AR + AI*AI
        END DO
      END DO
!   |phi_2t> = Sum_i u(2,i)exp(i eps_i t) |psi_i>
!   <phi_2t| phi_2t> = sum_ij exp(i (eps_i-eps_j)t <psi_j|d|phi_i>*u(2,i)*u(2,j)
      DR = 0.D0       
      DI = 0.D0
      DO i=1,3
        DO j=1,3
          phs = (HAM(i,i) - HAM(j,j))*t
          DR = DR + COSINE(phs)*UMT(2,i)*UMT(2,j)*DIP(i,j)
          DI = DI + SINE(phs)*UMT(2,i)*UMT(2,j)*DIP(i,j) 
        END DO
      END DO
      DIPOLE = SQR(DR*DR + DI*DI)*1.D3
      WRITE(*,'(10F12.4)') t,(PRD(K,2),K=1,NBS),DR,DI             
      WRITE(12,'(10F12.4)') t,(PRD(K,2),K=1,NBS),DIPOLE               
    END DO
    CLOSE(12)
    OPEN(12,file='plot_directions')
    WRITE(12,12)
    12 FORMAT('p for [col=2:5] "PLOT" u 1:col w lp')
    CLOSE(12)
    CALL system('gnuplot <plot_directions')
    STOP
  END SUBROUTINE HMODVR
END MODULE Applications

