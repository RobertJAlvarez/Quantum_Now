PROGRAM Main
  USE FortranFunctions, ONLY: DBL, ABSO, DIV, SQR
  USE ArrayFunctions
  USE Applications
!  USE Practice
  IMPLICIT NONE

  REAL(DBL) :: AA(4,4), BB(4,4)             !INVERSE4x4
  REAL(DBL) :: F(4), G(4)                   !
  REAL(DBL) :: H(2,2), O(2,2), E(2)!, V(2,2) !J2x2 & JAC2BY2GEN
  INTEGER :: i, j, l, whichOption, OPT

  DO
    WRITE(*,'(/,A)') 'DIAG options:'
    WRITE(*,'(A)') '1. STARK', '2. PIR', '3. HMO', '4. BOX', '5. DIAG'
    WRITE(*,'(A)') 'n. Exit the program'
    !READ (*,*) OPT
    OPT = 6

    IF (OPT == 1) THEN
      CALL STARKDVR()
    ELSE IF (OPT == 2) THEN
      CALL RINGDVR()
    ELSE IF (OPT == 3) THEN
      CALL HMODVR()
    ELSE IF (OPT == 4) THEN
      CALL BOXDVR()
!    ELSE IF (OPT == 5) THEN
!      CALL DIAGVR()
    ELSE
      EXIT
    END IF
  END DO

  CALL SOLVEGEN()

  whichOption = 0

  IF (whichOption == 1) THEN
    WRITE(*,*) '2 by 2(1) or general(2)?'
    READ(*,*) whichOption

    IF (whichOption == 1) THEN
      WRITE(*,*) 'H = ?'
      READ(*,*) H(1,1), H(1,2), H(2,2)
      H(2,1) = H(1,2)

      DO
        WRITE(*,*) 'O = ?'
        READ(*,*) O(1,1), O(1,2), O(2,2)
        O(2,1) = O(1,2)

        IF (O(1,1) < 0.0D0 .OR. O(2,2) < 0.0D0) CYCLE
        EXIT
      END DO

      !CALL JAC2BY2GEN(H, O, V, E)

      WRITE(*,'(/,A,2F12.6)') 'Eigenvalues:', E(1), E(2)

!      DO i=1, 2
!        WRITE(*,100) E(i), (O(j,i) j=1,2)
!        100 FORMAT(3F12.6)
!      END DO
!      WRITE(*,*) O(1,1)*O(1,2) + O(2,1)*O(2,2)
!      WRITE(*,*) O(1,1)*O(1,1) + O(2,1)*O(2,1)
!      WRITE(*,*) O(2,2)*O(2,2) + O(1,2)*O(1,2)
    ELSE IF (whichOption == 2) THEN
      WRITE(*,*) 'Enter H values:'
      READ(*,*) H(1,1), H(2,2), H(1,2)
      H(2,1) = H(1,2)

      WRITE(*,*) 'Enter O values:'
      READ(*,*) O(1,1), O(2,2), O(1,2)
      O(2,1) = O(1,2)

!      CALL JAC2BY2GEN(H, O, V, E)
      WRITE(*,'(/,2F12.6)') E(1), E(2)

!      CALL DIAGDVR()
    END IF
  END IF

!  WRITE(*,*) 'Invert matrix?(1=yes, 0=no)'
!  READ(*,*) whichOption

  IF (whichOption == 1) THEN
    WRITE(*,*) 'Do you want to invert a symmetric matrix(1) or a regular matrix(2)?'
    READ(*,*) whichOption

    IF (whichOption == 1) THEN
      WRITE(*,*) 'A(1,1),A(1,2),A(1,3),A(1,4)'
      READ(*,*) (AA(1,j),j=1,4)
      WRITE(*,*) 'A(2,2),A(2,3),A(2,4)'
      READ(*,*) (AA(2,j),j=2,4)
      WRITE(*,*) 'A(3,3),A(1,2),A(1,3),A(1,4)'
      READ(*,*) (AA(1,j),j=3,4)
      WRITE(*,*) 'A(4,4)'
      READ(*,*) AA(1,4)

      DO i=1, 4
        DO j=i, 4
          AA(j,i) = AA(i,j)
        END DO
      END DO
    ELSE IF (whichOption == 2) THEN
      WRITE(*,*) 'Write your 4x4 matrix:'
      DO i=1, 4
        READ(*,*) (AA(i,j), j=1,4)
      END DO
    ELSE
      WRITE(*,*) 'Do you want to invert a symmetric matrix(1) or a regular matrix(2)?'
    END IF

!    CALL INVERSE4x4(AA, BB)

    DO l=1, 4
      WRITE(*,*) 'F=?'
      READ(*,*) (F(i),i=1,4)
      DO i=1, 4
        G(i) =0.0D0
        DO j=1, 4
          G(i) = G(i) + BB(i,j)*F(j)
        END DO
        WRITE(*,*) G(i)
      END DO
    END DO
  END IF

  CONTAINS
 
  SUBROUTINE SOLVEGEN()
    IMPLICIT NONE
    REAL*8, DIMENSION(2,2) :: H, O, V
    REAL*8 :: E(2)!, F(4), G(4)
 
    INTEGER, PARAMETER :: seed = 86456
    INTEGER :: ND

    CALL SRAND(seed)
    WRITE(*,*) 'Dimension = ?'
    READ(*,*) ND

    BLOCK
      REAL*8, DIMENSION(ND,ND) :: HAM, OVR, EVC, TST, RST, UPD
      REAL*8 :: ERR
      INTEGER :: i,j, k, l, m, n, ITST, ITRY, in, ip
 
      DO i=1, ND
        DO j=i, ND
          HAM(i,j) = rand() - 0.5D0
          OVR(i,j) = 0.3D0*(rand() - 0.0D0)
          HAM(j,i) = HAM(i,j)
          OVR(j,i) = OVR(i,j)
        END DO
        HAM(i,i) = 5.0D0*HAM(i,i)
        OVR(i,i) = 1.0D0 + rand()
        WRITE(*,'(12F10.3)') (OVR(i,j),j=1,ND)
      END DO

      DO i=1, ND
        DO j=1, ND
          EVC(i,j) = 0.0D0
        END DO
        EVC(i,i) = 1.0D0
      END DO

      DO ITRY=1, 100
        DO i=1, ND
          DO j=i+1, ND
            H = 0.0D0
            O = 0.0D0
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

            CALL JAC2BY2GEN(H,O,V,E)

            UPD = 0.0D0
            DO k=1, ND
              UPD(k,k) = 1.0D0
            END DO

            UPD(i,i) = V(1,1)
            UPD(i,j) = V(1,2)
            UPD(j,i) = V(2,1)
            UPD(j,j) = V(2,2)

            ! RST1(l,k) = V1_ml HAM(m,n) V1_nk
            ! RST2(q,p) = V2_lq RST!(l,k) V2_kp = V1_ml V2_lq [Hmn V1_nk V2_kp]_np
            ! V2_ji V1i
            TST = 0.0D0
            DO in=1,ND
              DO ip=1,ND
                DO k=1, ND
                  TST(in,ip) = TST(in,ip) + EVC(in,k)*UPD(k,ip)
                END DO
              END DO
            END DO

            EVC = TST
            DO ITST=1,2
              IF (ITST == 1) TST = OVR
              IF (ITST == 2) TST = HAM
              DO k=1, ND
                DO l=1, ND
                  RST(l,k) = 0.0D0
                  DO m=1, ND
                    DO n=1, ND
                      RST(l,k) = RST(l,k) + EVC(m,k)*EVC(n,l)*TST(m,n)
                    END DO
                  END DO
                END DO
                WRITE(*,'(12F10.3)') (RST(l,k),l=1,ND)
              END DO
              ERR = 0.0D0
              DO m=1, ND
                DO n=m+1,ND
                  ERR = ERR + RST(m,n)*RST(n,m)
                END DO
              END DO
              WRITE(*,*) 'ERR:', SQRT(ERR)
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
    END BLOCK
  END SUBROUTINE SOLVEGEN
END PROGRAM Main
