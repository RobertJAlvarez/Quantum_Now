MODULE arrayFunctions
  USE fortranFunctions, ONLY: DBL, ABSO, SQR, DIV
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: print_mtx, print_EV, INVERSE, J2X2, JAC2BY2GEN, DIAGNxN, LEASTSQUARE

  CONTAINS

  !Author: Robert Alvarez
  !Last modification: August 4th, 2022
  !  Change writing format from ES to F
  SUBROUTINE print_mtx(X)   !Print matrix of any size
    IMPLICIT NONE

    REAL(DBL), INTENT(IN) :: X(:,:)
    CHARACTER(len=9) :: fmt_mt
    INTEGER :: i, n_r, n_c

    n_r = SIZE(X,1)
    n_c = SIZE(X,2)

    !WRITE(fmt_mt, '( "(",I2,"ES17.8E3)" )' ) n_c !Change len of fmt_mt to 13
    WRITE(fmt_mt, '( "(",I2,"F11.7)" )' ) n_c
    WRITE(*,fmt_mt) (X(i,:n_c), i=1,n_r)
  END SUBROUTINE print_mtx

  SUBROUTINE print_EV(E, O)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: E(:), O(:,:)
    CHARACTER(len=18) :: fmt_mt
    INTEGER :: i, j, NBS

    NBS = SIZE(E)
    WRITE(fmt_mt, '( "(F10.6,A1,",I2,"F11.7)" )' ) NBS

    WRITE(*,*) 'EIGENVALUES AND EIGENVECTORS:'
    DO i=1, NBS
      WRITE(*,fmt_mt) E(i),":",(O(j,i),j=1,NBS)
    END DO
  END SUBROUTINE print_EV

  !Author (Fortran 77): Dr. Mark Pederson
  !Date:  September 8th, 2021
  !Modifier: Robert Alvarez / July 9th, 2022
  ! Modifications:
  !  Convert Subroutine to function
  !  Generalization of 4x4 matrix to nxn matrix
  !  Improve matrix printing format
  FUNCTION INVERSE(AA) RESULT(B)
    IMPLICIT NONE

    REAL(DBL), INTENT(IN) :: AA(:,:)

    REAL(DBL), DIMENSION(SIZE(AA,1),SIZE(AA,1)) :: B, CC
    REAL(DBL), DIMENSION(SIZE(AA,1),SIZE(AA,1)*2) :: A
    REAL(DBL) :: TMAX, temp
    INTEGER :: i, j, k, n

    n = SIZE(AA,1)
    A = 0.D0

    !Create matrix A = [AA|I]
    DO i=1, n
      A(i,i+n) = 1.D0       !Create identity matrix
      A(i,1:n) = AA(i,1:n)  !Copy values of AA into first nxn spaces of A
    END DO

    !Print matrix A
    WRITE(*,'(/,A)') '[Matrix | Identity]:'
    CALL print_mtx(A)

    !Find invert matrix using Gaussian elimination
    DO i=1, n
      !Find largest value of column i
      k = i
      TMAX = ABSO(A(k,i))
      DO j=i+1, n
        IF (ABSO(A(j,i)) > TMAX) THEN
          k = j
          TMAX = ABSO(A(k,i))
        END IF
      END DO

      TMAX = A(k,i)
      IF (TMAX < 1.D-15) STOP 'A row is linearly dependent of one or more other rows'

      IF (k /= i) THEN !Swap rows
        DO j=i, 2*n
          temp = A(i,j)
          A(i,j) = A(k,j)
          A(k,j) = temp
        END DO
      END IF

      !Make leading one
      DO j=i,2*n
        A(i,j) = DIV(A(i,j),TMAX)
      END DO

      !Clear numbers above and below A(j,i)
      DO j=1, n
        IF (j /= i) THEN
          A(j,i:2*n) = A(j,i:2*n) - A(j,i)*A(i,i:2*n)
        END IF
      END DO
    END DO

    !Copy inverse matrix into B
    B = A(:,n+1:2*n)

    !Print inverse matrix
    WRITE(*,'(/,A)') 'Inverse matrix:'
    CALL print_mtx(B)

    !Multiplication of A and A inverse = identity matrix
    CC = 0.D0
    DO i=1, n
      DO j=1, n
        DO k=1, n
          CC(i,j) = CC(i,j) + AA(i,k)*B(k,j)
        END DO
      END DO
    END DO

    !Print identity matrix, CC
    WRITE(*,'(/,A)') 'A*A^(-1):'
    CALL print_mtx(CC)
  END FUNCTION INVERSE

  !Author (Fortran 77): Dr. Mark Pederson
  !Date:  September 15th, 2021
  SUBROUTINE J2X2(H, E, O)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: H(:,:)         ! Hamiltonian matrix
    REAL(DBL), INTENT(OUT) :: E(:), O(:,:)  ! Eigenvalues - orthogonal eigenvectors
    REAL(DBL) :: vMag, trc, dis !vector magnitude, trace, discriminant
    INTEGER :: i

    dis = SQR((H(1,1) - H(2,2))*(H(1,1) - H(2,2)) + 4.D0*H(1,2)*H(2,1))
    trc = H(1,1) + H(2,2)

    !Calculate eigenvalues
    E(1) = (0.5D0)*(trc+dis)
    E(2) = (0.5D0)*(trc-dis)

    !Calculate eigenvectors
    DO i=1, 2
      O(1,i) = -H(1,2)
      O(2,i) =  H(1,1) - E(i)
    END DO

    !Normalize eigenvectors
    DO i=1, 2
      vMag = SQR(O(1,i)*O(1,i) + O(2,i)*O(2,i))
      O(1,i) = DIV(O(1,i),vMag)
      O(2,i) = DIV(O(2,i),vMag)
    END DO
  END SUBROUTINE J2X2

  !Author (Fortran 77): Dr. Mark Pederson
  !Date:  September 22nd, 2021
  SUBROUTINE JAC2BY2GEN(H, O, V, E)
    IMPLICIT NONE
    REAL(DBL), INTENT(INOUT) :: H(:,:), O(:,:)
    REAL(DBL), INTENT(OUT) :: V(:,:), E(:) ! Eigenvalues

    REAL(DBL) :: T(2,2)    ! Eigenvectors
    REAL(DBL) :: D(2,2)    ! Product matrix
    REAL(DBL) :: TRC, RAD, A, B, C
    INTEGER :: i, j, k, l, iTry

    A = O(1,1)*O(2,2) - O(1,2)*O(2,1)
    IF (A < 1.D-15) STOP 'Non positive overlap matrix'

    B = -( H(1,1)*O(2,2) - O(1,2)*H(2,1) + O(1,1)*H(2,2) - H(1,2)*O(2,1) )
    C = H(1,1)*H(2,2) - H(1,2)*H(2,1)

    TRC = DIV(-B,A)
    RAD = DIV(SQR(B*B - 4.D0*A*C),A)

    E(1) = (0.5D0)*(TRC + RAD)
    E(2) = (0.5D0)*(TRC - RAD)

!    WRITE(*,'(/,A,2F12.6,/)') 'Eigenvalues:', E(1), E(2)

    !Calculate eigenvectors
    DO k=1, 2
      DO i=1, 2
        DO j=1, 2
          T(i,j) = H(i,j) - E(k)*O(i,j)
        END DO
      END DO
      V(1,k) = -T(k,2)
      V(2,k) =  T(k,1)
!      WRITE(*,*) T(1,1)*T(2,2) - T(1,2)*T(2,1)
    END DO

    ! <V_1 | V_2> = 0 ?
    DO iTry=1, 3
!      IF (iTry <= 2) THEN
!        WRITE(*,'(/,A)') 'Overlap matrix'
!      ELSE
!        WRITE(*,'(/,A)') 'Hamiltonian matrix'
!      END IF

      DO i=1, 2
        DO j=1, 2
          D(i,j) = 0.D0
          DO k=1, 2
            DO l=1, 2
              IF (iTry <= 2) THEN
                D(i,j) = D(i,j) + V(k,i)*V(l,j)*O(k,l)
              ELSE
                D(i,j) = D(i,j) + V(k,i)*V(l,j)*H(k,l)
              END IF
            END DO
          END DO
        END DO
!        WRITE(*,'(2F12.6)') (D(i,j),j=1,2)
      END DO

      IF (iTry == 1) THEN
        DO i=1, 2
          DO k=1, 2
            V(k,i) = DIV(V(k,i),SQR(D(i,i)))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE JAC2BY2GEN
!
!Subroutine sort and MergeIdx are only used only for DIAGNxN
!
  !Author: Robert Alvarez
  !Last modification: August 1st, 2022
  !  Delete lower boundary value and change it to 1
  !  Stop specifying the lower and upper boundary when slicing
  RECURSIVE SUBROUTINE SortIdx(idx, PRD)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: idx(:,:)
    REAL(DBL), INTENT(IN) :: PRD(:,:)
    INTEGER :: mid, high

    high = SIZE(idx,2)
    IF (1 < high) THEN
      mid = (high+1)/2
      CALL SortIdx(idx(:,:mid), PRD)
      CALL SortIdx(idx(:,mid+1:), PRD)
      idx(:,:) = MergeIdx(idx(:,:mid), idx(:,mid+1:), PRD)
    END IF
  END SUBROUTINE SortIdx

  !Author: Robert Alvarez
  !Last modification: August 1st, 2022
  !  Change > for >= inside DO WHILE to make it stable
  !  Upper boundary parameters deleted and taken from SIZE(mtx,2) instead
  FUNCTION MergeIdx(a, b, PRD)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(IN) :: a, b
    REAL(DBL), INTENT(IN) :: PRD(:,:)

    INTEGER :: MergeIdx(2,SIZE(a,2)+SIZE(b,2))
    INTEGER :: ai, a_high, bi, b_high, ci

    ai = 1
    a_high = SIZE(a,2)
    bi = 1
    b_high = SIZE(b,2)
    ci = 1

    DO WHILE (ai <= a_high .AND. bi <= b_high)
      IF (ABSO(PRD(a(1,ai),a(2,ai))) >= ABSO(PRD(b(1,bi),b(2,bi)))) THEN
        MergeIdx(:,ci) = a(:,ai)
        ai = ai + 1
      ELSE
        MergeIdx(:,ci) = b(:,bi)
        bi = bi + 1
      END IF
      ci = ci + 1
    END DO

    IF (ai > a_high) THEN
      MergeIdx(:,ci:) = b(:,bi:)
    ELSE
      MergeIdx(:,ci:) = a(:,ai:)
    END IF
  END FUNCTION MergeIdx

  !Author (Fortran 77): Dr. Mark Pederson
  !Date:  October 20th, 2021
  !Modifier: Robert Alvarez / August 4th, 2022
  ! Modifications:
  !  Change NBS, PRD, and SPC from subroutine parameters to local variables
  !  Change goto to break
  !  Use Merge Sort for index sorting
  !  Place index sorting and selection on BLOCK construction
  !  Reorder of calculation blocks
  SUBROUTINE DIAGNxN(HAM, UMT)
    IMPLICIT NONE
    REAL(DBL), INTENT(INOUT) :: HAM(:,:)  !Hamiltonian matrix
    REAL(DBL), INTENT(OUT) :: UMT(:,:)    !Unitary matrix

    REAL(DBL), DIMENSION(SIZE(HAM,1),SIZE(HAM,1)) :: PRD, SPC  !Product matrix
    REAL(DBL) :: H(2,2), E(2), O(2,2) !Hamiltonian, Eigenvalues, and Eigenvectors for J2X2
    REAL(DBL) :: ERROLD, ERRNW        !Error values

    INTEGER :: i, j, k, l
    INTEGER :: NBS, iTry, MXIT  !NBS = dimension of matrix

    NBS = SIZE(HAM,1)

    UMT = 0.D0
    DO i=1, NBS
      UMT(i,i) = 1.D0
    END DO

    PRD = HAM
    MXIT = NBS*NBS*2

    DO iTry = 1, MXIT
      ERROLD = 0.D0
      DO i=1, NBS
        DO j=i+1, NBS
          ERROLD = ERROLD + PRD(j,i)*PRD(j,i)
        END DO
      END DO

      bestIdxs: BLOCK
        LOGICAL :: useIdx(NBS)
        INTEGER :: idxAll(2,NBS*(NBS-1)/2)
        INTEGER :: idx(2,NBS/2)
        INTEGER :: n, idxSize

        n = 0
        DO i=1, NBS
          DO j=i+1, NBS
            IF (ABSO(PRD(j,i)) > 1.D-10) THEN   !Save indices of PRD entries whose absolute value is greater than 0
              n = n + 1
              idxAll(1,n) = i
              idxAll(2,n) = j
            END IF
          END DO
        END DO

!        WRITE(*,*) 'All indexes before sorting:'
!        DO i=1, n
!         WRITE(*,*) idxAll(1,i), idxAll(2,i), PRD(idxAll(1,i),idxAll(2,i))
!        END DO

        CALL SortIdx(idxAll(:,:n), PRD)  !Sort with respect to PRD values

!        WRITE(*,*) 'All indexes:'
!        DO i=1, n
!          WRITE(*,*) idxAll(1,i), idxAll(2,i), PRD(idxAll(1,i),idxAll(2,i))
!        END DO

        !Choose the indexes with the highest values at PRD without repeating indexes
        useIdx = .TRUE. !No indexes has been use
        idxSize = 1     !Keep track of position to add new non-repetitive idxs
        DO j=1, n
          IF (useIdx(idxAll(1,j)) .AND. useIdx(idxAll(2,j))) THEN !If none of this two indexes has been used
            useIdx(idxAll(:,j)) = .False. !Set both to false because they would be used
            idx(:,idxSize) = idxAll(:,j)  !Save both indexes in idx array
            idxSize = idxSize + 1         !Update position for next indexes
            IF (idxSize > NBS/2) EXIT     !Exit the loop if NBS (even) or NBS-1 (odd) indexes had been used
          END IF
        END DO

!        WRITE(*,*) 'Non repetitive indexes with highest values:'
!        idxSize = idxSize - 1
!        DO i=1, idxSize
!          WRITE(*,*) idx(1,i), idx(2,i), PRD(idx(1,i),idx(2,i))
!        END DO

        !Use best two indexes
        k = idx(1,1)
        l = idx(2,1)
      END BLOCK bestIdxs

      H(1,1) = PRD(k,k)
      H(1,2) = PRD(k,l)
      H(2,1) = PRD(l,k)
      H(2,2) = PRD(l,l)
      CALL J2X2(H, E, O)

      !CALL print_EV(E, O)

      !Update unitary matrix
      PRD = UMT
      DO i=1, NBS
        PRD(i,k) = UMT(i,k)*O(1,1) + UMT(i,l)*O(2,1)
        PRD(i,l) = UMT(i,k)*O(1,2) + UMT(i,l)*O(2,2)
      END DO
      UMT = PRD

      SPC = 0.D0
      DO i=1, NBS
        SPC(i,i) = 1.D0
      END DO
      SPC(k,k) = O(1,1)
      SPC(l,k) = O(2,1)
      SPC(l,l) = O(2,2)
      SPC(k,l) = O(1,2)
      DO i=1, NBS
        DO j=1, NBS
          SPC(j,i) = 0.D0
          DO k=1, NBS
            SPC(j,i) = SPC(j,i) + UMT(k,i)*HAM(k,j)
          END DO
        END DO
      END DO

      !Get new Hamiltonian matrix
      PRD = 0.D0
      DO i=1, NBS
        DO j=1, NBS
          DO k=1, NBS
            PRD(j,i) = PRD(j,i) + UMT(k,j)*SPC(k,i)
          END DO
        END DO
      END DO

      !CALL print_mtx(PRD)

      ERRNW = 0.D0
      DO i=1,NBS
        DO j=i+1, NBS
          ERRNW = ERRNW + PRD(j,i)*PRD(j,i)
        END DO
      END DO

      WRITE(*,'(I3,2G15.6)') iTry, ERRNW, ERROLD

      IF (ERRNW < 1.D-15) EXIT
    END DO

    IF (iTry == MXIT) WRITE(*,*) 'Warning: No Convergence'

    WRITE(*,'(A,I3,A)') 'Diagonalization took', iTry, ' iterations.'
    WRITE(*,'(A,ES14.6E2)') 'Diagonalization efficienty:', DIV(DBLE(iTry),DBLE(NBS*NBS))

    HAM = PRD
  END SUBROUTINE DIAGNxN

  SUBROUTINE LEASTSQUARE()
    IMPLICIT NONE
    REAL(DBL) :: F(100)    !Function to evaluate to
    REAL(DBL) :: DM(4,4)   !Design matrix
    REAL(DBL) :: DU(4,4)   !unitary matrix
    REAL(DBL) :: DI(4,4)   !Inverse matrix
    REAL(DBL) :: B(4)      !
    REAL(DBL) :: A(4)      !
    REAL(DBL) :: P(100,4)  !
    REAL(DBL) :: xMin, xMax, dx, x, pow, fit, dat
    INTEGER :: i, j, k

    WRITE(*,*) 'xMin, xMax=?'
    READ(*,*) xMin, xMax

    dx = DIV(xMax - xMin, 99.D0)
    x = xMin - dx

    DO i=1, 100
      x = x + dx
      F(i) = x*x*x*x
      pow = 1.D0

      DO j=1, 4
        P(i,j) = pow
        pow = pow*x
      END DO
    END DO

    !Evaluate design matrix
    B = 0.D0
    DO j=1, 4
      DO k=1, 4
        DM(k,j) = 0.D0
        DO i=1, 100
          DM(k,j) = DM(k,j) + P(i,k)*P(i,j)
        END DO
      END DO

      DO i=1, 100
        B(j) = B(j) + F(i)*P(i,j)
      END DO
    END DO

    CALL DIAGNxN(DM, DU)

    DI = INVERSE(DM)

    A = 0.D0
    DO j=1, 4
      DO k=1, 4
        A(j) = A(j) + DI(j,k)*B(k)
      END DO
    END DO

    WRITE(*,*) 'A=', A

    DO i=1, 10
      WRITE(*,*) 'x=?'
      READ(*,*) x

      fit = 0.D0
      pow = 1.D0

      DO j=1, 4
        fit = fit + A(j)*pow
        pow = pow*x
      END DO

      dat = x*x*x*x
      WRITE(*,*) x, dat, fit
    END DO
  END SUBROUTINE LEASTSQUARE
END MODULE arrayFunctions
