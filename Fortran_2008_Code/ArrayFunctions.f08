MODULE ArrayFunctions
  USE FortranFunctions, ONLY: DBL, ABSO, SQR, DIV
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE LEASTSQUARE()
    IMPLICIT NONE
    REAL(DBL) :: F(100)    !Function to evaluate to
    REAL(DBL) :: DM(4,4)   !Design matrix
    REAL(DBL) :: DI(4,4)   !Identity matrix
    REAL(DBL) :: B(4)      !
    REAL(DBL) :: A(4)      !
    REAL(DBL) :: P(100,4)  !
    REAL(DBL) :: xMin, xMax, dx, x, pow, fit, dat
    INTEGER :: i, j, k

    WRITE(*,*) 'xMin, xMax=?'
    READ(*,*) xMin, xMax

    dx = DIV(xMax - xMin, 99.0D0)
    x = xMin - dx

    DO i=1, 100
      x = x + dx
      F(i) = x*x*x*x
      pow = 1.0D0

      DO j=1, 4
        P(i,j) = pow
        pow = pow*x
      END DO
    END DO

    !Evaluate design matrix
    B = 0.0D0
    DO j=1, 4
      DO k=1, 4
        DM(k,j) = 0.0D0
        DO i=1, 100
          DM(k,j) = DM(k,j) + P(i,k)*P(i,j)
        END DO
      END DO

      DO i=1, 100
        B(j) = B(j) + F(i)*P(i,j)
      END DO
    END DO

!    CALL DIAGNxN(DM)

    STOP 'Before INVERSE4x4'

    CALL INVERSE4x4(DM, DI)

    A = 0.0D0
    DO j=1, 4
      DO k=1, 4
        A(j) = A(j) + DI(j,k)*B(k)
      END DO
    END DO

    WRITE(*,*) 'A=', A

    DO i=1, 10
      WRITE(*,*) 'x=?'
      READ(*,*) x

      fit = 0.0D0
      pow = 1.0D0

      DO j=1, 4
        fit = fit + A(j)*pow
        pow = pow*x
      END DO

      dat = x*x*x*x
      WRITE(*,*) x, dat, fit
    END DO
  END SUBROUTINE LEASTSQUARE

  SUBROUTINE INVERSE4x4(AA, BB)
    IMPLICIT NONE

    REAL(DBL), INTENT(IN) :: AA(:,:)   !AA holds matrix 1
    REAL(DBL), INTENT(OUT) :: BB(:,:)  !BB hold matrix 2

!    INTEGER :: n = SIZE(AA)
    REAL(DBL) :: A(4,8)                !A holds matrix 1 follow by identity matrix
    REAL(DBL) :: CC(4,4)               !CC holds identity matrix
    REAL(DBL) :: TMAX, temp
    INTEGER :: i, j, k

    A = 0.0D0

    !Create matrix 1
    DO i=1, 4
      A(i,i+4) = 1.0D0    !Create identity matrix

      !Set next columns of A to be the value of the previous column times x(i)
      DO j=1, 4
        A(i,j) = AA(i,j)
      END DO
    END DO

    !Output matrix 1 and identity matrix after it
    WRITE(*,*)
    WRITE(*,10) (A(i,1:8), i=1,4)
    10 FORMAT (8F12.6)

    !Find invert matrix
    DO i=1, 4
      !Find largest value of column i
      TMAX = ABSO(A(i,i))
      k = i
      DO j=i+1, 4
        IF (ABSO(A(j,i)) > TMAX) THEN
          k = j
          TMAX = ABSO(A(j,i))
        END IF
      END DO

      TMAX = A(k,i)
      IF (ABSO(TMAX) < 1.0D15) THEN
        WRITE(*,*) 'A row is linearly dependent of one or more other rows'
        STOP
      END IF

      !Swap row with highest value in column i
      IF (k /= i) THEN
        DO j=i, 8
          temp = A(i,j)
          A(i,j) = A(k,j)
          A(k,j) = temp
        END DO
      END IF

      !Normalize matrix
      DO j=i,8
        A(i,j) = A(i,j) / TMAX
      END DO
      !A(i,i:8) = A(i,i:8) / TMAX

      !Subtract value A(j,i) to every column at i
      DO j=1, 4
        temp = A(j,i)
        IF (j /= i) THEN
          A(j,i:8) = A(j,i:8) - temp*A(i,i:8)
        END IF
      END DO
    END DO

    !Copy matrix 2 to BB
    BB(1:4,1:4) = A(1:4,5:8)

    !Print matrix 2
    WRITE(*,*)
    WRITE(*,20) (BB(i,1:4), i=1, 4)
    20 FORMAT (4F12.6)

    !Dot matrix 1 with 2 to get CC, identity matrix
    CC = 0.0D0
    DO i=1, 4
      DO j=1, 4
        DO k=1, 4
          CC(i,j) = CC(i,j) + AA(i,k)*BB(k,j)
        END DO
      END DO
    END DO

    !Print identity matrix, CC
    WRITE(*,*)
    WRITE(*,20) (CC(i,1:4), i=1, 4)
  END SUBROUTINE INVERSE4x4

  SUBROUTINE J2X2(H, E, O)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: H(:,:)
    REAL(DBL), INTENT(OUT) :: E(:), O(:,:) !O is for orthogonal matrix
    REAL(DBL) :: dot, trc, rad
    INTEGER :: i

!    det(H) = H(1,1)*H(2,2) - H(1,2)*H(2,1)
!    An eigenvalue satisfies:
!    (H(1,1)-E)*(H(2,2)-E) -H(1,2)*H(2,1) = 0
!    E**2 - E*(H(1,1)+H(2,2)) - H(1,2)*H(2,1) = 0
!    E = (H(1,1) + H(2,2) +- SQRT((H(1,1) + H(2,2))**2 - 4*H(1,1)*H(2,2) + 4*H(1,2)*H(2,1))) / 2
!    E = (H(1,1) + H(2,2)) / 2 +- SQRT((H(1,1) + H(2,2))**2 + 4*H(1,2)*H(2,1)) / 2

    rad = SQR((H(1,1) - H(2,2))*(H(1,1) - H(2,2)) + 4.0D0*H(1,2)*H(2,1))
    trc = H(1,1) + H(2,2)
    E(1) = (0.5D0)*(trc+rad)
    E(2) = (0.5D0)*(trc-rad)

!The eigenvector O_n follows that (H-IE(n))O_n=0.
!So: (H(1,1)-E(n))O_n(1)+H(1,2)O_n(2)=0
!     H(2,1)O_n(1)+H(2,2)-E(n))O_n(2)=0
!But: H(2,1)=H(1,2)
!We can do: O_n(1)=-H(1,2) ; O_n(2)=(H(1,1)-E(n))
!This way the first equation is zero, and the second will have the form:
!-H(1,2)^2 + H(1,2)H(2,2)-E(n)(H(1,1)+H(2,2)) + E(n)^2 = 0
!Notice that this equation is equal to the determinant, which is equal to zero.
!So our solutions are valid, and any non zero scalar multiples of this O_n vectors.
!Notice that we can do the same using the second equation to obtain O_n(1) and O_n(2)
!And we will have that: O_n(1)=(H(2,2)-E(n)) and O_n(2)=-H(1,2), so
!       Left handed
    O(1,1) = -H(1,2)
    O(2,1) =  H(1,1) - E(1)
    O(1,2) =  H(2,2) - E(2)
    O(2,2) = -H(1,2)

    DO i=1, 2
      dot = SQR(O(1,i)*O(1,i) + O(2,i)*O(2,i))
      O(1,i) = DIV(O(1,i),dot)
      O(2,i) = DIV(O(2,i),dot)
    END DO
  END SUBROUTINE J2X2

  SUBROUTINE JAC2BY2GEN(H, O, V, E)
    IMPLICIT NONE
    REAL(DBL), INTENT(INOUT) :: H(:,:), O(:,:), V(:,:)
    REAL(DBL), INTENT(OUT) :: E(:) !Eigenvalues

    REAL(DBL) :: T(2,2)    !Eigenvectors
    REAL(DBL) :: D(2,2)    !Dot product matrix
    REAL(DBL) :: TRC, RAD, A, B, C
    INTEGER :: i, j, k, l, iTry

    A = O(1,1)*O(2,2) - O(1,2)*O(2,1)
    IF (A < 0.0D0) THEN
      STOP 'Non positive overlap matrix'
    END IF

    B = H(1,1)*O(2,2) - O(1,2)*H(2,1) + O(1,1)*H(2,2) - H(1,2)*O(2,1)
    B = -B
    C = H(1,1)*H(2,2) - H(1,2)*H(2,1)

    TRC = DIV(-B,A)
    RAD = DIV(SQR(B*B - 4.0D0*A*C),A)

    E(1) = (0.5D0)*(TRC + RAD)
    E(2) = (0.5D0)*(TRC - RAD)

    WRITE(*,'(/,A,2F12.6,/)') 'Eigenvalues:', E(1), E(2)

    !Calculate eigenvectors
    DO k=1, 2
      DO i=1, 2
        DO j=1, 2
          T(i,j) = H(i,j) - E(k)*O(i,j)
        END DO
      END DO
      V(1,k) = -T(k,2)
      V(2,k) =  T(k,1)
      WRITE(*,*) T(1,1)*T(2,2) - T(1,2)*T(2,1)
    END DO

    !C <V_1 | V_2> = 0 ?
    DO iTry=1, 3
      IF (iTry <= 2) THEN
        WRITE(*,'(/,A)') 'Overlap matrix'
      ELSE
        WRITE(*,'(/,A)') 'Hamiltonian matrix'
      END IF

      DO i=1, 2
        DO j=1, 2
          D(i,j) = 0.0D0
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
        WRITE(*,'(2F12.6)') (D(i,j),j=1,2)
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

  SUBROUTINE DIAGDVR()  !Diag driver
    IMPLICIT NONE
    INTEGER, PARAMETER :: NDH = 10
    REAL(DBL) :: HAM(NDH,NDH)  !Hamiltonian
    REAL(DBL) :: OVR(NDH,NDH)  !Overlap
    REAL(DBL) :: UMT(NDH,NDH)  !Unitary matrix
    REAL(DBL) :: SPC(NDH,NDH)  !Space
    REAL(DBL) :: PRD(NDH,NDH)  !Product

    REAL(DBL) :: G, X, P, TXR
    INTEGER :: i, j, NBS

    G = -1.0D0  !Ground
    X = -0.5D0  !Exited
    P =  0.1D0  !Perturbation
    TXR = 0.2D0 !Transfer

    OVR = 0.0D0
    HAM = 0.0D0

    DO i=1,NDH
      OVR(i,i) = 1.0D0
    END DO

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

    NBS = 7

    DO i=1, NBS
      DO j=i, NBS
        HAM(j,i) = HAM(i,j)
      END DO
      WRITE(*,'(10F7.2)') (HAM(i,j),j=1,NBS)
    END DO

    CALL DIAGNxN(NDH, NBS, HAM, OVR, UMT, PRD, SPC)

    WRITE(*,*) 'Updated Hamiltonian:'
    DO i=1,NBS
      WRITE(*,'(10F12.4)') (HAM(j,i), j=1,NBS)
    END DO
    STOP
  END SUBROUTINE DIAGDVR

  SUBROUTINE DIAGNxN(NDH, NBS, HAM, OVR, UMT, PRD, SPC)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDH, NBS
    REAL(DBL), INTENT(INOUT) :: HAM(:,:), OVR(:,:)
    REAL(DBL), INTENT(OUT) :: UMT(:,:), PRD(:,:), SPC(:,:)

    REAL(DBL) :: H(2,2), E(2), O(2,2)!, V(2,2), T(2,2), D(2,2)  !J2X2
    REAL(DBL) :: ERRPREV, ERRNW

    LOGICAL :: useIdx(NDH)
    INTEGER :: idxAll(2,NDH*(NDH-1)/2)
    INTEGER :: idx(2,NDH/2)
    INTEGER :: i, j, k, l, m, n, idxSize
    INTEGER :: iTry, MXIT

    IF (NBS > NDH) THEN
      WRITE(*,*) 'NDH must be larger than ', NBS
      STOP
    END IF

    UMT = 0.0D0
    DO i=1, NBS
      UMT(i,i) = 1.0D0
    END DO

    PRD = HAM
    MXIT = NBS*NBS*2

    DO iTry = 1, MXIT
      n = 0
      ERRPREV = 0.0D0
      DO i=1, NDH
        DO j=i+1, NDH
          ERRPREV = ERRPREV + PRD(i,j)*PRD(j,i)
          IF (ABSO(PRD(i,j)) .GT. 1.0D-10) THEN ! Save all indices that have values in PRD that are greater than 0
            n = n + 1
            idxAll(1,n) = i
            idxAll(2,n) = j
          END IF
        END DO
      END DO

!      WRITE(*,*) 'All indexes before sorting:'
!      DO i=1, n
!        WRITE(*,*) idxAll(1,i), idxAll(2,i), PRD(idxAll(1,i),idxAll(2,i))
!      END DO

      CALL sort(idxAll(:,1:n), PRD, NDH*(NDH-1)/2)  !Sort with respect to PRD values

!STOP 'Hi :)'
!      WRITE(*,*) 'All indexes:'
!      DO i=1, n
!        WRITE(*,*) idxAll(1,i), idxAll(2,i), PRD(idxAll(1,i),idxAll(2,i))
!      END DO

      !Choose the indexes with the highest values at PRD without repeating indexes
      useIdx = .TRUE. !No indexes has been use
      idxSize = 1     !Keep track of position to add new non-repetitive idxs
      DO j=1, n
        IF (useIdx(idxAll(1,j)) .AND. useIdx(idxAll(2,j))) THEN !If any of this two indexes has been used
          useIdx(idxAll(:,j)) = .False. !Set both to false because they would be used
          idx(:,idxSize) = idxAll(:,j)  !Save both indexes in idx array
          idxSize = idxSize + 1         !Update position for next indexes
          IF(idxSize > NDH/2) EXIT      !Exit the loop if NDH (even) or NDH-1 (odd) indexes had been used
        END IF
      END DO

!      WRITE(*,*) 'Non repetitive indexes with highest values:'
!      idxSize = idxSize - 1
!      DO i=1, idxSize
!        WRITE(*,*) idx(1,i), idx(2,i), PRD(idx(1,i),idx(2,i))
!      END DO
!      STOP 'Hi:)'

      !Use best two indexes
      k = idx(1,1)
      l = idx(2,1)

      H(1,1) = PRD(k,k)
      H(1,2) = PRD(k,l)
      H(2,1) = PRD(l,k)
      H(2,2) = PRD(l,l)
      CALL J2X2(H, E, O)

!      WRITE(*,*) 'E and O values:'
!      DO i=1, 2
!        WRITE(*,*) E(i), (O(i,j),j=1,2)
!      END DO

      SPC = 0.0D0
      DO i=1, NBS
        SPC(i,i) = 1.0D0
      END DO

      SPC(k,k) = O(1,1)
      SPC(l,k) = O(2,1)
      SPC(l,l) = O(2,2)
      SPC(k,l) = O(1,2)

      !Get new unitary matrix
      PRD = 0.0D0
      DO m=1, NBS
        IF (m /= k .AND. m /= l) THEN
          DO n=1, NBS
            PRD(n,m) = UMT(n,m)
          END DO
        END IF
      END DO

      DO n=1, NBS
        PRD(n,k) = PRD(n,k) + UMT(n,k)*O(1,1)
        PRD(n,k) = PRD(n,k) + UMT(n,l)*O(2,1)
        PRD(n,l) = PRD(n,l) + UMT(n,k)*O(1,2)
        PRD(n,l) = PRD(n,l) + UMT(n,l)*O(2,2)
      END DO

      UMT = PRD
      DO i=1, NBS
        DO k=1, NBS
          SPC(k,i) = 0.0D0
          DO l=1, NBS
            SPC(k,i) = SPC(k,i) + UMT(l,i)*HAM(l,k) !l,k better than k,l base on memory allocation
          END DO
        END DO
      END DO

      PRD = 0.0D0

      !Make new HAM
      DO i=1, NBS
        DO j=1, NBS
          DO k=1, NBS
            PRD(j,i) = PRD(j,i) + UMT(k,j)*SPC(k,i)
          END DO
        END DO
      END DO

      ERRNW = 0.0D0
      DO i=1,NBS
        DO j=i+1, NBS
          ERRNW = ERRNW + PRD(i,j)*PRD(j,i)
        END DO
      END DO

      WRITE(*,30) iTry, ERRNW, ERRPREV
      30 FORMAT(I3, 3G15.6)

      IF (ERRNW < 1.0D-12) EXIT
    END DO

    IF (iTry == MXIT) THEN
      WRITE(*,*) 'Warning: No Convergence'
    END IF

    WRITE(*,*) iTry, NBS, FLOAT(iTry)/(NBS*NBS), 'Diag Eff'

    HAM = PRD
  END SUBROUTINE DIAGNxN

  RECURSIVE SUBROUTINE sort(idx, PRD, high)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: idx(:,:)
    REAL(DBL), INTENT(IN) :: PRD(:,:)
    INTEGER, INTENT(IN) :: high
    INTEGER :: low, mid

    low = 1
    IF (low < high) THEN
      mid = low + (high-low)/2
      CALL sort(idx(:,low:mid), PRD, mid)
      CALL sort(idx(:,mid+1:high), PRD, high-mid)
      idx(:,low:high) = Merge(idx(:,low:mid), idx(:,mid+1:high), PRD, mid, high-mid)
    END IF
  END SUBROUTINE sort

  FUNCTION Merge(a, b, PRD, a_high, b_high)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: a, b
    REAL(DBL), INTENT(IN) :: PRD(:,:)
    INTEGER, INTENT(IN) :: a_high, b_high

    INTEGER :: Merge(2,a_high+b_high)
    INTEGER :: a_ptr, b_ptr, c_ptr

    a_ptr = 1
    b_ptr = 1
    c_ptr = 1

    DO WHILE (a_ptr <= a_high .AND. b_ptr <= b_high)
      IF (ABSO(PRD(a(1,a_ptr),a(2,a_ptr))) > ABSO(PRD(b(1,b_ptr),b(2,b_ptr)))) THEN
        Merge(:,c_ptr) = a(:,a_ptr)
        a_ptr = a_ptr + 1
      ELSE
        Merge(:,c_ptr) = b(:,b_ptr)
        b_ptr = b_ptr + 1
      END IF
      c_ptr = c_ptr + 1
    END DO

    IF (a_ptr > a_high) THEN
      Merge(:,c_ptr:) = b(:,b_ptr:b_high)
    ELSE
      Merge(:,c_ptr:) = a(:,a_ptr:a_high)
    END IF
  END FUNCTION Merge

  SUBROUTINE Class_DIAGNxN(NDH, NBS, HAM, OVR, UMT, PRD, SPC)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDH, NBS
    REAL(DBL), INTENT(INOUT) :: HAM(:,:), OVR(:,:) !OVR is not used yet
    REAL(DBL), INTENT(OUT) :: UMT(:,:), PRD(:,:), SPC(:,:)

    REAL(DBL) :: H(2,2), O(2,2), E(2)!, V(2,2), T(2,2), D(2,2)  !J2X2
    REAL(DBL) :: OFF_MAX, ERRPREV, ERRNW, FOM(1000)    !fom = Figure of marite

    INTEGER :: i, j, k, l, m, n, iTry
    INTEGER :: MXIT, nPair, idx(2,1000)
    LOGICAL :: update(1000), new
    INTEGER :: mPair, no

    IF (NDH > 1000) THEN
      WRITE(*,*) 'NDH must be smaller then 1000'
      STOP
    END IF

    IF (NBS > NDH) THEN
      WRITE(*,*) 'NDH must be larger than ', NBS
      STOP
    END IF

    UMT = 0.0D0
    DO i=1, NBS
      UMT(i,i) = 1.0D0
    END DO

    PRD = HAM
    MXIT = NBS*NBS*2

    DO iTry = 1, MXIT
      update = .TRUE.
      OFF_MAX = 0.0D0
      ERRPREV = 0.0D0
      nPair = 0

      DO i=1, NBS
        IF (update(i)) THEN
          new = .FALSE.
          DO j=1, NBS
!            IF (update(j)) THEN
            IF (j /= i) THEN
              ERRPREV = ERRPREV + PRD(i,j)*PRD(j,i)
              IF (ABSO(PRD(i,j)) >= OFF_MAX) THEN
                k = i
                l = j
                OFF_MAX = ABSO(PRD(i,j))
                new = .TRUE.
              END IF
            END IF
          END DO

          IF (new) THEN
            mPair = 0
            DO j=1, nPair
              no = 0
              IF (idx(1,j) == k) no=1
              IF (idx(2,j) == k) no=1
              IF (idx(1,j) == l) no=1
              IF (idx(2,j) == l) no=1
              IF (no == 0) THEN
                mPair = mPair + 1
                idx(1,mPair) = idx(1,j)
                idx(2,mPair) = idx(2,j)
                FOM(mPair) = FOM(j)
              ELSE
                update(idx(1,j)) = .TRUE.
                update(idx(2,j)) = .TRUE.
              END IF
            END DO

            nPair = mPair + 1
            DO j=nPair, 2, -1
              idx(1,j) = idx(1,j-1)
              idx(2,j) = idx(2,j-1)
              FOM(j) = FOM(j-1)
            END DO

            update(k) = .FALSE.
            update(l) = .FALSE.

            idx(1,1) = MIN(k,l)
            idx(2,1) = MAX(k,l)

            FOM(1) = OFF_MAX
          END IF
        END IF
      END DO

      DO i=1, nPair
        WRITE(*,*) idx(1,i), idx(2,i), PRD(idx(1,i),idx(2,i))
      END DO

      k = idx(1,1)
      l = idx(2,1)

      H(1,1) = PRD(k,k)
      H(1,2) = PRD(k,l)
      H(2,1) = PRD(l,k)
      H(2,2) = PRD(l,l)
      CALL J2X2(H, E, O)

!      WRITE(*,*) 'E and O values:'
!      DO i=1, 2
!        WRITE(*,*) E(i), (O(i,j),j=1,2)
!      END DO

      SPC = 0.0D0
      DO i=1, NBS
        SPC(i,i) = 1.0D0
      END DO

      SPC(k,k) = O(1,1)
      SPC(l,k) = O(2,1)
      SPC(l,l) = O(2,2)
      SPC(k,l) = O(1,2)

      !Get new unitary matrix
      PRD = 0.0D0
      DO m=1, NBS
        IF (m /= k .AND. m /= l) THEN
          DO n=1, NBS
            PRD(n,m) = UMT(n,m)
          END DO
        END IF
      END DO

      DO n=1, NBS
        PRD(n,k) = PRD(n,k) + UMT(n,k)*O(1,1)
        PRD(n,k) = PRD(n,k) + UMT(n,l)*O(2,1)
        PRD(n,l) = PRD(n,l) + UMT(n,k)*O(1,2)
        PRD(n,l) = PRD(n,l) + UMT(n,l)*O(2,2)
      END DO

      UMT = PRD
      DO i=1, NBS
        DO k=1, NBS
          SPC(k,i) = 0.0D0
          DO l=1, NBS
            SPC(k,i) = SPC(k,i) + UMT(l,i)*HAM(l,k) !l,k better than k,l base on memory allocation
          END DO
        END DO
      END DO

      PRD = 0.0D0
      !Make new HAM
      DO i=1, NBS
        DO j=1, NBS
          DO k=1, NBS
            PRD(j,i) = PRD(j,i) + UMT(k,j)*SPC(k,i)
          END DO
        END DO
      END DO

      ERRNW = 0.0D0
      DO i=1,NBS
        DO j=i+1, NBS
          ERRNW = ERRNW + PRD(i,j)*PRD(j,i)
        END DO
      END DO

      WRITE(*,30) iTry, ERRNW, ERRPREV, ERRPREV - OFF_MAX*OFF_MAX
      30 FORMAT(I3, 3G15.6)

      IF (ERRNW < 1.0D-12) EXIT
    END DO

    IF (iTry == MXIT) THEN
      WRITE(*,*) 'Warning: No Convergence'
    END IF

    WRITE(*,*) iTry, NBS, float(iTry)/float(NBS*NBS), 'Diag Eff'

    HAM = PRD
  END SUBROUTINE Class_DIAGNxN
END MODULE ArrayFunctions
