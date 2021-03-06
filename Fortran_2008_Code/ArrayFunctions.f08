MODULE ArrayFunctions
  USE FortranFunctions, ONLY: DBL, ABSO, SQR, DIV
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: INVERSE, J2x2, JAC2BY2GEN, DIAGNxN, LEASTSQUARE

  CONTAINS

  !Author (Fortran 77): Dr. Mark Pederson
  !Date:  September 8th, 2021
  !Modifier: Robert Alvarez / July 9th, 2022
  ! Modifications:
  !   Convert Subroutine to function
  !   Generalization of 4x4 matrix to nxn matrix
  !   Improve matrix printing format
  FUNCTION INVERSE(AA) RESULT(B)
    IMPLICIT NONE

    REAL(DBL), INTENT(IN) :: AA(:,:)

    REAL(DBL), DIMENSION(SIZE(AA,1),SIZE(AA,1)) :: B, CC
    REAL(DBL), DIMENSION(SIZE(AA,1),SIZE(AA,1)*2) :: A
    REAL(DBL) :: TMAX, temp
    CHARACTER(len=13) :: fmt_mt
    INTEGER :: i, j, k, n

    n = SIZE(AA,1)
    A = 0.D0

    !Create matrix A = [AA|I]
    DO i=1, n
      A(i,i+n) = 1.D0       !Create identity matrix
      A(i,1:n) = AA(i,1:n)  !Copy values of AA into first nxn spaces of A
    END DO

    !Print matrix A
    WRITE(*,*)
    WRITE(fmt_mt, '( "(",I2,"ES17.8E3))" )' ) 2*n
    WRITE(*,fmt_mt) (A(i,1:2*n), i=1,n)

    !Find invert matrix using Gaussian elimination
    DO i=1, n
      !Find largest value of column i
      TMAX = ABSO(A(i,i))
      k = i
      DO j=i+1, n
        IF (ABSO(A(j,i)) > TMAX) THEN
          k = j
          TMAX = ABSO(A(j,i))
        END IF
      END DO

      TMAX = A(k,i)
      IF (TMAX < 1.D-15) STOP 'A row is linearly dependent of one or more other rows'

      !Swap row with highest value in column i
      IF (k /= i) THEN
        DO j=i, 2*n
          temp = A(i,j)
          A(i,j) = A(k,j)
          A(k,j) = temp
        END DO
      END IF

      !Normalize matrix
      DO j=i,2*n
        A(i,j) = DIV(A(i,j),TMAX)
      END DO

      !Subtract value A(j,i) to every column at row i
      DO j=1, n
        temp = A(j,i)
        IF (j /= i) THEN
          A(j,i:2*n) = A(j,i:2*n) - temp*A(i,i:2*n)
        END IF
      END DO
    END DO

    !Copy inverse matrix into B
    B(1:n,1:n) = A(1:n,n+1:2*n)

    !Print inverse matrix
    WRITE(*,*)
    WRITE(fmt_mt, '( "(",I2,"ES17.8E3))" )' ) n
    WRITE(*,fmt_mt) (B(i,1:n), i=1, n)

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
    WRITE(*,*)
    WRITE(*,fmt_mt) (CC(i,1:4), i=1, n)
  END FUNCTION INVERSE

  !Author (Fortran 77): Dr. Mark Pederson
  !Date:  September 15th, 2021
  SUBROUTINE J2X2(H, E, O)
    IMPLICIT NONE
    REAL(DBL), INTENT(IN) :: H(:,:)
    REAL(DBL), INTENT(OUT) :: E(:), O(:,:) !O is for orthogonal matrix
    REAL(DBL) :: dot, trc, rad
    INTEGER :: i

!   if A = [[d, e], [f,g]]
!   det(A) = d*g - f*e
!   An eigenvalue satisfies:
!   (d-E)*(g-E) - e*f = 0 <=> E^2 - E*(d+g) + (d*g-f*e) = 0
!   From quadratic formula: b^2-4*a*c ->
!   (d+g)^2 - 4*1*(d*g-f*e) <=> (d-g)^2 + 4*e*f
!   E = answers of quadratic formula

    rad = SQR((H(1,1) - H(2,2))*(H(1,1) - H(2,2)) + 4.D0*H(1,2)*H(2,1))
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
    IF (A < 1.D-15) STOP 'Non positive overlap matrix'

    B = -( H(1,1)*O(2,2) - O(1,2)*H(2,1) + O(1,1)*H(2,2) - H(1,2)*O(2,1) )
    C = H(1,1)*H(2,2) - H(1,2)*H(2,1)

    TRC = DIV(-B,A)
    RAD = DIV(SQR(B*B - 4.D0*A*C),A)

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

    ! <V_1 | V_2> = 0 ?
    DO iTry=1, 3
      IF (iTry <= 2) THEN
        WRITE(*,'(/,A)') 'Overlap matrix'
      ELSE
        WRITE(*,'(/,A)') 'Hamiltonian matrix'
      END IF

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
!
!Subroutine sort and MergeIdx are only used for DIAGNxN
!
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
      idx(:,low:high) = MergeIdx(idx(:,low:mid), idx(:,mid+1:high), PRD, mid, high-mid)
    END IF
  END SUBROUTINE sort

  FUNCTION MergeIdx(a, b, PRD, a_high, b_high)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(IN) :: a, b
    REAL(DBL), INTENT(IN) :: PRD(:,:)
    INTEGER, INTENT(IN) :: a_high, b_high

    INTEGER :: MergeIdx(2,a_high+b_high), ai, bi, ci

    ai = 1
    bi = 1
    ci = 1

    DO WHILE (ai <= a_high .AND. bi <= b_high)
      IF (ABSO(PRD(a(1,ai),a(2,ai))) > ABSO(PRD(b(1,bi),b(2,bi)))) THEN
        MergeIdx(:,ci) = a(:,ai)
        ai = ai + 1
      ELSE
        MergeIdx(:,ci) = b(:,bi)
        bi = bi + 1
      END IF
      ci = ci + 1
    END DO

    IF (ai > a_high) THEN
      MergeIdx(:,ci:) = b(:,bi:b_high)
    ELSE
      MergeIdx(:,ci:) = a(:,ai:a_high)
    END IF
  END FUNCTION MergeIdx

  SUBROUTINE DIAGNxN(NDH, NBS, HAM, UMT, PRD)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDH, NBS
    REAL(DBL), INTENT(INOUT) :: HAM(:,:)
    REAL(DBL), INTENT(OUT) :: UMT(:,:), PRD(:,:)

    REAL(DBL) :: SPC(NDH,NDH)
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

    UMT = 0.D0
    DO i=1, NBS
      UMT(i,i) = 1.D0
    END DO

    PRD = HAM
    MXIT = NBS*NBS*2

    DO iTry = 1, MXIT
      n = 0
      ERRPREV = 0.D0
      DO i=1, NDH
        DO j=i+1, NDH
          ERRPREV = ERRPREV + PRD(i,j)*PRD(j,i)
          IF (ABSO(PRD(i,j)) > 1.D-10) THEN ! Save all indices that have values in PRD that are greater than 0
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

!      WRITE(*,*) 'All indexes:'
!      DO i=1, n
!        WRITE(*,*) idxAll(1,i), idxAll(2,i), PRD(idxAll(1,i),idxAll(2,i))
!      END DO

      !Choose the indexes with the highest values at PRD without repeating indexes
      useIdx = .TRUE. !No indexes has been use
      idxSize = 1     !Keep track of position to add new non-repetitive idxs
      DO j=1, n
        IF (useIdx(idxAll(1,j)) .AND. useIdx(idxAll(2,j))) THEN !If none of this two indexes has been used
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

      SPC = 0.D0
      DO i=1, NBS
        SPC(i,i) = 1.D0
      END DO

      SPC(k,k) = O(1,1)
      SPC(l,k) = O(2,1)
      SPC(l,l) = O(2,2)
      SPC(k,l) = O(1,2)

      !Get new unitary matrix
      PRD = 0.D0
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
          SPC(k,i) = 0.D0
          DO l=1, NBS
            SPC(k,i) = SPC(k,i) + UMT(l,i)*HAM(l,k) !l,k better than k,l base on memory allocation
          END DO
        END DO
      END DO

      PRD = 0.D0

      !Make new HAM
      DO i=1, NBS
        DO j=1, NBS
          DO k=1, NBS
            PRD(j,i) = PRD(j,i) + UMT(k,j)*SPC(k,i)
          END DO
        END DO
      END DO

      ERRNW = 0.D0
      DO i=1,NBS
        DO j=i+1, NBS
          ERRNW = ERRNW + PRD(i,j)*PRD(j,i)
        END DO
      END DO

      WRITE(*,30) iTry, ERRNW, ERRPREV
      30 FORMAT(I3, 3G15.6)

      IF (ERRNW < 1.D-12) EXIT
    END DO

    IF (iTry == MXIT) THEN
      WRITE(*,*) 'Warning: No Convergence'
    END IF

    WRITE(*,*) iTry, NBS, REAL(iTry)/REAL(NBS*NBS), 'Diag Eff'

    HAM = PRD
  END SUBROUTINE DIAGNxN

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

!    CALL DIAGNxN(DM)

    STOP 'Before INVERSE'

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
END MODULE ArrayFunctions
