      PROGRAM MAIN
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AA(4,4),BB(4,4)
      DIMENSION H(2,2), O(2,2), V(2,2), E(2)
      DIMENSION F(4), G(4)
      
      DO
       PRINT'(/,A)','DIAG options:'
       PRINT'(A)','1. DIAG','2. STARK','3. PIR','4. HMO','5. BOX'
       PRINT'(A)', 'n. exit the program'
       READ*,NOPT

       IF(NOPT.EQ.1) THEN
        CALL DIAGDVR
       ELSE IF(NOPT.EQ.2) THEN
        CALL STARKDVR
       ELSE IF(NOPT.EQ.3) THEN
        CALL RINGDVR
       ELSE IF (NOPT .EQ. 4) THEN
        CALL HMODVR
       ELSE IF (NOPT .EQ. 5) THEN
         CALL BOXDVR
       ELSE
        EXIT
       END IF
      END DO
c     PRINT*,'H=?'
C     
C      PRINT*,'O=?'
C      READ*, O(1,1), O(1,2), O(2,2)
C      O(2,1)=O(1,2)
C
C      CALL JAC2BY2(H, O, V, E)
C      WRITE(*,'(/,2F12.6)') E(1), E(2)
      
C      DO I=1,2
C       PRINT 100,E(I),(O(J,I),J=1,2)
C      END DO
C     PRINT*, O(1,1)*O(1,2)+O(2,1)*O(2,2)
C     PRINT*, O(1,1)*O(1,1)+O(2,1)*O(2,1)
C     PRINT*, O(2,2)*O(2,2)+O(1,2)*O(1,2)
C100  FORMAT(3F12.6)

C     CALL INVERSE4X4(AA,BB)
C     
C     DO L=1,4
C      PRINT*,"F=?"
C      READ*,(F(I),I=1,4)
C      DO I=1,4
C       G(I)=0.0D0
C       DO J=1,4
C        G(I)=G(I)+BB(I,J)*F(J)
C       END DO
C       PRINT*,G(I)
C      END DO
C     END DO
      END PROGRAM

      SUBROUTINE JAC2BY2(H, O, V, E)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(2,2), O(2,2), V(2,2), T(2,2), E(2)
      DIMENSION D(2,2)  !Dot product matrix

      A = O(1,1)*O(2,2) - O(1,2)*O(2,1)
      IF (A .LT. 0.0D0) THEN
       WRITE(*,*) 'Non positive overlap matrix'
       STOP
      END IF
      B = H(1,1)*O(2,2) - O(1,2)*H(2,1) 
     &  + O(1,1)*H(2,2) - H(1,2)*O(2,1)
      B = -B
      C = H(1,1)*H(2,2) - H(1,2)*H(2,1)

      TRC = DIV(-B,A)
      RAD = DIV(SQR(B*B - 4.0D0*A*C),A)

      E(1) = DIV((TRC + RAD),2.0D0)
      E(2) = DIV((TRC - RAD),2.0D0)

C     WRITE(*,'(/,A,2F12.6,/)') 'Eigenvalues:', E(1), E(2)

      !Calculate eigenvectors
      DO k=1, 2
       DO i=1, 2
        DO j=1, 2
         T(i,j) = H(i,j) - E(k)*O(i,j)
        END DO
       END DO
       V(1,k) = -T(k,2)
       V(2,k) =  T(k,1)
C      WRITE(*,*) T(1,1)*T(2,2) - T(1,2)*T(2,1)
      END DO

      !C <V_1 | V_2> = 0 ?
      DO iTry=1, 3
       IF (iTry .LE. 2) THEN
C       WRITE(*,'(/,A)') 'Overlap matrix'
       ELSE
C       WRITE(*,'(/,A)') 'Hamiltonian matrix'
       END IF
       DO i=1, 2
        DO j=1, 2
         D(i,j) = 0.0D0
         DO k=1, 2
          DO l=1, 2
           IF (iTry .LE. 2) THEN
            D(i,j) = D(i,j) + V(k,i)*V(l,j)*O(k,l)
           ELSE
            D(i,j) = D(i,j) + V(k,i)*V(l,j)*H(k,l)
           END IF
          END DO
         END DO
        END DO
C       WRITE(*,'(2F12.6)') (D(i,j),j=1,2)
       END DO
       IF (iTry .EQ. 1) THEN
        DO i=1, 2
         DO k=1, 2
          V(k,i) = DIV(V(k,i),SQR(D(i,i)))
         END DO
        END DO
       END IF
      END DO
      END SUBROUTINE JAC2BY2

      SUBROUTINE LEASTSQUARE
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(100),DM(4,4),DI(4,4),B(4),A(4)
      DIMENSION P(100, 4)
      PRINT*,"XMIN, XMAX=?"
      READ*, XMIN, XMAX
      DX=DIV(XMAX-XMIN,99.0D0) 
      X=XMIN-DX
      DO I=1,100
      X=X+DX
      F(I)=X*X*X*X
      POW=1.0D0
       DO J=1,4
       P(I,J)=POW
       POW=POW*X
       END DO
      END DO
C     EVALUATE DESIGN MATRIX
      DO J=1,4
      DO K=1,4
      DM(K,J)=0.0D0
        DO I=1,100
        DM(K,J)=DM(K,J)+P(I,K)*P(I,J)
        END DO
      END DO
        B(J)=0
        DO I=1,100
         B(J)=B(J)+F(I)*P(I,J)
        END DO
      END DO
      
      CALL INVERSE4X4(DM,DI)
      
      DO J=1,4
       A(J)=0.0D0
       DO K=1,4
        A(J)=A(J)+DI(J,K)*B(K)
       END DO
      END DO
      PRINT*,"A=",A
      
      DO I=1,10
       PRINT*,"X=?"
       READ*,X
       FIT=0.0D0
       POW=1.0D0
       DO J=1,4
        FIT=FIT+A(J)*POW
        POW=POW*X
       END DO
       DAT=X*X*X*X
       PRINT*,X,DAT,FIT
      END DO
      STOP
      END SUBROUTINE LEASTSQUARE

      SUBROUTINE INVERSE4x4(AA,BB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(4),A(4,8)
      DIMENSION AA(4,4),BB(4,4),CC(4,4)
C     A contains an original matrix and a Unitary matrix at the right.
      DX=1.0D0
      X(1)=DX
      DO I=2,4
      X(I)=X(I-1)+DX
      END DO

      A=0.0D0
      PRINT*,X
      PRINT*, 'Original Matrix:'
      DO I=1,4
      A(I,I+4)=1.0D0 !Unitary matrix in the right part of A
      A(I,1)=1.0D0   !1.0 on the first column of A
      DO J=1,4
      A(I,J)=AA(I,J)
      END DO
      PRINT 10,(A(I,J),J=1,8)
      END DO
  10  FORMAT (8F10.4)

C     Saving the original left matrix to AA
      DO I=1,4
      DO J=1,4
      AA(I,J)=A(I,J)
      END DO
      END DO

C     We need to find the greater value of each column and put that element's row on top
C     Then go to the next column and do the same but starting on the next row, ignoring 
c     The line we already moved
      DO I=1,4
       TMAX=ABSO(A(I,I))
       K=I
       DO J=I+1,4
        IF(ABSO(A(J,I)).GT.TMAX)THEN
         K=J
         TMAX=ABSO(A(J,I))
        END IF
       END DO

C     We have to check for TMAX to not be equal to zero, if that happens it means
C     one row is linearly dependent of other lines.
       TMAX=A(K,I)
       IF(TMAX.EQ.0.0D0) THEN
       PRINT*, 'A row is linearly dependant of one or more other rows'
       STOP
       END IF

       IF(K.NE.I) THEN
        DO J=I,8
         SAV=A(I,J)
         A(I,J)=A(K,J)
         A(K,J)=SAV
        END DO
       END IF
       DO K=I,8 
        A(I,K)=DIV(A(I,K),TMAX)
       END DO
       DO J=1,4
        SAV=A(J,I)
        IF(J.NE.I)THEN
         DO K=I,8
          A (J,K)=A(J,K)-SAV*A(I,K)
         END DO
        END IF
       END DO
       PRINT*,K,TMAX
       DO L=1,4
        PRINT 10,(A(L,J),J=1,8)
       END DO
      END DO
      
C     Save the inversed matrix as BB      
      DO I=1,4
       DO J=1,4
        BB(I,J)=A(I,J+4)
       END DO
      END DO
C Lets check there's no row equal to zero, if theres a row equal to zero
      
C     If we did this right, C should come out as the unitary matrix.
      CC=0.D0      
      PRINT*,"IS THIS ONE?"
      DO I=1,4
      DO J=1,4
      DO K=1,4
      CC(I,J)=CC(I,J)+AA(I,K)*BB(K,J)
      END DO
      END DO
      PRINT 10,(CC(I,J),J=1,4)
      END DO
      RETURN
      END SUBROUTINE INVERSE4x4
     
      SUBROUTINE J2X2(H,E,O)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(2,2),E(2),O(2,2)
C     DET(H)=H(1,1)*H(2,2)-H(1,2)*H(2,1)
C     An eigenvalue satisfies; 
C     (H(1,1)-E)*(H(2,2)-E)-H(1,2)*H(2,1)=0
C     E^2 -E*(H(1,1)+H(2,2))+(H(1,1)*H(2,2)-H(1,2)*H(2,1))=0
C     E = (1/2)[H(1,1)+H(2,2) +- SQR[(H(1,1)-H(2,2))^2 + H(1,2)*H(2,1)]]
      RAD=SQR((H(1,1)-H(2,2))*(H(1,1)-H(2,2))+4.0D0*H(1,2)*H(2,1))
      TRC=     H(1,1)+H(2,2) !Trace=Sum of the diagonal of a square matrix.
      E(1)=(0.5D0)*(TRC+RAD)
      E(2)=(0.5D0)*(TRC-RAD) !Its easy to see that E(1)+E(2)=TRC
C     PRINT*,RAD,TRC,E,"HI"
C     PRINT*,H
C     The eigenvector O_n follows that (H-IE(n))O_n=0.
C     so: (H(1,1)-E(n))O_n(1)+H(1,2)O_n(2)=0
C         H(2,1)O_n(1)+(H(2,2)-E(n))O_n(2)=0
C     But: H(2,1)=H(1,2)
C     We can do: O_n(1)=-H(1,2) ; O_n(2)=(H(1,1)-E(n))
C     This way the first equation is zero, and the second will have the form:
C     -H(1,2)^2 + H(1,1)H(2,2)-E(n)(H(1,1)+H(2,2)) + E(n)^2 = 0
C     Notice that this equation is equal to the determinant, which is equal to zero.
C     So our solutions are valid, and any non zero scalar multiples of this O_n vectors.
C     Notice that we can do the same using the second equation to obtain O_n(1) and O_n(2)
C     And we will have that:
C     O_n(1)=(H(2,2)-E(n)) and O_n(2)=-H(1,2), so
C     Because O_n are vertical vectors we can form a matrix O as:
C     O=(O_1(1),O_2(1))
C       (O_1(2),O_2(2))
C     And use this last two equations to get O_n(1) and O_n(2) for n=2
C     And the  previous two equations to get O_n(1) and O_n(2) for n=1
C     This way: (O(BASIS, EIGENVECTOR INDEX))
      O(1,1)=-H(1,2)
      O(2,1)= H(1,1)-E(1)
      O(1,2)= H(2,2)-E(2)
      O(2,2)=-H(1,2) 
C     This vectors are not normalized, so lets normalize them:
      DO I=1,2
       DOT=O(1,I)*O(1,I)+O(2,I)*O(2,I)  
       DOT=SQR(DOT)                     
       O(1,I)=DIV(O(1,I),DOT)       
       O(2,I)=DIV(O(2,I),DOT)
      ENDDO 
      RETURN
      END SUBROUTINE J2x2
      
      FUNCTION ABSO(B)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(B.LT.0.0D0) THEN  
      ABSO=-B
      ELSE
      ABSO=B
      ENDIF
      RETURN
      END FUNCTION ABSO
      
      FUNCTION SQR(B) !This function finds the square root of B
      IMPLICIT REAL*8 (A-H,O-Z)
      SQR=1.1D0 ! SQR needs a first value to start the iterations
      SQA=-0.1D0 ! SQA, SQB & SQC are dummy variables to keep track
      SQB=-0.2D0 ! of how SQR evolves through the iterations.
      SQC=-0.3D0
C      S =~ (B)^0.5
C  E + S =  (B)^0.5 
C (E + S)^2 = B
C  E^2 + 2SE + S^2 = B, But E is very small, so E^2=0:
C        E = DIV(B - S^2,2S), this way:
C  S =~ B^0.5 = S + E -> S =~ S + E, this updates S to a new value.
C So by repeating this multiple times we will get B^0.5

      DO I=1,100
       EPS=DIV(B-SQR*SQR,2.0D0*SQR) !Get the value of E according to the current S
       SQR=SQR+EPS !Update S to S+E
       SQA=SQB !Update the values of our dummys
       SQB=SQC
       SQC=SQR
       IF(SQC.EQ.SQB)THEN !If the last two values are the same then the code stops.
        RETURN
       ELSE IF(SQC.EQ.SQA)THEN !If the Ith and Ith-2 value are the same, then the
        SQR=DIV(SQA+SQB,2.0D0)!code is cycling, so we stop it and average between this
        RETURN                !values
       END IF
C      PRINT*,I,SQR
      END DO
C      STOP
      END FUNCTION SQR

      FUNCTION DIV(AA,BB) !This function returns A/B
      IMPLICIT REAL*8 (A-H,O-Z)
      A = AA
      B = BB
      DIV = AA*RECIPROCAL(BB) 
      RETURN
      END FUNCTION DIV

      FUNCTION RECIPROCAL(ZZ)
      IMPLICIT REAL*8(A-H,O-Z)
      Z=ZZ
C      DO ITEST=1,30
C      Z=Z/2.
C stop if 1/0
      IF(Z.EQ.0.0D0)STOP '1/0 = Inf'

C ALLOW NEGATIVE NUMBERS
      FAC=1.0D0
      IF(Z.LT.0.0D0)THEN
       FAC=-1.0D0
       Z=ABSO(Z)
      END IF

C FROM CLASS
      IF(Z.GT.1.0D0)THEN
       X=Z-1.0D0
       CALL GT1DIVIDE(X,REC)
C IF 1/1, return 1
      ELSE IF(Z.EQ.1.0D0)THEN
       REC=1.0D0
C IF 0<Z<1, use new subroutine
      ELSE IF(Z.GT.0.1D0.AND.Z.LT.1.0D0)THEN
       X=1.0D0-Z
       CALL LT1DIVIDE(X,REC)
      ELSE IF(Z.GT.0.AND. Z.LE.0.1)THEN
      CALL TNYDIVIDE(Z,REC)
      END IF
      REC=FAC*REC
C      END DO
      END FUNCTION RECIPROCAL

C DIVIDER (1/(1+x))   0<x<1
      FUNCTION DIVIDER(X) 
      IMPLICIT REAL*8(A-H,O-Z)

      TWOTHIRDS = 0.666666666666667D0
C      PRINT*,"X = ?"
C      READ*,  X
      IF(X.LT.0.0D0) STOP'X LT 0'
      IF(X.GT.1.0D0) STOP'X GT 1'

      IF(X.GT.0.5D0) THEN
      Y = TWOTHIRDS*(X-0.5D0)
      REC = DIVIDE(Y)
      REC = REC*TWOTHIRDS
      ELSE 
      REC = DIVIDE(X)
      ENDIF
      DIVIDER = REC
      END FUNCTION DIVIDER

C DIVIDE (1/(1+x))   0<x<0.5
      FUNCTION DIVIDE(X)
      IMPLICIT REAL*8(A-H,O-Z)
      N=50
      IF(X.LT.0.0D0) STOP'X LT 0   DIVIDE'
      IF(X.GT.0.5D0) STOP'X GT 0.5 DIVIDE'
       REC = 0.0D0
         P = 1.0D0

       DO i=1,N
        REC = REC + P
          P = -P*X
       ENDDO

      DIVIDE = REC
      END FUNCTION DIVIDE

      SUBROUTINE GT1DIVIDE(X,REC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ONED2N(0:64),ARE(0:64),TWON(0:64) 
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      
      IF(FIRST) THEN
      ONED2N(0)=1.0D0
      DO I=1,64
       ONED2N(I)=ONED2N(I-1)*0.5D0
c       PRINT*,'2**',I,ONED2N(I)
      ENDDO
      
      ARE(0)=0.5D0
      DO I=1,64
       ARE(I)=ONED2N(I)*DIVIDE(ONED2N(I)) 
C      PRINT*,I,ARE(I)
      ENDDO
      
      TWON(0)=1.0D0
      DO I=1,64
       TWON(I)=TWON(I-1)*2
C      PRINT*, I, TWON(I)
      END DO
      END IF
      FIRST=.TRUE. !Temporary Change
      
      IF(X.LT.1.0D0) THEN
      REC = DIVIDER(X)
      ELSE
      DO I=0,63 
       IF(X.GE.TWON(I).AND.X.LT.TWON(I+1)) THEN
       Y = (X-TWON(I))*ARE(I)
C       PRINT'(I3,6F12.6)',I,X,TWON(I),X-TWON(I),ARE(I-1)
C     $,1.0D0/(1.0D0+Y)  
        REC = DIVIDER(Y)
C       PRINT'(3F12.6)', X,Y,REC
        REC = ARE(I)*REC 
       END IF
      END DO
      END IF
      END SUBROUTINE GT1DIVIDE

C new subroutine for 1/Z, where 0<Z<1
C 1/(1-x)=(1+x+x**2+....x**N)
      SUBROUTINE LT1DIVIDE(X,REC)
      IMPLICIT REAL*8(A-H,O-Z)
      N=1000
      IF(X.LT.0.0D0) STOP'X LT 0 LT1DIVIDE'
      IF(X.GT.1.0D0) STOP'X GT 1 LT1DIVIDE'
       REC = 0.0D0
         P = 1.0D0

       DO i=1,N
        REC = REC + P
        P=P*X
       END DO
       RETURN 
      END SUBROUTINE LT1DIVIDE 

      SUBROUTINE TNYDIVIDE(Z,REC)
       IMPLICIT REAL*8 (A-H,O-Z)
C 1/Z = 1/(0.5**N + (Z-0.5**N)
C     = (1./0.5**N) 1/(1+ (Z-0.5**N)/(0.5**N))
C     = 2**N [1/(1 + (2^N(Z-0.5^N)))]
       PP=1.0D0
       PR=1.0D0
       DO I=0,64
        IF(Z.GE.PR.and.Z.LT.PR*2.0D0)THEN
         X=(Z-PR)*PP
C         PRINT*,'X:',X
         REC=PP*DIVIDER(X)
        END IF
       PR=PR*0.5D0
       PP=PP*2.0D0
       END DO 
      RETURN
      END SUBROUTINE TNYDIVIDE
      
      SUBROUTINE DIAGDVR
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NDH=12)
      DIMENSION HAM(NDH,NDH),OVR(NDH,NDH),UMT(NDH,NDH),SPC(NDH,NDH)
      DIMENSION PRD(NDH,NDH)

      G  = -1.0D0 
      X  = -0.5D0
      P  =  0.1D0
      TXR=  0.2D0
      HAM=  0.0D0
      OVR=  0.0D0
      DO I=1,NDH
       OVR(I,I)=1.0D0
      END DO
      HAM(1,1)=G
      HAM(2,2)=X
      HAM(3,3)=X
      HAM(4,4)=G
      HAM(5,5)=X
      HAM(6,6)=X
      HAM(7,7)=G
      HAM(1,2)=P
      HAM(2,3)=TXR
      HAM(3,4)=P
      HAM(4,5)=P
      HAM(5,6)=TXR
      HAM(6,7)=P
      NBS=7
C      HAM(1,7)=0.73
C      HAM(7,1)=0.73
      DO I=1,NBS
      DO J=I,NBS
       HAM(J,I)=HAM(I,J)
      END DO
       PRINT'(10F7.2)',(HAM(I,J),J=1,NBS)
      END DO
C      READ*, NDUMMY
C     HAM(2,1)=P
C     HAM(4,3)=P
C     HAM(2,3)=TXR
C     HAM(3,2)=TXR

      CALL DIAGnxn(NDH,NBS,HAM,OVR,UMT,PRD,SPC)
      PRINT*,"UPDATED HAM"
      DO I=1,NBS
       PRINT"(10F12.4)",(HAM(J,I),J=1,NBS)
      END DO

      STOP
      END SUBROUTINE DIAGDVR

      SUBROUTINE DIAGnxn(NDH,NBS,HAM,OVR,UMT,PRD,SPC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION HAM(NDH,NDH),OVR(NDH,NDH),UMT(NDH,NDH),PRD(NDH,NDH)
      DIMENSION H(2,2), O(2,2), V(2,2), T(2,2), E(2)
      DIMENSION D(2,2)  !Dot product matrix
      DIMENSION SPC(NDH,NDH)
      INTEGER :: PROCS(2,1000)
      LOGICAL   UPDATE(1000),FOUND
      DIMENSION INDEX(2,1000)
      DIMENSION FOM(1000)
      IF(NDH.GT.1000)THEN
       PRINT*, 'NDH MUST BE LESS THAN 1000'
       STOP
      END IF

      UPDATE=.TRUE.
C     Can we get rid of SPC?
      
      IF(NBS.GE.NDH) THEN 
      PRINT*, 'NDH MUST BE LARGER THAN', NBS
      STOP
      END IF
      UMT=0.0D0
      DO I=1,NBS
       UMT(I,I)=1.0D0
      END DO

      PRD=HAM

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     G=-1.0D0 
C     X=-0.5D0
C     P=0.1D0
C     TXR=0.2D0
C     HAM=0.0D0
C     OVR=0.0D0
C     DO I=1,4
C      OVR(I,I)=1.0D0
C     END DO
C     HAM(1,1)=G
C     HAM(4,4)=G
C     HAM(2,2)=X
C     HAM(3,3)=X
C     HAM(1,2)=P
C
C     HAM(2,1)=P
C     HAM(3,4)=P
C     HAM(4,3)=P
C     HAM(2,3)=TXR
C     HAM(3,2)=TXR
C      PRD=HAM
      
      MXIT=NBS*NBS*2
      DO 1000 ITRY=1,MXIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Alex 10-28-2021
      FOM=0.0D0
      INDEX=0
      KNT=0
      DO I=1,NBS
       DO J=I+1,NBS
C ignore zero elements
        IF(ABS(PRD(J,I)).GT.0.0D0)THEN
         KNT=KNT+1
         FOM(KNT)=ABS(PRD(J,I))
         INDEX(1,KNT)=I
         INDEX(2,KNT)=J
        END IF
       END DO
      END DO
C Sort FOM
C todo: use NlogN sort instead:
C add counter
      ICOUNT=0
      DO I=2,KNT
       DO J=1,I-1
        IF(FOM(I).GT.FOM(J))THEN
         ICOUNT=ICOUNT+1
         SAV=FOM(I)
         FOM(I)=FOM(J)
         FOM(J)=SAV
         DO M=1,2
          ISAV=INDEX(M,I)
          INDEX(M,I)=INDEX(M,J)
          INDEX(M,J)=ISAV
         END DO
        END IF
       END DO
      END DO
C      DO I=1,KNT
C       PRINT*,FOM(I),(INDEX(M,I),M=1,2)
C      END DO
C NOW WE LOOP TO GIVE PROCESSORS THEIR JOBS
C when processor picks the next pair, we can turn off already chosen indices with update:
       NPROC=2
       PROCS=0.0D0
       DO IPROC=1,NPROC
        FOUND=.FALSE.
        DO I=1,KNT 
         K=INDEX(1,I)
         J=INDEX(2,I)
         IF(UPDATE(K).AND.UPDATE(J).AND.FOM(I).GT.0.0D0)THEN
          UPDATE(K)=.FALSE.
          UPDATE(J)=.FALSE.
          PROCS(1,IPROC)=K
          PROCS(2,IPROC)=J
C          PRD(K,J)=0.0D0 ! JUST TO TEST IT GETS THE 0.1 next try
C          PRD(J,K)=0.0D0 ! JUST TO TEST IT GETS THE 0.1 next try
          GOTO 555
         END IF
        END DO
 555    CONTINUE
C send processor data:
       END DO
C just to avoid printing while not actually diagonalizing
C       IF(PROCS(1,1).EQ.0)THEN 
C         GOTO 1000
C       END IF
CCCCCCCCCCCCCC
C      PRINT*,'ITRY',ITRY
C       DO IPROC=1,NPROC
C       PRINT*,'PROC, INDICES:',IPROC,(PROCS(M,IPROC),M=1,2)
C       END DO
C       IF(PROCS(1,IPROC).EQ.0)DONT SEND WORK
       UPDATE=.TRUE.
C       GOTO 1000
C       STOP 'ALEX'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
c      DO I=1,4
c       PRINT"(4F12.4)",(PRD(I,J),J=1,4)
c      END DO
      ERROLD=0.0D0
      OFF_NRM=0.0D0
      DO I=1,NBS
      DO J=I+1,NBS
       ERROLD=ERROLD+PRD(I,J)*PRD(I,J)
      END DO
      END DO
      DO 544 IPROC=1,NPROC
      K=PROCS(1,IPROC)
      L=PROCS(2,IPROC)
C     PRINT*,'IPROC=',IPROC,K,L
C      DO I=1,NBS
C      PRINT '(7F12.6)',(PRD(J,I),J=1,NBS)
C      END DO
C      OFF_NRM=OFF_NRM+PRD(K,L)*PRD(K,L)
C      IF(OFF_NRM.LT.1.0D-20)THEN
C       PRINT*,'OFF_NRM IS ZERO',OFF_NRM,K,L
C      END IF
       
C      OFF_MAX=0.0D0 
C      ERROR=0.0D0
C      NPAIR=0
C      UPDATE=.TRUE.
C      KPAIR=0
C      DO 400 IOUT=1,1
C      UPDATE=.TRUE.
C      DO IPAIR=1,NPAIR
C       UPDATE(INDEX(1,IPAIR))=.FALSE.
C       UPDATE(INDEX(2,IPAIR))=.FALSE.
C      END DO
C       OFF_MAX=0.0D0
C      DO I=1,NBS
C       IF(UPDATE(I)) THEN
C        NEW=0
C        DO J=1,NBS!I+1,NBS
CC         IF(UPDATE(J).AND.J.NE.I)THEN
C         IF(J.NE.I) THEN
C          ERROR=ERROR+PRD(I,J)*PRD(J,I)
C          IF(ABS(PRD(I,J)).GE.OFF_MAX)THEN
C           NEW=NEW+1
C           K=I
C           L=J
C           OFF_MAX=ABS(PRD(I,J))
C          END IF
C         END IF
C        END DO
C        IF(NEW.NE.0) THEN
C         MPAIR=0
C         DO IPAIR=KPAIR+1,NPAIR
C
C                                 NO=0
C          IF(INDEX(1,IPAIR).EQ.K)NO=1
C          IF(INDEX(2,IPAIR).EQ.K)NO=1
C          IF(INDEX(1,IPAIR).EQ.L)NO=1
C          IF(INDEX(2,IPAIR).EQ.L)NO=1
C          IF(NO.EQ.0) THEN
C           MPAIR=MPAIR+1
C           INDEX(1,MPAIR)=INDEX(1,IPAIR)
C           INDEX(2,MPAIR)=INDEX(2,IPAIR)
C               FOM(MPAIR)=    FOM(IPAIR)
C          ELSE 
C           UPDATE(INDEX(1,IPAIR))=.TRUE.
C           UPDATE(INDEX(2,IPAIR))=.TRUE.
C          END IF
C
C         END DO
C         NPAIR=MPAIR+1
C         
C         DO IPAIR=NPAIR,2,-1
C          
C          INDEX(1,IPAIR)=INDEX(1,IPAIR-1)
C          INDEX(2,IPAIR)=INDEX(2,IPAIR-1)
C              FOM(IPAIR)=    FOM(IPAIR-1)
C         END DO
C         UPDATE(K)=.FALSE.
C         UPDATE(L)=.FALSE.
C         INDEX(1,1)=MIN(K,L)
C         INDEX(2,1)=MAX(K,L)
C             FOM(1)=OFF_MAX
CC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CC            UPDATE=.TRUE.
CC             NPAIR=1
CC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CC              DO IPAIR=      1,NPAIR
CC              DO JPAIR=IPAIR+1,NPAIR
CC               IF(FOM(JPAIR).GT.FOM(IPAIR)) THEN
CC               FSAVE=FOM(JPAIR)
CC               FOM(JPAIR)=FOM(IPAIR)
CC               FOM(IPAIR)=FSAVE
CC               DO J=1,2 
CC                ISAVE=INDEX(J,JPAIR)
CC                INDEX(J,JPAIR)=INDEX(J,IPAIR)
CC                INDEX(J,IPAIR)=ISAVE
CC               END DO
CC               END IF
CC              END DO
CC              END DO
C        END IF
C       END IF
C      END DO
C 400  CONTINUE
C      DO JPAIR=1,NPAIR
C       PRINT*, INDEX(1,JPAIR),INDEX(2,JPAIR),
C     $     PRD(INDEX(1,JPAIR),INDEX(2,JPAIR))
C      END DO

C     STOP
C      K=INDEX(1,1)
C      L=INDEX(2,1)
C     PRINT*, K,L,OFF_MAX
      H(1,1)=PRD(K,K)
      H(1,2)=PRD(K,L)
      H(2,1)=PRD(L,K)
      H(2,2)=PRD(L,L)
      OFF_NRM=OFF_NRM+H(1,2)*H(1,2)
      CALL J2X2(H,E,O)
      
C     PRINT*, "E AND O VALUES:"
C     DO I=1,2
c      PRINT*,E(I),(O(I,J),J=1,2)
C     END DO
      
      SPC=0.0D0
      DO I=1,NBS
       SPC(I,I)=1.0D0
      END DO
      
C store O for each proc?
      SPC(K,K)=O(1,1)
      SPC(L,K)=O(2,1)
      SPC(L,L)=O(2,2)
      SPC(K,L)=O(1,2)
CCC   Stop using PRD between here***************
      PRD=0.0D0
      DO M=1,NBS
       IF(M.NE.K.AND.M.NE.L) THEN
        DO N=1,NBS
         PRD(N,M)=UMT(N,M)
        END DO
       END IF
      END DO
C     
       DO N=1,NBS
        PRD(N,K)=PRD(N,K)+UMT(N,K)*O(1,1)
        PRD(N,K)=PRD(N,K)+UMT(N,L)*O(2,1)
        PRD(N,L)=PRD(N,L)+UMT(N,K)*O(1,2)
        PRD(N,L)=PRD(N,L)+UMT(N,L)*O(2,2)
       END DO

CCC   And here*****************
      UMT=PRD
      DO I=1,NBS
      DO K=1,NBS
       SPC(K,I)=0.0D0
      DO L=1,NBS
       SPC(K,I)=SPC(K,I)+UMT(L,I)*HAM(L,K)
      END DO
      END DO
      END DO

C     UMT=PRD
      PRD=0.0D0
      DO I=1,NBS
      DO J=1,NBS
      DO K=1,NBS
       PRD(J,I)=PRD(J,I)+UMT(K,J)*SPC(K,I)
      END DO
      END DO
      END DO
      
c      PRINT*,"UPDATED HAM"
      ERRNW=0.0D0
      DO I=1,NBS
       DO J=I+1,NBS
        ERRNW=ERRNW+PRD(I,J)*PRD(J,I)
       END DO
c       PRINT"(4F12.4)",(PRD(J,I),J=1,4)
      END DO

      PRINT 197,ITRY,ERRNW,ERROR,ERROR-OFF_MAX*OFF_MAX
 197  FORMAT(I3,3G15.6)
      IF(ERRNW.LT.1.0D-12) GO TO 1010
 544  CONTINUE
      ERRNEW=0.0D0
      DO I=1,NBS
      DO J=I+1,NBS
       ERRNEW=ERRNEW+PRD(I,J)*PRD(I,J)
      END DO
      END DO
      PRINT*,'tmp',ERROLD,ERRNEW,ERROLD-OFF_NRM
1000  CONTINUE
      PRINT*, 'WARNING: NO CONVERGENCE',MXIT
1010  CONTINUE
      PRINT*, ITRY,NBS,Float(ITRY)/(NBS*NBS), 'Diag Eff'
      HAM=PRD
      RETURN
      END SUBROUTINE DIAGnxn

      SUBROUTINE STARKDVR
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NDH=10)
      DIMENSION HAM(NDH,NDH),OVR(NDH,NDH),UMT(NDH,NDH),SPC(NDH,NDH)
      DIMENSION PRD(NDH,NDH),DIP(NDH,NDH)
      PRINT*,'WELCOME TO STARK DRIVER, EFIELD=?'
      READ*,EFIELD
      
C     PRINT*,
      G=-1.0D0 
      X=-0.5D0
      P=0.1D0
      TXR=0.2D0
      HAM=0.0D0
      OVR=0.0D0
      DIP=0.0D0
      DIP(1,3)=0.01
      DIP(2,3)=0.03
      DO I=1,3
      DO J=I+1,3
      DIP(J,I)=DIP(I,J)
      END DO
      END DO
      DO I=1,NDH
       OVR(I,I)=1.0D0
      END DO
      HAM(1,1)=-0.5
      HAM(2,2)=-0.125
      HAM(3,3)=-0.125
      HAM(1,3)=DIP(1,3)*EFIELD
      HAM(3,1)=HAM(1,3)
      HAM(2,3)=DIP(2,3)*EFIELD
      HAM(3,2)=HAM(2,3)
      NBS=3
      CALL DIAGnxn(NDH,NBS,HAM,OVR,UMT,PRD,SPC)
      PRINT*,"UPDATED HAM"
      DO I=1,NBS
       PRINT"(10F12.4)",(HAM(J,I),J=1,NBS)  !LAMDA_J
      END DO
      PRINT*,'EIGENVALUES AND EIGENVECTORS:'
      DO I=1,NBS
      PRINT "(10F12.4)",HAM(I,I),(UMT(J,I),J=1,NBS)
C PROVE THAT THE INVERSE OF UMT IS THE TRANPOSE OF UMT
      END DO
      DO J=1,NBS
      DO K=1,NBS
      OVR(K,J)=UMT(J,K)
      END DO
      END DO 
      PRINT*,'2s Wavefunction in terms of new eigenstates'
      PRINT*,(OVR(2,K),K=1,3)
C AT TIME=0, OCCUPY THE 2S FUNCTION:
C |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(K,J)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
      TAU=4*8.0D0*ATAN(1.0D0)/abs(ham(1,1))
      DO M=0,1000
      t=M*(TAU/1000)
      OPEN(12,FILE='PLOT')
      DO K=1,3
      DO L=1,3  
          AR=0.0D0
          AI=0.0d0
          DO I=1,3
          AR=AR+cos(HAM(i,i)*t)*umt(k,i)*umt(l,i) 
          AI=AI+sin(HAM(i,i)*t)*umt(k,i)*umt(l,i) 
          END DO
          PRD(K,L)=AR*AR+AI*AI
      END DO
      END DO
c |phi_2t> = Sum_i u(2,i)exp(i eps_i t) |psi_i>
c <phi_2t| phi_2t> = sum_ij exp(i (eps_i-eps_j)t <psi_j|d|phi_i>*u(2,i)*u(2,j)
      DR=0.0D0       
      DI=0.0D0
         DO i=1,3
         DO j=1,3
         phs=(ham(i,i)-ham(j,j))*t
          DR=DR+cos(phs)*umt(2,i)*umt(2,j)*dip(i,j)
          DI=DI+sin(phs)*umt(2,i)*umt(2,j)*dip(i,j) 
         END DO
         END DO
      DIPOLE=SQRT(DR*DR+DI*DI)*1000
      PRINT "(10F12.4)",t,(PRD(K,2),K=1,NBS),DR,DI             
      WRITE(12,"(10F12.4)")t,(PRD(K,2),K=1,NBS),DIPOLE               
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
C         write(12,12)
C12   FORMAT('p for [col=2:5] "PLOT" u 1:col w lp')    
      CLOSE(12)
      CALL SYSTEM('chmod +x plot_directions')
      CALL SYSTEM('./plot_directions')
      CALL SYSTEM('eog output.png')
      STOP
      END SUBROUTINE STARKDVR

      SUBROUTINE RINGDVR
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NDH=12)
      DIMENSION HAM(NDH,NDH),OVR(NDH,NDH),UMT(NDH,NDH),SPC(NDH,NDH)
      DIMENSION PRD(NDH,NDH),DIP(NDH,NDH)
      PRINT*, 'HOW LARGE IS YOUR RING?'
      READ* ,  RINGSZ
      PRINT*, 'ELECTRIC(1) OR MAGNETIC(2)?'
      READ*, iAP
      PRINT*, 'HOW BIG THE FIELD?'
      READ*, E

      DIP=0.0D0
      DO I=1,NDH
       OVR(I,I)=1.0D0
      END DO
         
      NBS=0
      aPi=3.14159265359D0
      HAM=0.0D0
        aM=-5.0d0
      DO M=-5,5 
       NBS=NBS+1

       IF(NBS.GT.1) THEN
        DIP(NBS-1,NBS)=DIV(E*RINGSZ,2.0D0)
        DIP(NBS,NBS-1)=DIP(NBS-1,NBS)
C       HAM(NBS-1,NBS)=DIV(E*RINGSZ,2.0D0)
C       HAM(NBS,NBS-1)=HAM(NBS-1,NBS)
       END IF       
       HAM(NBS,NBS)=DIV(aM*aM,RINGSZ*RINGSZ)
        aM=aM+1.0d0
      END DO

      IF(iAP.EQ.2)THEN
       DIP=0.0D0
       NBS=0
C      E=-1*DIV(2.0D0,RINGSZ*RINGSZ)
       DO M=-5,5
        NBS=NBS+1
        DIP(NBS,NBS)=M*E 
       END DO
      END IF
      
      DO I=1,NBS
      PRINT ('(11F7.2)'), (HAM(I,J),J=1,NBS)
      END DO
      
      READ*, DUMMY
      DIP=DIP+HAM
      CALL DIAGnxn(NDH,NBS,HAM,OVR,UMT,PRD,SPC)
      PRINT*,"UPDATED HAM"
      DO I=1,NBS
       PRINT"(10F12.4)",(HAM(J,I),J=1,NBS)  !LAMDA_J
      END DO
      PRINT*,'EIGENVALUES AND EIGENVECTORS:'
      DO I=1,NBS
      PRINT "(10F12.4)",HAM(I,I),(UMT(J,I),J=1,NBS)
C PROVE THAT THE INVERSE OF UMT IS THE TRANPOSE OF UMT
      END DO
      DO J=1,NBS
      DO K=1,NBS
      OVR(K,J)=UMT(J,K)
      END DO
      END DO 
      PRINT*,'2s Wavefunction in terms of new eigenstates'
      PRINT*,(OVR(2,K),K=1,NBS)
C AT TIME=0, OCCUPY THE 2S FUNCTION:
C |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(K,J)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
      TAU=4*8.0D0*ATAN(1.0D0)/abs(ham(1,1))
      OPEN(12,FILE='PLOT')
      DO M=-100,100000
      IF(M.EQ.0)THEN
       HAM=DIP
       OVR=0.0D0
       DO I=1,NBS
        OVR(I,I)=1.0D0
       END DO
       CALL DIAGnxn(NDH,NBS,HAM,OVR,UMT,PRD,SPC)
       DO I=1,NBS
        PRINT*, HAM(I,I)
       END DO 
C      STOP
      END IF
      t=M*(TAU/50)
      DO K=1,NBS
      DO L=1,NBS
          AR=0.0D0
          AI=0.0d0
          DO I=1,NBS
          AR=AR+cos(HAM(i,i)*t)*umt(k,i)*umt(l,i) 
          AI=AI+sin(HAM(i,i)*t)*umt(k,i)*umt(l,i) 
          END DO
          PRD(K,L)=AR*AR+AI*AI
      END DO
      END DO
c |phi_2t> = Sum_i u(2,i)exp(i eps_i t) |psi_i>
c <phi_2t| phi_2t> = sum_ij exp(i (eps_i-eps_j)t <psi_j|d|phi_i>*u(2,i)*u(2,j)
      DR=0.0D0       
      DI=0.0D0
         DO i=1,NBS
         DO j=1,NBS
         phs=(ham(i,i)-ham(j,j))*t
          DR=DR+cos(phs)*umt(2,i)*umt(2,j)*dip(i,j)
          DI=DI+sin(phs)*umt(2,i)*umt(2,j)*dip(i,j) 
         END DO
         END DO
      DIPOLE=SQRT(DR*DR+DI*DI)*100
      PRINT "(12F12.4)",t,(PRD(K,2),K=1,NBS),DR,DI             
C     WRITE(12,"(3F12.4)")t,PRD(2,2),PRD(6,2)!,K=1,NBS)!,DIPOLE               
      WRITE(12,"(12F12.4)")t,(PRD(K,2),K=1,NBS)!,DIPOLE               
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
C         write(12,12)
C12   FORMAT('p for [col=2:5] "PLOT" u 1:col w lp')    
      CLOSE(12)
      CALL SYSTEM('chmod +x plot_directions')
      CALL SYSTEM('./plot_directions')
      CALL SYSTEM('eog output.png')
      STOP
      END SUBROUTINE RINGDVR

      SUBROUTINE BOXDVR
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NDH=10)
      DIMENSION HAM(NDH,NDH),OVR(NDH,NDH),UMT(NDH,NDH),SPC(NDH,NDH)
      DIMENSION PRD(NDH,NDH),DIP(NDH,NDH)
      PRINT*, 'WELCOME TO BOX DRIVER, HOW LARGE IS YOUR BOX?'
      READ* ,  BOXSZ
C      PRINT*, DIV(1.0D0,2.0D0)
C      STOP
      aPi=3.14159265359
      Alpha=DIV(2*aPi,BOXSZ)
       A=REAL(M*M)*aPi
C     PRINT*,
      DO I=1,NDH
       OVR(I,I)=1.0D0
      END DO

      NBS=0
      HAM=0.0D0
      DO M=-4,4
       NBS=NBS+1
       TWOM=M*2.0D0*aPi
       PRINT*,M,TWOM
       HAM(NBS,NBS)=DIV(TWOM,BOXSZ)*DIV(TWOM,BOXSZ)
      END DO
      PRINT*,"INITIAL HAM"
      DO I=1,NBS
       PRINT"(10F12.4)",(HAM(J,I),J=1,NBS)  !LAMDA_J
      END DO
      CALL DIAGnxn(NDH,NBS,HAM,OVR,UMT,PRD,SPC)
C      STOP
      PRINT*,"UPDATED HAM"
      DO I=1,NBS
       PRINT"(10F12.4)",(HAM(J,I),J=1,NBS)  !LAMDA_J

      END DO
      PRINT*,'EIGENVALUES AND EIGENVECTORS:'
      DO I=1,NBS
      PRINT "(10F12.4)",HAM(I,I),(UMT(J,I),J=1,NBS)
C PROVE THAT THE INVERSE OF UMT IS THE TRANPOSE OF UMT
      END DO
      DO J=1,NBS
      DO K=1,NBS
      OVR(K,J)=UMT(J,K)
      END DO
      END DO
      PRINT*,'2s Wavefunction in terms of new eigenstates'
      PRINT*,(OVR(2,K),K=1,3)
C AT TIME=0, OCCUPY THE 2S FUNCTION:
C |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(K,J)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
      TAU=4*8.0D0*ATAN(1.0D0)/abs(ham(1,1))
      TAU=1.0D0/abs(ham(1,1)/4.0D0)
      DO M=0,1000
      t=M*(TAU/1000)
      OPEN(12,FILE='PLOT')
      DO K=1,3
      DO L=1,3
          AR=0.0D0
          AI=0.0d0
          DO I=1,3
          AR=AR+cos(HAM(i,i)*t)*umt(k,i)*umt(l,i)
          AI=AI+sin(HAM(i,i)*t)*umt(k,i)*umt(l,i)
      END DO
          PRD(K,L)=AR*AR+AI*AI
      END DO
      END DO
c |phi_2t> = Sum_i u(2,i)exp(i eps_i t) |psi_i>
c <phi_2t| phi_2t> = sum_ij exp(i (eps_i-eps_j)t <psi_j|d|phi_i>*u(2,i)*u(2,j)
      DR=0.0D0
      DI=0.0D0
         DO i=1,3
         DO j=1,3
         phs=(ham(i,i)-ham(j,j))*t
          DR=DR+cos(phs)*umt(2,i)*umt(2,j)*dip(i,j)
          DI=DI+sin(phs)*umt(2,i)*umt(2,j)*dip(i,j)
         END DO
         END DO
      DIPOLE=SQRT(DR*DR+DI*DI)*1000
      PRINT "(10F12.4)",t,(PRD(K,2),K=1,NBS),DR,DI
      WRITE(12,"(10F12.4)")t,(PRD(K,2),K=1,NBS),DIPOLE
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
C         write(12,12)
C12   FORMAT('p for [col=2:5] "PLOT" u 1:col w lp')
      CLOSE(12)
      CALL SYSTEM('chmod +x plot_directions')
      CALL SYSTEM('./plot_directions')
      CALL SYSTEM('xdg-open output.png')
      STOP
      END SUBROUTINE BOXDVR

      SUBROUTINE HMODVR
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NDH=10)
      DIMENSION HAM(NDH,NDH),OVR(NDH,NDH),UMT(NDH,NDH),SPC(NDH,NDH)
      DIMENSION PRD(NDH,NDH),DIP(NDH,NDH)
      print*,'NRLMOL:',0.36005*1.6e-19/6.626E-34
      read*,dummy
      rkmx=1.0D30
      rkmn=0.0D0
      emax=1.0D30
      emin=0.0    
      PRINT*,'WELCOME TO HARMONIC OSCILLATOR DRIVER, EFIELD=?'
      PRINT*,'HCl Diatomic in and Electric Field' 
      PRINT*,'Clorine atom has infinite mass'
      PRINT*,'IR frequency of HCl molecule is: 8.88*10^13 Hz' 
      energy=8.88E13
      print*,'energy:',energy
      PRINT*,'Please googlge this to determine if this is correct'
      PRINT*,'The mass of the Hydrogen is 1.67*10^{-27} kg'
      PRINT*,'Strength of Electric field?'
      hmass=1.67e-27!*(35.)/(36.)
      read*,EField
      DO II=1,30
      PRINT*,'Guess the spring constant in (Newtons/Meter)'
      PRINT*,'Guess should be between:',rkmn,' and ',rkmx
      read*,guess 
      DIP(1,2)=0.01
      DIP(2,1)=DIP(1,2)
      OVR=0.0D0
      OVR=0.0D0
      DO I=1,2
      DO J=I+1,2
      DIP(J,I)=DIP(I,J)
      END DO
      END DO
      DO I=1,NDH
       OVR(I,I)=1.0D13
      END DO
      HAM=0.0D0 
      print '(4G15.6)',guess, hmass, guess/hmass,sqrt(guess/hmass)
      scale=sqrt(guess/hmass)
      HAM(1,1)= 0.5*sqrt(guess/hmass)/(8.0*ATAN(1.0D0))
      HAM(2,2)= 1.5*sqrt(guess/hmass)/(8.0*ATAN(1.0D0))
      HAM(1,2)=0.0!DIP(1,2) 
      HAM(2,1)=0.0!DIP(2,1) 
      NBS=2
            DO I=1,NBS
            PRINT*,(HAM(I,J),J=1,NBS)
            END DO
            DO I=1,NBS
            DO J=1,NBS
            HAM(I,J)=HAM(I,J)/scale                
            END DO
            END DO
      CALL DIAGnxn(NDH,NBS,HAM,OVR,UMT,PRD,SPC)
      DO I=1,NBS
      HAM(I,I)=HAM(I,I)*SCALE
      END DO
      HAM(3,3)=MAX(HAM(2,2),HAM(1,1))-MIN(HAM(2,2),HAM(1,1))
      print"('frequency:',2G12.4)",HAM(3,3),energy
      if(HAM(3,3).GT.energy)then
           print*,'Too Big!'
           if(HAM(3,3).LT.EMAX)THEN
           EMAX=HAM(3,3)
           RKMX=guess
           end if
      end if
      if(HAM(3,3).LT.energy)then
          print*,'Too SMALL!'
           if(HAM(3,3).GT.EMIN)THEN
           EMIN=HAM(3,3)
           RKMN=guess
           end if
      end if
      END DO 
      PRINT*,"UPDATED HAM"
      DO I=1,NBS
       PRINT"(10F12.4)",(HAM(J,I),J=1,NBS)  !LAMDA_J
      END DO
      PRINT*,'EIGENVALUES AND EIGENVECTORS:'
      DO I=1,NBS
      PRINT "(10F12.4)",HAM(I,I),(UMT(J,I),J=1,NBS)
C PROVE THAT THE INVERSE OF UMT IS THE TRANPOSE OF UMT
      END DO
      DO J=1,NBS
      DO K=1,NBS
      OVR(K,J)=UMT(J,K)
      END DO
      END DO 
      PRINT*,'2s Wavefunction in terms of new eigenstates'
      PRINT*,(OVR(2,K),K=1,3)
C AT TIME=0, OCCUPY THE 2S FUNCTION:
C |PHI_K (t) > = SUM_J exp(ie(J) t)* OVR(K,J)|PSI_J> = SUM_JL exp(iEjt)ovr(j,k)*umt(k,l)|PHI_l> !E_j=LAMDA_J
      TAU=4*8.0D0*ATAN(1.0D0)/abs(ham(1,1))
      DO M=0,1000
      t=M*(TAU/1000)
      OPEN(12,FILE='PLOT')
      DO K=1,3
      DO L=1,3  
          AR=0.0D0
          AI=0.0d0
          DO I=1,3
          AR=AR+cos(HAM(i,i)*t)*umt(k,i)*umt(l,i) 
          AI=AI+sin(HAM(i,i)*t)*umt(k,i)*umt(l,i) 
          END DO
          PRD(K,L)=AR*AR+AI*AI
      END DO
      END DO
c |phi_2t> = Sum_i u(2,i)exp(i eps_i t) |psi_i>
c <phi_2t| phi_2t> = sum_ij exp(i (eps_i-eps_j)t <psi_j|d|phi_i>*u(2,i)*u(2,j)
      DR=0.0D0       
      DI=0.0D0
         DO i=1,3
         DO j=1,3
         phs=(ham(i,i)-ham(j,j))*t
          DR=DR+cos(phs)*umt(2,i)*umt(2,j)*dip(i,j)
          DI=DI+sin(phs)*umt(2,i)*umt(2,j)*dip(i,j) 
         END DO
         END DO
      DIPOLE=SQRT(DR*DR+DI*DI)*1000
      PRINT "(10F12.4)",t,(PRD(K,2),K=1,NBS),DR,DI             
      WRITE(12,"(10F12.4)")t,(PRD(K,2),K=1,NBS),DIPOLE               
      END DO
      CLOSE(12)
          open(12,file='plot_directions')
          write(12,12)
 12   FORMAT('p for [col=2:5] "PLOT" u 1:col w lp')
          close(12)
      call system('gnuplot <plot_directions')
      STOP
      END SUBROUTINE HMODVR
