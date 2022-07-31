      PROGRAM temp_f
      IMPLICIT REAL*8 (A-H,O-Z)
      CALL DIAGDVR
      END PROGRAM
      
      SUBROUTINE J2X2(H,E,O)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(2,2),E(2),O(2,2)
      RAD=SQRT((H(1,1)-H(2,2))*(H(1,1)-H(2,2))+4.D0*H(1,2)*H(2,1))
      TRC=     H(1,1)+H(2,2) !Trace=Sum of the diagonal of a square matrix.
      E(1)=(0.5D0)*(TRC+RAD)
      E(2)=(0.5D0)*(TRC-RAD) !Its easy to see that E(1)+E(2)=TRC
      
      O(1,1)=-H(1,2)
      O(2,1)= H(1,1)-E(1)
      O(1,2)= H(2,2)-E(2)
      O(2,2)=-H(1,2) 
      
      DO I=1,2
       DOT=O(1,I)*O(1,I)+O(2,I)*O(2,I)  
       DOT=SQRT(DOT)                     
       O(1,I)=O(1,I)/DOT
       O(2,I)=O(2,I)/DOT
      ENDDO 
      RETURN
      END SUBROUTINE J2X2

      SUBROUTINE DIAGDVR
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NDH=12)
      DIMENSION HAM(NDH,NDH),UMT(NDH,NDH)
      DIMENSION PRD(NDH,NDH)

      G  = -1.D0 
      X  = -0.5D0
      P  =  0.1D0
      TXR=  0.2D0
      HAM=  0.D0
      
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
      
      DO I=1,NBS
       DO J=I,NBS
        HAM(J,I)=HAM(I,J)
       END DO
       PRINT'(10F7.2)',(HAM(I,J),J=1,NBS)
      END DO

      CALL DIAGnxn(NDH,NBS,HAM,UMT,PRD)
      PRINT*,"UPDATED HAM"
      DO I=1,NBS
       WRITE(*,'(10F11.6)') (HAM(j,i), j=1,NBS)
      END DO
      END SUBROUTINE DIAGDVR

      SUBROUTINE DIAGnxn(NDH,NBS,HAM,UMT,PRD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION HAM(NDH,NDH),UMT(NDH,NDH),PRD(NDH,NDH)
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
      
      IF(NBS.GE.NDH) THEN 
       PRINT*, 'NDH MUST BE LARGER THAN', NBS
       STOP
      END IF
      UMT=0.D0
      DO I=1,NBS
       UMT(I,I)=1.D0
      END DO

      PRD=HAM
      
      MXIT=NBS*NBS*2
      DO 1000 ITRY=1,MXIT
      FOM=0.D0
      INDEX=0
      KNT=0
      DO I=1,NBS
       DO J=I+1,NBS
        IF(ABS(PRD(J,I)).GT.0.D0)THEN
         KNT=KNT+1
         FOM(KNT)=ABS(PRD(J,I))
         INDEX(1,KNT)=I
         INDEX(2,KNT)=J
        END IF
       END DO
      END DO
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
      NPROC=2
      PROCS=0.D0
      DO IPROC=1,NPROC
       FOUND=.FALSE.
       DO I=1,KNT 
        K=INDEX(1,I)
        J=INDEX(2,I)
        IF(UPDATE(K).AND.UPDATE(J).AND.FOM(I).GT.0.D0)THEN
         UPDATE(K)=.FALSE.
         UPDATE(J)=.FALSE.
         PROCS(1,IPROC)=K
         PROCS(2,IPROC)=J
         GOTO 555
        END IF
       END DO
 555   CONTINUE
      END DO
      UPDATE=.TRUE.
      
      ERROLD=0.D0
      OFF_NRM=0.D0
      DO I=1,NBS
      DO J=I+1,NBS
       ERROLD=ERROLD+PRD(I,J)*PRD(I,J)
      END DO
      END DO
      
      WRITE(*,*) 'Non repetitive indexes with highest values:'
      DO i=1, NPROC
        WRITE(*,*) procs(1,i), procs(2,i), PRD(procs(1,i),procs(2,i))
      END DO
      
      K=PROCS(1,1)
      L=PROCS(2,1)
      
      H(1,1)=PRD(K,K)
      H(1,2)=PRD(K,L)
      H(2,1)=PRD(L,K)
      H(2,2)=PRD(L,L)
      OFF_NRM=OFF_NRM+H(1,2)*H(1,2)
      CALL J2X2(H,E,O)
      
      WRITE(*,*) 'E and O values:'
      DO i=1, 2
        WRITE(*,'(3ES20.12E2)') E(i), (O(i,j),j=1,2)
      END DO
      
      SPC=0.D0
      DO I=1,NBS
       SPC(I,I)=1.D0
      END DO
      
      SPC(K,K)=O(1,1)
      SPC(L,K)=O(2,1)
      SPC(L,L)=O(2,2)
      SPC(K,L)=O(1,2)
      
      PRD=0.D0
      DO i=1,NBS
       IF(i.NE.K .AND. i.NE.L) THEN
        DO j=1,NBS
         PRD(j,i)=UMT(j,i)
        END DO
       END IF
      END DO
C     
      DO i=1,NBS
       PRD(i,K)=PRD(i,K)+UMT(i,K)*O(1,1)
       PRD(i,K)=PRD(i,K)+UMT(i,L)*O(2,1)
       PRD(i,L)=PRD(i,L)+UMT(i,K)*O(1,2)
       PRD(i,L)=PRD(i,L)+UMT(i,L)*O(2,2)
      END DO

      UMT=PRD
      DO I=1,NBS
      DO K=1,NBS
       SPC(K,I)=0.D0
      DO L=1,NBS
       SPC(K,I)=SPC(K,I)+UMT(L,I)*HAM(L,K)
      END DO
      END DO
      END DO

      PRD=0.D0
      DO I=1,NBS
      DO J=1,NBS
      DO K=1,NBS
       PRD(J,I)=PRD(J,I)+UMT(K,J)*SPC(K,I)
      END DO
      END DO
      END DO
      
      ERRNW=0.D0
      DO I=1,NBS
       DO J=I+1,NBS
        ERRNW=ERRNW+PRD(I,J)*PRD(J,I)
       END DO
      END DO

      DO I=1,NBS
       WRITE(*,'(7F10.5)') (PRD(j,i), j=1,NBS)
      END DO

      PRINT 197,ITRY,ERRNW,ERROLD
 197  FORMAT(I3,3G15.6)
      IF(ERRNW.LT.1.D-12) GO TO 1010
1000  CONTINUE
      PRINT*, 'WARNING: NO CONVERGENCE',MXIT
1010  CONTINUE
      PRINT*, ITRY,NBS,Float(ITRY)/(NBS*NBS), 'Diag Eff'
      HAM=PRD
      RETURN
      END SUBROUTINE DIAGnxn

