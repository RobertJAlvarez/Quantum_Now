PROGRAM main
  USE applications, ONLY: STARKDVR, RINGDVR, HMODVR, BOXDVR
  IMPLICIT NONE

  INTEGER :: OPT

  DO
    WRITE(*,'(/,A)') 'DIAG options:'
    WRITE(*,'(A)') '1. STARK', '2. RING', '3. HMO', '4. BOX'
    WRITE(*,'(A)') 'n. Exit the program'
    READ (*,*) OPT

    IF (OPT == 1) THEN
      CALL STARKDVR()
    ELSE IF (OPT == 2) THEN
      CALL RINGDVR()
    ELSE IF (OPT == 3) THEN
      CALL HMODVR()
    ELSE IF (OPT == 4) THEN
      CALL BOXDVR()
    ELSE
      EXIT
    END IF
  END DO
END PROGRAM main

