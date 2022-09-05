PROGRAM main
  USE applications, ONLY: STARKDVR, RINGDVR, BOXDVR, HMODVR, write_plot_instructions, open_plot
  IMPLICIT NONE

  INTEGER :: input

  DO
    WRITE(*,'(/,A)') 'DIAG options:'
    WRITE(*,'(A)') '1. STARKDVR', '2. RINGDVR', '3. HMODVR', '4. BOXDVR'
    WRITE(*,'(A)') 'n. Exit the program'
    READ (*,*) input

    SELECT CASE (input)
    CASE (1)
      CALL STARKDVR()
    CASE (2)
      CALL RINGDVR()
    CASE (3)
      CALL HMODVR()
    CASE (4)
      CALL BOXDVR()
    CASE DEFAULT
      WRITE(*,*) 'Have a nice day:)'
      EXIT
    END SELECT

    CALL write_plot_instructions()
    CALL open_plot()
  END DO
END PROGRAM main

