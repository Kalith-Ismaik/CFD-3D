!WRITE OUTPUT FILES FOR SURTHER ANALYSIS.
!THIS SUBROUTINE IS WRITTEN BY KALITH M ISMAIL, NOVEMBER-2022, NUS SINGAPORE.

    SUBROUTINE Swrite()
        INCLUDE 'common.h'
        CHARACTER*5 chtime
        REAL          :: II,JJ,KK

        ITIME=ANINT(TIME)
        IF (ITIME.LT.1000) Then
           WRITE(chtime, FMT='(I3.3)') ITIME
        ELSEIF (ITIME.LT.10000) Then
           WRITE(chtime, FMT='(I4.4)') ITIME
        ELSEIF (ITIME.LT.100000) Then
           WRITE(chtime, FMT='(I5.5)') ITIME
        ELSE
           PRINT *, "Simulation is too long... Time=",Time," ps."
           STOP
        ENDIF

        OPEN(UNIT=969,FILE='./'//DRNAME(1:LDR)//'/time'//trim(chtime)//'.d')
        REWIND 969

        REWIND 969
        I_LOOP: DO I=2, IN+1
            J_LOOP: DO J=2, IN+1
                K_LOOP: DO K=2, IN+1
                    II = (I-1)*XI
                    JJ = (J-1)*YI
                    KK = (K-1)*ZI
                    WRITE(969,110) II, JJ, KK, UU(I,J,K), VV(I,J,K), WW(I,J,K)
                ENDDO K_LOOP
            ENDDO J_LOOP
        ENDDO I_LOOP

110     FORMAT(6(E10.4,1x))

        CLOSE (UNIT = 969)

        RETURN
    END SUBROUTINE Swrite