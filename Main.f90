!THIS PROGRAM IS WRITTEN BY KALITH M ISMAIL, NOVEMBER-2022, NUS SINGAPORE.

    PROGRAM DEM

        INCLUDE 'common.h'

        INTEGER       :: S

        CALL CFDinit()

        CALL CPU_TIME(timestart) ! For simulation time calculation

        OPEN(UNIT = 200,FILE = 'c.out')

        NSTEP = 1000

        S=0
        Main_loop: DO S=1,NSTEP

            TIME=TIME+DELTA

            CALL Navier()

            IF(MOD(S,10).EQ.0) CALL Swrite()

        END DO Main_loop

        CALL CPU_TIME(timeend) ! For simulation time calculation

        print *,"Total CPU time =",timeend-timestart

        WRITE(200,*) VV

        CLOSE (UNIT=200)

        STOP
         
    END PROGRAM DEM