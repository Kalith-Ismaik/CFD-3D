!THIS SUBROUTINE SETS THE BOUNDARY CONDITION FOR CFD SIMULATION.
!THIS SUBROUTINE IS WRITTEN BY KALITH M ISMAIL, DECEMBER-2022, NUS SINGAPORE.

    SUBROUTINE CBound()

        INCLUDE 'common.h'

        B_loop: DO I=2,IN
           
           IF(LIDX.EQ.1) THEN

                UU(1,1:IN+1,1:IN+1)    = 0
                UU(IN-1,1:IN+1,1:IN+1) = 0
                VV(1,1:IN+1,1:IN+1)    = 0
                VV(IN-1,1:IN+1,1:IN+1) = 0
                WW(1,1:IN+1,1:IN+1)    = 0
                WW(IN-1,1:IN+1,1:IN+1) = 0
              
           ENDIF
           
           IF(LIDZ.EQ.1) THEN

                UU(1:IN+1,1:IN+1,1)    = 0
                UU(1:IN+1,1:IN+1,IN-1) = 0
                VV(1:IN+1,1:IN+1,1)    = 0
                VV(1:IN+1,1:IN+1,IN-1) = 0
                WW(1:IN+1,1:IN+1,1)    = 0
                WW(1:IN+1,1:IN+1,IN-1) = 0
               
           ENDIF
           
           IF(LIDY.EQ.1) THEN

                UU(1:IN+1,1,1:IN+1)    = 0
                UU(1:IN+1,IN-1,1:IN+1) = 0
                VV(1:IN+1,1,1:IN+1)    = 0
                VV(1:IN+1,IN-1,1:IN+1) = 0
                WW(1:IN+1,1,1:IN+1)    = 0
                WW(1:IN+1,IN-1,1:IN+1) = 0
              
           ENDIF

           !WW(1:IN+1,IN-1,1:IN+1) = 1.0

        END DO B_loop
        
        RETURN

    END SUBROUTINE CBound

    SUBROUTINE PBound()

        INCLUDE 'common.h'

        P_loop: DO I=2,IN
           
           IF(LIDX.EQ.1) THEN

                PN(1,1:IN+1,1:IN+1)    = PN(2,1:IN+1,1:IN+1)
                PN(IN-1,1:IN+1,1:IN+1) = PN(IN,1:IN+1,1:IN+1)
              
           ENDIF
           
           IF(LIDZ.EQ.1) THEN

                PN(1:IN+1,1:IN+1,1)    = PN(1:IN+1,1:IN+1,2)
                PN(1:IN+1,1:IN+1,IN-1) = PN(1:IN+1,1:IN+1,IN)
               
           ENDIF
           
           IF(LIDY.EQ.1) THEN

                PN(1:IN+1,1,1:IN+1)    = PN(1:IN+1,2,1:IN+1)
                PN(1:IN+1,IN-1,1:IN+1) = PN(1:IN+1,IN,1:IN+1)
              
           ENDIF

           !PN(1:IN+1,IN-1,1:IN+1) = 1.0

        END DO P_loop
        
        RETURN

    END SUBROUTINE PBound