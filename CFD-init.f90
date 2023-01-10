!THIS SUBROUTINE SETS THE INITIAL PARAMETERS FOR CFD SIMULATION.
!THIS SUBROUTINE IS WRITTEN BY KALITH M ISMAIL, DECEMBER-2022, NUS SINGAPORE.

    SUBROUTINE CFDinit()

        INCLUDE 'common.h'

        RHO = 1.0
        MU  = 0.1

        DELTA = 0.1
        TIME = 0.0

        XL = 10.0
        YL = 10.0
        ZL = 10.0

        CF(1) = 1.0
        CF(2) = 1.0
        CF(3) = 1.0

        XI  = XL/IN
        YI  = YL/IN
        ZI  = ZL/IN

        DXY2 = (XI**2)*(YI**2)
        DXZ2 = (XI**2)*(ZI**2)
        DYZ2 = (YI**2)*(ZI**2)

        DXYZ2 = (XI**2)*(YI**2)*(ZI**2)
        Dco2  = DXY2+DXZ2+DYZ2

        PRINT *, DXY2
        PRINT *, DXZ2
        PRINT *, DYZ2
        PRINT *, DXYZ2
        PRINT *, Dco2

        UU(1:IN2,1:IN2,1:IN2) = 0.0
        VV(1:IN2,1:IN2,1:IN2) = 0.0
        WW(1:IN2,1:IN2,1:IN2) = 0.0

        DC(1:IN2,1:IN2,1:IN2) = 0.0
        PN(1:IN2,1:IN2,1:IN2) = 1.0

        LIDX = 1
        LIDY = 1
        LIDZ = 1

        RETURN

    END SUBROUTINE CFDinit