!NAVIER STOKES ALGORITHM TO INTEGRATE FLUID FLOW.
!THIS SUBROUTINE IS WRITTEN BY KALITH M ISMAIL, DECEMBER-2022, NUS SINGAPORE.

    SUBROUTINE Navier()

        INCLUDE 'common.h'

        REAL, DIMENSION(IN2,IN2,IN2) :: UI,VI,WI

        CALL Poisson()
        !PN(1:IN2,1:IN2,1:IN2) = 1.0

        UI = UU
        VI = VV
        WI = WW

        X_LOOP: DO I = 2, IN+1

            Y_LOOP: DO J = 2, IN+1

                Z_LOOP: DO K=2, IN+1

                    UI(I,J,K) = UU(I,J,K) - UU(I,J,K)*(DELTA/XI)*(UU(I,J,K) - UU(I-1,J,K)) &
                                - VV(I,J,K)*(DELTA/YI)*(UU(I,J,K) - UU(I,J-1,K)) &
                                - WW(I,J,K)*(DELTA/ZI)*(UU(I,J,K) - UU(I,J,K-1)) &
                                - (DELTA/(2*XI*RHO))*(PN(I+1,J,K)-PN(I-1,J,K)) &
                                - MU*((DELTA/XI**2)*(UU(I-1,J,K)+UU(I+1,J,K)-2*UU(I,J,K)) &
                                + (DELTA/YI**2)*(UU(I,J-1,K)+UU(I,J+1,K)-2*UU(I,J,K)) &
                                + (DELTA/ZI**2)*(UU(I,J,K-1)+UU(I,J,K+1)-2*UU(I,J,K))) + CF(1)*DELTA

                    VI(I,J,K) = VV(I,J,K) - UU(I,J,K)*(DELTA/XI)*(VV(I,J,K) - VV(I-1,J,K)) &
                                - VV(I,J,K)*(DELTA/YI)*(VV(I,J,K) - VV(I,J-1,K)) &
                                - WW(I,J,K)*(DELTA/ZI)*(WW(I,J,K) - WW(I,J,K-1)) &
                                - (DELTA/(2*XI*RHO))*(PN(I,J+1,K)-PN(I,J-1,K)) &
                                - MU*((DELTA/XI**2)*(VV(I-1,J,K)+VV(I+1,J,K)-2*VV(I,J,K)) &
                                + (DELTA/YI**2)*(VV(I,J-1,K)+VV(I,J+1,K)-2*VV(I,J,K)) &
                                + (DELTA/ZI**2)*(VV(I,J,K-1)+VV(I,J,K+1)-2*VV(I,J,K))) + CF(2)*DELTA

                    WI(I,J,K) = WW(I,J,K) - UU(I,J,K)*(DELTA/XI)*(WW(I,J,K) - WW(I-1,J,K)) &
                                - VV(I,J,K)*(DELTA/YI)*(WW(I,J,K) - WW(I,J-1,K)) &
                                - WW(I,J,K)*(DELTA/ZI)*(WW(I,J,K) - WW(I,J,K-1)) &
                                - (DELTA/(2*XI*RHO))*(PN(I,J,K+1)-PN(I,J,K-1)) &
                                - MU*((DELTA/XI**2)*(WW(I-1,J,K)+WW(I+1,J,K)-2*WW(I,J,K)) &
                                + (DELTA/YI**2)*(WW(I,J-1,K)+WW(I,J+1,K)-2*WW(I,J,K)) &
                                + (DELTA/ZI**2)*(WW(I,J,K-1)+WW(I,J,K+1)-2*WW(I,J,K))) + CF(3)*DELTA

                ENDDO Z_LOOP

            ENDDO Y_LOOP

        ENDDO X_LOOP

        IF(KBOUND.EQ.0) THEN

            I_LOOP: DO I=2, IN+1

                J_LOOP: DO J = 2, IN+1

                    K_LOOP: DO K=2, IN+1

                        UI(IN2,J,K) = UU(IN2,J,K) - UU(IN2,J,K)*(DELTA/XI)*(UU(IN2,J,K) - UU(IN2-1,J,K)) &
                                    - VV(IN2,J,K)*(DELTA/YI)*(UU(IN2,J,K) - UU(IN2,J-1,K)) &
                                    - WW(IN2,J,K)*(DELTA/ZI)*(UU(IN2,J,K) - UU(IN2,J,K-1)) &
                                    - (DELTA/(2*XI*RHO))*(PN(1,J,K)-PN(IN2-1,J,K)) &
                                    - MU*((DELTA/XI**2)*(UU(IN2-1,J,K)+UU(1,J,K)-2*UU(IN2,J,K)) &
                                    + (DELTA/YI**2)*(UU(IN2,J-1,K)+UU(IN2,J+1,K)-2*UU(IN2,J,K)) &
                                    + (DELTA/ZI**2)*(UU(IN2,J,K-1)+UU(IN2,J,K+1)-2*UU(IN2,J,K))) + CF(1)*DELTA

                        VI(I,IN2,K) = VV(I,IN2,K) - UU(I,IN2,K)*(DELTA/XI)*(VV(I,IN2,K) - VV(I-1,IN2,K)) &
                                    - VV(I,IN2,K)*(DELTA/YI)*(VV(I,IN2,K) - VV(I,IN2-1,K)) &
                                    - WW(I,IN2,K)*(DELTA/ZI)*(WW(I,IN2,K) - WW(I,IN2,K-1)) &
                                    - (DELTA/(2*XI*RHO))*(PN(I,1,K)-PN(I,IN2-1,K)) &
                                    - MU*((DELTA/XI**2)*(VV(I-1,IN2,K)+VV(I+1,IN2,K)-2*VV(I,IN2,K)) &
                                    + (DELTA/YI**2)*(VV(I,IN2-1,K)+VV(I,1,K)-2*VV(I,IN2,K)) &
                                    + (DELTA/ZI**2)*(VV(I,IN2,K-1)+VV(I,IN2,K+1)-2*VV(I,IN2,K))) + CF(2)*DELTA

                        WI(I,J,IN2) = WW(I,J,IN2) - UU(I,J,IN2)*(DELTA/XI)*(WW(I,J,IN2) - WW(I-1,J,IN2)) &
                                    - VV(I,J,IN2)*(DELTA/YI)*(WW(I,J,IN2) - WW(I,J-1,IN2)) &
                                    - WW(I,J,IN2)*(DELTA/ZI)*(WW(I,J,IN2) - WW(I,J,IN2-1)) &
                                    - (DELTA/(2*XI*RHO))*(PN(I,J,1)-PN(I,J,IN2-1)) &
                                    - MU*((DELTA/XI**2)*(WW(I-1,J,IN2)+WW(I+1,J,IN2)-2*WW(I,J,IN2)) &
                                    + (DELTA/YI**2)*(WW(I,J-1,IN2)+WW(I,J+1,IN2)-2*WW(I,J,IN2)) &
                                    + (DELTA/ZI**2)*(WW(I,J,IN2-1)+WW(I,J,1)-2*WW(I,J,IN2))) + CF(3)*DELTA

                        UI(1,J,K) = UU(1,J,K) - UU(1,J,K)*(DELTA/XI)*(UU(1,J,K) - UU(IN2,J,K)) &
                                    - VV(1,J,K)*(DELTA/YI)*(UU(1,J,K) - UU(1,J-1,K)) &
                                    - WW(1,J,K)*(DELTA/ZI)*(UU(1,J,K) - UU(1,J,K-1)) &
                                    - (DELTA/(2*XI*RHO))*(PN(2,J,K)-PN(IN2,J,K)) &
                                    - MU*((DELTA/XI**2)*(UU(2,J,K)+UU(IN2,J,K)-2*UU(1,J,K)) &
                                    + (DELTA/YI**2)*(UU(1,J-1,K)+UU(1,J+1,K)-2*UU(1,J,K)) &
                                    + (DELTA/ZI**2)*(UU(1,J,K-1)+UU(1,J,K+1)-2*UU(1,J,K))) + CF(1)*DELTA

                        VI(I,1,K) = VV(I,1,K) - UU(I,1,K)*(DELTA/XI)*(VV(I,1,K) - VV(I-1,1,K)) &
                                    - VV(I,1,K)*(DELTA/YI)*(VV(I,1,K) - VV(I,IN2,K)) &
                                    - WW(I,1,K)*(DELTA/ZI)*(WW(I,1,K) - WW(I,1,K-1)) &
                                    - (DELTA/(2*XI*RHO))*(PN(I,2,K)-PN(I,IN2,K)) &
                                    - MU*((DELTA/XI**2)*(VV(I-1,1,K)+VV(I+1,1,K)-2*VV(I,1,K)) &
                                    + (DELTA/YI**2)*(VV(I,IN2,K)+VV(I,2,K)-2*VV(I,1,K)) &
                                    + (DELTA/ZI**2)*(VV(I,1,K-1)+VV(I,1,K+1)-2*VV(I,1,K))) + CF(2)*DELTA

                        WI(I,J,1) = WW(I,J,1) - UU(I,J,1)*(DELTA/XI)*(WW(I,J,1) - WW(I-1,J,1)) &
                                    - VV(I,J,1)*(DELTA/YI)*(WW(I,J,1) - WW(I,J-1,1)) &
                                    - WW(I,J,1)*(DELTA/ZI)*(WW(I,J,1) - WW(I,J,IN2)) &
                                    - (DELTA/(2*XI*RHO))*(PN(I,J,2)-PN(I,J,IN2)) &
                                    - MU*((DELTA/XI**2)*(WW(I-1,J,1)+WW(I+1,J,1)-2*WW(I,J,1)) &
                                    + (DELTA/YI**2)*(WW(I,J-1,1)+WW(I,J+1,1)-2*WW(I,J,1)) &
                                    + (DELTA/ZI**2)*(WW(I,J,IN2)+WW(I,J,2)-2*WW(I,J,1))) + CF(3)*DELTA

                    ENDDO K_LOOP

                ENDDO J_LOOP

            ENDDO I_LOOP

        ENDIF

        UU = UI
        VV = VI
        WW = WI

        CALL CBound()

        RETURN

    END SUBROUTINE Navier