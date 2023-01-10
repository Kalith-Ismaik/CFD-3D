!POISSON PRESSURE CALCULATION TO ESTIMATE THE FLUID FLOW.
!THIS SUBROUTINE IS WRITTEN BY KALITH M ISMAIL, DECEMBER-2022, NUS SINGAPORE.

    SUBROUTINE Poisson()

        INCLUDE 'common.h'

        REAL, DIMENSION(IN2,IN2,IN2) :: P

        CALL Buildup()

        !DC(1:IN2,1:IN2,1:IN2) = 1.0

        P = PN

        X_LOOP: DO I = 2, IN+1

            Y_LOOP: DO J = 2, IN+1

                Z_LOOP: DO K=2, IN+1

                    P(I,J,K) =  (DYZ2*(P(I+1,J,K)+P(I-1,J,K)) + DXZ2*(P(I,J+1,K)+P(I,J-1,K)) &
                                + DXY2*(P(I,J,K+1)+P(I,J,K-1))) / (2*Dco2)
                    P(I,J,K) = P(I,J,K) - DC(I,J,K) 

                ENDDO Z_LOOP

            ENDDO Y_LOOP

        ENDDO X_LOOP

        IF(LIDX.EQ.2) THEN

            J1_LOOP: DO J=2, IN+1

                K1_LOOP: DO K=2, IN+1

                    P(IN2,J,K) =  (DYZ2*(P(1,J,K)+P(IN2-1,J,K)) + DXZ2*(P(IN2,J+1,K)+P(IN2,J-1,K)) &
                                + DXY2*(P(IN2,J,K+1)+P(IN2,J,K-1))) / (2*Dco2)
                    P(IN2,J,K) = P(IN2,J,K) - DC(IN2,J,K) 

                    P(1,J,K) =  (DYZ2*(P(2,J,K)+P(IN2,J,K)) + DXZ2*(P(1,J+1,K)+P(1,J-1,K)) &
                                + DXY2*(P(1,J,K+1)+P(1,J,K-1))) / (2*Dco2)
                    P(1,J,K) = P(1,J,K) - DC(1,J,K) 

                ENDDO K1_LOOP

            ENDDO J1_LOOP

        ENDIF

        IF(LIDY.EQ.2) THEN

            I2_LOOP: DO I=2, IN+1

                K2_LOOP: DO K=2, IN+1

                    P(I,IN2,K) =  (DYZ2*(P(I+1,IN2,K)+P(I-1,IN2,K)) + DXZ2*(P(I,1,K)+P(I,IN2-1,K)) &
                                    + DXY2*(P(I,IN2,K+1)+P(I,IN2,K-1))) / (2*Dco2)
                    P(I,IN2,K) = P(I,IN2,K) - DC(I,IN2,K) 

                    P(I,1,K) =  (DYZ2*(P(I+1,1,K)+P(I-1,1,K)) + DXZ2*(P(I,2,K)+P(I,IN2,K)) &
                                    + DXY2*(P(I,1,K+1)+P(I,1,K-1))) / (2*Dco2)
                    P(I,1,K) = P(I,1,K) - DC(I,1,K) 

                ENDDO K2_LOOP

            ENDDO I2_LOOP

        ENDIF

        IF(LIDZ.EQ.2) THEN

            I3_LOOP: DO I=2, IN+1

                J3_LOOP: DO J=2, IN+1

                    P(I,J,IN2) =  (DYZ2*(P(I+1,J,IN2)+P(I-1,J,IN2)) + DXZ2*(P(I,J+1,IN2)+P(I,J-1,IN2)) &
                                    + DXY2*(P(I,J,1)+P(I,J,IN2-1))) / (2*Dco2)
                    P(I,J,IN2) =  P(I,J,IN2) - DC(I,J,IN2)

                    P(I,J,1)   =  (DYZ2*(P(I+1,J,1)+P(I-1,J,1)) + DXZ2*(P(I,J+1,1)+P(I,J-1,1)) &
                                    + DXY2*(P(I,2,K+1)+P(I,IN2,K-1))) / (2*Dco2)
                    P(I,J,1)   =  P(I,J,K) - DC(I,J,K) 

                ENDDO J3_LOOP
        
            ENDDO I3_LOOP

        ENDIF

        PN = P

        CALL PBound()

        RETURN

    END SUBROUTINE Poisson