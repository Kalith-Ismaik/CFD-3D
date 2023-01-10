!COMPUTATIONAL CELL DISCRETZIATION TO ESTIMATE THE FLUID FLOW.
!THIS SUBROUTINE IS WRITTEN BY KALITH M ISMAIL, DECEMBER-2022, NUS SINGAPORE.

    SUBROUTINE Buildup()

        INCLUDE 'common.h'

        REAL :: DIFco
        REAL, DIMENSION(IN2,IN2,IN2) :: D

        D = DC

        X_LOOP: DO I = 2, IN+1

            Y_LOOP: DO J = 2, IN+1

                Z_LOOP: DO K=2, IN+1

                    DIFco = ( RHO * DXYZ2 ) / (2 * Dco2)

                    D(I,J,K) =  DIFco * ( (1/DELTA) * ( (UU(I+1,J,K) - UU(I-1,J,K)/2*XI &
                                + (VV(I,J+1,K) - VV(I,J-1,K))/2*YI + (WW(I,J,K+1) - WW(I,J,K-1))/2*ZI) &
                                - (UU(I+1,J,K) - UU(I-1,J,K)/2*XI)**2 - (VV(I,J+1,K) - VV(I,J-1,K)/2*YI)**2 &
                                - (WW(I,J,K+1) - WW(I,J,K-1))/2*ZI)**2 &
                                - 2*((WW(I,J,K+1) - WW(I,J,K-1))/2*ZI)*((UU(I+1,J,K) - UU(I-1,J,K)/2*XI)) &
                                - 2*((WW(I,J,K+1) - WW(I,J,K-1))/2*ZI)*((VV(I,J+1,K) - VV(I,J-1,K)/2*YI)) &
                                - 2*((UU(I+1,J,K) - UU(I-1,J,K))/2*XI)*((VV(I,J+1,K) - VV(I,J-1,K)/2*YI)))

                ENDDO Z_LOOP

            ENDDO Y_LOOP

        ENDDO X_LOOP
           
        IF(LIDX.EQ.1) THEN

            J1_LOOP: DO J = 2, IN

                K1_LOOP: DO K = 2,IN

                    DIFco = ( RHO * DXYZ2 ) / (2 * Dco2)

                    D(IN2,J,K) =  DIFco * ( (1/DELTA) * ( (UU(1,J,K) - UU(IN2-1,J,K)/2*XI &
                                + (VV(IN2,J+1,K) - VV(IN2,J-1,K))/2*YI + (WW(IN2,J,K+1) - WW(I,J,K-1))/2*ZI) &
                                - (UU(1,J,K) - UU(IN2-1,J,K)/2*XI)**2 - (VV(IN2,J+1,K) - VV(IN2,J-1,K)/2*YI)**2 &
                                - (WW(IN2,J,K+1) - WW(IN2,J,K-1))/2*ZI)**2 &
                                - 2*((WW(IN2,J,K+1) - WW(IN2,J,K-1))/2*ZI)*((UU(1,J,K) - UU(IN2-1,J,K)/2*XI)) &
                                - 2*((WW(IN2,J,K+1) - WW(IN2,J,K-1))/2*ZI)*((VV(IN2,J+1,K) - VV(IN2,J-1,K)/2*YI)) &
                                - 2*((UU(1,J,K) - UU(IN2-1,J,K))/2*XI)*((VV(IN2,J+1,K) - VV(IN2,J-1,K)/2*YI)))

                    D(1,J,K)   =  DIFco * ( (1/DELTA) * ( (UU(2,J,K) - UU(IN2,J,K)/2*XI &
                                + (VV(1,J+1,K) - VV(1,J-1,K))/2*YI + (WW(1,J,K+1) - WW(1,J,K-1))/2*ZI) &
                                - (UU(2,J,K) - UU(IN2,J,K)/2*XI)**2 - (VV(1,J+1,K) - VV(1,J-1,K)/2*YI)**2 &
                                - (WW(1,J,K+1) - WW(1,J,K-1))/2*ZI)**2 &
                                - 2*((WW(1,J,K+1) - WW(1,J,K-1))/2*ZI)*((UU(2,J,K) - UU(IN2,J,K)/2*XI)) &
                                - 2*((WW(1,J,K+1) - WW(1,J,K-1))/2*ZI)*((VV(1,J+1,K) - VV(1,J-1,K)/2*YI)) &
                                - 2*((UU(I+1,J,K) - UU(I-1,J,K))/2*XI)*((VV(1,J+1,K) - VV(1,J-1,K)/2*YI)))

                ENDDO K1_LOOP

            ENDDO J1_LOOP

        ENDIF
           
        IF(LIDZ.EQ.1) THEN

            I2_LOOP: DO I = 2, IN

                J2_LOOP: DO J = 2,IN

                    DIFco = ( RHO * DXYZ2 ) / (2 * Dco2)

                    D(I,J,IN2) =  DIFco * ( (1/DELTA) * ( (UU(I+1,J,IN2) - UU(I-1,J,IN2)/2*XI &
                                + (VV(I,J+1,IN2) - VV(I,J-1,IN2))/2*YI + (WW(I,J,1) - WW(I,J,IN2-1))/2*ZI) &
                                - (UU(I+1,J,IN2) - UU(I-1,J,IN2)/2*XI)**2 - (VV(I,J+1,IN2) - VV(I,J-1,IN2)/2*YI)**2 &
                                - (WW(I,J,1) - WW(I,J,IN2-1))/2*ZI)**2 &
                                - 2*((WW(I,J,1) - WW(I,J,IN2-1))/2*ZI)*((UU(I+1,J,IN2) - UU(I-1,J,IN2)/2*XI)) &
                                - 2*((WW(I,J,1) - WW(I,J,IN2-1))/2*ZI)*((VV(I,J+1,IN2) - VV(I,J-1,IN2)/2*YI)) &
                                - 2*((UU(I+1,J,IN2) - UU(I-1,J,IN2))/2*XI)*((VV(I,J+1,IN2) - VV(I,J-1,IN2)/2*YI)))

                    D(I,J,1) =  DIFco * ( (1/DELTA) * ( (UU(I+1,J,1) - UU(I-1,J,1)/2*XI &
                                + (VV(I,J+1,1) - VV(I,J-1,1))/2*YI + (WW(I,J,2) - WW(I,J,IN2))/2*ZI) &
                                - (UU(I+1,J,1) - UU(I-1,J,1)/2*XI)**2 - (VV(I,J+1,1) - VV(I,J-1,1)/2*YI)**2 &
                                - (WW(I,J,2) - WW(I,J,IN2))/2*ZI)**2 &
                                - 2*((WW(I,J,2) - WW(I,J,IN2))/2*ZI)*((UU(I+1,J,1) - UU(I-1,J,1)/2*XI)) &
                                - 2*((WW(I,J,2) - WW(I,J,IN2))/2*ZI)*((VV(I,J+1,1) - VV(I,J-1,1)/2*YI)) &
                                - 2*((UU(I+1,J,1) - UU(I-1,J,1))/2*XI)*((VV(I,J+1,1) - VV(I,J-1,1)/2*YI)))
                    
                ENDDO J2_LOOP

            ENDDO I2_LOOP
               
        ENDIF
           
        IF(LIDY.EQ.1) THEN

            I3_LOOP: DO I=2,IN

                K3_LOOP: DO K=2,IN

                    DIFco = ( RHO * DXYZ2 ) / (2 * Dco2)

                    D(I,IN2,K) =  DIFco * ( (1/DELTA) * ( (UU(I+1,IN2,K) - UU(I-1,IN2,K)/2*XI &
                                    + (VV(I,1,K) - VV(I,IN2-1,K))/2*YI + (WW(I,IN2,K+1) - WW(I,IN2,K-1))/2*ZI) &
                                    - (UU(I+1,IN2,K) - UU(I-1,IN2,K)/2*XI)**2 - (VV(I,1,K) - VV(I,IN2-1,K)/2*YI)**2 &
                                    - (WW(I,IN2,K+1) - WW(I,IN2,K-1))/2*ZI)**2 &
                                    - 2*((WW(I,IN2,K+1) - WW(I,IN2,K-1))/2*ZI)*((UU(I+1,IN2,K) - UU(I-1,IN2,K)/2*XI)) &
                                    - 2*((WW(I,IN2,K+1) - WW(I,IN2,K-1))/2*ZI)*((VV(I,1,K) - VV(I,IN2-1,K)/2*YI)) &
                                    - 2*((UU(I+1,IN2,K) - UU(I-1,IN2,K))/2*XI)*((VV(I,1,K) - VV(I,IN2-1,K)/2*YI)))

                    D(I,1,K) =  DIFco * ( (1/DELTA) * ( (UU(I+1,1,K) - UU(I-1,1,K)/2*XI &
                                    + (VV(I,2,K) - VV(I,IN2,K))/2*YI + (WW(I,1,K+1) - WW(I,1,K-1))/2*ZI) &
                                    - (UU(I+1,1,K) - UU(I-1,1,K)/2*XI)**2 - (VV(I,2,K) - VV(I,IN2,K)/2*YI)**2 &
                                    - (WW(I,1,K+1) - WW(I,1,K-1))/2*ZI)**2 &
                                    - 2*((WW(I,IN2,K+1) - WW(I,IN2,K-1))/2*ZI)*((UU(I+1,IN2,K) - UU(I-1,IN2,K)/2*XI)) &
                                    - 2*((WW(I,IN2,K+1) - WW(I,IN2,K-1))/2*ZI)*((VV(I,1,K) - VV(I,IN2-1,K)/2*YI)) &
                                    - 2*((UU(I+1,IN2,K) - UU(I-1,IN2,K))/2*XI)*((VV(I,1,K) - VV(I,IN2-1,K)/2*YI)))
            
                ENDDO K3_LOOP

            ENDDO I3_LOOP

        ENDIF

        DC = D

        RETURN

    END SUBROUTINE Buildup