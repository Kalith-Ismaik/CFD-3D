      IMPLICIT REAL*8(A-H,O-Z)
      character*14 drname
      LOGICAL FullList
      INCLUDE 'parameters.h'

!     Sizes of the computational cell
      COMMON/CELL/ XL,YL,ZL

!     Sizes of the CFD computational cell
      COMMON/NELL/ XI,YI,ZI

!     Different constants of system
      COMMON/CON/ DXYZ2,Dco2,DXY2,DXZ2,DYZ2

!     CFD velocities
      COMMON/CXX/ UU(IN2,IN2,IN2),VV(IN2,IN2,IN2),WW(IN2,IN2,IN2),CF(ID)

!     CFD Pressure
      COMMON/PRR/DC(IN2,IN2,IN2),PN(IN2,IN2,IN2)

!     CFD parameters
      COMMON/RVS/RHO,MU

!     INPUT
      COMMON/INP/LIDX,LIDY,LIDZ
      
!     T-factor
      COMMON/TIM/TIME,DELTA
