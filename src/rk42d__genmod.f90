        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct  5 14:16:05 2016
        MODULE RK42D__genmod
          INTERFACE 
            SUBROUTINE RK42D(TSTEP_SEC)
              USE PARAMETER_MODULE, ONLY :                              &
     &          GRIDX1,                                                 &
     &          GRIDY1,                                                 &
     &          GRIDX2,                                                 &
     &          GRIDY2,                                                 &
     &          DP,                                                     &
     &          PI,                                                     &
     &          FORWARD,                                                &
     &          X_INPUT,                                                &
     &          Y_INPUT,                                                &
     &          U,                                                      &
     &          V,                                                      &
     &          NX,                                                     &
     &          NY,                                                     &
     &          NX_INPUT,                                               &
     &          NY_INPUT
              REAL(KIND=8), INTENT(IN) :: TSTEP_SEC
            END SUBROUTINE RK42D
          END INTERFACE 
        END MODULE RK42D__genmod
