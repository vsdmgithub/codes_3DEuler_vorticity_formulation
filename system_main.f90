! --------------------------------------------------------------
! -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
! CODE BY:
! --------   |         |   ---------        /\        |\      |
! |          |         |  |                /  \       | \     |
! |          |         |  |               /    \      |  \    |
! --------   |         |  |   ------|    /------\     |   \   |
!         |  |         |  |         |   /        \    |    \  |
!         |  |         |  |         |  /          \   |     \ |
! ---------   ----------  ----------  /            \  |      \|
! --------------------------------------------------------------

! ##################
! MODULE: system_main
! LAST MODIFIED: 3 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MAIN RUN  MODULE TO RUN THE TIME EVOLUTION FOR 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_main
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This is the main module. Consists of subroutines
! 1. Pre-analysis
! 2. Time-evolution
! 3. Inter-analysis
! 4. Post-analysis
! All the other modules are sub-modules to this.
! Then nesting of all major sub-modules is as follows

! MAIN-RUN MODULE
!   |
!   ∟ ---> OUTPUT
!   |
!   ∟ ---> SOLVER
!   |
!   ∟ ---> FUNCTIONS
!   |
!   ∟ ---> INITIAL CONDITION
!   |
!   ∟ ---> VARIABLES
!
! There are other modules, which are subsidary modules.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_functions
  USE system_advectionsolver
  USE system_vorticitysolver
  USE system_output
  USE system_analysis

  IMPLICIT NONE
  ! _________________________
  ! LOCAL VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION:: temp_data
  CONTAINS

  SUBROUTINE pre_analysis
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! Loop of time steps, where at each step the spectral velocities
  ! are updated through any of the algoritm. Meanwhile, analysis and
  ! outputs are printed respectively.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !       T    I    M     E              S    T    E    P              C   H    E   C   K
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( ( dt .LE. dt_max ) ) THEN

      check_status = 1

      CALL allocate_velocity
      ! Allocates velocity arrays for the system to start initialisation

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  I  N  I  T  I  A  L        C  O  N  D  I  T  I  O  N
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL normalized_initial_condition
      ! Calls initial condition, then computes energy, then normalizes.

      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
      !      S  O  L  V  E  R     T  Y  P  E
      ! ----------------------------------------------------------------
      !      'ad'- ADVECTION TYPE SOLVER
      !      'vo'- VORTICITY TYPE SOLVER
              solver_type = 'vo'
      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
      !      S  O  L  V  E  R     A  L  G  O  R  I  T  H  M
      ! ----------------------------------------------------------------
      !      'ab'- ADAMBASHFORTH PRED & CORRECTOR ALG
      !      'rk'- RUNGA KUTTA 4TH ORDER ALG
      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
              solver_alg  = 'rk'
      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      CALL allocate_solver

      IF ( run_code .EQ. 'y' ) THEN

        CALL create_output_directories
        ! Creates folders to save files, open files in them to write data.

        CALL allocate_vorticity_moments
        ! Allocates arrays for the moments of Vorticity

        CALL allocate_strain_tensor
        ! Allocates arrays for the strain tensor

        CALL allocate_bck_strain_opr
        ! Allocates arrays for the background strain tensor calculation

      END IF

    ELSE

      check_status =   0

      CALL print_error_timestep()

    END IF

  END

  SUBROUTINE time_evolution
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! Loop of time steps, where at each step the spectral velocities
  ! are updated through any of the algoritm. Meanwhile, analysis and
  ! outputs are printed respectively.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    !             S        T         A         R       T
    ! 8888888888888888888888888888888888888888888888888888888888888888

    !  ---------------------------------------------------------------
    !             T   I   M   E       L  O  O  P
    ! _________________________________________________________________
    DO t_step = 0, t_step_total

      CALL inter_analysis

      IF (debug_error .EQ. 1) THEN
        EXIT  ! Meaning some error in computation.
      END IF

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L  G
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF ( solver_type .EQ. 'ad' ) THEN
        IF ( solver_alg .EQ. 'ab') THEN
          CALL advectionsolver_AB4_algorithm
        ELSE
          CALL advectionsolver_RK4_algorithm
        END IF

      ELSE

        IF ( solver_alg .EQ. 'ab') THEN
          CALL vorticitysolver_AB4_algorithm
        ELSE
          CALL vorticitysolver_RK4_algorithm
        END IF

      END IF
      ! Updates v_x,v_y,v_z for next time step

    END DO
    ! _________________________________________________________________

    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    !                    E     N     D
    ! 8888888888888888888888888888888888888888888888888888888888888888

	END

  SUBROUTINE inter_analysis
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This does all the inter_analysis, making calls to write output during the evolution, debug and statistics part.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    CALL step_to_time_convert(t_step,time_now,dt)
    ! Converts the 't_step' to actual time 'time_now'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  N  A  L  Y  S  I  S       C   A   L   C  .
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    CALL compute_spectral_data

    CALL write_temporal_data

    ! CALL write_test_data

    IF (MOD(t_step,t_step_save) .EQ. 0) THEN

      CALL write_spectral_data

      CALL compute_vorticity_moments
      ! WIll compute the moments and write it from there.

      CALL compute_strain_tensor

      CALL compute_vortex_stretching

      CALL compute_vorticity_dot_moments

      CALL compute_bck_strain_tensor

      CALL compute_bck_vortex_stretching

      CALL compute_vorticity_dot_loc_moments

    END IF

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  D  E  B  U  G             F  O  R          N  a   N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF (MOD(t_step,t_step_debug) .EQ. 0) THEN

      CALL compute_energy_spectral_data

      CALL compute_compressibility

      CALL print_running_status

    END IF

   END

  SUBROUTINE post_analysis
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This does all the post analysis, making calls to write output after the evolution, debug and statistics part.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    CALL deallocate_solver

    ! CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
    ! Making sure, 'v' and 'u' are upto same evolution step

    ! CALL write_spectral_velocity

    ! CALL write_velocity

    CALL deallocate_velocity

    CALL deallocate_operators

    state_sim = 1
    ! Stating that the simulation has ended.

  END

  SUBROUTINE allocate_solver
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This allocates the corresponding solver
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    IF ( solver_type .EQ. 'ad' ) THEN

      CALL allocate_advectionsolver

      IF ( solver_alg .EQ. 'ab') THEN
        CALL allocate_advectionsolver_AB4
      ELSE
        CALL allocate_advectionsolver_RK4
      END IF

    ELSE

      CALL allocate_vorticitysolver

    IF ( solver_alg .EQ. 'ab') THEN
      CALL allocate_vorticitysolver_AB4
    ELSE
      CALL allocate_vorticitysolver_RK4
    END IF

    END IF
    ! Allocates arrays for solving

  END

  SUBROUTINE deallocate_solver
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This deallocates the corresponding solver
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    IF ( solver_type .EQ. 'ad' ) THEN

      CALL deallocate_advectionsolver

      IF ( solver_alg .EQ. 'ab') THEN
        CALL deallocate_advectionsolver_AB4
      ELSE
        CALL deallocate_advectionsolver_RK4
      END IF

    ELSE

      CALL deallocate_vorticitysolver

      IF ( solver_alg .EQ. 'ab') THEN
        CALL deallocate_vorticitysolver_AB4
      ELSE
        CALL deallocate_vorticitysolver_RK4
      END IF

    END IF
    ! Deallocates arrays for solving

  END

 END MODULE system_main
