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

! #########################
! MODULE: system_basicfunctions
! LAST MODIFIED: 21 JUNE 2021
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! BASIC FUNCTIONS MODULE FOR 3D EULER ANALYSIS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_basicfunctions
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major system basic functions for the code to run.
! THese are standard functions, more advanced are done in advanced functions module
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_initialcondition
  USE system_basicoutput

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE normalized_initial_condition
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS TO :
  ! Get a normalized initial condition
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  I  N  I  T  I  A  L        C  O  N  D  I  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL init_initcondn
    ! Calls the subroutine to get a initial condition
    ! REF-> <<< system_initialcondition >>>

    CALL compute_energy_spectral_data
    ! REF-> <<< system_initialcondition >>>

    CALL perform_debug

    IF ( ( NaN_count .EQ. 0 ) .AND. ( k_dot_v_error .EQ. 0) ) THEN

      check_status = 1

      CALL compute_spectral_data
      ! Gets the energy,enstrophy from spectral space

      CALL fft_c2r_vec( v_x, v_y, v_z, u_x, u_y, u_z )
      ! FFT spectral to real velocity

    END IF

  END

  SUBROUTINE compute_spectral_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This calculates energy, spectral shell wise. It goes through each
  ! spectral mode and puts the energy in the corresponding shell.
  ! This gives the ENERGY SPECTRUM.  When the time is right, it saves
  ! them too.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    spectral_energy     = zero
    spectral_energy_avg = zero
    ! Reset the array

    energy              = zero
    enstrophy           = zero
    ! Reset the variables

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y     S  P  E  C  T  R  U  M     C  A   L   C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Gets the energy in that particular mode (i_x,i_y,i_z) into 'en_temp'
    ! Keeps adding energy to that particular shell (|k| fixed), from all solid angles
    DO i_x =    1 , kTru_x
    DO i_y = - kTru_y, kTru_y
    DO i_z = - kTru_z, kTru_z
    IF ( k_2 ( i_x, i_y, i_z ) .LT. k_G_2 ) THEN
      energy_mode                                   = CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                                      CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                                      CDABS( v_z( i_x, i_y, i_z ) ) ** two
      spectral_energy( shell_no ( i_x, i_y, i_z ) ) = spectral_energy( shell_no ( i_x, i_y, i_z ) ) + energy_mode
      enstrophy                                     = enstrophy + k_2( i_x, i_y, i_z) * energy_mode
    END IF
    END DO
    END DO
    END DO

    i_x    =   0
    DO i_y = - kTru_y, kTru_y
    DO i_z = - kTru_z, -1
    IF ( k_2 ( i_x, i_y, i_z ) .LT. k_G_2 ) THEN
      energy_mode                                   = CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                                      CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                                      CDABS( v_z( i_x, i_y, i_z ) ) ** two
      spectral_energy( shell_no ( i_x, i_y, i_z ) ) = spectral_energy( shell_no ( i_x, i_y, i_z ) ) + energy_mode
      enstrophy                                     = enstrophy + k_2( i_x, i_y, i_z) * energy_mode
    END IF
    END DO
    END DO
    i_z    = 0
    DO i_y = 0, kTru_y
    IF ( k_2 ( i_x, i_y, i_z ) .LT. k_G_2 ) THEN
      energy_mode                                   = CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                                      CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                                      CDABS( v_z( i_x, i_y, i_z ) ) ** two
      spectral_energy( shell_no ( i_x, i_y, i_z ) ) = spectral_energy( shell_no ( i_x, i_y, i_z ) ) + energy_mode
      enstrophy                                     = enstrophy + k_2( i_x, i_y, i_z) * energy_mode
    END IF
    END DO

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  S  H  E  L  L      A  V  E  R  A  G  I  N  G
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    spectral_energy_avg( 1 )            = qtr * ( thr * spectral_energy( 1 ) + spectral_energy( 2 ) )
    spectral_energy_avg( max_wave_no ) = qtr * ( thr * spectral_energy( max_wave_no ) + &
                                                        spectral_energy( max_wave_no - 1 ) )

    DO k_no                             = 2, max_wave_no - 1

        spectral_energy_avg( k_no )     = qtr * ( spectral_energy( k_no - 1 ) + spectral_energy( k_no + 1 ) ) + &
                                           hf * ( spectral_energy( k_no ) )

    END DO

    energy = SUM( spectral_energy( : ) )
    ! Computes the net energy

  END

  SUBROUTINE compute_energy
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the Kinetic energy in real space
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    energy = hf * SUM( u_x ** two + u_y ** two + u_z ** two ) / N3

  END

  SUBROUTINE compute_vorticity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get vorticity field
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    w_vx = i * ( k_y * v_z - k_z * v_y )
    w_vy = i * ( k_z * v_x - k_x * v_z )
    w_vz = i * ( k_x * v_y - k_y * v_x )
    ! Spectral Vorticity

    CALL fft_c2r_vec( w_vx, w_vy, w_vz, w_ux, w_uy, w_uz )
    ! Real Vorticity

  END

  SUBROUTINE compute_enstrophy
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the enstrophy in real space
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    enstrophy = hf * SUM( w_ux ** two + w_uy ** two + w_uz ** two ) / N3

  END

  SUBROUTINE perform_debug
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check incompressibility criterion and Nan in data
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    CALL check_nan

    CALL compute_compressibility

  END

  SUBROUTINE compute_compressibility
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check incompressibility condition. Sums over all residues
  ! of incompressibility and prints it. Of order 10^(-12).
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    k_dot_v_norm   = zero

    DO i_x         = kMin_x + 1, kMax_x - 1
    DO i_y         = kMin_y , kMax_y
    DO i_z         = kMin_z , kMax_z
      k_dot_v_norm = k_dot_v_norm + CDABS( k_x( i_x, i_y, i_z ) * v_x( i_x, i_y, i_z ) + &
                                           k_y( i_x, i_y, i_z ) * v_y( i_x, i_y, i_z ) + &
                                           k_z( i_x, i_y, i_z ) * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO
    END DO

    i_x            = kMin_x
    DO i_y         = kMin_y , kMax_y
    DO i_z         = kMin_z , kMax_z
      k_dot_v_norm = k_dot_v_norm + hf* CDABS( k_x( i_x, i_y, i_z ) * v_x( i_x, i_y, i_z ) + &
                                               k_y( i_x, i_y, i_z ) * v_y( i_x, i_y, i_z ) + &
                                               k_z( i_x, i_y, i_z ) * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO

    i_x            = kMax_x
    DO i_y         = kMin_y , kMax_y
    DO i_z         = kMin_z , kMax_z
      k_dot_v_norm = k_dot_v_norm + hf* CDABS( k_x( i_x, i_y, i_z ) * v_x( i_x, i_y, i_z ) + &
                                               k_y( i_x, i_y, i_z ) * v_y( i_x, i_y, i_z ) + &
                                               k_z( i_x, i_y, i_z ) * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO

    k_dot_v_norm = DSQRT( k_dot_v_norm )

    IF (k_dot_v_norm .GT. tol_float ) THEN

      k_dot_v_error = 1

      debug_error   = 1
      ! This will jump out of evolution loop, if caught during that.

      CALL print_error_incomp
      ! This will prompt error when checked
      ! REF-> <<< system_output >>>

    END IF
  END

  SUBROUTINE check_nan
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check if there is any NaN in the data.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    NaN_count  = 0

    DO i_x         = kMin_x + 1, kMax_x - 1
    DO i_y         = kMin_y , kMax_y
    DO i_z         = kMin_z , kMax_z
      IF ( v_x( i_x, i_y, i_z ) .NE. v_x( i_x, i_y, i_z ) ) THEN
        NaN_count = NaN_count + 1
      END IF
    END DO
    END DO
    END DO

    i_x            = kMin_x
    DO i_y         = kMin_y , kMax_y
    DO i_z         = kMin_z , kMax_z
      IF ( v_x( i_x, i_y, i_z ) .NE. v_x( i_x, i_y, i_z ) ) THEN
        NaN_count = NaN_count + 1
      END IF
    END DO
    END DO

    i_x            = kMax_x
    DO i_y         = kMin_y , kMax_y
    DO i_z         = kMin_z , kMax_z
      IF ( v_x( i_x, i_y, i_z ) .NE. v_x( i_x, i_y, i_z ) ) THEN
        NaN_count = NaN_count + 1
      END IF
    END DO
    END DO

    IF (NaN_count .NE. 0) THEN

      debug_error = 1
      ! This will jump out of evolution loop, if caught during that.

      CALL print_error_nan
      ! This will prompt error when checked
      ! REF-> <<< system_output >>>

    END IF

  END

END MODULE system_basicfunctions
