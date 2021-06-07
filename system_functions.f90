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
! MODULE: system_functions
! LAST MODIFIED: 3 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! RUN AND OUTPUT MODULE TO RUN THE TIME EVOLUTION FOR 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_functions
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major system_functions for the Code to run.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_initialcondition
  USE system_variables
  USE system_fftw

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
    ! Calls the subroutine to get a initial condition with norm_factor=1

    CALL compute_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_initial / energy )
    ! Normalizing the norm_factor, so that we get energy='en_initial'

    CALL init_initcondn
    ! Calls it again to get normalized velocities.

    CALL compute_spectral_data
    ! Gets the energy from spectral space

    CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
    ! FFT spectral to real velocity

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
    ! Reset the variable

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y     S  P  E  C  T  R  U  M     C  A   L   C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Gets the energy in that particular mode (i_x,i_y,i_z) into 'en_temp'
    ! Keeps adding energy to that particular shell (|k| fixed), from all solid angles
    DO i_x =    1 , k_G
    DO i_y = - k_G, k_G
    DO i_z = - k_G, k_G
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
    DO i_y = - k_G, k_G
    DO i_z = - k_G, -1
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
    DO i_y = 0, k_G
      energy_mode                                   = CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                                      CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                                      CDABS( v_z( i_x, i_y, i_z ) ) ** two
      spectral_energy( shell_no ( i_x, i_y, i_z ) ) = spectral_energy( shell_no ( i_x, i_y, i_z ) ) + energy_mode
      enstrophy                                     = enstrophy + k_2( i_x, i_y, i_z) * energy_mode
    END DO

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  S  H  E  L  L      A  V  E  R  A  G  I  N  G
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    spectral_energy_avg( 1 )            = qtr * ( thr * spectral_energy( 1 ) + spectral_energy( 2 ) )
    spectral_energy_avg( max_shell_no ) = qtr * ( thr * spectral_energy( max_shell_no ) + &
                                                        spectral_energy( max_shell_no - 1 ) )
    DO k_no                             = 2, max_shell_no - 1
        spectral_energy_avg( k_no )     = qtr * ( spectral_energy( k_no - 1 ) + spectral_energy( k_no + 1 ) ) + &
                                           hf * ( spectral_energy( k_no ) )
    END DO

    energy = SUM( spectral_energy( : ) )

  END

  SUBROUTINE compute_energy_spectral_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check the presence of NaN in your spectral velocity data (v(k)),
  ! and also the L2 norm or the Kinetic energy.
  ! NOTE: Count certain modes once, certain modes half (owing to 1/2 factor)
  ! in the first loop i_x=0 plane is left. later it is considered
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    energy      = zero

    DO i_x      =  1, Nh - 1
    DO i_y      = -Nh, Nh - 1
    DO i_z      = -Nh, Nh - 1
      energy    = energy + CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                           CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                           CDABS( v_z( i_x, i_y, i_z ) ) ** two

    IF ( v_x( i_x, i_y, i_z ) .NE. v_x( i_x, i_y, i_z ) ) THEN
      NaN_count = NaN_count + 1
    END IF

    END DO
    END DO
    END DO

    i_x         =   0
    DO i_y      = - Nh, Nh - 1
    DO i_z      = - Nh, Nh - 1
      energy    = energy + hf * ( CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_z( i_x, i_y, i_z ) ) ** two )
    END DO
    END DO

    i_x         =   Nh
    DO i_y      = - Nh, Nh - 1
    DO i_z      = - Nh, Nh - 1
      energy    = energy + hf * ( CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_z( i_x, i_y, i_z ) ) ** two )
    END DO
    END DO

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

  SUBROUTINE compute_compressibility
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check incompressibility condition. Sums over all residues
  ! of incompressibility and prints it. Of order 10^(-12).
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    k_dot_v_norm   = zero

    DO i_x         = 0, Nh
    DO i_y         = -Nh, Nh - 1
    DO i_z         = -Nh, Nh - 1
      k_dot_v_norm = k_dot_v_norm + CDABS( i_x * v_x( i_x, i_y, i_z ) + &
                                           i_y * v_y( i_x, i_y, i_z ) + &
                                           i_z * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO
    END DO

    i_x            = 0
    DO i_y         = - Nh, Nh - 1
    DO i_z         = - Nh, Nh - 1
      k_dot_v_norm = k_dot_v_norm + hf* CDABS( i_x * v_x( i_x, i_y, i_z ) + &
                                               i_y * v_y( i_x, i_y, i_z ) + &
                                               i_z * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO

    i_x            = Nh
    DO i_y         = - Nh, Nh - 1
    DO i_z         = - Nh, Nh - 1
      k_dot_v_norm = k_dot_v_norm + hf* CDABS( i_x * v_x( i_x, i_y, i_z ) + &
                                               i_y * v_y( i_x, i_y, i_z ) + &
                                               i_z * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO

    k_dot_v_norm = DSQRT( k_dot_v_norm )

    IF (k_dot_v_norm .GT. tol_float ) THEN

      k_dot_v_error = 1

    END IF
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

    CALL fft_c2r( w_vx, w_vy, w_vz, N, Nh, w_ux, w_uy, w_uz )
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


END MODULE system_functions
