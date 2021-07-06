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
! MODULE: system_basicvariables
! LAST MODIFIED: 2 June 2021
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! BASIC VARIABLES AND ARRAYS FOR 3D EULER EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_basicvariables
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the global variables and arrays for the simulation space are declared and given values
! here, wheras temporary variables (IF necessary) are declared within the SUBROUTINEs
! Further, each variable is classified based on where its purpose suits apt.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_constants
  USE system_auxilaries
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  IMPLICIT NONE
  ! _________________________
  ! REAL SPACE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND =4) ::N_x,N_y,N_z
  INTEGER(KIND =4) ::N_min,N_max
  INTEGER(KIND =4) ::kMax_x,kMax_y,kMax_z
  INTEGER(KIND =4) ::kMin_x,kMin_y,kMin_z
  INTEGER(KIND =4) ::k_G,k_no
  INTEGER(KIND =4) ::kTru_x,kTru_y,kTru_z
  INTEGER(KIND =4) ::i_x,i_y,i_z
  INTEGER(KIND =4) ::j_x,j_y,j_z
  INTEGER(KIND =4) ::max_shell_no
  INTEGER(KIND =4) ::max_wave_no
  INTEGER(KIND =4) ::tot_active_modes
  INTEGER(KIND =4) ::tot_modes
  INTEGER(KIND =4),DIMENSION(3)::N_ar
  DOUBLE PRECISION ::N3,vol,dxdydz
  DOUBLE PRECISION ::k_G_2
  DOUBLE PRECISION ::trunc_const
  DOUBLE PRECISION ::L_x,L_y,L_z
  DOUBLE PRECISION ::dx,dy,dz
  DOUBLE PRECISION ::K_scale_x,K_scale_y,K_scale_z
  ! _________________________
  ! TIME VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION ::dt,dt_max
  DOUBLE PRECISION ::time_total,time_now,time_grid
  DOUBLE PRECISION ::time_save
  INTEGER (KIND=4) ::t_step,t_step_total,t_step_save,no_of_saves
  INTEGER (KIND=4) ::t_step_PVD_save,no_of_PVD_saves
  INTEGER (KIND=4) ::CFL_min,CFL_system
  ! _________________________
  ! FLUID VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION ::initial_circ
  INTEGER(KIND=4)  ::k_integral
  ! _________________________
  ! FUNCTION VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION ::energy
  DOUBLE PRECISION ::energy_initial,v_rms_1D
  DOUBLE PRECISION ::k_dot_v_norm
  DOUBLE PRECISION ::enstrophy
  DOUBLE PRECISION ::norm_factor
  DOUBLE PRECISION ::energy_mode
  ! _________________________
  ! CHARACTERS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  CHARACTER( LEN = 20) :: N_char
  CHARACTER( LEN = 20) :: IC_type
  CHARACTER( LEN = 3 ) :: run_code
  CHARACTER( LEN = 3 ) :: test_code
  CHARACTER( LEN = 3 ) :: solver_type
  CHARACTER( LEN = 3 ) :: solver_alg
  ! _________________________
  ! DEBUGGING VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  ::state_sim,check_status,debug_error
  INTEGER(KIND=4)  ::NaN_count,k_dot_v_error
  INTEGER(KIND=4)  ::no_of_debug,t_step_debug
  ! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::axis_x,axis_y,axis_z
  ! Arrays that store the grid
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::u_x,u_y,u_z
  ! Real velocity matrix  (updated after every time step)
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::w_ux,w_uy,w_uz
  ! Real vorticity
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::proj_xx,proj_yy,proj_zz
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::proj_xy,proj_yz,proj_zx
  ! projection operators matrix \mathbb{P}_{ij}=\delta_{ij}-\frac{k_ik_j}{k^2}}
  ! _________________________________________
  ! FOURIER SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4),DIMENSION(:,:,:),ALLOCATABLE  ::shell_no
  ! Every grid has its modulus, |k|
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::k_x,k_y,k_z,k_2,truncator
  ! wavenumber,truncator matrix
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   ::v_x,v_y,v_z
  ! Spectral velocity matrix (will be updated after every time step)
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   ::w_vx,w_vy,w_vz
  ! Spectral vorticity
  ! _________________________________________
  ! FOURIER (SHELL AVG) - SPECTRAL ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::spectral_energy,spectral_energy_avg
  ! Spectral data Integrated over spherical shell of radius k_mod
  INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE      ::count_modes_shell
  ! This counts no of modes that have (k-1)<|\vec{k}|<k

  !TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

  CONTAINS

  SUBROUTINE read_input
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Read simulation parameters from a file 'input_file'
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER( LEN = 60 )::input_file

    input_file  = 'parameters.dat'
    ! This file contains all major Input parameters to be fed from outside file

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! R  E  A  D  I  N  G       I  N  P  U  T       F  I  L  E
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    OPEN( UNIT = 1001, FILE = TRIM( ADJUSTL(input_file) ) )

    READ( 1001, f_i4,  ADVANCE ='yes')
    READ( 1001, f_i4,  ADVANCE ='yes')  N_x
    READ( 1001, f_i4,  ADVANCE ='yes')  N_y
    READ( 1001, f_i4,  ADVANCE ='yes')  N_z
    ! Resolution

    ! Total time to simulate
    READ( 1001, f_d8p4,  ADVANCE ='yes')
    READ( 1001, f_d8p4,  ADVANCE ='yes')  time_total
    ! Total time to simulate

    READ( 1001, f_d8p4, ADVANCE ='yes')
    READ( 1001, f_i4,   ADVANCE ='yes')  no_of_saves
    ! No of saves

    READ( 1001, f_d8p4, ADVANCE ='yes')
    READ( 1001, f_i4,   ADVANCE ='yes')  no_of_PVD_saves
    ! No of PVD saves

    CLOSE(1001)
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  END

  SUBROUTINE init_global_variables
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       Initialize all the variables used in the code and
  ! are explained too. Only variables are initialized here.
  ! Arrays are initialized in  another SUBROUTINE.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER( LEN = 10 ):: x_char,y_char,z_char
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! S P A C E    A N D     T I M E
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    N_ar            = (/ N_x, N_y, N_z /)
    N_min           = MINVAL( N_ar )
    N_max           = MAXVAL( N_ar )
    N3              = N_x * N_y * N_z
    ! Resolution of the cuboid

  	L_x             = ( N_x / N_min ) * two_pi
  	L_y             = ( N_y / N_min ) * two_pi
  	L_z             = ( N_z / N_min ) * two_pi
    ! Length of periodic cuboid

  	dx              = L_x / N_x
  	dy              = L_y / N_y
  	dz              = L_z / N_z
    ! Grid spacing - made to be equal in all direction

  	K_scale_x       = two_pi / L_x
  	K_scale_y       = two_pi / L_y
  	K_scale_z       = two_pi / L_z
    ! Scaling factor to be multiplied to the index to get the actual wavenumber

  	norm_factor     = DBLE( N_x * N_y * N_z )
    ! No of poINts IN real space

  	! ---------------------------
  	! MAX & MIN WAVENUMBERS ALLOWED
  	! ---------------------------
  	kMax_x          = + INT( N_x / 2 )
  	kMin_x          = + 0
  	kMax_y          = + INT( N_y / 2 ) - 1
  	kMin_y          = - INT( N_y / 2 )
  	kMax_z          = + INT( N_z / 2 ) - 1
  	kMin_z          = - INT( N_z / 2 )
    ! Maximum and Minimum wavemubers in the First Brillouin zone.

    WRITE( x_char, f_i4 ) N_x
    WRITE( y_char, f_i4 ) N_y
    WRITE( z_char, f_i4 ) N_z
    WRITE( N_char, "(A12)" ) TRIM( ADJUSTL( 'L' ) ) // TRIM( ADJUSTL ( x_char ) ) // &
     TRIM( ADJUSTL( 'W' ) ) // TRIM( ADJUSTL ( y_char ) ) // TRIM( ADJUSTL( 'H' ) ) // TRIM( ADJUSTL ( z_char ) )
    ! Converting resolution value to character

    dxdydz          = dx * dy * dz
    ! Grid volume

    vol             = L_x * L_y * L_z
    ! Volume of domain

    trunc_const     = thr ! Simply N/3 truncation const.
    ! Ratio of N to the truncation wavenumber

    k_G             = FLOOR( DBLE( N_min ) / trunc_const ) - 1
    k_G_2           = DBLE( k_G * k_G )
    ! Truncation wavenumber shell radius. For |k|>k_G, modes are truncated to remove aliasing error.

    kTru_x          = CEILING( DBLE( k_G ) / K_scale_x )
    kTru_y          = CEILING( DBLE( k_G ) / K_scale_y )
    kTru_z          = CEILING( DBLE( k_G ) / K_scale_z )
    ! Appropriate index limits along each wavenumber direction , beyond which it will be definitely truncated.

    max_shell_no    = CEILING( DSQRT( kMax_x ** two + kMin_y ** two + kMin_z ** two ) )
    ! Maximum '|k|' that can possibly reach. Even in the untruncated case.

    max_wave_no     = k_G + 1
    ! size of cube enclosing the truncation sphere. This is the size of sphere that encloses k_G cube.

    norm_factor     = one
    ! Normalization factor for energy - later changed so that initial energy is obtained.

    energy_initial  = one
    ! Initial energy of the system

    v_rms_1D        = DSQRT( two * energy_initial / thr )
    ! RMS Velocity

    time_grid       = dx / v_rms_1D
    ! Time scale for particle to cross a grid

    CFL_min         = 5
    ! - Courant-Friedrichs-Lewy (CFL) condition - CFL no is inverse of the above ratio
    ! No of steps (minimum) that should take to cross a grid

    dt_max          = time_grid / DBLE( CFL_min )
    ! Maximum value of time step

    CALL find_CFL_timestep( time_grid, dt_max, dt )
    ! Finds a time smaller than 'dt_max' in terms of '0.0..0p' , where p being '1' or '5'
    ! REF-> <<< system_auxilaries >>>

    CFL_system      = FLOOR( time_grid / dt )

    CALL time_to_step_convert( time_total, t_step_total, dt )
    ! returns the no of time_steps (\delta t) in a given time
    ! REF-> <<< system_auxilaries >>>

    t_step_save     = t_step_total / no_of_saves
    ! Determines how many time steps after the save has to be made.

    t_step_PVD_save = t_step_total / no_of_PVD_saves
    ! Determines how many time steps after the PVD save has to be made.

    CALL step_to_time_convert( t_step_save, time_save, dt )
    ! Determines the saving time intervals
    ! REF-> <<< system_auxilaries >>>

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! A U X I L A R Y
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    no_of_debug     = 2
    ! No of times that the program looks for any 'NaN' while marching forward in time.

    t_step_debug    = t_step_total / no_of_debug
    ! Determines how many time steps after the checking has to be made

    state_sim       = 0
    ! MeanINg it is INitializINg , '1' means fINal, will be changed after the time_marchINg is DOne.

    check_status    = 0
    ! Status that is checked before the time evolution

    NaN_count       = 0
    ! No of NaN IN v_x

    k_dot_v_error   = 0
    ! Incompressibility error

    debug_error     = 0
    ! Error found in debug, either NaN_count or k_dot_v_error
	END

  SUBROUTINE init_global_arrays
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !   This defines few arrays that do not evolve. k^2 matrix,k matrix, projection matrix,
  ! shell no matrix, count_modes_shell matrix.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::diff_ceiling,k_mod,k_mod_2
    DOUBLE PRECISION::kx,ky,kz

    CALL allocate_operators
    ! Allocates the arrays declared here.

    shell_no          = 0
    count_modes_shell = 0
    tot_modes         = 0
    tot_active_modes  = 0

    !  +++++++++++++++++++++++++++++++++
    !  A  X  I  S      G   R  I   D   S
    !  +++++++++++++++++++++++++++++++++
    DO i_x          = 0, N_x-1

      axis_x( i_x ) = DBLE( i_x ) * dx

    END DO
    DO i_y          = 0, N_y-1

      axis_y( i_y ) = DBLE( i_y ) * dy

    END DO
    DO i_z          = 0, N_z-1

      axis_z( i_z ) = DBLE( i_z ) * dz

    END DO

    DO j_x = kMin_x, kMax_x
  	DO j_y = kMin_y, kMax_y
  	DO j_z = kMin_z, kMax_z

      kx                    = K_scale_x * DBLE( j_x )
      ky                    = K_scale_y * DBLE( j_y )
      kz                    = K_scale_z * DBLE( j_z )

      k_x( j_x, j_y, j_z )  = kx
      k_y( j_x, j_y, j_z )  = ky
      k_z( j_x, j_y, j_z )  = kz
      ! Just the k component matrix storing its grid points.

      k_mod_2               = kx ** two + ky ** two + kz ** two
      k_2( j_x, j_y, j_z )  = k_mod_2
      k_mod                 = DSQRT( k_mod_2 )
      ! Square of distance to origin

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  P  R  O  J  E  C  T  I  O  N             M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Projection matrix \mathbb{P}_{ij}=\delta_{ij}-\frac{k_ik_j}{k^2}}

      IF ( k_mod_2 .GT. tol ) THEN
        ! Checking IF k^2 is not too low, to cause NaN (this will happen only for (0,0,0))
        proj_xx( j_x, j_y, j_z ) = one - ( kx * kx ) / k_mod_2
        proj_yy( j_x, j_y, j_z ) = one - ( ky * ky ) / k_mod_2
        proj_zz( j_x, j_y, j_z ) = one - ( kz * kz ) / k_mod_2
        proj_xy( j_x, j_y, j_z ) = - ( kx * ky ) / k_mod_2
        proj_yz( j_x, j_y, j_z ) = - ( ky * kz ) / k_mod_2
        proj_zx( j_x, j_y, j_z ) = - ( kz * kx ) / k_mod_2
      ELSE
        proj_xx( j_x, j_y, j_z ) = + twothird
        proj_yy( j_x, j_y, j_z ) = + twothird
        proj_zz( j_x, j_y, j_z ) = + twothird
        proj_xy( j_x, j_y, j_z ) = - onethird
        proj_yz( j_x, j_y, j_z ) = - onethird
        proj_zx( j_x, j_y, j_z ) = - onethird
      END IF

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  T  R  U  N  C  A  T  I  O  N  ,   S  H  E  L  L           M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Truncation mask matrix (multiply this with any spectral matrix to DO the truncation)
      IF ( k_mod_2 .LT. k_G_2 ) THEN
        truncator( j_x, j_y, j_z )                     =  one
        ! Spherical truncation filter matrix

        tot_active_modes                               =  tot_active_modes + 1
      ELSE

        truncator( j_x, j_y, j_z )                     =  zero
        ! Outised the truncation sphere.

      END IF

      diff_ceiling                                   =  DBLE( CEILING( k_mod ) ) - k_mod
      ! figuring the decimal part of the |k|, so that it will be alloted to shell k or k+1

      IF (diff_ceiling .GE. hf) THEN
        shell_no( j_x, j_y, j_z )                    =  FLOOR( k_mod )
      ELSE
        shell_no( j_x, j_y, j_z )                    =  CEILING( k_mod )
      END IF

      count_modes_shell( shell_no( j_x, j_y, j_z ) ) =  count_modes_shell( shell_no( j_x, j_y, j_z ) ) + 1
      ! counts no of grid poINts that belong to a particular shell, it should go as ~s^2

    END DO
    END DO
    END DO

    k_2( 0, 0, 0 ) = one
    ! Just to make sure , when something is divided by k_2 , to avoid NaN. Numerator would be zero anyways,
    ! in most of such cases. So nothing wrong here.

    ! Total no of active modes inside the truncation sphere.

    tot_modes = SUM( count_modes_shell )
    ! Total no of modes present

  END

  SUBROUTINE allocate_operators
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate arrays which are constants basically, k2,truncator, shell no etc.,
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y         A  L  L  O  C  A  T  I  O  N .
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( axis_x( 0 : N_x - 1 ) )
    ALLOCATE( axis_y( 0 : N_y - 1 ) )
    ALLOCATE( axis_z( 0 : N_z - 1 ) )
    ALLOCATE( k_2( kMin_x       : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( truncator( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( k_x( kMin_x       : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( k_y( kMin_x       : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( k_z( kMin_x       : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( proj_xx( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( proj_yy( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( proj_zz( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( proj_xy( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( proj_yz( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( proj_zx( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( shell_no( kMin_x  : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( count_modes_shell( 0 : max_shell_no ) )

  END

  SUBROUTINE allocate_velocity
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays related to velocity(real and spectral) and it spectrum
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N    -   V  E  L  O  C  I  T  Y
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( u_x( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( u_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( u_z( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( v_x( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( v_y( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( v_z( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( spectral_energy( 0     : max_wave_no ) )
    ALLOCATE( spectral_energy_avg( 0 : max_wave_no ) )

	END

  SUBROUTINE allocate_vorticity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate arrays related to vorticity(real and spectral)
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N    -   V  O  R  T  I  C  I  T  Y
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( w_ux( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( w_uy( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( w_uz( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( w_vx( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( w_vy( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( w_vz( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )

  END

	SUBROUTINE deallocate_velocity
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE( v_x, v_y, v_z )
		DEALLOCATE( u_x, u_y, u_z )
		DEALLOCATE( spectral_energy)
		DEALLOCATE( spectral_energy_avg)

	END

	SUBROUTINE deallocate_vorticity
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE( w_ux, w_uy, w_uz )
		DEALLOCATE( w_vx, w_vy, w_vz )

	END

  SUBROUTINE deallocate_operators
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE( axis_x, axis_y, axis_z )
    DEALLOCATE( k_2, k_x, k_y, k_z )
    DEALLOCATE( truncator )
    DEALLOCATE( proj_xx, proj_yy, proj_zz )
    DEALLOCATE( proj_xy, proj_yz, proj_zx )
    DEALLOCATE( shell_no, count_modes_shell )

	END

END MODULE system_basicvariables
