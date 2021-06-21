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
! MODULE: system_initialcondition
! LAST MODIFIED: 21 June 2021
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! INITIAL CONDITION FOR 3D EULER EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_initialcondition
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Different initial conditions are coded here
! 1. Exponentially decaying spectrum
! 2. Thermalized spectrum
! 3. Kolmogorov like spectrum
! 4. Taylor-Green Initial Condition
! 5. Kida Peltz Initial Condition
! 6. ABC Flow initial condition
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicvariables
  USE system_fftw
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  IMPLICIT  NONE

  CONTAINS

  SUBROUTINE init_initcondn
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Initialize initial condition for velocity, Choose one.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! Initializing the initial velocity (spectral) and projecting it so that the flow is incompressible.

    CALL IC_exp_decaying_spectrum(energy_initial)
    ! Generic randomized initial condition, with energy mainly in integral scale (spectrally)

    ! CALL IC_Kolmogorov_spectrum(energy_initial)
    ! Generic initial condition, with energy mainly in inertial range with a k^-(5/3) spectrum.

    ! CALL IC_perfect_thermalized_spectrum(energy_initial)
    ! Create its own thermalized spectrum by equiparition, (no permanence of large eddies in this case)

    ! CALL IC_TG(energy_initial)
    ! TAYLOR_GREEN Initial condition - lots of symmetries (although solver is not using them)

    ! CALL IC_KP(energy_initial)
    ! KIDA-PELTZ Initial condition - lots of symmetries

    ! CALL IC_ABC(energy_initial)
    ! Arnold-Beltrami-Childress Initial condition

    ! CALL IC_vortex_sheet(energy_initial)
    ! Creates a vortex sheets at z = +pi/2, -pi/2, pointing along y direction.
    ! With a background field from IC_exp_decaying_spectrum

    ! CALL IC_vortex_tube(energy_initial)
    ! Creates a vortex tube at z = 0, along z direction.
    ! With a background field from IC_exp_decaying_spectrum

    ! CALL IC_from_file_spectral
    ! Read from file.
    ! *****Check whether file is available already.

    ! CALL IC_from_file_real
    ! Read from file.
    ! *****Check whether file is available already.

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! UNCOMMENT TO TRUNCATE IF NEEDED - (most of I.C are already truncated)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! v_x = truncator * v_x
    ! v_y = truncator * v_y
    ! v_z = truncator * v_z

  END

  SUBROUTINE IC_exp_decaying_spectrum(energy_input)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! A typical initial condition IN which only few large modes are excited.
  ! CALL this to give 3 components of COMPLEX spectral velocity at spectral grid 'i_x,i_y,i_z'
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)  ::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION             ::phi,theta
    DOUBLE PRECISION,DIMENSION(3)::ph
    DOUBLE PRECISION             ::V_k_mod,norm_const,k_ratio
    DOUBLE COMPLEX,DIMENSION(3)  ::V_k
    INTEGER(KIND=4)              ::integral_exponent

    IC_type = 'EXP-DEC'

    CALL init_random_seed
    ! Randomizes seed for random numbers (in 'auxilary_functions' module )
    ! REF-> <<< system_auxilaries >>>

    k_integral        = 2
    ! Integral scale wavenumber

    integral_exponent = 2
    ! The power in the spectrum E(k) for k<k_integral
    ! Generally either 2 or 4 or 6. But has to be a even number

    CALL normalization_exponential_spectrum( integral_exponent, k_integral, norm_const)
    ! Returns the 'norm_const', so that theoretically net energy is O(1).
    ! Additionally normalization to any energy can be done with norm_factor
    ! REF-> <<< system_auxilaries >>>

    DO i_x = 0, Nh
    DO i_y = -Nh, Nh - 1
    DO i_z = -Nh, Nh - 1
    IF ( k_2( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(theta)
      CALL RANDOM_NUMBER(ph)

      phi     = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      theta   = DACOS( one - two * theta )! Polar angle of \hat{u}_k vector
      ph      = two_pi * ph ! Phases of \hat{u}_k components

      k_ratio = DSQRT( k_2( i_x, i_y, i_z) ) / DBLE(k_integral)

      V_k_mod = norm_factor * norm_const * k_ratio**( hf * integral_exponent - 1 ) &
                * DEXP( - qtr * integral_exponent * ( k_ratio ** two ) )

      V_k(1)  = V_k_mod * DSIN( theta ) * DCOS( phi ) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
      V_k(2)  = V_k_mod * DSIN( theta ) * DSIN( phi ) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
      V_k(3)  = V_k_mod * DCOS( theta ) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
      ! 3 COMPLEX values for spectral velocity

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  V  E  L  O  C  I  T  Y          P  R  O  J  E  C  T  O  R
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Assigning the velocity along with projection so that it is Incompressible -- u(k).k=0
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  N  O  T  E  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! The following steps ensure, that except for i_x=0 or Nh plane, all the other planes are given Initial velocity
      ! But those planes require special attention, becoz, their Inversion lies IN the same plane,
      ! so conjugates must be placed accordingly.

      IF (((i_x .NE. 0) .AND. (i_x .NE. Nh)) .OR. (i_z .LT. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z ) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z ) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z ) * V_k(3)

      IF (((i_x .EQ. 0) .OR. (i_x .EQ. Nh)) .AND. ((i_y .NE. -Nh) .AND. (i_z .NE. -Nh))) THEN
        v_x( i_x, - i_y, - i_z ) = DCONJG( v_x( i_x, i_y, i_z ) )
        v_y( i_x, - i_y, - i_z ) = DCONJG( v_y( i_x, i_y, i_z ) )
        v_z( i_x, - i_y, - i_z ) = DCONJG( v_z( i_x, i_y, i_z ) )
      END IF

      ELSE IF ((i_z .EQ. 0) .AND. (i_y .LE. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z) * V_k(3)

        v_x( i_x, - i_y, i_z )   = DCONJG( v_x ( i_x, i_y, i_z ) )
        v_y( i_x, - i_y, i_z )   = DCONJG( v_y ( i_x, i_y, i_z ) )
        v_z( i_x, - i_y, i_z )   = DCONJG( v_z ( i_x, i_y, i_z ) )
      ELSE IF ((i_y .EQ. -Nh) .AND. (i_z .GT. 0)) THEN
        v_x(i_x,i_y,i_z)         = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z) * V_k(3)
        v_y(i_x,i_y,i_z)         = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z) * V_k(3)
        v_z(i_x,i_y,i_z)         = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z) * V_k(3)
      END IF
      ! ----------------------------------------------------------------------------------------

    END IF
    END DO
    END DO
    END DO

    ! Making sure, that the average velocity is zero.
    v_x( 0, 0, 0 )    =     zero
    v_y( 0, 0, 0 )    =     zero
    v_z( 0, 0, 0 )    =     zero

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

  END

  SUBROUTINE IC_Kolmogorov_spectrum(energy_input)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! A typical initial condition with kolmogorov spectrum model, referred from Pope's Turbulence.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)  ::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION             ::phi,theta
    DOUBLE PRECISION,DIMENSION(3)::ph
    DOUBLE PRECISION             ::V_k_mod,norm_const,k_ratio
    DOUBLE PRECISION             ::factor_integral,factor_dissipation
    DOUBLE PRECISION             ::eleven_by_twelve
    DOUBLE COMPLEX,DIMENSION(3)  ::V_k
    INTEGER(KIND=4)              ::integral_exponent
    INTEGER(KIND=4)              ::k_kol

    IC_type = 'KOL-INE'

    CALL init_random_seed
    ! Randomizes seed for random numbers (in 'auxilary_functions' module )
    ! REF-> <<< system_auxilaries >>>

    k_integral        = FLOOR( DLOG( DBLE( N ) ) / DLOG( 4.0D0 ) ) - 1
    ! Integral scale wavenumber

    k_kol             = FLOOR ( DBLE( N ) / 4.0D0 ) - 1
    ! End of kolmogorov spectrum

    integral_exponent = 2
    ! The power in the spectrum E(k) for k<k_integral
    ! Generally either 2 or 4 or 6. But has to be a even number

    eleven_by_twelve  = 11.0D0 /  12.0D0

    CALL normalization_kolmogorov_spectrum( k_integral, k_kol, norm_const)
    ! Returns the 'norm_const', so that theoretically net energy is O(1).
    ! Additionally normalization to any energy can be done with norm_factor
    ! REF-> <<< system_auxilaries >>>

    DO i_x = 0, Nh
    DO i_y = -Nh, Nh - 1
    DO i_z = -Nh, Nh - 1
    IF ( k_2( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(theta)
      CALL RANDOM_NUMBER(ph)

      phi                = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      theta              = DACOS( one - two * theta )! Polar angle of \hat{u}_k vector
      ph                 = two_pi * ph ! Phases of \hat{u}_k components

      k_ratio            = DSQRT( k_2( i_x, i_y, i_z) ) / DBLE(k_integral)
      factor_integral    = one
      factor_dissipation = one

      IF ( k_ratio .LT. 1 ) THEN
        CALL kolmogorov_spectrum_integralscale_subpart(k_ratio,integral_exponent,factor_integral)
        ! REF-> <<< system_auxilaries >>>
      END IF

      k_ratio            = DSQRT( k_2( i_x, i_y, i_z) ) / DBLE(k_kol)

      IF ( k_ratio .GT. qtr ) THEN
        CALL kolmogorov_spectrum_dissipationscale_subpart(k_ratio,factor_dissipation)
        ! REF-> <<< system_auxilaries >>>
      END IF

      V_k_mod            = norm_factor * norm_const * ( k_2( i_x, i_y, i_z ) ** ( - eleven_by_twelve ) ) &
                          * factor_integral * factor_dissipation

      V_k(1)             = V_k_mod * DSIN( theta ) * DCOS( phi ) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
      V_k(2)             = V_k_mod * DSIN( theta ) * DSIN( phi ) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
      V_k(3)             = V_k_mod * DCOS( theta ) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
      ! 3 COMPLEX values for spectral velocity

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  V  E  L  O  C  I  T  Y          P  R  O  J  E  C  T  O  R
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Assigning the velocity along with projection so that it is INcompressible -- u(k).k=0
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  N  O  T  E  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! The following steps ensure, that except for i_x=0 or Nh plane, all the other planes are given INitial velocity
      ! But those planes require special attention, becoz, their INversion lies IN the same plane,
      ! so conjugates must be placed accordINgly.

      IF (((i_x .NE. 0) .AND. (i_x .NE. Nh)) .OR. (i_z .LT. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z ) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z ) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z ) * V_k(3)

      IF (((i_x .EQ. 0) .OR. (i_x .EQ. Nh)) .AND. ((i_y .NE. -Nh) .AND. (i_z .NE. -Nh))) THEN
        v_x( i_x, - i_y, - i_z ) = DCONJG( v_x( i_x, i_y, i_z ) )
        v_y( i_x, - i_y, - i_z ) = DCONJG( v_y( i_x, i_y, i_z ) )
        v_z( i_x, - i_y, - i_z ) = DCONJG( v_z( i_x, i_y, i_z ) )
      END IF

      ELSE IF ((i_z .EQ. 0) .AND. (i_y .LE. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z) * V_k(3)

        v_x( i_x, - i_y, i_z )   = DCONJG( v_x ( i_x, i_y, i_z ) )
        v_y( i_x, - i_y, i_z )   = DCONJG( v_y ( i_x, i_y, i_z ) )
        v_z( i_x, - i_y, i_z )   = DCONJG( v_z ( i_x, i_y, i_z ) )
      ELSE IF ((i_y .EQ. -Nh) .AND. (i_z .GT. 0)) THEN
        v_x(i_x,i_y,i_z)         = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z) * V_k(3)
        v_y(i_x,i_y,i_z)         = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z) * V_k(3)
        v_z(i_x,i_y,i_z)         = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z) * V_k(3)
      END IF
      ! ----------------------------------------------------------------------------------------

    END IF
    END DO
    END DO
    END DO

    ! Making sure, that the average velocity is zero.
    v_x( 0, 0, 0 )    =     zero
    v_y( 0, 0, 0 )    =     zero
    v_z( 0, 0, 0 )    =     zero

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

  END

  SUBROUTINE IC_perfect_thermalized_spectrum(energy_input)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Initialize initial condition for thermalized spectrum
  ! where equipartition is given.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)  ::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::phi,theta,norm_const
    DOUBLE PRECISION,DIMENSION(3)::ph
    DOUBLE COMPLEX,DIMENSION(3)::V_k

    IC_type = 'THER-K2'

    CALL init_random_seed
    ! Randomizes seed for random numbers (in 'auxilary_functions' module )
    ! REF-> <<< system_auxilaries >>>

    norm_const = DSQRT( thr / ( two_pi * DBLE( k_G ** 3 ) ) )

    DO i_x = 0, Nh
    DO i_y = -Nh, Nh - 1
    DO i_z = -Nh, Nh - 1
    IF ( k_2( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(theta)
      CALL RANDOM_NUMBER(ph)

      phi     = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      theta   = DACOS( one - two * theta )! Polar angle of \hat{u}_k vector
      ph      = two_pi * ph ! Phases of \hat{u}_k components

      V_k(1)  = norm_factor * norm_const * DSIN( theta ) * DCOS( phi ) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
      V_k(2)  = norm_factor * norm_const * DSIN( theta ) * DSIN( phi ) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
      V_k(3)  = norm_factor * norm_const * DCOS( theta ) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
      ! 3 COMPLEX values for spectral velocity

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  V  E  L  O  C  I  T  Y          P  R  O  J  E  C  T  O  R
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Assigning the velocity along with projection so that it is INcompressible -- u(k).k=0
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  N  O  T  E  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! The following steps ensure, that except for i_x=0 or Nh plane, all the other planes are given INitial velocity
      ! But those planes require special attention, becoz, their INversion lies IN the same plane,
       ! so conjugates must be placed accordINgly.
      IF (((i_x .NE. 0) .AND. (i_x .NE. Nh)) .OR. (i_z .LT. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z ) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z ) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z ) * V_k(3)

        IF (((i_x .EQ. 0) .OR. (i_x .EQ. Nh)) .AND. ((i_y .NE. -Nh) .AND. (i_z .NE. -Nh))) THEN
          v_x( i_x, - i_y, - i_z ) = DCONJG( v_x( i_x, i_y, i_z ) )
          v_y( i_x, - i_y, - i_z ) = DCONJG( v_y( i_x, i_y, i_z ) )
          v_z( i_x, - i_y, - i_z ) = DCONJG( v_z( i_x, i_y, i_z ) )
        END IF
      ELSE IF ((i_z .EQ. 0) .AND. (i_y .LE. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z) * V_k(3)

        v_x( i_x, - i_y, i_z )   = DCONJG( v_x ( i_x, i_y, i_z ) )
        v_y( i_x, - i_y, i_z )   = DCONJG( v_y ( i_x, i_y, i_z ) )
        v_z( i_x, - i_y, i_z )   = DCONJG( v_z ( i_x, i_y, i_z ) )
      ELSE IF ((i_y .EQ. -Nh) .AND. (i_z .GT. 0)) THEN
        v_x(i_x,i_y,i_z)         = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z) * V_k(3)
        v_y(i_x,i_y,i_z)         = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z) * V_k(3)
        v_z(i_x,i_y,i_z)         = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z ) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z) * V_k(3)
      END IF
      ! ----------------------------------------------------------------------------------------

    END IF
    END DO
    END DO
    END DO

    ! Making sure, that the average velocity is zero.
    v_x( 0, 0, 0 )    =     zero
    v_y( 0, 0, 0 )    =     zero
    v_z( 0, 0, 0 )    =     zero

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

  END

  SUBROUTINE IC_TG(energy_input)
  ! This is TAYLOR_GREEN Vortex initial condition which has a lot of symmetries
    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX::v0

    IC_type = 'TAY-GRE'

    v0            = i / 8.0D0
    !-------- 'x' velocity--------------
    v_x(+1,+1,+1) = - v0
    v_x(+1,+1,-1) = - v0
    v_x(+1,-1,+1) = - v0
    v_x(+1,-1,-1) = - v0
    !-------- 'y' velocity--------------
    v_y(+1,+1,+1) = + v0
    v_y(+1,+1,-1) = + v0
    v_y(+1,-1,+1) = - v0
    v_y(+1,-1,-1) = - v0

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v0           = norm_factor * v0
    !-------- 'x' velocity--------------
    v_x(+1,+1,+1) = - v0
    v_x(+1,+1,-1) = - v0
    v_x(+1,-1,+1) = - v0
    v_x(+1,-1,-1) = - v0
    !-------- 'y' velocity--------------
    v_y(+1,+1,+1) = + v0
    v_y(+1,+1,-1) = + v0
    v_y(+1,-1,+1) = - v0
    v_y(+1,-1,-1) = - v0
  END

  SUBROUTINE IC_KP(energy_input)
  ! This is Kida Peltz Vortex INitial condition which has a lot of symmetries
    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX::v1

    IC_type = 'KID-PEL'

    v1           = - i / 8.0D0

    !-------- 'x' velocity--------------
    v_x(+1,+3,+1) = + v1
    v_x(+1,+3,-1) = + v1
    v_x(+1,-3,+1) = + v1
    v_x(+1,-3,-1) = + v1
    v_x(+1,+1,+3) = - v1
    v_x(+1,+1,-3) = - v1
    v_x(+1,-1,+3) = - v1
    v_x(+1,-1,-3) = - v1

    !-------- 'y' velocity--------------
    v_y(+1,+1,+3) = + v1
    v_y(+1,-1,-3) = - v1
    v_y(+1,+1,-3) = + v1
    v_y(+1,-1,+3) = - v1
    v_y(+3,+1,+1) = - v1
    v_y(+3,-1,-1) = + v1
    v_y(+3,+1,-1) = - v1
    v_y(+3,-1,+1) = + v1

    !-------- 'z' velocity--------------
    v_z(+3,+1,+1) = + v1
    v_z(+3,-1,+1) = + v1
    v_z(+3,-1,-1) = - v1
    v_z(+3,+1,-1) = - v1
    v_z(+1,+3,+1) = - v1
    v_z(+1,-3,+1) = - v1
    v_z(+1,-3,-1) = + v1
    v_z(+1,+3,-1) = + v1

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v1           = norm_factor * v1

    !-------- 'x' velocity--------------
    v_x(+1,+3,+1) = + v1
    v_x(+1,+3,-1) = + v1
    v_x(+1,-3,+1) = + v1
    v_x(+1,-3,-1) = + v1
    v_x(+1,+1,+3) = - v1
    v_x(+1,+1,-3) = - v1
    v_x(+1,-1,+3) = - v1
    v_x(+1,-1,-3) = - v1

    !-------- 'y' velocity--------------
    v_y(+1,+1,+3) = + v1
    v_y(+1,-1,-3) = - v1
    v_y(+1,+1,-3) = + v1
    v_y(+1,-1,+3) = - v1
    v_y(+3,+1,+1) = - v1
    v_y(+3,-1,-1) = + v1
    v_y(+3,+1,-1) = - v1
    v_y(+3,-1,+1) = + v1

    !-------- 'z' velocity--------------
    v_z(+3,+1,+1) = + v1
    v_z(+3,-1,+1) = + v1
    v_z(+3,-1,-1) = - v1
    v_z(+3,+1,-1) = - v1
    v_z(+1,+3,+1) = - v1
    v_z(+1,-3,+1) = - v1
    v_z(+1,-3,-1) = + v1
    v_z(+1,+3,-1) = + v1

  END

  SUBROUTINE IC_ABC(energy_input)
  ! Call this to implement Arnold Beltrami Childress Initial condition

    IMPLICIT NONE
    ! _________________________
    ! TRANSFER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX  ::v0
    DOUBLE PRECISION:: A,B,C

    IC_type = 'ABC-FLW'

    ! ----------------------
    ! ABC PARAMETERS
    ! ----------------------
    A = one
    B = one
    C = one
    ! ----------------------
    v0           = one
    !--------- 'x' velocity --------------
    v_x(0,0,1)  = - v0 * A * i
    v_x(0,1,0)  = + v0 * C
    !--------- 'y' velocity --------------
    v_y(1,0,0)  = - v0 * B * i
    v_y(0,0,1)  = + v0 * A
    !--------- 'z' velocity --------------
    v_z(1,0,0)  = + v0 * B
    v_z(0,1,0)  = - v0 * C * i

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v0           = norm_factor * v0
    !--------- 'x' velocity --------------
    v_x(0,0,1)  = - v0 * A * i
    v_x(0,1,0)  = + v0 * C
    !--------- 'y' velocity --------------
    v_y(1,0,0)  = - v0 * B * i
    v_y(0,0,1)  = + v0 * A
    !--------- 'z' velocity --------------
    v_z(1,0,0)  = + v0 * B
    v_z(0,1,0)  = - v0 * C * i

  END

  SUBROUTINE IC_vortex_sheet(energy_input)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! An initial condition with vortex sheet imposed with a background field.
  ! Ratio of energy split between sheet and background is adjustable
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::u0,smooth_pm
    DOUBLE PRECISION::energy_sheet,energy_ratio
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_sheet_x

    ALLOCATE(u_sheet_x(0:N-1,0:N-1,0:N-1))

    u0           = one
    ! Normalizing parameter

    smooth_pm    = 0.5D0
    ! How thick the sheet is, smaller the parameter thicker it is

    energy_ratio = 0.02D0
    ! Percentage of energy in Background field

    DO i_x = 0, N - 1
    DO i_y = 0, N - 1
    DO i_z = 0, Nh - 1

      u_sheet_x( i_x, i_y, i_z ) = u0 * DTANH( - smooth_pm * hf * DBLE( i_z - ( N / 4 ) ) )

    END DO
    DO i_z = Nh, N - 1

      u_sheet_x( i_x, i_y, i_z ) = u0 * DTANH( smooth_pm * hf * DBLE( i_z - 3 * ( N / 4 ) ) )

    END DO
    END DO
    END DO

    energy_sheet = hf * SUM( u_sheet_x ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_sheet )
    u_sheet_x    = u0 * u_sheet_x
    ! Normalization of sheet

    CALL IC_exp_decaying_spectrum( energy_ratio * energy_input )
    ! Gets a background flow for remaining energy

    CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
    ! Getting the real velocity to add the sheet

    u_x = u_x + u_sheet_x

    CALL fft_r2c_scalar( u_x, N, Nh, v_x )
    ! FFT spectral to real velocity

    IC_type = 'VOR-SHT'

    DEALLOCATE(u_sheet_x)

  END

  SUBROUTINE IC_vortex_tube(energy_input)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! An initial condition with vortex sheet imposed with a background field.
  ! Ratio of energy split between sheet and background is adjustable
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::u0,smooth_pm
    DOUBLE PRECISION::tube_x0,tube_y0
    DOUBLE PRECISION::tube_x,tube_y
    DOUBLE PRECISION::energy_tube,energy_ratio
    DOUBLE PRECISION::u_ang,radius
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_tube_x,u_tube_y

    ALLOCATE(u_tube_x(0:N-1,0:N-1,0:N-1))
    ALLOCATE(u_tube_y(0:N-1,0:N-1,0:N-1))

    u0           = one
    ! Normalizing parameter

    smooth_pm    = 0.5D0
    ! How thick the sheet is, smaller the parameter thicker it is

    energy_ratio = 0.02D0
    ! Percentage of energy in Background field

    tube_x0      = DBLE( Nh )
    tube_y0      = DBLE( Nh )
    ! Center of the tube

    DO i_x = 0, N - 1
    DO i_y = 0, N - 1

      tube_x  = DBLE( i_x ) - tube_x0
      tube_y  = DBLE( i_y ) - tube_y0
      radius  = DSQRT( tube_x ** two + tube_y ** two )
      u_ang   = u0 * smooth_pm * DEXP( - hf * ( smooth_pm * radius ) ** two )
      u_tube_x( i_x, i_y, : ) = - u_ang * tube_y
      u_tube_y( i_x, i_y, : ) = + u_ang * tube_x

    END DO
    END DO

    energy_tube  = hf * SUM( u_tube_x ** two + u_tube_y ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_tube )
    u_tube_x     = u0 * u_tube_x
    u_tube_y     = u0 * u_tube_y
    ! Normalization of tube

    CALL IC_exp_decaying_spectrum( energy_ratio * energy_input )
    ! Gets a background flow for remaining energy

    CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
    ! Getting the real velocity to add the tube velocity

    u_x = u_x + u_tube_x
    u_y = u_y + u_tube_y

    CALL fft_r2c_scalar( u_x, N, Nh, v_x )
    CALL fft_r2c_scalar( u_y, N, Nh, v_y )
    ! FFT spectral to real velocity

    CALL compute_projected_velocity
    ! Projects the velocity to remove some incompressibility

    IC_type = 'VOR-TUB'

    DEALLOCATE(u_tube_x,u_tube_y)

  END

  SUBROUTINE IC_from_file_spectral
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Read a file for the velocity data in spectral space.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::real_part,imag_part
    CHARACTER(LEN=80)::IC_file_name

    IC_type = 'INP-SPE'

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! V   E  L  O  C  I  T  Y       I  N  P  U  T     F  I  L  E
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IC_file_name  = 'spectral_velocity_' // TRIM( ADJUSTL( N_char ) ) // '_input.dat'
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 43, FILE = IC_file_name )

    ! '+++++++++++++++++++++'
    ! 'I.C FROM FILE'
    ! '+++++++++++++++++++++'

    DO i_x =    0, Nh
    DO i_y = - Nh, Nh - 1
    DO i_z = - Nh, Nh - 1

      READ(43,f_c32p17,ADVANCE='NO') real_part, imag_part
      v_x( i_x, i_y, i_z ) = DCMPLX( real_part, imag_part )
      READ(43,f_c32p17,ADVANCE='NO') real_part, imag_part
      v_y( i_x, i_y, i_z ) = DCMPLX( real_part, imag_part )
      READ(43,f_c32p17,ADVANCE='YES') real_part, imag_part
      v_z( i_x, i_y, i_z ) = DCMPLX( real_part, imag_part )

    END DO
    END DO
    END DO

    CLOSE(43)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE IC_from_file_real
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Read a file for the velocity data in real space.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN=80)::IC_file_name

    IC_type = 'INP-REA'

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! V   E  L  O  C  I  T  Y       I  N  P  U  T     F  I  L  E
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IC_file_name  = 'velocity_' // TRIM( ADJUSTL( N_char ) ) // '_input.dat'
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 44, FILE = IC_file_name )

    ! '+++++++++++++++++++++'
    ! 'I.C FROM FILE'
    ! '+++++++++++++++++++++'

    DO i_x = 0 , N - 1
    DO i_y = 0 , N - 1
    DO i_z = 0 , N - 1

      READ(44,f_d32p17,ADVANCE = 'NO')  u_x( i_x, i_y, i_z)
      READ(44,f_d32p17,ADVANCE = 'NO')  u_y( i_x, i_y, i_z)
      READ(44,f_d32p17,ADVANCE = 'YES') u_z( i_x, i_y, i_z)

    END DO
    END DO
    END DO

    CLOSE(44)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    CALL fft_r2c( u_x, u_y, u_z, N, Nh, v_x, v_y, v_z )
    ! FFT spectral to real velocity

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

  SUBROUTINE compute_projected_velocity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to project the spectral velocity, to make it incompressible
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   ::v_P_x,v_P_y,v_P_z
    ALLOCATE(v_P_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_P_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_P_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

    v_P_x = v_x
    v_P_y = v_y
    v_P_z = v_z

    v_x   = proj_xx * v_P_x + proj_xy * v_P_y + proj_zx * v_P_z
    v_y   = proj_xy * v_P_x + proj_yy * v_P_y + proj_yz * v_P_z
    v_z   = proj_zx * v_P_x + proj_yz * v_P_y + proj_zz * v_P_z

    DEALLOCATE(v_P_x,v_P_y,v_P_z)

  END

END MODULE system_initialcondition
