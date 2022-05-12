! <f
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
! LAST MODIFIED: 22 JAN 2022
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! INITIAL CONDITION FOR 3D EULER EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_initialcondition
! <f
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Different initial conditions are coded here
! 1. Exponentially decaying spectrum
! 2. Thermalized spectrum
! 3. Kolmogorov like spectrum
! 4. Vortex sheets
! 5. Vortex sheet with Taylor Green flow
! 6. Vortex sheet with ABC Flow
! 7. Vortex sheet with Vortex tube
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! </f>
! <f
  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicvariables
  USE system_fftw_adv
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  IMPLICIT  NONE

! </f>
  CONTAINS

  SUBROUTINE init_initcondn
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Initialize initial condition for velocity, Choose one.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    CALL init_fft_size( N_x, N_y, N_z )
    ! Initializing the size of domain for FFT's in simulation - One time procedure.
    ! REF-> <<< system_fftw_adv >>>

    ! Initializing the initial velocity (spectral) and projecting it so that the flow is incompressible.

    ! CALL IC_exp_decaying_spectrum(energy_initial)
    ! Generic randomized initial condition, with energy mainly in integral scale (spectrally)

    ! CALL IC_Kolmogorov_spectrum(energy_initial)
    ! Generic initial condition, with energy mainly in inertial range with a k^-(5/3) spectrum.

    ! CALL IC_perfect_thermalized_spectrum(energy_initial)
    ! Create its own thermalized spectrum by equiparition, (no permanence of large eddies in this case)

    CALL IC_vortex_sheet(energy_initial)
    ! Creates a vortex sheets at z = +pi/2, -pi/2, pointing along y direction.
    ! With a background field from IC_exp_decaying_spectrum

    ! CALL IC_vortex_sheet_with_TG(energy_initial)
    ! Creates a vortex sheets at z = +pi/2, -pi/2, pointing along y direction.
    ! With a background field from Taylor Green

    ! CALL IC_vortex_sheet_with_ABC(energy_initial)
    ! Creates a vortex sheets at z = +pi/2, -pi/2, pointing along y direction.
    ! With a background field from ABC flow

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
    v_x = truncator * v_x
    v_y = truncator * v_y
    v_z = truncator * v_z

  END
	! </f>

  SUBROUTINE IC_exp_decaying_spectrum(energy_input)
! <f
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

    DO i_x = kMin_x, kMax_x
  	DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

    IF ( k_2( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(theta)
      CALL RANDOM_NUMBER(ph)

      phi     = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      theta   = DACOS( one - two * theta )! Polar angle of \hat{u}_k vector
      ph      = two_pi * ph ! Phases of \hat{u}_k components

      k_ratio = DSQRT( k_2( i_x, i_y, i_z) ) / DBLE( k_integral )

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

      IF (((i_x .NE. 0) .AND. (i_x .NE. kMax_x)) .OR. (i_z .LT. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z ) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z ) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z ) * V_k(3)

      IF (((i_x .EQ. 0) .OR. (i_x .EQ. kMax_x)) .AND. ((i_y .NE. kMin_y) .AND. (i_z .NE. kMin_z))) THEN
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
      ELSE IF ((i_y .EQ. kMin_y) .AND. (i_z .GT. 0)) THEN
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
	! </f>

  SUBROUTINE IC_Kolmogorov_spectrum(energy_input)
! <f
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

    k_integral        = FLOOR( DLOG( DBLE( N_max ) ) / DLOG( 4.0D0 ) ) - 1
    ! Integral scale wavenumber

    k_kol             = FLOOR ( DBLE( N_max ) / 4.0D0 ) - 1
    ! End of kolmogorov spectrum

    integral_exponent = 2
    ! The power in the spectrum E(k) for k<k_integral
    ! Generally either 2 or 4 or 6. But has to be a even number

    eleven_by_twelve  = 11.0D0 /  12.0D0

    CALL normalization_kolmogorov_spectrum( k_integral, k_kol, norm_const)
    ! Returns the 'norm_const', so that theoretically net energy is O(1).
    ! Additionally normalization to any energy can be done with norm_factor
    ! REF-> <<< system_auxilaries >>>

    DO i_x = kMin_x, kMax_x
  	DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

    IF ( k_2( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(theta)
      CALL RANDOM_NUMBER(ph)

      phi                = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      theta              = DACOS( one - two * theta )! Polar angle of \hat{u}_k vector
      ph                 = two_pi * ph ! Phases of \hat{u}_k components

      k_ratio            = DSQRT( k_2( i_x, i_y, i_z) ) / DBLE( k_integral )
      factor_integral    = one
      factor_dissipation = one

      IF ( k_ratio .LT. 1 ) THEN
        CALL kolmogorov_spectrum_integralscale_subpart(k_ratio,integral_exponent,factor_integral)
        ! REF-> <<< system_auxilaries >>>
      END IF

      k_ratio            = DSQRT( k_2( i_x, i_y, i_z) ) / DBLE( k_kol )

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

      IF (((i_x .NE. 0) .AND. (i_x .NE. kMax_x)) .OR. (i_z .LT. 0)) THEN
        v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zx( i_x, i_y, i_z ) * V_k(3)
        v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z) * V_k(2) + &
                                   proj_yz( i_x, i_y, i_z ) * V_k(3)
        v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z) * V_k(2) + &
                                   proj_zz( i_x, i_y, i_z ) * V_k(3)

      IF (((i_x .EQ. 0) .OR. (i_x .EQ. kMax_x)) .AND. ((i_y .NE. kMin_y) .AND. (i_z .NE. kMin_z))) THEN
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
      ELSE IF ((i_y .EQ. kMin_y) .AND. (i_z .GT. 0)) THEN
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
	! </f>

  SUBROUTINE IC_perfect_thermalized_spectrum(energy_input)
! <f
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

    DO i_x = kMin_x, kMax_x
  	DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

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

      IF (((i_x .NE. 0) .AND. (i_x .NE. kMax_x)) .OR. (i_z .LT. 0)) THEN
       v_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * V_k(1) + proj_xy( i_x, i_y, i_z) * V_k(2) + &
                                  proj_zx( i_x, i_y, i_z ) * V_k(3)
       v_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * V_k(1) + proj_yy( i_x, i_y, i_z) * V_k(2) + &
                                  proj_yz( i_x, i_y, i_z ) * V_k(3)
       v_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * V_k(1) + proj_yz( i_x, i_y, i_z) * V_k(2) + &
                                  proj_zz( i_x, i_y, i_z ) * V_k(3)

      IF (((i_x .EQ. 0) .OR. (i_x .EQ. kMax_x)) .AND. ((i_y .NE. kMin_y) .AND. (i_z .NE. kMin_z))) THEN
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
      ELSE IF ((i_y .EQ. kMin_y) .AND. (i_z .GT. 0)) THEN
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
	! </f>

  SUBROUTINE IC_vortex_sheet(energy_input)
! <f
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
    INTEGER(KIND=4) ::i_x0,i_y0,i_z0
    INTEGER(KIND=4) ::i_x1,i_x3
    DOUBLE PRECISION::u0,smooth_pm,c_factor
    DOUBLE PRECISION::energy_sheet,energy_ratio
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_sheet_y

    ALLOCATE( u_sheet_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

    u0                 = one
    ! Normalizing parameter

    smooth_pm          = 0.5D0
    ! How thick the sheet is, smaller the parameter thicker it is, has to be less than 1

    c_factor           = smooth_pm * two_pi / thr
    ! TO KEEP UP THE NOMENCLATURE FOR THIS STUDY.
    ! With this factor => c_factor * i_x = smooth_pm * k_G * x = k_0 * x

    energy_ratio       = 0.01D0
    ! Percentage of energy in Background field

    ! i_x0             = INT( N_x / 8)
    i_x0               = 0
    i_y0               = 0
    i_z0               = 0
    i_x1               = 1 * INT( N_x / 4 )
    i_x3               = 3 * INT( N_x / 4 )

    DO i_z = 0, N_z - 1
  	DO i_y = 0, N_y - 1
  	DO i_x = 0, N_x - 1

      u_sheet_y( i_x, i_y, i_z ) = one + DTANH( - c_factor * DBLE( i_x - i_x1 ) ) &
                                       + DTANH( + c_factor * DBLE( i_x - i_x3 ) )

    END DO
    END DO
    END DO

    energy_sheet = hf * SUM( u_sheet_y ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_sheet )
    u_sheet_y    = u0 * u_sheet_y
    ! Normalization of sheet

    CALL IC_exp_decaying_spectrum( energy_ratio * energy_input )
    ! Gets a background flow for remaining energy

    CALL fft_c2r_vec( v_x, v_y, v_z, u_x, u_y, u_z )
    ! Getting the real velocity to add the sheet

    u_y = u_y + u_sheet_y

    CALL fft_r2c( u_y, v_y )
    ! FFT spectral to real velocity

    IC_type = 'VX-LE'

    DEALLOCATE(u_sheet_y)

  END
	! </f>

  SUBROUTINE IC_vortex_sheet_with_TG(energy_input)
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! An initial condition with vortex sheet imposed with a TAYLOR GREEN FLOW as background.
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
    DOUBLE PRECISION::u0,smooth_pm,c_factor
    INTEGER(KIND=4) ::i_x0,i_y0,i_z0
    INTEGER(KIND=4) ::i_x1,i_x3
    DOUBLE PRECISION::energy_sheet,energy_ratio,energy_TG
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_sheet_y

    ALLOCATE( u_sheet_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

    u0           = one
    ! Normalizing parameter

    smooth_pm    = 0.8D0
    ! How thick the sheet is, smaller the parameter thicker it is, has to be less than 1

    c_factor = smooth_pm * two_pi / thr
    ! TO KEEP UP THE NOMENCLATURE FOR THIS STUDY.
    ! With this factor => c_factor * i_x = smooth_pm * k_G * x = k_0 * x

    energy_ratio = 0.005D0
    ! Percentage of energy in Background field

    ! i_x0 =     INT( N_x /  8)
    i_x0 = 0
    i_y0 = 0
    i_z0 = 0
    i_x1 = 1 * INT( N_x / 4 )
    i_x3 = 3 * INT( N_x / 4 )

    DO i_z = 0, N_z - 1
  	DO i_y = 0, N_y - 1
  	DO i_x = 0, N_x - 1

      u_sheet_y( i_x, i_y, i_z ) = two + DTANH( - c_factor * DBLE( i_x - i_x1 ) ) &
                                       + DTANH( + c_factor * DBLE( i_x - i_x3 ) )

      u_x( i_x, i_y, i_z )       = + DCOS( DBLE( i_x - i_x0 ) * dx )&
                                   * DSIN( DBLE( i_y - i_y0 ) * dy )&
                                   * DCOS( DBLE( i_z - i_z0 ) * dz )
      u_y( i_x, i_y, i_z )       = - DSIN( DBLE( i_x - i_x0 ) * dx )&
                                   * DCOS( DBLE( i_y - i_y0 ) * dy )&
                                   * DCOS( DBLE( i_z - i_z0 ) * dz )

    END DO
    END DO
    END DO

    energy_sheet = hf * SUM( u_sheet_y ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_sheet )
    u_sheet_y    = u0 * u_sheet_y
    ! Normalization of sheet

    energy_TG    = hf * SUM( ( u_x ** two ) + ( u_y ** two ) ) / N3
    u0           = DSQRT( energy_ratio * energy_input / energy_TG )
    u_x          = u0 * u_x
    u_y          = u0 * u_y
    ! Normalisation of TG flow

    u_y = u_y + u_sheet_y
    ! Combining both flows (superposition)

    IC_type = 'VOR-STG'
    ! vortex sheet with Taylor Green

    CALL fft_r2c_vec( u_x, u_y, u_z, v_x, v_y, v_z )
    ! Getting spectral velocity

    DEALLOCATE(u_sheet_y)

  END
	! </f>

  SUBROUTINE IC_vortex_sheet_with_ABC(energy_input)
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! An initial condition with vortex sheet imposed with a ABC FLOW as background.
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
    DOUBLE PRECISION::u0,smooth_pm,c_factor
    DOUBLE PRECISION::A_f,B_f,C_f
    INTEGER(KIND=4) ::i_x0,i_y0,i_z0
    INTEGER(KIND=4) ::i_x1,i_x3
    DOUBLE PRECISION::energy_sheet,energy_ratio,energy_ABC
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_sheet_y

    ALLOCATE( u_sheet_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

    u0           = one
    ! Normalizing parameter

    smooth_pm    = 0.250D0
    ! How thick the sheet is, smaller the parameter thicker it is, has to be less than 1

    c_factor = smooth_pm * two_pi / thr
    ! TO KEEP UP THE NOMENCLATURE FOR THIS STUDY.
    ! With this factor => c_factor * i_x = smooth_pm * k_G * x = k_0 * x

    energy_ratio = 0.05D0
    ! Percentage of energy in Background field

    ! i_x0 =     INT( N_x /  8)
    i_x0 = 0
    i_y0 = 0
    i_z0 = 0
    i_x1 = 1 * INT( N_x / 4 )
    i_x3 = 3 * INT( N_x / 4 )

    A_f = one
    B_f = one
    C_f = one
    ! ABC Flow parameters

    DO i_z = 0, N_z - 1
  	DO i_y = 0, N_y - 1
  	DO i_x = 0, N_x - 1

      u_sheet_y( i_x, i_y, i_z ) = one + DTANH( - c_factor * DBLE( i_x - i_x1 ) ) &
                                       + DTANH( + c_factor * DBLE( i_x - i_x3 ) )

      u_x( i_x, i_y, i_z )       = A_f * DSIN( DBLE( i_z * dz ) )  + C_f * DCOS( DBLE( i_y * dy ) )

      u_y( i_x, i_y, i_z )       = B_f * DSIN( DBLE( i_x * dx ) )  + A_f * DCOS( DBLE( i_z * dz ) )

      u_z( i_x, i_y, i_z )       = C_f * DSIN( DBLE( i_y * dy ) )  + B_f * DCOS( DBLE( i_x * dx ) )

    END DO
    END DO
    END DO

    energy_sheet = hf * SUM( u_sheet_y ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_sheet )
    u_sheet_y    = u0 * u_sheet_y
    ! Normalization of sheet

    energy_ABC   = hf * SUM( ( u_x ** two ) + ( u_y ** two ) + ( u_z ** two ) ) / N3
    u0           = DSQRT( energy_ratio * energy_input / energy_ABC )
    u_x          = u0 * u_x
    u_y          = u0 * u_y
    u_z          = u0 * u_z
    ! Normalisation of ABC flow

    u_y = u_y + u_sheet_y
    ! Combining both flows (superposition)

    IC_type = 'VOR-ABC'
    ! vortex sheet with Taylor Green

    ! v_x = c_zero
    ! v_y = c_zero
    ! v_z = c_zero

    CALL fft_r2c_vec( u_x, u_y, u_z, v_x, v_y, v_z )
    ! Getting spectral velocity

    v_x(0,0,0) = c_zero
    v_y(0,0,0) = c_zero
    v_z(0,0,0) = c_zero
    ! Setting the mean velocity to be zero.

    DEALLOCATE(u_sheet_y)

  END
	! </f>

  SUBROUTINE IC_vortex_tube(energy_input)
! <f
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
    DOUBLE PRECISION::tube_y0,tube_z0
    DOUBLE PRECISION::tube_y,tube_z
    DOUBLE PRECISION::energy_tube,energy_ratio
    DOUBLE PRECISION::u_ang,radius,arg
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_tube_y,u_tube_z

    ALLOCATE(u_tube_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE(u_tube_z( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

    u0           = one
    ! Normalizing parameter

    smooth_pm    = 0.04D0
    ! How thick the sheet is, smaller the parameter thicker it is

    energy_ratio = 0.01D0
    ! Percentage of energy in Background field

    tube_y0      = DBLE( N_y / 2 )
    tube_z0      = DBLE( N_z / 2 )
    ! Center of the tube

    DO i_y = 0, N_y - 1
    DO i_z = 0, N_z - 1

      tube_y                  = DBLE( i_y ) - tube_y0
      tube_z                  = DBLE( i_z ) - tube_z0
      radius                  = DSQRT( tube_y ** two + tube_z ** two )
      arg                     = radius * smooth_pm * two_pi / thr
      u_ang                   = arg * DEXP( - hf * ( arg ** two ) )
      u_tube_y( :, i_y, i_z ) = - u_ang * tube_z / radius
      u_tube_z( :, i_y, i_z ) = + u_ang * tube_y / radius

    END DO
    END DO

    energy_tube  = hf * SUM( u_tube_y ** two + u_tube_z ** two ) / N3
    u0           = DSQRT( ( one - energy_ratio ) * energy_input / energy_tube )
    u_tube_y     = u0 * u_tube_y
    u_tube_z     = u0 * u_tube_z
    ! Normalization of tube

    CALL IC_exp_decaying_spectrum( energy_ratio * energy_input )
    ! Gets a background flow for remaining energy

    CALL fft_c2r_vec( v_x, v_y, v_z, u_x, u_y, u_z )
    ! Getting the real velocity to add the tube velocity

    u_y = u_y + u_tube_y
    u_z = u_z + u_tube_z

    CALL fft_r2c( u_y, v_y )
    CALL fft_r2c( u_z, v_z )
    ! FFT spectral to real velocity

    CALL compute_projected_velocity
    ! Projects the velocity to remove some compressibility

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_initial / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

    IC_type = 'VOR-TUB'

    DEALLOCATE( u_tube_y, u_tube_z )

  END
	! </f>

  SUBROUTINE IC_from_file_spectral
! <f
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

    DO i_x = kMin_x, kMax_x
    DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

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
	! </f>

  SUBROUTINE IC_from_file_real
! <f
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

    DO i_x = 0 , N_x - 1
    DO i_y = 0 , N_y - 1
    DO i_z = 0 , N_z - 1

      READ(44,f_d32p17,ADVANCE = 'NO')  u_x( i_x, i_y, i_z)
      READ(44,f_d32p17,ADVANCE = 'NO')  u_y( i_x, i_y, i_z)
      READ(44,f_d32p17,ADVANCE = 'YES') u_z( i_x, i_y, i_z)

    END DO
    END DO
    END DO

    CLOSE(44)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    CALL fft_r2c_vec( u_x, u_y, u_z, v_x, v_y, v_z )
    ! FFT spectral to real velocity

  END
	! </f>

  SUBROUTINE compute_energy_spectral_data
! <f
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

    DO i_x      =  kMin_x + 1, kMax_x - 1
    DO i_y      =  kMin_y    , kMax_y
    DO i_z      =  kMin_z    , kMax_z
      energy    = energy + CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                           CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                           CDABS( v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO
    END DO

    i_x         =   kMin_x
    DO i_y      =   kMin_y , kMax_y
    DO i_z      =   kMin_z , kMax_z
      energy    = energy + hf * ( CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_z( i_x, i_y, i_z ) ) ** two )
    END DO
    END DO

    i_x         =   kMax_x
    DO i_y      =   kMin_y , kMax_y
    DO i_z      =   kMin_z , kMax_z
      energy    = energy + hf * ( CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_z( i_x, i_y, i_z ) ) ** two )
    END DO
    END DO

  END
	! </f>

  SUBROUTINE compute_projected_velocity
! <f
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
    ALLOCATE( v_P_x( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( v_P_y( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( v_P_z( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )

    v_P_x = v_x
    v_P_y = v_y
    v_P_z = v_z

    v_x   = proj_xx * v_P_x + proj_xy * v_P_y + proj_zx * v_P_z
    v_y   = proj_xy * v_P_x + proj_yy * v_P_y + proj_yz * v_P_z
    v_z   = proj_zx * v_P_x + proj_yz * v_P_y + proj_zz * v_P_z

    DEALLOCATE(v_P_x,v_P_y,v_P_z)

  END
	! </f>

END MODULE system_initialcondition
