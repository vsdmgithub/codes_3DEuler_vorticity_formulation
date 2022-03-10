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
! MODULE: system_advvariables
! LAST MODIFIED: 22 JAN 2022
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED VARIABLES TO DO ANALYSIS IN 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advvariables
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major advanced functions involving analysis.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicvariables
  USE system_basicfunctions

  IMPLICIT NONE
  ! _________________________
  ! VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION :: loc_stretching, vx_stretching_max, loc_vx_stretching_max
  DOUBLE PRECISION :: energy_filter,energy_filter_spectral
  DOUBLE PRECISION :: th_norm,purging_alpha,time_purging
  DOUBLE PRECISION :: th_coeff_exponent
  DOUBLE PRECISION :: purg_beta
  INTEGER(KIND=4)  :: t_step_purging,k_P
  ! _________________________
  ! ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::str_xx,str_yy,str_zz
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::str_xy,str_yz,str_zx
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::w_mod_2
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::vx_stretching

  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_str_xx,bck_str_yy,bck_str_zz
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_str_xy,bck_str_yz,bck_str_zx
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::th_coeff
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_str_opr
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_vx_stretching

  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::tr_wave_filter

  CONTAINS

  SUBROUTINE allocate_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the strain tensor array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( w_mod_2( 0       : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( str_xx( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( str_yy( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( str_zz( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( str_xy( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( str_yz( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( str_zx( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ! ALLOCATE( vx_stretching( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

  END

  SUBROUTINE allocate_bck_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the bckground strain operator array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::arg, r_local

    ! r_local = 3.5D0 * dx
    r_local = L_x / DBLE(k_G)
    ! Size of the region to integrate and consider as local

    purging_alpha = 0.8D0
    ! Exponent in deciding the time to purge based on kG

    time_purging = DBLE(k_G) ** (-purging_alpha)

    CALL time_to_step_convert( time_purging , t_step_purging, dt )
    ! REF-> <<< system_auxilaries >>>
print*,time_purging,t_step_purging,dt

    t_step_purging = 5

    WRITE(*,"(A40)")'======================================'
    WRITE(*,"(A20,I4)")'Purging T-Step =',t_step_purging
    WRITE(*,"(A40)")'======================================'

    th_coeff_exponent = 2.0D0
    ! This determines, how steep the erf function goes from 0 to 1

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( bck_str_opr( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( bck_str_xx( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( bck_str_yy( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( bck_str_zz( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( bck_str_xy( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( bck_str_yz( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( bck_str_zx( 0        : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( th_coeff( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ! ALLOCATE( bck_vx_stretching( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )


    DO i_x = kMin_x, kMax_x
  	DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

      arg                          = r_local * DSQRT( k_2( i_x, i_y, i_z ) )

      bck_str_opr( i_x, i_y, i_z ) = thr * ( DSIN( arg ) - arg * DCOS( arg ) ) /  (arg ** thr)

    END DO
    END DO
    END DO

    bck_str_opr( 0, 0, 0 ) = one
    ! Zero mode which is background

  END

  SUBROUTINE allocate_tr_wave_filter
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the wavenumber filter to extract the modes
  ! with truncation waves.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION             ::threshold_angle
    DOUBLE PRECISION             ::threshold_magn,k_mod
    ! DOUBLE PRECISION,DIMENSION(3)::tr_wave_dir

    ! tr_wave_dir     = (/ one, zero, zero /)
    ! Unit normal of truncation wave (expected)

    threshold_angle = ( two_pi / 360.0D0 ) * 15
    ! 75deg angle
    ! threshold_angle = 5.0D0 * two_pi / 24.0D0

    threshold_angle = DCOS( threshold_angle )
    ! Greater than this, will be taken into the tr_wave_filter

    purg_beta       = 0.4D0
    threshold_magn  = k_G - k_G ** ( purg_beta )
    k_P             = FLOOR( threshold_magn )
    ! Threshold (lower) for wavenumber filter
    ! print*,k_P,k_G

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( tr_wave_filter( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )

    DO i_x = kMin_x, kMax_x
  	DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

      k_mod = DSQRT( k_2( i_x, i_y, i_z ) )
      IF ( ( k_mod .GT. DBLE( k_P ) ) .AND. ( k_x( i_x, i_y, i_z ) / k_mod  .GT. threshold_angle ) ) THEN
      ! Checking the angle between \k and \k_G less than threshold

      tr_wave_filter( i_x, i_y, i_z ) = one * truncator( i_x, i_y, i_z )
      END IF

    END DO
    END DO
    END DO

    ! print*,sum(tr_wave_filter),sum(truncator)
  END

	SUBROUTINE deallocate_strain_tensor
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays related to strain
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE( w_mod_2 )
		DEALLOCATE( str_xx, str_yy, str_zz )
		DEALLOCATE( str_xy, str_yz, str_zx )
		! DEALLOCATE( vx_stretching )

	END

	SUBROUTINE deallocate_bck_strain_tensor
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays related to background strain
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE( bck_str_xx, bck_str_yy, bck_str_zz)
    DEALLOCATE( bck_str_xy, bck_str_yz, bck_str_zx)
		DEALLOCATE( th_coeff )
		! DEALLOCATE( bck_vx_stretching )
		DEALLOCATE( bck_str_opr )

	END

  SUBROUTINE deallocate_tr_wave_filter
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays related to truncation wave filter
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE( tr_wave_filter )

	END

END MODULE system_advvariables
