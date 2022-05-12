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
  DOUBLE PRECISION :: energy_thermalised
  DOUBLE PRECISION :: th_norm,purging_alpha,time_purging
  DOUBLE PRECISION :: th_threshold
  DOUBLE PRECISION :: purg_beta,k_P_2
  INTEGER(KIND=4)  :: t_step_purging,k_P
  ! _________________________
  ! ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::w_mod_2
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::w_loc_mod_2
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::w_fil_mod_2
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::w_k_mod_2
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::w_ux_loc,w_uy_loc,w_uz_loc
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::th_filter
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::loc_opr

  CONTAINS

  SUBROUTINE allocate_local_vorticity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the local vorticity operator array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::arg, r_local

    r_local       = L_x / DBLE(k_G)
    ! Size of the region to integrate and consider as local

    k_P           = k_G - 10
    k_P_2         = DBLE( k_P * k_P )
    ! Purging wavenumber

    ! purging_alpha = 0.8D0
    ! ! Exponent in deciding the time to purge based on kG
    !
    ! time_purging  = DBLE(k_G) ** (-purging_alpha)
    !
    ! CALL time_to_step_convert( time_purging , t_step_purging, dt )
    ! REF-> <<< system_auxilaries >>>

    ! t_step_purging =20

    ! WRITE(*,"(A40)")'======================================'
    ! WRITE(*,"(A20,I4)")'Purging T-Step =',t_step_purging
    ! WRITE(*,"(A40)")'======================================'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( loc_opr(   kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( th_filter( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( w_k_mod_2( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( w_ux_loc( 0         : N_x - 1, 0     : N_y - 1, 0     : N_z - 1 ) )
    ALLOCATE( w_uy_loc( 0         : N_x - 1, 0     : N_y - 1, 0     : N_z - 1 ) )
    ALLOCATE( w_uz_loc( 0         : N_x - 1, 0     : N_y - 1, 0     : N_z - 1 ) )
    ALLOCATE( w_mod_2( 0       : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( w_loc_mod_2( 0       : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
    ALLOCATE( w_fil_mod_2( 0       : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )

    th_filter = zero

    DO i_x = kMin_x, kMax_x
  	DO i_y = kMin_y, kMax_y
  	DO i_z = kMin_z, kMax_z

      arg                      = r_local * DSQRT( k_2( i_x, i_y, i_z ) )

      loc_opr( i_x, i_y, i_z ) = one - ( thr * ( DSIN( arg ) - arg * DCOS( arg ) ) / (arg ** thr) )

      IF ( k_2 ( i_x, i_y, i_z ) .GT. k_P_2 ) THEN

        th_filter( i_x, i_y, i_z ) = one

      END IF

    END DO
    END DO
    END DO

    loc_opr( 0, 0, 0 ) = zero
    ! Zero mode which is background

  END

	SUBROUTINE deallocate_local_vorticity
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays related to background strain
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE( w_ux_loc, w_uy_loc, w_uz_loc )
		DEALLOCATE( loc_opr )
    DEALLOCATE( th_filter )
    DEALLOCATE( w_mod_2 )
    DEALLOCATE( w_loc_mod_2 )
    DEALLOCATE( w_fil_mod_2 )
		DEALLOCATE( w_k_mod_2 )

	END

END MODULE system_advvariables
