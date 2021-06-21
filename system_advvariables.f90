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
! LAST MODIFIED: 21 JUNE 2021
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
  INTEGER(KIND=4)  :: m_ind, m_ind_max, m_ind_min
  DOUBLE PRECISION :: circulation
  DOUBLE PRECISION :: loc_stretching, vx_stretching_max, loc_vx_stretching_max
  ! _________________________
  ! ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::m_arr,m2_inv_arr,alpha_arr,alpha_by_2m_arr
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vx_D_moment,vx_O_moment
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vx_dot_moment
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::loc_vx_dot_moment
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::str_xx,str_yy,str_zz,str_xy,str_yz,str_zx
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_str_xx,bck_str_yy,bck_str_zz,bck_str_xy,bck_str_yz,bck_str_zx
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_str_opr,w_mod_2
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_vx_stretching,vx_stretching

  CONTAINS

  SUBROUTINE allocate_vorticity_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the moments array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    m_ind_min = 1
    m_ind_max = 9

    CALL compute_vorticity
    ! REF-> <<< system_basicfunctions >>>

    CALL compute_enstrophy
    ! REF-> <<< system_basicfunctions >>>

    circulation = DSQRT( two * enstrophy )
    ! Calculates a time scale fron initial condition

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( m_arr( m_ind_min            : m_ind_max ) )
    ALLOCATE( m2_inv_arr( m_ind_min       : m_ind_max ) )
    ALLOCATE( alpha_arr( m_ind_min        : m_ind_max ) )
    ALLOCATE( alpha_by_2m_arr( m_ind_min  : m_ind_max ) )
    ALLOCATE( vx_O_moment( m_ind_min : m_ind_max ) )
    ALLOCATE( vx_D_moment( m_ind_min : m_ind_max ) )
    ALLOCATE( w_mod_2( 0 : N - 1, 0 : N - 1, 0 : N - 1 ) )

    m_arr =  (/1.0D0, 1.25D0, 1.5D0, 5.0D0 / 3.0D0, 2.0D0, 2.5D0, 3.0D0, 4.0D0, 6.0D0 /)
    ! List of moment exponents.

    DO m_ind = m_ind_min, m_ind_max

      ! m_arr( m_ind )           = DBLE( m_ind )
      m2_inv_arr( m_ind )      = one / ( two * m_arr( m_ind ) )
      alpha_arr( m_ind )       = two * m_arr( m_ind ) / ( two * two * m_arr( m_ind ) - thr )
      alpha_by_2m_arr( m_ind ) = one / ( two * two * m_arr( m_ind ) - thr )

    END DO

  END

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
    ALLOCATE(str_xx(0:N-1,0:N-1,0:N-1),str_yy(0:N-1,0:N-1,0:N-1),str_zz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(str_zx(0:N-1,0:N-1,0:N-1),str_xy(0:N-1,0:N-1,0:N-1),str_yz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(vx_stretching(0:N-1,0:N-1,0:N-1))
    ALLOCATE(vx_dot_moment( m_ind_min : m_ind_max ) )

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

    r_local = 3.5D0 * dx
    ! Size of the region to integrate and consider as local

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(bck_str_opr(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(bck_str_xx(0:N-1,0:N-1,0:N-1),bck_str_yy(0:N-1,0:N-1,0:N-1),bck_str_zz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(bck_str_zx(0:N-1,0:N-1,0:N-1),bck_str_xy(0:N-1,0:N-1,0:N-1),bck_str_yz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(bck_vx_stretching(0:N-1,0:N-1,0:N-1))
    ALLOCATE(loc_vx_dot_moment( m_ind_min : m_ind_max ) )

    DO i_x = 0, Nh
    DO i_y = -Nh, Nh - 1
    DO i_z = -Nh, Nh - 1

      arg                 = r_local * DSQRT( k_2( i_x, i_y, i_z ) )

      bck_str_opr( i_x, i_y, i_z ) = thr * ( DSIN( arg ) - arg * DCOS( arg ) ) /  (arg ** thr)

    END DO
    END DO
    END DO

    bck_str_opr( 0, 0, 0 ) = one
    ! Zero mode which is background

  END

	SUBROUTINE deallocate_vorticity_moments
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays related to vorticity moments
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(w_mod_2)
    DEALLOCATE(m_arr,m2_inv_arr)
    DEALLOCATE(alpha_arr,alpha_by_2m_arr)
    DEALLOCATE(vx_O_moment,vx_D_moment)

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
		DEALLOCATE(str_xx,str_yy,str_zz)
		DEALLOCATE(str_xy,str_yz,str_zx)
		DEALLOCATE(vx_stretching)
    DEALLOCATE(vx_dot_moment)

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
    DEALLOCATE(bck_str_xx,bck_str_yy,bck_str_zz)
    DEALLOCATE(bck_str_xy,bck_str_yz,bck_str_zx)
		DEALLOCATE(bck_vx_stretching)
		DEALLOCATE(bck_str_opr)
    DEALLOCATE(loc_vx_dot_moment)

	END

END MODULE system_advvariables
