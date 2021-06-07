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
! MODULE: system_analysis
! LAST MODIFIED: 3 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! RUN AND OUTPUT MODULE TO RUN THE TIME EVOLUTION FOR 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_analysis
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major system_analysis for the Code to run.
! The standard analysis are done in system_functions module
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_variables
  USE system_fftw
  USE system_output
  USE system_functions

  IMPLICIT NONE
  ! _________________________
  ! MOMENTS VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  :: m_ind, m_ind_max, m_ind_min
  DOUBLE PRECISION :: circulation
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::m_array,m2_inv_array,alpha_array,alpha_by_2m_array
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vorticity_D_moment,vorticity_O_moment

  CONTAINS

  SUBROUTINE allocate_vorticity_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the moments array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    m_ind_min = 1
    m_ind_max = 7

    CALL compute_vorticity

    CALL compute_enstrophy

    circulation = DSQRT( two * enstrophy )

    ALLOCATE( m_array( m_ind_min           : m_ind_max ) )
    ALLOCATE( m2_inv_array( m_ind_min      : m_ind_max ) )
    ALLOCATE( alpha_array( m_ind_min       : m_ind_max ) )
    ALLOCATE( alpha_by_2m_array( m_ind_min : m_ind_max ) )
    ALLOCATE( vorticity_O_moment( m_ind_min : m_ind_max ) )
    ALLOCATE( vorticity_D_moment( m_ind_min : m_ind_max ) )

    DO m_ind = m_ind_min, m_ind_max

      m_array( m_ind ) = DBLE( m_ind )
      m2_inv_array( m_ind ) = one / ( two * m_array( m_ind ) )
      alpha_array( m_ind ) = two * m_array( m_ind ) / ( 4.0D0 * m_array( m_ind ) - thr )
      alpha_by_2m_array( m_ind ) = one / ( 4.0D0 * m_array( m_ind ) - thr )

    END DO

    file_name = TRIM( ADJUSTL( file_address ) ) // 'vorticity_moments_O.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(unit = 8888, file = file_name )

    file_name = TRIM( ADJUSTL( file_address ) ) // 'vorticity_moments_D.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(unit = 8889, file = file_name )

  END

  SUBROUTINE compute_vorticity_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the L-p moments and Rescaled Dimensionless moments of
  ! the vorticity field
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   V  O  R  T  I  C  I  T  Y      M  O  M  E  N  T  S
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO m_ind = m_ind_min, m_ind_max

      vorticity_O_moment( m_ind ) = SUM( ( w_ux ** two + w_uy ** two + w_uz ** two ) ** m_array( m_ind ) )
      vorticity_O_moment( m_ind ) = vorticity_O_moment( m_ind ) / N3
      vorticity_D_moment( m_ind ) = vorticity_O_moment( m_ind ) ** alpha_by_2m_array( m_ind )
      vorticity_D_moment( m_ind ) = vorticity_D_moment( m_ind ) / ( circulation ** alpha_array( m_ind ) )
      vorticity_O_moment( m_ind ) = vorticity_O_moment( m_ind ) ** m2_inv_array( m_ind )

    END DO

    CALL write_vorticity_moments

  END

  SUBROUTINE write_vorticity_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to write the moments
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(8888,f_d8p4,ADVANCE   ='no')  time_now
    WRITE(8889,f_d8p4,ADVANCE   ='no')  time_now

    DO m_ind = m_ind_min, m_ind_max - 1

      WRITE(8888,f_d32p17,ADVANCE ='no') vorticity_O_moment( m_ind )
      WRITE(8889,f_d32p17,ADVANCE ='no') vorticity_D_moment( m_ind )

    END DO

    WRITE(8888,f_d32p17,ADVANCE ='yes') vorticity_O_moment( m_ind )
    WRITE(8889,f_d32p17,ADVANCE ='yes') vorticity_D_moment( m_ind)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IF ( t_step .EQ. t_step_total ) THEN

      CLOSE(8888)
      CLOSE(8889)

    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    END

END MODULE system_analysis
