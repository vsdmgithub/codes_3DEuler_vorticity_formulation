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
! MODULE: system_advoutput
! LAST MODIFIED: 21 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED OUTPUT MODULE - RELATED TO ADVANCED FUNCTION MODULE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advoutput
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all output calls from advanced functions module.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_advvariables
  USE system_basicoutput

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE write_moment_exponents
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to write the list of exponents
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_mom ) ) &
    // 'moments_exponent_list.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(unit = 8880, file = file_name )

    DO m_ind = m_ind_min, m_ind_max

      WRITE(8880,f_d8p4,ADVANCE ='yes') m_arr( m_ind )
      ! SAVING THE MOMENT EXPONENTS

    END DO

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CLOSE(8880)

  END

  SUBROUTINE write_vorticity_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to write the vorticity moments at a given time
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    IF ( t_step .EQ. 0 ) THEN

      CALL write_moment_exponents

      file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_mom ) ) &
      // 'VX_mom_O.dat'

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 8888, file = file_name )

      file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_mom ) ) &
      // 'VX_mom_D.dat'

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 8889, file = file_name )

    END IF

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(8888,f_d8p4,ADVANCE   ='no')  time_now
    WRITE(8889,f_d8p4,ADVANCE   ='no')  time_now

    DO m_ind = m_ind_min, m_ind_max - 1

      WRITE(8888,f_d32p17,ADVANCE ='no') vx_O_moment( m_ind )
      WRITE(8889,f_d32p17,ADVANCE ='no') vx_D_moment( m_ind )

    END DO

    WRITE(8888,f_d32p17,ADVANCE ='yes') vx_O_moment( m_ind )
    WRITE(8889,f_d32p17,ADVANCE ='yes') vx_D_moment( m_ind )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IF ( t_step .EQ. t_step_total ) THEN

      CLOSE(8888)
      CLOSE(8889)

    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_vorticity_dot_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to write the moments for the vorticity_dot
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    IF ( t_step .EQ. 0 ) THEN

      file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_mom ) ) &
      // 'VX_dot_mom.dat'

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 7777, file = file_name )

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE(7777,f_d8p4,ADVANCE   ='no')  time_now

    END IF

    DO m_ind = m_ind_min, m_ind_max - 1

      WRITE(7777,f_d32p17,ADVANCE ='no') vx_dot_moment( m_ind )

    END DO

    WRITE(7777,f_d32p17,ADVANCE ='yes') vx_dot_moment( m_ind )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IF ( t_step .EQ. t_step_total ) THEN

      CLOSE(7777)

    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_vorticity_dot_loc_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to write the moments for the vorticity_dot
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    IF ( t_step .EQ. 0 ) THEN

      file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_mom ) ) &
      // 'VX_dot_loc_mom.dat'

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 7779, file = file_name )

    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(7779,f_d8p4,ADVANCE   ='no')  time_now

    DO m_ind = m_ind_min, m_ind_max - 1

      WRITE(7779,f_d32p17,ADVANCE ='no') loc_vx_dot_moment( m_ind )

    END DO

    WRITE(7779,f_d32p17,ADVANCE ='yes') loc_vx_dot_moment( m_ind )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IF ( t_step .EQ. t_step_total ) THEN

      CLOSE(7779)

    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_vorticity_section()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a real space data given for a particular section
  ! into a .dat file named d_nam
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'vty_sec_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 455, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       i_y = Nh - 1
    DO i_x = 0, N - 1
    DO i_z = 0, N - 1
         WRITE(455,f_d32p17,ADVANCE='no')  w_ux(i_x, i_y, i_z)
         WRITE(455,f_d32p17,ADVANCE='no')  w_uy(i_x, i_y, i_z)
         WRITE(455,f_d32p17,ADVANCE='yes') w_uz(i_x, i_y, i_z)
    END DO
    END DO

    CLOSE(455)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_vx_stretching_section()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a real space data given for a particular section
  ! into a .dat file named d_nam
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'VX_dot_sec_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 788, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! vx_stretching_max = MAXVAL( DABS( vx_stretching ) )
    ! For Normalization by the infinite norm of the vx stretching factor.

       i_y = Nh - 1
    DO i_x = 0, N - 1
    DO i_z = 0, N - 1
         WRITE(788,f_d32p17,ADVANCE='yes') vx_stretching(i_x,i_y,i_z)
    END DO
    END DO

    CLOSE(788)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_loc_vx_stretching_section()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a real space data given for a particular section
  ! into a .dat file named d_nam
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_2D ) ) &
                // 'VX_dot_loc_sec_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 789, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! loc_vx_stretching_max = MAXVAL( DABS( vx_stretching - bck_vx_stretching ) )
    ! For Normalization by the infinite norm of the vx stretching factor.

       i_y = Nh - 1
    DO i_x = 0, N - 1
    DO i_z = 0, N - 1
         loc_stretching = vx_stretching(i_x,i_y,i_z) - bck_vx_stretching(i_x,i_y,i_z)
         WRITE(789,f_d32p17,ADVANCE='yes') loc_stretching
    END DO
    END DO

    CLOSE(789)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

END MODULE system_advoutput
