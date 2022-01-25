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
! LAST MODIFIED: 22 JAN 2022
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

  SUBROUTINE write_vx_dot_section()
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

    OPEN( UNIT = 456, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       i_z = 0
    DO i_x = 0, N_x - 1
    DO i_y = 0, N_y - 1
         WRITE(456,f_d32p17,ADVANCE='no')  hf * w_uy( i_x, i_y, i_z ) + str_zx( i_x, i_y, i_z )
         WRITE(456,f_d32p17,ADVANCE='no')  hf * w_ux( i_x, i_y, i_z ) + str_yz( i_x, i_y, i_z )
         WRITE(456,f_d32p17,ADVANCE='yes') str_zz( i_x, i_y, i_z )
    END DO
    END DO

    CLOSE(456)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_vx_section()
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
                // 'VX_sec_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 455, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       i_z = 0
    DO i_x = 0, N_x - 1
    DO i_y = 0, N_y - 1
      WRITE(455,f_d32p17,ADVANCE='no')  w_ux( i_x, i_y, i_z )
      WRITE(455,f_d32p17,ADVANCE='no')  w_uy( i_x, i_y, i_z )
      WRITE(455,f_d32p17,ADVANCE='yes') w_uz( i_x, i_y, i_z )
    END DO
    END DO

    CLOSE(455)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_strain_section()
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
                // 'STR_sec_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 459, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       i_z = 0
    DO i_x = 0, N_x - 1
    DO i_y = 0, N_y - 1
         WRITE(459,f_d32p17,ADVANCE='no')  str_xx( i_x, i_y, i_z )
         WRITE(459,f_d32p17,ADVANCE='no')  str_yy( i_x, i_y, i_z )
         WRITE(459,f_d32p17,ADVANCE='no')  str_zz( i_x, i_y, i_z )
         WRITE(459,f_d32p17,ADVANCE='no')  str_xy( i_x, i_y, i_z )
         WRITE(459,f_d32p17,ADVANCE='no')  str_yz( i_x, i_y, i_z )
         WRITE(459,f_d32p17,ADVANCE='yes') str_zx( i_x, i_y, i_z )
    END DO
    END DO

    CLOSE(459)
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
                // 'VX_alp_sec_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 788, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! vx_stretching_max = MAXVAL( DABS( vx_stretching ) )
    ! For Normalization by the infinite norm of the vx stretching factor.

       i_z = 0
    DO i_x = 0, N_x - 1
    DO i_y = 0, N_y - 1
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

       i_z = 0
    DO i_x = 0, N_x - 1
    DO i_y = 0, N_y - 1
      loc_stretching = vx_stretching(i_x,i_y,i_z) - bck_vx_stretching(i_x,i_y,i_z)
      WRITE(789,f_d32p17,ADVANCE='no')  vx_stretching( i_x, i_y, i_z )
      WRITE(789,f_d32p17,ADVANCE='yes') loc_stretching
    END DO
    END DO

    CLOSE(789)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_energy_filter
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! writes the energy in the filter
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   E N E R G Y    V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      file_name = TRIM( ADJUSTL( file_address ) ) // 'energy_filter.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 8009, file = file_name )
      ! File where energy vs time will be written. With additional data
    END IF

    WRITE(8009,f_d8p4,ADVANCE   ='no')  time_now
    WRITE(8009,f_d32p17,ADVANCE ='no')  spectral_energy( k_G )
    WRITE(8009,f_d32p17,ADVANCE ='no')  energy_filter
    WRITE(8009,f_d32p17,ADVANCE ='yes') energy_filter_spectral

    IF ( t_step .EQ. t_step_total ) THEN

      CLOSE(8009)

    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

END MODULE system_advoutput
