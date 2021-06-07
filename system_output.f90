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
! MODULE: system_output
! LAST MODIFIED: 3 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! OUTPUT MODULE, WHERE DATA IS WRITTEN IN FILES
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_output
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major system_output for the Code to run.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_variables

  IMPLICIT NONE
  ! _________________________
  ! OUTPUT VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  CHARACTER(LEN=140)::file_name
  CHARACTER(LEN=40)::file_time
  CHARACTER(LEN=40)::path_dir
  CHARACTER(LEN=40)::type_sim
  CHARACTER(LEN=60)::name_sim
  CHARACTER(LEN=100)::file_address
  CHARACTER(LEN=40)::sub_dir_3D,sub_dir_sp

  CONTAINS

  SUBROUTINE create_output_directories
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Name the folders, create them, open files to write system_output .
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE

    path_dir    =   '../euler_data/'
    ! path of the main directory relative to this file.

    sub_dir_3D  =   '3D_data/'
    ! Sub directory name to store 3D data - large file sizes.

    sub_dir_sp  =   'spectral_data/'
    ! Sub directory name to store spectral data

    type_sim    =   'classic_N' // TRIM( ADJUSTL( N_char ) ) // '/'
    ! type of simulation, the data is storing

    CALL get_simulation_name(name_sim)
    ! Creating dated and timed name for the simulation for this particular type

    ! name_sim    =   'test_sim'
    ! Use this to give CUSTOM SIMULATION NAME

    file_address =   TRIM( ADJUSTL( path_dir ) ) // TRIM( ADJUSTL( type_sim ) ) //  &
                     TRIM( ADJUSTL( name_sim ) ) // '/'
    ! Address should be added to all file names, if needed sub-dir can be declared later and appended to
    ! this address

    ! CALL SYSTEM('mkdir ' // TRIM( ADJUSTL( path_dir ) ) )

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL( path_dir ) ) // TRIM( ADJUSTL( type_sim ) ) )

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) )

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) // TRIM( ADJUSTL( sub_dir_sp ) ) )

    CALL SYSTEM('mkdir '// TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_3D ) ) )

    ! Command to create the main directory and sub directories (name_sim) in the desired path
    ! If exists already, it won't be an error

    CALL write_simulation_details
    ! Writes the parameters used in the simulation

	END

  SUBROUTINE write_test_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the test_data of any new verification/validation conducted in the code
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN=80)::test_file_name

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   E N E R G Y    V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      test_file_name = TRIM( ADJUSTL( file_address ) ) // 'test_data.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 9009, file = test_file_name )
      ! File where energy vs time will be written. With additional data
    END IF

    WRITE(9009,f_d8p4,advance   ='no')   time_now
    WRITE(9009,f_d32p17,advance ='yes')  u_x( 2, 2, 2)

    IF ( t_step .EQ. t_step_total ) THEN
      CLOSE(9009)
    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_simulation_details
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !   Write the details of the simulation
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE

    file_name = TRIM(ADJUSTL(file_address))//'system_details'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  I  M  U  L  A  T  I  O  N     D  E  T  A  I  L  S
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(UNIT =233,FILE=TRIM( ADJUSTL( file_name ) ) // '.dat')

    WRITE(233,"(A50)")TRIM(ADJUSTL('-------------------------------------------------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('------3D EULER EQUATION (INCOMPRESSIBLE)--------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('_______________________________________________________'))
    WRITE(233,*)
    WRITE(233,"(A50)")TRIM(ADJUSTL('--------------------------------------------------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('-----------PARAMETERS OF SIMULATION------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('--------------------------------------------------------'))
    WRITE(233,"(A20,A2,I8)")     'Resolution    ','= ',N
    WRITE(233,"(A20,A2,I8)")     'Trunc. Mode  ','= ',k_G
    WRITE(233,"(A20,A2,ES8.2)")  'Time step   ','= ',dt
    WRITE(233,"(A20,A2,I8)")     'Total time steps   ','= ',t_step_total
    WRITE(233,"(A20,A2,F8.4)")   'Total time ','= ',time_total
    WRITE(233,"(A20,A2,I8)")     'No of saves   ','= ',no_of_saves
    WRITE(233,"(A20,A2,F8.4)")   'Initial energy ','= ',energy
    WRITE(233,"(A20,A2,F8.2)")   'Initial enstrophy ','= ',enstrophy
    WRITE(233,*)
    WRITE(233,"(A50)")TRIM(ADJUSTL('_______________________________________________________'))

    CLOSE(233)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_temporal_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in time, for every timestep
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   E N E R G Y    V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      file_name = TRIM( ADJUSTL( file_address ) ) // 'energy_vs_time.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 4004, file = file_name )
      ! File where energy vs time will be written. With additional data
    END IF

    WRITE(4004,f_d8p4,advance   ='no')  time_now
    WRITE(4004,f_d32p17,advance ='no')  energy
    WRITE(4004,f_d32p17,advance ='yes') enstrophy

    IF ( t_step .EQ. t_step_total ) THEN
      CLOSE(4004)
    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_spectral_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Write the spectral energy into a file with time index
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_sp ) ) &
                // 'spectral_energy_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y      S  P  E  C  T  R  U  M
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 1001, FILE = file_name )
    DO k_no = 1 , k_G

      WRITE(1001,f_i8,advance  ='no')       k_no
      WRITE(1001,f_d32p17,advance ='no')    spectral_energy(k_no)
      WRITE(1001,f_d32p17,advance ='yes')   spectral_energy_avg(k_no)

    END DO
    CLOSE(1001)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_spectral_velocity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a spectral velocity, in such a way that it can be used
  ! for input in the next simulation using 'IC_from_file' subroutine.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL ( sub_dir_3D )) &
                // 'spectral_velocity_' //TRIM( ADJUSTL( N_char ) ) // '_input.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! File where energy vs time will be written. With additional data

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( unit = 73, file = file_name )
    DO i_x =   0, Nh
    DO i_y = -Nh, Nh - 1
    DO i_z = -Nh, Nh - 1

      WRITE(73,f_c32p17,ADVANCE ='no')    DREAL(v_x(i_x,i_y,i_z)),DIMAG(v_x(i_x,i_y,i_z))
      WRITE(73,f_c32p17,ADVANCE ='no')    DREAL(v_y(i_x,i_y,i_z)),DIMAG(v_y(i_x,i_y,i_z))
      WRITE(73,f_c32p17,ADVANCE ='yes')   DREAL(v_z(i_x,i_y,i_z)),DIMAG(v_z(i_x,i_y,i_z))

    END DO
    END DO
    END DO

    CLOSE(73)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_velocity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a real velocity. To read and plot velocity functions
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL ( sub_dir_3D )) // 'velocity_' &
              //TRIM( ADJUSTL( N_char ) ) // '_t_' // TRIM( ADJUSTL ( file_time ) ) // '.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! File where energy vs time will be written. With additional data

    OPEN( unit = 74, file = file_name )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO i_x = 0 , N - 1
    DO i_y = 0 , N - 1
    DO i_z = 0 , N - 1

      WRITE(74,f_d32p17,ADVANCE ='NO')    u_x(i_x,i_y,i_z)
      WRITE(74,f_d32p17,ADVANCE ='NO')    u_y(i_x,i_y,i_z)
      WRITE(74,f_d32p17,ADVANCE ='YES')   u_z(i_x,i_y,i_z)

    END DO
    END DO
    END DO

    CLOSE(74)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE print_error_timestep()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! print error if time step chosen is large
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,'(A60)')       TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
    WRITE(*,'(A40)')       TRIM( ADJUSTL( 'ERROR: TIME STEP TOO LARGE') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '-----------------------------------------------------------------') )
    WRITE(*,'(A40,F10.6)') TRIM( ADJUSTL( ' RESET THE TIME STEP (AT MAX) AS :') ),dt_max
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '_________________________________________________________________') )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE print_running_status
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! print the running status of the program, when called
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    IF ( t_step .EQ. 0 ) THEN
      PRINT*,'-----------------------------------------------------------'
      PRINT*,'TIME  |   ENERGY   |   ENSTROPHY    |   INCOMPRESSIBILITY '
      PRINT*,'-----------------------------------------------------------'
    END IF

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,'(F6.3,A3,F8.4,A3,F12.4,A7,E12.4)') time_now,'   ',energy,'   ',enstrophy,'     ',k_dot_v_norm
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IF (NaN_count .NE. 0) then
      PRINT*,' '
      PRINT*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      PRINT*,' '
      PRINT*,"NaN ENCOUNTERED BEFORE T = ",time_now
      PRINT*,' '
      PRINT*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      PRINT*,"-------------------SIMULATION STOPPED -------------------"
      PRINT*,' '
      debug_error = 1
      ! IF any NaN is encountered, the loop is exited (in time_evolution subroutine) without any further continuation.
    END IF

    IF ( k_dot_v_error .NE. 0 ) then
      PRINT*,' '
      PRINT*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      PRINT*,' '
      PRINT*,"INCOMPRESSIBILITY LOST BEYOND TOLERANCE BEFORE T = ",time_now
      PRINT*,' '
      PRINT*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      PRINT*,"-------------------SIMULATION STOPPED -------------------"
      PRINT*,' '
      debug_error = 1
      ! IF the incompressibility norm is beyond the tolerance of 'tol_float'. This error arises
    END IF

    IF ( t_step .EQ. t_step_total ) THEN
      PRINT*,'-----------------------------------------------------------'
    END IF

  END

END MODULE system_output
