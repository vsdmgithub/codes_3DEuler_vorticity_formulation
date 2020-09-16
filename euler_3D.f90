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
! PROGRAM THAT SOLVES EULER EQUATION
! LAST MODIFIED: 29 JUNE 2020
! ##################

! ##################
! LIST OF MODULES:
! ------------------------------------
!   1. run_and_output
!   2. solver
!   3. variables_and_arrays
!   4. FFT_mod
!   5. STAT_mod
!   6. VTR_mod
! _____________________________________
! ##################

program euler_3D
  
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! All the work is done in the modules. Calling a few would finish
    ! the code. 
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	USE run_and_output
	IMPLICIT NONE
    CHARACTER(LEN=20)::c
    DOUBLE PRECISION::time_start,time_end
    CALL CPU_TIME(time_start)
    
    CALL init_parameters
    ! We get all the variables, and arrays ready to be allocated
    
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y         A  L  L  O  C  A  T  I  O  N . 
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(axis(0:N-1))
    ALLOCATE(k_2(0:Nh,-Nh:Nh-1,-Nh:Nh-1),truncator(0:Nh,-Nh:Nh-1,-Nh:Nh-1),purger(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(k_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),k_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),k_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(proj_xx(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_yy(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_zz(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(proj_xy(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_yz(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_zx(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(shell_no(0:k_G,-k_G:k_G,-k_G:k_G))
    ALLOCATE(count_modes_shell(0:k_G+1))
    ALLOCATE(en_sp(0:max_shell_no))     
    ALLOCATE(en_sp_av(0:max_shell_no))     
    
    CALL init_arrays
    ! Now arrays can be declared values.
    
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  I  M  U  L  A  T  I  O  N       -        R   U   N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    c='y'
    IF (c .EQ. 'y') then
        PRINT*,'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        PRINT*,'   S  I  M  U  L  A  T  I  O  N        S  T  A  R  T  S '
        PRINT*,'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
        CALL time_forward_march
    END IF
    IF ( state_sim .EQ. 1) then
        PRINT*,'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
        PRINT*,'   S  I  M  U  L  A  T  I  O  N        E  N  D  S '
        PRINT*,'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
    END IF
    CALL CPU_TIME(time_end)
    write(*,'(A30,F20.4)'),'Total Run time= ',time_end-time_start
END program euler_3D
