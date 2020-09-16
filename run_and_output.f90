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
! MODULE: run_and_output
! LAST MODIFIED: 29 JUNE 2020
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! RUN AND OUTPUT MODULE TO RUN, STUDY, PRINT OUTPUT FOR 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE run_and_output
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This module runs the Euler code, using solver module and other modules.
! 1. Main subroutine is 'time_forward_march' which moves through time steps
! 2. Subroutine 'spectral_data' collects spectral data and writes a shell wise energy.
! 3. Subroutine 'vortex_stretching' analyzes the way vorticity is stretched at every grid.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    USE solver
    USE paraview
    USE statistics
    ! List of modules refered and used.
    
    IMPLICIT NONE
    TYPE(VTR_file_handle)::fd_en          ! creating dataypes to store paraview files for 3d viewing   
    CONTAINS
    SUBROUTINE time_forward_march
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! Loop of time steps, where at each step the spectral velocities
    ! are updated through any of the algoritm. Meanwhile, analysis and
    ! outputs are printed respectively.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        
        ALLOCATE(u_x(0:N-1,0:N-1,0:N-1),u_y(0:N-1,0:N-1,0:N-1),u_z(0:N-1,0:N-1,0:N-1))
        ALLOCATE(v_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  I  N  I  T  I  A  L        C  O  N  D  I  T  I  O  N
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL init_initcondn
        ! Calls the subroutine to get a initial condition with norm_factor=1
        
        CALL compute_spectral_energy
        ! Gets the energy from spectral space
        
        norm_factor=DSQRT(en_initial/norm_k)
        ! Normalizing the norm_factor, so that we get energy='en_initial'
        
        CALL init_initcondn
        ! Calls it again to get normalized velocities.

         CALL compute_spectral_energy
        ! Gets the energy from spectral space   
        
        CALL write_parameters
        ! Writes the parameters

        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        OPEN(unit=4,file=TRIM(ADJUSTL(file_dir))//TRIM(ADJUSTL(file_name_en))//'.dat')
        ! File where energy vs time will be written. With additional data
        
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
        !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
        !     E E U U L L E E R R    E E V V O O L L U U T T I I O O N N
        !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
        ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
        ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        !             S        T         A         R       T
        ! 8888888888888888888888888888888888888888888888888888888888888888
 
        PRINT*,'-----------------------------------------------------------'
        PRINT*,'TIME  |  ENERGY  |  INCOMPRESSIBILITY  |  HELICITY'
        PRINT*,'-----------------------------------------------------------'
        DO t_step=0,t_step_total

            CALL step_to_time_convert(t_step,time_now)
            ! Converts the 't_step' to actual time 'time_now'

            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  A  N  A  L  Y  S  I  S       C   A   L   C  .
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            CALL spectral_data  ! ALSO PRINTS AT SAVE TIMES
            
            IF (viscosity_status .EQ. 0) THEN  
                CALL find_k_th
            END IF
            
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  E  N  E  R  G  Y      V  S      T  I  M  E
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                WRITE(4,f_d8p4,advance='no'),time_now
                WRITE(4,f_d32p17,advance='no'),SUM(en_sp)
                WRITE(4,f_d32p17,advance='no'),Z_total
            IF (viscosity_status .EQ. 0) THEN  
                WRITE(4,f_d32p17,advance='no'),en_th
                WRITE(4,f_i4,advance='no'),k_th
            ELSE
                WRITE(4,f_d32p17,advance='no'),SUM(en_sp(k_kol:))
                WRITE(4,f_d32p17,advance='no'),SUM(en_sp(:k_kol-1))
            END IF

            
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  P  U  R  G  I  N  G         T  I  M  E
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            IF ((purging_status .EQ. 1) .AND. (MOD(t_step,t_step_purg) .EQ. 0)) then
            ! Its time to purge !
                v_x=purger*v_x
                v_y=purger*v_y
                v_z=purger*v_z
                purging_count=purging_count+1
            END IF
            IF (purging_status .EQ. 1 ) THEN
                WRITE(4,f_d32p17,advance='no'),SUM(en_sp(k_P:))
                WRITE(4,f_d32p17,advance='no'),flux_k_P(1)
                WRITE(4,f_d32p17,advance='no'),flux_k_P(2)
            END IF
            WRITE(4,f_i4,advance='yes'),purging_count
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ALLOCATE(v_x_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(conv_v_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),conv_v_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),conv_v_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(conv_u_x(0:N-1,0:N-1,0:N-1),conv_u_y(0:N-1,0:N-1,0:N-1),conv_u_z(0:N-1,0:N-1,0:N-1))
            ALLOCATE(grad_x(0:N-1,0:N-1,0:N-1),grad_y(0:N-1,0:N-1,0:N-1),grad_z(0:N-1,0:N-1,0:N-1))
            ALLOCATE(dv1_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(dv3_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(dv1_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(dv3_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(dv1_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(dv3_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            ALLOCATE(v_x_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L  G
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            CALL evolve_algorithm_rk4
            ! Updates v_x,v_y,v_z for next time step
            
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  D  E  A  L  L  O  C  A  T  I  O  N
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            DEALLOCATE(v_x_dot,v_y_dot,v_z_dot)
            DEALLOCATE(conv_v_x,conv_v_y,conv_v_z)
            DEALLOCATE(conv_u_x,conv_u_y,conv_u_z)
            DEALLOCATE(grad_x,grad_y,grad_z)
            DEALLOCATE(dv1_x,dv2_x)
            DEALLOCATE(dv3_x,dv4_x)
            DEALLOCATE(dv1_y,dv2_y)
            DEALLOCATE(dv3_y,dv4_y)
            DEALLOCATE(dv1_z,dv2_z)
            DEALLOCATE(dv3_z,dv4_z)
            DEALLOCATE(v_x_temp,v_y_temp,v_z_temp)
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  D  E  B  U  G             F  O  R          N  a   N
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            IF (MOD(t_step,t_step_debug) .EQ. 0) then
                CALL compute_spectral_energy
                CALL compute_real_energy
                CALL compute_helicity
                WRITE(*,'(F6.3,A3,F8.4,A3,F12.8,A8,F12.8)'),time_now,'   ',norm_k,'   ',k_dot_v,'         ',hely
                CALL compute_compressibility
                IF (count_nan .NE. 0) then
                    PRINT*,"NaN encountered before t=",time_now
                    exit
                    ! IF any NaN is encountered, the loop is exited without any further continuation.
               END IF
            END IF

        END DO
        PRINT*,'-----------------------------------------------------------'
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
        !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
        !     E E U U L L E E R R    E E V V O O L L U U T T I I O O N N
        !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
        ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
        ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        !                    E     N     D 
        ! 8888888888888888888888888888888888888888888888888888888888888888
        IF (pvd_status .EQ. 1) then
            CALL VTR_collect_file(FD=fd_en)
        END IF
        ! Standard routine to collect the pvd file
        
        CLOSE(4)
        
        state_sim=1
        ! Stating that the simulation has ended.
	END
    SUBROUTINE spectral_data
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! This calculates energy, spectral shell wise. It goes through each
    ! spectral mode and puts the energy in the corresponding shell.
    ! This gives the ENERGY SPECTRUM.  When the time is right, it saves
    ! them too.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        Z_total=0.0D0
        en_sp=0.0D0
        ! Reset the array
        
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  E  N  E  R  G  Y     S  P  E  C  T  R  U  M     C  A   L   C.
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DO j_z=-k_G,k_G
        DO j_y=-k_G,k_G
        DO j_x=1,k_G
            CALL mode_energy
            ! Gets the energy in that particular mode (j_x,j_y,j_z) into 'en_temp'
            en_sp(shell_no(j_x,j_y,j_z))=en_sp(shell_no(j_x,j_y,j_z))+en_mode
            Z_total=Z_total+k_2(j_x,j_y,j_z)*en_mode ! Enstrophy net
        END DO
        END DO
        END DO
        j_x=0
        DO j_z=-k_G,-1
        DO j_y=-k_G,k_G
            CALL mode_energy             ! Gets the energy in that particular mode (j_x,j_y,j_z) into 'en_temp'
            en_sp(shell_no(j_x,j_y,j_z))=en_sp(shell_no(j_x,j_y,j_z))+en_mode
            ! Keeps adding energy to that particular shell (|k| fixed), from all solid angles
            Z_total=Z_total+k_2(j_x,j_y,j_z)*en_mode ! Enstrophy net
        END DO
        END DO
        j_z=0
        DO j_y=0,k_G
            CALL mode_energy             ! Gets the energy in that particular mode (j_x,j_y,j_z) into 'en_temp'
            en_sp(shell_no(j_x,j_y,j_z))=en_sp(shell_no(j_x,j_y,j_z))+en_mode
            ! Keeps adding energy to that particular shell (|k| fixed), from all solid angles
            Z_total=Z_total+k_2(j_x,j_y,j_z)*en_mode ! Enstrophy net
        END DO
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  F  L  U  X       A  T      P  U  R  G  I  N  G      B  D  Y
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF (purging_status .EQ. 1) THEN
             IF (t_step .EQ. 0) THEN
                  flux_k_P(1)=SUM(en_sp(:k_P-1))
             END IF
             flux_k_P(2)=(flux_k_P(1)-SUM(en_sp(:k_P-1)))/dt
             flux_k_P(1)=SUM(en_sp(:k_P-1))
        END IF
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  S  H  E  L  L      A  V  E  R  A  G  I  N  G  
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        en_sp_av(1)=en_sp(1)
        en_sp_av(max_shell_no)=en_sp(max_shell_no)
        DO s=2,max_shell_no-1
            en_sp_av(s)=0.25D0*(en_sp(s-1)+en_sp(s+1))+0.5D0*(en_sp(s))
        END DO
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  E  N  E  R  G  Y        O  U  T  P  U  T
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF (MOD(t_step,t_step_save) .EQ. 0) then
            CALL write_spectral_energy
            ! Writes the energy spectrum.
        END IF
    END
   SUBROUTINE mode_energy
    ! CALL this to get the energy of that particular mode
        IMPLICIT NONE
        en_mode=CDABS(v_x(j_x,j_y,j_z))**two+CDABS(v_y(j_x,j_y,j_z))**two+CDABS(v_z(j_x,j_y,j_z))**two
    END
    SUBROUTINE find_k_th
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! NOTE: Only for Truncated/Purged systems.
    ! As thermalization begins, the bottom of the trench in the
    ! energy spectrum is denoted as the thermalization wavenumber.
    ! Trend is that starts near k_G and decreases all the way to
    ! the inertial range, which is a completely thermalized state.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        
        IF (en_sp(k_P) .GT. tol) then
            act_therm=1
        END IF ! Activation of last shell

        IF (act_therm .EQ. 1) THEN
            k_th=MINLOC(en_sp(k_int:k_P), DIM=1)
            en_th=SUM(en_sp(k_th:k_P))
            ! Allowing the 'k_th' to reach till the integral scale 'k_int'
        ELSE
            k_th=k_P
            en_th=0.0D0
        END IF
    END
    SUBROUTINE write_spectral_energy
    ! CALL this to WRITE the spectral energy into a file with time index
        IMPLICIT NONE
        file_name='spectral_energy_t_'
        WRITE (file_time,f_d8p4),time_now
        ! Writes 'time_now' as a CHARACTER
        OPEN(unit=1,file=TRIM(ADJUSTL(file_dir))//TRIM(ADJUSTL(file_name))//TRIM(ADJUSTL(file_time))//'.dat')
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  P  R  I  N   T          O  U  T  P  U  T   -   SPECTRUM
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DO s1=1,k_G-1
             WRITE(1,f_d12p2,advance='no'),DFLOAT(s1)
             WRITE(1,f_d32p17,advance='no'),en_sp(s1)
             WRITE(1,f_d32p17,advance='yes'),en_sp_av(s1)
        END DO
        CLOSE(1)
    END
   END MODULE run_and_output
