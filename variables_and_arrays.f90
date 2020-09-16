! --------------------------------------------------------------
! -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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
! MODULE: variables_and_arrays 
! LAST MODIFIED: 29 JUNE 2020
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! VARIABLES AND ARRAYS FOR 3D EULER EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE variables_and_arrays
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the global variables and arrays for the simulation space are declared and given values
! here, wheras temp_vorary variables (IF necessary) are declared withIN the SUBROUTINEs
! Further, each variable is classified as space, time, fluid, viscosity, constants, output,
! and etc. based on where its purpose suits apt.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	IMPLICIT  NONE
    ! _________________________
    ! SPACE VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::lth,dx,N3,vol,dx_3
	INTEGER (KIND=4)::N,Nh,k_G,k_G_2
    INTEGER (KIND=4)::i_x,i_y,i_z
    INTEGER(KIND=4)::j_x,j_y,j_z
    ! _________________________
    ! TIME VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::dt,time_total,time_save,time_now
	INTEGER (KIND=4)::t_step,t_step_total,t_step_save
    INTEGER (KIND=4)::t_step_save_pvd,t_step_debug,t_step_purg
    INTEGER (KIND=4)::no_of_saves,no_of_pvd_saves,no_of_debug
    ! _________________________
    ! PURGING VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER(KIND=4)::k_P,purging_count
    DOUBLE PRECISION::k_P_2,time_purg
    DOUBLE PRECISION::purg_alpha,purg_beta
    DOUBLE PRECISION,DIMENSION(2)::flux_k_P
    ! _________________________
    ! FLUID VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::viscosity,dissip_rate
    INTEGER(KIND=4)::k_int,k_kol
    ! _________________________
    ! CONSTANTS
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX,parameter::i=DCMPLX(0.0D0,1.0D0),c0=DCMPLX(0.0D0,0.0D0)
	DOUBLE PRECISION,parameter::two_pi=ATAN(1.0D0)*8.0D0
    DOUBLE PRECISION,parameter::hf=0.5D0,two=2.0D0,thr=3.0D0,six=6.0D0
    DOUBLE PRECISION,parameter::tol=0.000001D0
    INTEGER (KIND=4)::viscosity_status,purging_status,pvd_status
    ! _________________________
    ! ENERGY CALC VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::s,s0,s1,max_shell_no
    INTEGER(KIND=4)::k_th,count_nan,act_therm
    DOUBLE PRECISION::norm_k,norm_x,en_th,en_initial,hely
    DOUBLE PRECISION::en_mode,norm_factor,Z_total,k_dot_v
    ! _________________________
    ! OUTPUT VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER (KIND=4)::state_sim
    CHARACTER(LEN=60)::file_dir,file_dir_pvd
    CHARACTER(LEN=60)::file_name,file_name_en,file_name_param,file_name_input
    CHARACTER(LEN=60)::file_name_pvd_es,file_name_pvd_1,file_time
    CHARACTER(LEN=60)::N_char,alpha_char,beta_char
    ! _________________________
    ! FORMATS FOR I/O
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN=20),PARAMETER::f_i2='(I2)',f_i4='(I4)',f_i8='(I8)',f_i16='(I16)'
    CHARACTER(LEN=20),PARAMETER::f_d8p4='(F8.4)',f_d12p6='(F12.6)',f_d16p8='(F16.8)',f_d5p2='(F5.2)'
    CHARACTER(LEN=20),PARAMETER::f_d12p2='(F12.2)',f_d32p17='(F32.17)',f_d12p8='(F12.8)'
    CHARACTER(LEN=20),PARAMETER::f_e5p2='(ES6.2)',f_e10p4='(ES12.4)'
    ! _________________________________________
    ! SPACE, VELOCITY, PSEUDO-SPECTRAL ARRAYS
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX,DIMENSION(3)::temp_v
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::axis                                        ! Arrays that store the grid
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x,v_y,v_z                               ! Spectral velocity matrix (will be updated after every time step)
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::u_x,u_y,u_z                             ! Real velocity matrix (used for gettINg convection term IN fourier space)
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::w_ux,w_uy,w_uz                           ! Real vorticity
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::k_x,k_y,k_z,k_2,truncator,purger        ! wavenumber,truncator matrix
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::proj_xx,proj_yy,proj_zz
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::proj_xy,proj_yz,proj_zx                 ! projection operators matrix \mathbb{P}_{ij}=\delta_{ij}-\frac{k_ik_j}{k^2}}
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_dot,v_y_dot,v_z_dot                   ! spectral velocity derivative matrix 
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::conv_v_x,conv_v_y,conv_v_z                ! Spectral convection term matrix
	DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::conv_u_x,conv_u_y,conv_u_z              ! Real convection term matrix
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::grad_x,grad_y,grad_z                    ! Calculates gradient IN real space from FFT 
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::dv1_x,dv2_x,dv3_x,dv4_x,dv1_y,dv2_y
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::dv3_y,dv4_y,dv1_z,dv2_z,dv3_z,dv4_z        ! RK4 INtermediate matrices
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_temp,v_y_temp,v_z_temp           ! temporary matrices to store velocities durINg RK4 algorithm
    ! _________________________________________
    ! SPECTRAL ENERGY ARRAYS
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::en_sp,en_sp_av                               ! Spectral data INtegrated over spherical shell of radius k_mod
    INTEGER(KIND=4),DIMENSION(:,:,:),ALLOCATABLE::shell_no                                  ! Every grid has its modulus, |k|
    INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE::count_modes_shell                             ! This counts no of modes that have (k-1)<|\vec{k}|<k
    !TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
	contaINs
	SUBROUTINE init_parameters
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !       INitialize all the variables, either IN the code, or READ from a file 'file_name_INput'
    ! variables are explaINed too.Only variables are INitialized here. Arrays are INitialized IN
    ! another SUBROUTINE.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
        IMPLICIT  NONE
        
		file_name_input='input.dat'
        ! This file contaINs all major INput parameters to be fed
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! R  E  A  D  I  N  G       I  N  P  U  T       F  I  L  E
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        OPEN(unit=10,file=TRIM(ADJUSTL(file_name_input)))
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_i4,ADVANCE='yes'),N 	! No of collocation poINts IN physical space IN one DIMENSION
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_d8p4,ADVANCE='yes'),dt  ! Time step for marchINg forward IN time
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_d8p4,ADVANCE='yes'),time_total ! Total time to simulate
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_i4,ADVANCE='yes'),no_of_saves ! No of saves
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_i4,ADVANCE='yes'),purging_status ! PurgINg ON or OFF
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_d8p4,ADVANCE='yes'),purg_alpha ! PurgINg alpha for time to purge 0.<alpha<1.2 
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_d8p4,ADVANCE='yes'),purg_beta ! PurgINg beta for k_P wavenumber cutoff
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_i4,ADVANCE='yes'),viscosity_status ! Viscosity ON or OFF
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_d12p8,ADVANCE='yes'),viscosity ! Viscosity value
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_i4,ADVANCE='yes'),pvd_status ! Paraview files storage ON or OFF
        READ(10,f_d8p4,ADVANCE='yes')
        READ(10,f_i4,ADVANCE='yes'),no_of_pvd_saves ! No of such Paraview files saved
        CLOSE(10)
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! S  P  A  C  E       A  N  D         T  I  M  E  
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        lth=two_pi
		! Length of periodic cube
        
        Nh=N/2
        ! Maximum wavemuber IN the First BrillouIN zone.
        
		dx=lth/DBLE(N)
        dx_3=dx**thr
		! Grid distance, same for 'y' and 'z' directions too.
        
        N3=DBLE(N*N*N)
        vol=(lth**thr)
        ! No of poINts IN real space
        
        k_G=FLOOR(DBLE(N)/thr)
        ! Truncation wavenumber
        
        k_th=k_G-1
        ! initial thermalization wavenumber
        
        k_int=3
        ! integral wavenumber, should be same if different resolutions are
        ! compared.

        k_G_2=DBLE(k_G*k_G)
        
        max_shell_no=FLOOR(DSQRT(thr)*DBLE(k_G))
        ! Truncation wavenumber shell radius. For k_i>k_G, modes are truncated to remove dealiasINg error.
        ! size of cube enclosINg the truncation sphere. This is the size of sphere that enCLOSEs k_G cube.

        no_of_debug=2
        ! No of times that the program looks for any 'NaN' while marchINg forward IN time.
        
        CALL time_to_step_convert(time_total,t_step_total)
        ! returns the no of time_steps (\delta t) IN a given time
        
    	t_step_save=t_step_total/no_of_saves
        ! DetermINes how many time steps after the save has to be made.
        
        t_step_debug=t_step_total/no_of_debug
        ! DetermINes how many time steps after the checkINg has to be made

        t_step_save_pvd=t_step_total/no_of_pvd_saves
        ! DetermINes how many time steps after the checkINg has to be made

        CALL step_to_time_convert(t_step_total,time_total)
        CALL step_to_time_convert(t_step_save,time_save)
        ! Converts the $\delta_t$ to time

        norm_factor=1.0D0
        en_initial=1.0D0
        ! Normalization factor to obtaIN net energy =1.0

        WRITE (N_char,f_i8),N
        ! convertINg value to CHARACTER

        file_dir='data_'//TRIM(ADJUSTL(N_char))//'/'
        file_dir_pvd='data_'//TRIM(ADJUSTL(N_char))//'_3D/'
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! P  U  R  G  I  N  G         V  A  R  I  B  L  E  S
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        IF (purging_status .EQ. 1) THEN
        
            k_P=k_G-FLOOR(DBLE(k_G)**purg_beta)
            ! PurgINg wavenumber
            
            time_purg=DBLE(k_G)**(-purg_alpha)
            ! PurgINg time

            CALL time_to_step_convert(time_purg,t_step_purg)
            ! returns the no of time_steps (\delta t) IN a given time
            
            k_P_2=DBLE(k_P*k_P)
            ! Good programmINg

            purging_count=0
            ! Purging count -  no of times purged

            WRITE (beta_char,f_d5p2),purg_beta
            WRITE (alpha_char,f_d5p2),purg_alpha
            ! convertINg values to CHARACTERs
            
            file_dir='data_'//TRIM(ADJUSTL(N_char))//'_a'//TRIM(ADJUSTL(alpha_char))//'b'// &
            TRIM(ADJUSTL(beta_char))//'_/'
            ! SpecIFiyINg which directory to save data
            
        END IF
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! V  I  S  C  O  S  I  T  Y           V  A  R  I  A  B  L  E  S
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        IF (viscosity_status .EQ. 1) THEN
        
            dissip_rate=0.2D0
            ! GuessINg a tentative dissipation rate on which the k_kolmogorov c
            ! can be calculated and tested
            
            k_kol=FLOOR((dissip_rate/viscosity**thr)**0.25)
            ! kolmogorov viscous scale '1/\eta'

            file_dir='data_viscous'//TRIM(ADJUSTL(N_char))//'/'
    
        END IF
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! A   U   X   I   L   A   R  Y   
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        state_sim=0
        ! MeanINg it is INitializINg  , '1' means fINal, will be changed after the time_marchINg is DOne.

        count_nan=0
        ! No of NaN IN v_x
        
        file_name_param='parameters'
        ! SpecIFiyINg file name to save the parameters of the simulation, File names are given IN the SUBROUTINEs later.

        file_name_en='energy_vs_time'
        ! SPecifying file name where energy evolution with other data will be stored with time.
        
        CALL SYSTEM('mkdir '//TRIM(ADJUSTL(file_dir)))
        ! Creating the file directory, where all the .dat files will be stored

        IF (pvd_status .EQ. 1) then
            CALL SYSTEM('mkdir '//TRIM(ADJUSTL(file_dir_pvd)))
            file_name_pvd_1='quaternion_3D'
        END IF
        
	END
    SUBROUTINE init_arrays
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !   This defines few arrays that do not evolve. k^2 matrix,k matrix, projection matrix,
    ! shell no matrix, count_modes_shell matrix.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT  NONE
        DOUBLE PRECISION::diff_ceiling
        
        shell_no=0
        count_modes_shell=0
        !  +++++++++++++++++++++++++++++++++
        !  A  X  I  S      G   R  I   D   S
        !  +++++++++++++++++++++++++++++++++
        DO i_y=0,N-1
            axis(i_y)=DBLE(i_y)*dx ! Location of grid points along principal axis
        END DO
        
        DO i_x=0,Nh
        DO i_z=-Nh,Nh-1
        DO i_y=-Nh,Nh-1
            k_x(i_x,i_y,i_z)=DBLE(i_x)
            k_y(i_x,i_y,i_z)=DBLE(i_y)
            k_z(i_x,i_y,i_z)=DBLE(i_z)
            ! Just the k component matrix storINg its grid poINts.
            
            k_2(i_x,i_y,i_z)=DBLE(i_x**2+i_y**2+i_z**2)
            ! Square of distance to origIN
            
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  P  R  O  J  E  C  T  I  O  N             M  A  T  R  I  X
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ! Projection matrix \mathbb{P}_{ij}=\delta_{ij}-\frac{k_ik_j}{k^2}}
            IF (k_2(i_x,i_y,i_z).GT.tol) THEN
            ! CheckINg IF k^2 is not too low, to cause NaN (this will happen only for (0,0,0))
                proj_xx(i_x,i_y,i_z)=1.0D0-DBLE(i_x*i_x)/k_2(i_x,i_y,i_z)
                proj_yy(i_x,i_y,i_z)=1.0D0-DBLE(i_y*i_y)/k_2(i_x,i_y,i_z)
                proj_zz(i_x,i_y,i_z)=1.0D0-DBLE(i_z*i_z)/k_2(i_x,i_y,i_z)
                proj_xy(i_x,i_y,i_z)=-DBLE(i_x*i_y)/k_2(i_x,i_y,i_z)
                proj_yz(i_x,i_y,i_z)=-DBLE(i_y*i_z)/k_2(i_x,i_y,i_z)
                proj_zx(i_x,i_y,i_z)=-DBLE(i_z*i_x)/k_2(i_x,i_y,i_z)
            ELSE
                proj_xx(i_x,i_y,i_z)=DBLE(two/thr)
                proj_yy(i_x,i_y,i_z)=DBLE(two/thr)
                proj_zz(i_x,i_y,i_z)=DBLE(two/thr)
                proj_xy(i_x,i_y,i_z)=DBLE(-1.0D0/thr)
                proj_yz(i_x,i_y,i_z)=DBLE(-1.0D0/thr)
                proj_zx(i_x,i_y,i_z)=DBLE(-1.0D0/thr)
            END IF
            
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  T  R  U  N  C  A  T  I  O  N  ,   S  H  E  L  L           M  A  T  R  I  X
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ! Truncation mask matrix (multiply this with any spectral matrix to DO the truncation)
            IF (k_2(i_x,i_y,i_z) .LT. k_G_2)THEN 

                truncator(i_x,i_y,i_z)=1.0D0
                ! Spherical truncation filter matrix

                diff_ceiling=-DSQRT(k_2(i_x,i_y,i_z))+DBLE(CEILING(DSQRT(k_2(i_x,i_y,i_z))))
                ! figuring the decimal part of the |k|, so that it will be alloted to shell k or k+1
                
                IF (diff_ceiling .GE. hf) THEN
                    shell_no(i_x,i_y,i_z)=FLOOR(DSQRT(k_2(i_x,i_y,i_z)))
                ELSE
                    shell_no(i_x,i_y,i_z)=CEILING(DSQRT(k_2(i_x,i_y,i_z)))
                END IF
                
                count_modes_shell(shell_no(i_x,i_y,i_z))=count_modes_shell(shell_no(i_x,i_y,i_z))+1
                ! counts no of grid poINts that belong to a particular shell, it should go as ~s^2
            ELSE
                truncator(i_x,i_y,i_z)=0.0D0
                ! Outised the truncation sphere.
            END IF
        END DO
        END DO
        END DO
            IF (purging_status .EQ. 1) THEN
                CALL init_purger
            END IF
    END
    SUBROUTINE init_initcondn
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !   Initialize a initial condition for velocity. In spectral space
    !   We have exp. decaying I.C, thermalized I.C, Taylor-Green I.C
    ! Kido-Peltz I.C. It has to be made incompressible too.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
        IMPLICIT  NONE
        DO i_x=0,Nh
        DO i_z=-Nh,Nh-1
        DO i_y=-Nh,Nh-1
            ! INitializINg the INitial velocity (spectral) and projectINg it so that the flow is INcompressible.
            CALL IC_exp_decaying_spectrum
!            CALL IC_thermalized_spectrum

            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  V  E  L  O  C  I  T  Y          P  R  O  J  E  C  T  O  R
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ! Assigning the velocity along with projection so that it is INcompressible -- u(k).k=0
            ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  N  O  T  E  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
            ! The following steps ensure, that except for i_x=0 or Nh plane, all the other planes are given INitial velocity
            ! But those planes require special attention, becoz, their INversion lies IN the same plane, so conjugates must be placed accordINgly.
            IF (((i_x .NE. 0) .AND. (i_x .NE. Nh)) .OR. (i_z .LT. 0)) THEN
                v_x(i_x,i_y,i_z)=proj_xx(i_x,i_y,i_z)*temp_v(1)+proj_xy(i_x,i_y,i_z)*temp_v(2)+proj_zx(i_x,i_y,i_z)*temp_v(3)
                v_y(i_x,i_y,i_z)=proj_xy(i_x,i_y,i_z)*temp_v(1)+proj_yy(i_x,i_y,i_z)*temp_v(2)+proj_yz(i_x,i_y,i_z)*temp_v(3)
                v_z(i_x,i_y,i_z)=proj_zx(i_x,i_y,i_z)*temp_v(1)+proj_yz(i_x,i_y,i_z)*temp_v(2)+proj_zz(i_x,i_y,i_z)*temp_v(3)
                    IF (((i_x .EQ. 0) .OR. (i_x .EQ. Nh)) .AND. ((i_y .NE. -Nh) .AND. (i_z .NE. -Nh))) THEN
                        v_x(i_x,-i_y,-i_z)=DCONJG(v_x(i_x,i_y,i_z))
                        v_y(i_x,-i_y,-i_z)=DCONJG(v_y(i_x,i_y,i_z))
                        v_z(i_x,-i_y,-i_z)=DCONJG(v_z(i_x,i_y,i_z))
                    END IF
                ELSE IF ((i_z .EQ. 0) .AND. (i_y .LE. 0)) THEN
                    v_x(i_x,i_y,i_z)=proj_xx(i_x,i_y,i_z)*temp_v(1)+proj_xy(i_x,i_y,i_z)*temp_v(2)+proj_zx(i_x,i_y,i_z)*temp_v(3)
                    v_y(i_x,i_y,i_z)=proj_xy(i_x,i_y,i_z)*temp_v(1)+proj_yy(i_x,i_y,i_z)*temp_v(2)+proj_yz(i_x,i_y,i_z)*temp_v(3)
                    v_z(i_x,i_y,i_z)=proj_zx(i_x,i_y,i_z)*temp_v(1)+proj_yz(i_x,i_y,i_z)*temp_v(2)+proj_zz(i_x,i_y,i_z)*temp_v(3)
                    v_x(i_x,-i_y,i_z)=DCONJG(v_x(i_x,i_y,i_z))
                    v_y(i_x,-i_y,i_z)=DCONJG(v_y(i_x,i_y,i_z))
                    v_z(i_x,-i_y,i_z)=DCONJG(v_z(i_x,i_y,i_z))
                ELSE IF ((i_y .EQ. -Nh) .AND. (i_z .GT. 0)) THEN
                    v_x(i_x,i_y,i_z)=proj_xx(i_x,i_y,i_z)*temp_v(1)+proj_xy(i_x,i_y,i_z)*temp_v(2)+proj_zx(i_x,i_y,i_z)*temp_v(3)
                    v_y(i_x,i_y,i_z)=proj_xy(i_x,i_y,i_z)*temp_v(1)+proj_yy(i_x,i_y,i_z)*temp_v(2)+proj_yz(i_x,i_y,i_z)*temp_v(3)
                    v_z(i_x,i_y,i_z)=proj_zx(i_x,i_y,i_z)*temp_v(1)+proj_yz(i_x,i_y,i_z)*temp_v(2)+proj_zz(i_x,i_y,i_z)*temp_v(3)
            END IF
            ! ----------------------------------------------------------------------------------------
            ! Similar steps has to be followed IN energy calculation too.
        END DO
        END DO
        END DO
!        CALL IC_TG
!        CALL IC_KP
    END
    SUBROUTINE write_parameters
    ! --------------------------------------------------------------------------------------------
    ! PRINTS ALL THE PARAMETERS OF THE SIMULATION IN THE FILE file_name_param
    ! --------------------------------------------------------------------------------------------
     IMPLICIT  NONE
     OPEN(unit=23,file=TRIM(ADJUSTL(file_dir))//TRIM(ADJUSTL(file_name_param))//'.dat')
        WRITE(23,"(A80)")'------------------------------------------'
        WRITE(23,"(A80)")'P  A  R  A  M  E  T  E  R  S  of  Simulation'
        WRITE(23,"(A80)")'------------------------------------------'
        WRITE(23,"(A30,I4)")'N                  = ',N
        WRITE(23,"(A30,I6)")'k_Galerkin.Trunc   = ',k_G
        WRITE(23,"(A30,I6)")'k_Integr           = ',k_int
        WRITE(23,"(A30,ES10.3)")'Time step      = ',dt
        WRITE(23,"(A30,I8)")'Total time steps   = ',t_step_total
        WRITE(23,"(A30,F4.2)")'Total time       = ',time_total
        WRITE(23,"(A60,I2)")'PURGING STATUS     = ',purging_status
        WRITE(23,"(A30,F20.16)")'Purging time   = ',time_purg
        WRITE(23,"(A30,I6)")'k_Purging          = ',k_P
        WRITE(23,"(A60,I2)")'VISCOSITY STATUS   = ',viscosity_status
        WRITE(23,"(A30,F20.16)")'Viscosity      = ',viscosity
        WRITE(23,"(A30,I4)")'No of saves        = ',no_of_saves
        WRITE(23,"(A30,F20.16)")'Initial energy = ',norm_k 
    CLOSE(23)
    END 
    SUBROUTINE init_purger
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !   This decalares the purging matrix. 
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE
!    INTEGER(KIND=8)::purging_tolerance=0.6D0
    ! this is a decimated purging version
! -----------------------------------------------------------------
! THIS PRODUCES NEW RANDOM NUMBERS EVERYTIME CALLED (uniform dist.)
! -----------------------------------------------------------------
!    INTEGER, DIMENSION(8) :: seed_values
!    INTEGER(KIND=4)::seed_size 
!    ! Declare an assumed shape, dynamic array
!    INTEGER, DIMENSION(:),ALLOCATABLE :: seed
!    ! gfortran SUBROUTINE to return date and time INFORMATion 
!    ! from the real time system clock. Works DOwn to milliseconds 
!    ! and stores the eight return values IN array values.
!    CALL date_and_time(VALUES=seed_values)
!    ! restart the state of the pseuDOranDOm number generator
!    ! k = mINimum size of seed (12 on my system)
!    seed_size=20
!    CALL random_seed(size=seed_size)
!    ! ALLOCATE memory to seed
!    ALLOCATE(seed(seed_size))
!    ! assign INFORMATion IN values to seed
!    seed(:) = seed_values(:)
!    ! seed the ranDOm number generator
!    CALL random_seed(put=seed)    
! -----------------------------------------------------------------------
        DO i_x=0,Nh
        DO i_z=-Nh,Nh-1
        DO i_y=-Nh,Nh-1
!             CALL RANDOM_NUMBER(rd)
             IF (k_2(i_x,i_y,i_z) .LE. k_P_2) THEN
                purger(i_x,i_y,i_z)=1.0D0
!            ELSEIF (rd .GT. purgINg_tolerance) THEN
!                purger(i_x,i_y,i_z)=1.0D0
            ELSE
                purger(i_x,i_y,i_z)=0.0D0
            END IF
        END DO
        END DO
        END DO
    END 
    SUBROUTINE IC_exp_decaying_spectrum
    ! A typical initial condition IN which only few large modes are excited.
    ! CALL this to give 3 components of COMPLEX spectral velocity at spectral grid 'i_x,i_y,i_z'
        IMPLICIT  NONE
        DOUBLE PRECISION::phi,theta
        DOUBLE PRECISION,DIMENSION(3)::ph   
        DOUBLE PRECISION::pre_factor,A
        CALL RANDOM_NUMBER(phi)
        CALL RANDOM_NUMBER(theta)
        CALL RANDOM_NUMBER(ph)
        phi=two_pi*phi ! Azimuthal angle of \hat{u}_k vector
        theta=DACOS(1.0D0-two*theta)! Polar angle of \hat{u}_k vector
        ph=two_pi*ph ! Phases of \hat{u}_k components
        A=norm_factor*DSQRT(8.0D0/(thr*DSQRT(two_pi*hf))) ! Normalization factor IN the contINous k limit, but make sure energy is of order O(1)
        pre_factor=DSQRT((k_2(i_x,i_y,i_z))/(two*two_pi))/DBLE(k_INt*k_INt)*(DEXP(-k_2(i_x,i_y,i_z)/DBLE(two*(k_INt**two))))
        ! pre_factor is chosen such the peak of INitial energy happenINg at |k|=k0 is 1.0
        temp_v(1)=A*pre_factor*DSIN(theta)*DCOS(phi)*DCMPLX(DCOS(ph(1)),DSIN(ph(1)))
        temp_v(2)=A*pre_factor*DSIN(theta)*DSIN(phi)*DCMPLX(DCOS(ph(2)),DSIN(ph(2)))
        temp_v(3)=A*pre_factor*DCOS(theta)*DCMPLX(DCOS(ph(3)),DSIN(ph(3)))
        ! 3 COMPLEX values for spectral velocity
    END
    SUBROUTINE IC_thermalized_spectrum
    ! This is a thermalized INitial condition, where all the modes have same velocity.
    ! CALL this to give 3 components of COMPLEX spectral velocity at spectral grid 'i_x,i_y,i_z'
        IMPLICIT  NONE
        DOUBLE PRECISION::phi,theta
        DOUBLE PRECISION,DIMENSION(3)::ph   
        DOUBLE PRECISION::A
        CALL RANDOM_NUMBER(phi)
        CALL RANDOM_NUMBER(theta)
        CALL RANDOM_NUMBER(ph)
        phi=two_pi*phi ! Azimuthal angle of \hat{u}_k vector
        theta=DACOS(1.0D0-two*theta)! Polar angle of \hat{u}_k vector
        ph=two_pi*ph ! Phases of \hat{u}_k 
        A=DBLE(1.0/10000.0)
        ! pre_factor is chosen such the peak of INitial energy happenINg at |k|=k0 is 1.0
        temp_v(1)=norm_factor*A*SIN(theta)*COS(phi)*CMPLX(COS(ph(1)),SIN(ph(1)))
        temp_v(2)=norm_factor*A*SIN(theta)*SIN(phi)*CMPLX(COS(ph(2)),SIN(ph(2)))
        temp_v(3)=norm_factor*A*COS(theta)*CMPLX(COS(ph(3)),SIN(ph(3)))
        ! 3 COMPLEX values for spectral velocity
    END
    SUBROUTINE IC_TG
    ! This is TAYLOR_GREEN Vortex INitial condition which has a lot of symmetries
    IMPLICIT  NONE
    DOUBLE COMPLEX::v1
    v1=i/8.0D0
    v_x=DCMPLX(0.0D0,0.0D0)
    v_y=DCMPLX(0.0D0,0.0D0)
    v_z=DCMPLX(0.0D0,0.0D0)
    v_x(1,1,1)=-v1
    v_x(1,1,-1)=-v1
    v_x(1,-1,1)=-v1
    v_x(1,-1,-1)=-v1
    v_y(1,1,1)=v1
    v_y(1,1,-1)=v1
    v_y(1,-1,1)=-v1
    v_y(1,-1,-1)=-v1
    END
    SUBROUTINE IC_KP
    ! This is Kida Peltz Vortex INitial condition which has a lot of symmetries
    IMPLICIT  NONE
    DOUBLE COMPLEX::v1
    v1=norm_factor*i/8.0D0
    v_x=DCMPLX(0.0D0,0.0D0)
    v_y=DCMPLX(0.0D0,0.0D0)
    v_z=DCMPLX(0.0D0,0.0D0)
    !-------- 'x' velocity--------------
    v_x(1,-1,-3)=v1
    v_x(1,1,-3)=v1
    v_x(1,-3,-1)=-v1
    v_x(1,3,-1)=-v1
    v_x(1,-3,1)=-v1
    v_x(1,3,1)=-v1
    v_x(1,-1,3)=v1
    v_x(1,1,3)=v1
    !-------- 'y' velocity--------------
    v_y(1,-1,-3)=v1
    v_y(1,1,-3)=-v1
    v_y(3,-1,-1)=-v1
    v_y(3,1,-1)=v1
    v_y(3,-1,1)=-v1
    v_y(3,1,1)=v1
    v_y(1,-1,3)=v1
    v_y(1,1,3)=-v1
    !--------- 'z' velocity --------------
    v_z(1,-3,-1)=-v1
    v_z(3,-1,-1)=v1
    v_z(3,1,-1)=v1
    v_z(1,3,-1)=-v1
    v_z(1,-3,1)=v1
    v_z(3,-1,-1)=-v1
    v_z(3,1,1)=-v1
    v_z(1,3,1)=v1
    END
    SUBROUTINE step_to_time_convert(step,time)
	! CALL this to convert time step INto actual time of simulation
		IMPLICIT  NONE
		INTEGER (KIND=4),INTENT(IN)::step
		DOUBLE PRECISION,INTENT(OUT)::time
		time=DBLE(step)*dt
	END
    SUBROUTINE time_to_step_convert(time,step)
	! CALL this to convert time step INto actual time of simulation
		IMPLICIT  NONE
		INTEGER (KIND=4),INTENT(OUT)::step
		DOUBLE PRECISION,INTENT(IN)::time
        step=FLOOR(time/dt)
	END
END MODULE variables_and_arrays
