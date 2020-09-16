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
! MODULE: solver
! LAST MODIFIED: 29 JUNE 2020
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER MODULE TO SOLVE 3D EULER EQUATIONS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE solver
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Takes the spectral velocity and moves it one step forward in euler equation using the given algorithm, uses FFTW.
! This has RK4, RK2, AB2, AB4 algorithm. Additionally subroutines 'compute_spectral_energy', 'compute_real_energy', 'compute_compressibility'
! gives energy in spectral, real space and gives the output of \vec{k}\cdot \vec{v}.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	USE variables_and_arrays
	USE fft
	IMPLICIT NONE
    CONTAINS
	SUBROUTINE evolve_algorithm_rk4
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
    ! Alg: - Runga kutta 4th order
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        ! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'v(k)'
        v_x_temp=truncator*v_x
		v_y_temp=truncator*v_y
        v_z_temp=truncator*v_z
        CALL time_derivative(dv1_x,dv1_y,dv1_z) ! This call provides \vec{dv} for the existing \vec{v}
		v_x=v_x_temp+hf*dv1_x
		v_y=v_y_temp+hf*dv1_y
        v_z=v_z_temp+hf*dv1_z
        CALL time_derivative(dv2_x,dv2_y,dv2_z)
		v_x=v_x_temp+hf*dv2_x
		v_y=v_y_temp+hf*dv2_y
        v_z=v_z_temp+hf*dv2_z
        CALL time_derivative(dv3_x,dv3_y,dv3_z)
		v_x=v_x_temp+dv3_x
		v_y=v_y_temp+dv3_y
        v_z=v_z_temp+dv3_z
        CALL time_derivative(dv4_x,dv4_y,dv4_z)
        ! Final increment for 'v(k)'
		v_x=v_x_temp+(dv1_x+two*dv2_x+two*dv3_x+dv4_x)/six
		v_y=v_y_temp+(dv1_y+two*dv2_y+two*dv3_y+dv4_y)/six
        v_z=v_z_temp+(dv1_z+two*dv2_z+two*dv3_z+dv4_z)/six
    END 
	SUBROUTINE evolve_algorithm_rk2
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL this to USE RK2 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
    ! Alg: - Runga kutta 2th order
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        ! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'v(k)'
        v_x_temp=truncator*v_x
		v_y_temp=truncator*v_y
        v_z_temp=truncator*v_z
        CALL time_derivative(dv1_x,dv1_y,dv1_z) ! This call provides \vec{dv} for the existing \vec{v}
		v_x=v_x_temp+hf*dv1_x
		v_y=v_y_temp+hf*dv1_y
        v_z=v_z_temp+hf*dv1_z
        CALL time_derivative(dv2_x,dv2_y,dv2_z)
        ! Final increment for 'v(k)'
		v_x=v_x_temp+(dv1_x+dv2_x)/two
		v_y=v_y_temp+(dv1_y+dv2_y)/two
        v_z=v_z_temp+(dv1_z+dv2_z)/two
    END 
    SUBROUTINE time_derivative(dv_x,dv_y,dv_z)
        ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! ------------
        ! CALL this to get the time derivative matrix for matrix 'v(k)'
        ! This is the EULER EQUATION implemented for numerical computation
        ! spectral space.
        ! -------------
        ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        DOUBLE complex,dimension(0:Nh,-Nh:Nh-1,-Nh:Nh-1),intent(out)::dv_x,dv_y,dv_z
        CALL convection_term
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   3   D  -   E   U   L   E   R           E   Q   N.
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! Get the convection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.
        ! Viscous term is also optional
        IF (viscosity_status .EQ. 1) THEN
            dv_x=-dt*(truncator*(proj_xx*conv_v_x+proj_xy*conv_v_y+proj_zx*conv_v_z)+viscosity*k_2*v_x)
            dv_y=-dt*(truncator*(proj_xy*conv_v_x+proj_yy*conv_v_y+proj_yz*conv_v_z)+viscosity*k_2*v_y)
            dv_z=-dt*(truncator*(proj_zx*conv_v_x+proj_yz*conv_v_y+proj_zz*conv_v_z)+viscosity*k_2*v_z)
        ELSE
            dv_x=-dt*truncator*(proj_xx*conv_v_x+proj_xy*conv_v_y+proj_zx*conv_v_z)
            dv_y=-dt*truncator*(proj_xy*conv_v_x+proj_yy*conv_v_y+proj_yz*conv_v_z)
            dv_z=-dt*truncator*(proj_zx*conv_v_x+proj_yz*conv_v_y+proj_zz*conv_v_z)
        END IF
	END
    SUBROUTINE convection_term
        ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! ------------
        ! CALL this to give spectral convection term using v_k.
        ! 1. First v --> u  FFT is done
        ! 2. Next i*k*v --> du/dx  FFT is done
        ! 3. Next u.du/dx --> Fourier(u.du/dx)
        ! -------------
        ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   C   O   N   V   E   C   T   I   O   N       T   E   R   M
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! First FFT spectral to real velocity
        CALL fft_c2r(v_x,v_y,v_z,N,Nh,u_x,u_y,u_z)
        ! u_x gradient in real space
        CALL fft_c2r(i*k_x*v_x,i*k_y*v_x,i*k_z*v_x,N,Nh,grad_x,grad_y,grad_z)
        ! u.Nabla(u) term in x direction
        conv_u_x=(u_x*grad_x+u_y*grad_y+u_z*grad_z)
        ! u_y gradient in real space
        CALL fft_c2r(i*k_x*v_y,i*k_y*v_y,i*k_z*v_y,N,Nh,grad_x,grad_y,grad_z)
        ! u.Nabla(u) term in z direction
        conv_u_y=(u_x*grad_x+u_y*grad_y+u_z*grad_z)
        ! u_y gradient in real space
        CALL fft_c2r(i*k_x*v_z,i*k_y*v_z,i*k_z*v_z,N,Nh,grad_x,grad_y,grad_z)
        ! u.Nabla(u) term in z direction
        conv_u_z=(u_x*grad_x+u_y*grad_y+u_z*grad_z)
        ! Calculate the convection term in spectral space by doing iFFT
        CALL fft_r2c(conv_u_x,conv_u_y,conv_u_z,N,Nh,conv_v_x,conv_v_y,conv_v_z)
    END
    SUBROUTINE compute_spectral_energy
        ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! ------------
        ! CALL this to check the presence of NaN in your spectral velocity data (v(k)),
        ! and also the L2 norm or the Kinetic energy.
        ! NOTE: Count certain modes once, certain modes half (owing to 1/2 factor)
        ! in the first loop i_x=0 plane is left. later it is considered 
        ! -------------
        ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        norm_k=0.0D0
        DO i_z=-Nh,Nh-1
        DO i_y=-Nh,Nh-1
        DO i_x=1,Nh-1
        norm_k=norm_k+CDABS(v_x(i_x,i_y,i_z))**two+CDABS(v_y(i_x,i_y,i_z))**two+CDABS(v_z(i_x,i_y,i_z))**two
        IF (v_x(i_x,i_y,i_z).NE.v_x(i_x,i_y,i_z)) then
            count_nan=count_nan+1
        END IF
        END DO
        END DO
        END DO
        i_x=0
        DO i_z=-Nh,Nh-1
        DO i_y=-Nh,Nh-1
        norm_k=norm_k+hf*(CDABS(v_x(i_x,i_y,i_z))**two+CDABS(v_y(i_x,i_y,i_z))**two+CDABS(v_z(i_x,i_y,i_z))**two)
        END DO
        END DO
        i_x=Nh
        DO i_z=-Nh,Nh-1
        DO i_y=-Nh,Nh-1
        norm_k=norm_k+hf*(CDABS(v_x(i_x,i_y,i_z))**two+CDABS(v_y(i_x,i_y,i_z))**two+CDABS(v_z(i_x,i_y,i_z))**two)
        END DO
        END DO
    END
    SUBROUTINE compute_real_energy
        ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! ------------
        ! CALL this to get the Kinetic energy in real space
        ! -------------
        ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        norm_x=hf*SUM(u_x**two+u_y**two+u_z**two)/N3
    END
    SUBROUTINE compute_helicity
        ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! ------------
        ! CALL this find helicity of he given configuration
        ! -------------
        ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        ALLOCATE(w_ux(0:N-1,0:N-1,0:N-1),w_uy(0:N-1,0:N-1,0:N-1),w_uz(0:N-1,0:N-1,0:N-1))
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  V  O  R  T  I  C  I  T  Y       C  A  L  C.
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL fft_c2r(i*(k_y*v_z-k_z*v_y),i*(k_z*v_x-k_x*v_z),i*(k_x*v_y-k_y*v_x),N,Nh,w_ux,w_uy,w_uz)
        ! Spectral to Real Vorticity
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  H  E  L  I  C  I  T  Y         C  A  L  C.
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        hely=SUM(u_x*w_ux+u_y*w_uy+u_z*w_uz)
        hely=hely/N3
        DEALLOCATE(w_ux,w_uy,w_uz)
    END
    SUBROUTINE compute_compressibility
        ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! ------------
        ! CALL this to check incompressibility condition. Sums over all residues
        ! of incompressibility and prints it. Of order 10^(-12).
        ! -------------
        ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        k_dot_v=0.0D0
        DO i_z=-Nh,Nh-1
        DO i_y=-Nh,Nh-1
        DO i_x=0,Nh
            k_dot_v=k_dot_v+DABS(i_x*DREAL(v_x(i_x,i_y,i_z))+i_y*DREAL(v_y(i_x,i_y,i_z))+i_z*DREAL(v_z(i_x,i_y,i_z)))
            k_dot_v=k_dot_v+DABS(i_x*DIMAG(v_x(i_x,i_y,i_z))+i_y*DIMAG(v_y(i_x,i_y,i_z))+i_z*DIMAG(v_z(i_x,i_y,i_z)))
        END DO
        END DO
        END DO
    END
 END MODULE solver
