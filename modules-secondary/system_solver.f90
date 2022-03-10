! <f Stamp
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
! MODULE: system_solver
! LAST MODIFIED: 22 JAN 2022
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER MODULE TO SOLVE 3D EULER EQUATIONS IN VORTICITY FORMULATION
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! </f>

MODULE system_solver
! <f Info
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Takes the spectral vorticity and moves it one step forward in euler equation using the given algorithm, uses FFTW.
! This has RK4, AB4 algorithm.
! The spectral equation is
! dw_i(k)/dt = w_j . \grad_j u_i + u_j . \grad_k w_i
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! </f>

! <f Glob Dec
	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_advfunctions

	IMPLICIT NONE
	! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::adv_w_ux,adv_w_uy,adv_w_uz
  ! Real advection term matrix
	DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::sth_u_x,sth_u_y,sth_u_z
	! Real stretching term matrix

  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::grad_x,grad_y,grad_z
  ! Calculates gradient in real space from FFT

	! _________________________________________
  ! FOURIER SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::adv_w_vx,adv_w_vy,adv_w_vz
	! Spectral advection term matrix
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::sth_v_x,sth_v_y,sth_v_z
	! Spectral stretching term matrix

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::w_vx_dot_m1,w_vy_dot_m1,w_vz_dot_m1
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::w_vx_dot_m2,w_vy_dot_m2,w_vz_dot_m2
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::w_vx_dot_m3,w_vy_dot_m3,w_vz_dot_m3
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::w_vx_dot,w_vy_dot,w_vz_dot
	! Spectral derivatives for AB4

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::w_vx_pred,w_vy_pred,w_vz_pred
  ! Spectral velocity predictor for AB4 matrix

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::dw1_x,dw2_x,dw3_x,dw4_x,dw1_y,dw2_y
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::dw3_y,dw4_y,dw1_z,dw2_z,dw3_z,dw4_z
  ! Intermediate matrices for RK4 algorithm

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::w_vx_temp,w_vy_temp,w_vz_temp
  ! temporary matrices to store velocities during RK4 algorithm

! </f>
	CONTAINS

	SUBROUTINE allocate_solver_main
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate main arrays used in solving
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		ALLOCATE( adv_w_ux( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
		ALLOCATE( adv_w_uy( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
		ALLOCATE( adv_w_uz( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
		ALLOCATE( sth_u_x( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
		ALLOCATE( sth_u_y( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
		ALLOCATE( sth_u_z( 0 : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
		ALLOCATE( grad_x( 0  : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
		ALLOCATE( grad_y( 0  : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
		ALLOCATE( grad_z( 0  : N_x - 1, 0 : N_y - 1, 0 : N_z - 1 ) )
		ALLOCATE( adv_w_vx( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( adv_w_vy( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( adv_w_vz( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( sth_v_x( kMin_x  : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( sth_v_y( kMin_x  : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( sth_v_z( kMin_x  : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )

	END
	! </f>
! --------------------------------------------------------------

	SUBROUTINE allocate_solver_RK4
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays for RK4 algoritm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(dw1_x( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw2_x( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw3_x( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw4_x( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw1_y( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw2_y( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw3_y( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw4_y( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw1_z( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw2_z( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw3_z( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(dw4_z( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(w_vx_temp( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(w_vy_temp( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE(w_vz_temp( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )

	END
	! </f>

	SUBROUTINE solver_RK4_algorithm
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
	! Alg: - Runga kutta 4th order
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'v(k)'
		w_vx_temp = w_vx
		w_vy_temp = w_vy
		w_vz_temp = w_vz
		CALL time_increment_RK1() ! This call provides \vec{dv} for the existing \vec{v}
		w_vx      = w_vx_temp + hf * dw1_x
		w_vy      = w_vy_temp + hf * dw1_y
		w_vz      = w_vz_temp + hf * dw1_z
		CALL time_increment_RK2()
		w_vx      = w_vx_temp + hf * dw2_x
		w_vy      = w_vy_temp + hf * dw2_y
		w_vz      = w_vz_temp + hf * dw2_z
		CALL time_increment_RK3()
		w_vx      = w_vx_temp + dw3_x
		w_vy      = w_vy_temp + dw3_y
		w_vz      = w_vz_temp + dw3_z
		CALL time_increment_RK4()
		! Final increment for 'v(k)'
		w_vx      = w_vx_temp + ( dw1_x + two * dw2_x + two * dw3_x + dw4_x ) / six
		w_vy      = w_vy_temp + ( dw1_y + two * dw2_y + two * dw3_y + dw4_y ) / six
		w_vz      = w_vz_temp + ( dw1_z + two * dw2_z + two * dw3_z + dw4_z ) / six

	END
	! </f>

	SUBROUTINE time_increment_RK1()
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'w(k)'
	! This is the EULER EQUATION implemented for numerical computation
	! spectral space in vorticity formulation.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		CALL advection
		CALL stretching
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R        E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla w' and the stretching term 'w\cdot \Nabla u ' in spectral space
		dw1_x = dt * truncator * ( sth_v_x - adv_w_vx )
		dw1_y = dt * truncator * ( sth_v_y - adv_w_vy )
		dw1_z = dt * truncator * ( sth_v_z - adv_w_vz )

	END
	! </f>

	SUBROUTINE time_increment_RK2()
	! <f
		IMPLICIT NONE
		CALL advection
		CALL stretching
		dw2_x = dt * truncator * ( sth_v_x - adv_w_vx )
		dw2_y = dt * truncator * ( sth_v_y - adv_w_vy )
		dw2_z = dt * truncator * ( sth_v_z - adv_w_vz )
	END
	! </f>

	SUBROUTINE time_increment_RK3()
	! <f
		IMPLICIT NONE
		CALL advection
		CALL stretching
		dw3_x = dt * truncator * ( sth_v_x - adv_w_vx )
		dw3_y = dt * truncator * ( sth_v_y - adv_w_vy )
		dw3_z = dt * truncator * ( sth_v_z - adv_w_vz )
	END
	! </f>

	SUBROUTINE time_increment_RK4()
	! <f
		IMPLICIT NONE
		CALL advection
		CALL stretching
		dw4_x = dt * truncator * ( sth_v_x - adv_w_vx )
		dw4_y = dt * truncator * ( sth_v_y - adv_w_vy )
		dw4_z = dt * truncator * ( sth_v_z - adv_w_vz )
	END
	! </f>

	SUBROUTINE advection
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to give spectral advection term using v_k.
	! 1. First w 			--> v  velocity from vorticity
	! 2. Next v 			--> u  FFT is done
	! 3. Next i*k*w 	--> dw/dx  FFT is done
	! 4. Next u.dw/dx --> Fourier(u.dw/dx)
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   A   D   V   E   C   T   I   O   N       T   E   R   M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		! First getting the spectral velocity from spectral vorticity
		CALL compute_velocity

		! FFT spectral to real velocity
		CALL fft_c2r_vec( v_x, v_y, v_z, u_x, u_y, u_z )

		! w_ux gradient in real space
		CALL fft_c2r_vec( i * k_x * w_vx, i * k_y * w_vx, i * k_z * w_vx, grad_x, grad_y, grad_z )

		! u.Nabla(w) term in x direction
		adv_w_ux = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! w_uy gradient in real space
		CALL fft_c2r_vec( i * k_x * w_vy, i * k_y * w_vy, i * k_z * w_vy, grad_x, grad_y, grad_z )

		! u.Nabla(w) term in z direction
		adv_w_uy = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! w_vz gradient in real space
		CALL fft_c2r_vec( i * k_x * w_vz, i * k_y * w_vz, i * k_z * w_vz, grad_x, grad_y, grad_z )

		! u.Nabla(w) term in z direction
		adv_w_uz = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! Calculate the advection term in spectral space by doing iFFT
		CALL fft_r2c_vec( adv_w_ux, adv_w_uy, adv_w_uz, adv_w_vx, adv_w_vy, adv_w_vz )

  END
	! </f>

	SUBROUTINE stretching
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to compute spectral vortex stretching term using w_k.
	! Here, we use the Strain tensor to calculate the stretching.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  S  T  R  E  T  C  H  I  N  G       T  E  R  M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		! IF ( t_step .GT. 2 ) THEN
		!
		! 	CALL compute_filtered_strain_tensor
	  !   ! REF-> <<< system_advfunctions >>>
		!
		! ELSE
		!
		! 	CALL compute_strain_tensor
	  !   ! REF-> <<< system_advfunctions >>>
		!
		! END IF

		CALL compute_strain_tensor
    ! REF-> <<< system_advfunctions >>>

		! First getting the real vorticity
		CALL fft_c2r_vec( w_vx, w_vy, w_vz, w_ux, w_uy, w_uz )

		sth_u_x = str_xx * w_ux + str_xy * w_uy + str_zx * w_uz
		sth_u_y = str_xy * w_ux + str_yy * w_uy + str_yz * w_uz
		sth_u_z = str_zx * w_ux + str_yz * w_uy + str_zz * w_uz

		! Calculate the stretching term in spectral space by doing iFFT
		CALL fft_r2c_vec( sth_u_x, sth_u_y, sth_u_z, sth_v_x, sth_v_y, sth_v_z )

  END
	! </f>

	SUBROUTINE deallocate_solver_RK4
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE( dw1_x, dw2_x )
		DEALLOCATE( dw3_x, dw4_x )
		DEALLOCATE( dw1_y, dw2_y )
		DEALLOCATE( dw3_y, dw4_y )
		DEALLOCATE( dw1_z, dw2_z )
		DEALLOCATE( dw3_z, dw4_z )
		DEALLOCATE( w_vx_temp, w_vy_temp, w_vz_temp)

	END
	! </f>
! --------------------------------------------------------------
	SUBROUTINE allocate_solver_AB4
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays for AB4 algorithm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		ALLOCATE( w_vx_dot( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vy_dot( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vz_dot( kMin_x    : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vx_dot_m1( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vy_dot_m1( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vz_dot_m1( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vx_dot_m2( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vy_dot_m2( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vz_dot_m2( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vx_dot_m3( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vy_dot_m3( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vz_dot_m3( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vx_pred( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vy_pred( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
		ALLOCATE( w_vz_pred( kMin_x   : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )

	END
	! </f>

	SUBROUTINE solver_AB4_algorithm
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE AB4 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
	! Alg: - Adam Bashforth 4th Order algorithm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		IF ( t_step .EQ. 0 ) THEN

			! Using initial condition to get -3 step value
			CALL time_derivative_AB()
			w_vx_dot_m3 = w_vx_dot
			w_vy_dot_m3 = w_vy_dot
			w_vz_dot_m3 = w_vz_dot

			CALL allocate_solver_RK4
			! Allocating RK4 for first 3 steps

			CALL solver_RK4_algorithm

		ELSE IF ( t_step .EQ. 1 ) THEN

			! Filling up for -2 step value
			CALL time_derivative_AB()
			w_vx_dot_m2 = w_vx_dot
			w_vy_dot_m2 = w_vy_dot
			w_vz_dot_m2 = w_vz_dot

			CALL solver_RK4_algorithm

		ELSE IF ( t_step .EQ. 2 ) THEN

			! Filling up for -1 step value
			CALL time_derivative_AB()
			w_vx_dot_m1 = w_vx_dot
			w_vy_dot_m1 = w_vy_dot
			w_vz_dot_m1 = w_vz_dot

			CALL solver_RK4_algorithm

			CALL deallocate_solver_RK4
			! No need for RK4 anymore

		ELSE

			CALL time_derivative_AB()

			w_vx_pred = w_vx + dt * ( - 9.0D0 * w_vx_dot_m3 + 37.0D0 * w_vx_dot_m2 - 59.0D0 * w_vx_dot_m1 + 55.0D0 * w_vx_dot ) / 24.0D0
			w_vy_pred = w_vy + dt * ( - 9.0D0 * w_vy_dot_m3 + 37.0D0 * w_vy_dot_m2 - 59.0D0 * w_vy_dot_m1 + 55.0D0 * w_vy_dot ) / 24.0D0
			w_vz_pred = w_vz + dt * ( - 9.0D0 * w_vz_dot_m3 + 37.0D0 * w_vz_dot_m2 - 59.0D0 * w_vz_dot_m1 + 55.0D0 * w_vz_dot ) / 24.0D0

			CALL time_derivative_AB_pred()

			w_vx      = w_vx + dt * ( w_vx_dot_m2 - 5.0D0 * w_vx_dot_m1 + 19.0D0 * w_vx_dot + 9.0D0 * w_vx_dot_m3 ) / 24.0D0
			w_vy      = w_vy + dt * ( w_vy_dot_m2 - 5.0D0 * w_vy_dot_m1 + 19.0D0 * w_vy_dot + 9.0D0 * w_vy_dot_m3 ) / 24.0D0
			w_vz      = w_vz + dt * ( w_vz_dot_m2 - 5.0D0 * w_vz_dot_m1 + 19.0D0 * w_vz_dot + 9.0D0 * w_vz_dot_m3 ) / 24.0D0
			! Predicted 'v_dot' is stored in 'v_dot_m3' - to save space :)

			! Shiting the known velocities for next step
			w_vx_dot_m3 = w_vx_dot_m2
			w_vx_dot_m2 = w_vx_dot_m1
			w_vx_dot_m1 = w_vx_dot

			w_vy_dot_m3 = w_vy_dot_m2
			w_vy_dot_m2 = w_vy_dot_m1
			w_vy_dot_m1 = w_vy_dot

			w_vz_dot_m3 = w_vz_dot_m2
			w_vz_dot_m2 = w_vz_dot_m1
			w_vz_dot_m1 = w_vz_dot

		END IF

	END
	! </f>

	SUBROUTINE time_derivative_AB()
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the EULER EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   A   D   V   E   C   T   I   O   N       T   E   R   M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		! First getting the spectral velocity from spectral vorticity
		v_x = i * ( k_y * w_vz - k_z * w_vy ) / k_2
		v_y = i * ( k_z * w_vx - k_x * w_vz ) / k_2
		v_z = i * ( k_x * w_vy - k_y * w_vx ) / k_2

		! FFT spectral to real velocity
		CALL fft_c2r_vec( v_x, v_y, v_z, u_x, u_y, u_z )

		! w_ux gradient in real space
		CALL fft_c2r_vec( i * k_x * w_vx, i * k_y * w_vx, i * k_z * w_vx, grad_x, grad_y, grad_z )

		! u.Nabla(w) term in x direction
		adv_w_ux = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! w_uy gradient in real space
		CALL fft_c2r_vec( i * k_x * w_vy, i * k_y * w_vy, i * k_z * w_vy, grad_x, grad_y, grad_z )

		! u.Nabla(w) term in z direction
		adv_w_uy = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! w_vz gradient in real space
		CALL fft_c2r_vec( i * k_x * w_vz, i * k_y * w_vz, i * k_z * w_vz, grad_x, grad_y, grad_z )

		! u.Nabla(w) term in z direction
		adv_w_uz = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! Calculate the advection term in spectral space by doing iFFT
		CALL fft_r2c_vec( adv_w_ux, adv_w_uy, adv_w_uz, adv_w_vx, adv_w_vy, adv_w_vz )

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  S  T  R  E  T  C  H  I  N  G       T  E  R  M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		CALL compute_strain_tensor
    ! REF-> <<< system_advfunctions >>>

		! First getting the real vorticity
		CALL fft_c2r_vec( w_vx, w_vy, w_vz, w_ux, w_uy, w_uz )

		sth_u_x = str_xx * w_ux + str_xy * w_uy + str_zx * w_uz
		sth_u_y = str_xy * w_ux + str_yy * w_uy + str_yz * w_uz
		sth_u_z = str_zx * w_ux + str_yz * w_uy + str_zz * w_uz

		! Calculate the stretching term in spectral space by doing iFFT
		CALL fft_r2c_vec( sth_u_x, sth_u_y, sth_u_z, sth_v_x, sth_v_y, sth_v_z )

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla w' and the stretching term 'w\cdot \Nabla u ' in spectral space

		w_vx_dot = truncator * ( sth_v_x - adv_w_vx )
		w_vy_dot = truncator * ( sth_v_y - adv_w_vy )
		w_vz_dot = truncator * ( sth_v_z - adv_w_vz )

	END
	! </f>

	SUBROUTINE time_derivative_AB_pred()
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the EULER EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   A   D   V   E   C   T   I   O   N       T   E   R   M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		! First getting the spectral velocity from spectral vorticity
		v_x = i * ( k_y * w_vz_pred - k_z * w_vy_pred ) / k_2
		v_y = i * ( k_z * w_vx_pred - k_x * w_vz_pred ) / k_2
		v_z = i * ( k_x * w_vy_pred - k_y * w_vx_pred ) / k_2

		! FFT spectral to real velocity
		CALL fft_c2r_vec( v_x, v_y, v_z, u_x, u_y, u_z )

		! w_ux gradient in real space
		CALL fft_c2r_vec( i * k_x * w_vx_pred, i * k_y * w_vx_pred, i * k_z * w_vx_pred, grad_x, grad_y, grad_z )

		! u.Nabla(w) term in x direction
		adv_w_ux = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! w_uy gradient in real space
		CALL fft_c2r_vec( i * k_x * w_vy_pred, i * k_y * w_vy_pred, i * k_z * w_vy_pred, grad_x, grad_y, grad_z )

		! u.Nabla(w) term in z direction
		adv_w_uy = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! w_vz gradient in real space
		CALL fft_c2r_vec( i * k_x * w_vz_pred, i * k_y * w_vz_pred, i * k_z * w_vz_pred, grad_x, grad_y, grad_z )

		! u.Nabla(w) term in z direction
		adv_w_uz = ( u_x * grad_x + u_y * grad_y + u_z * grad_z )

		! Calculate the advection term in spectral space by doing iFFT
		CALL fft_r2c_vec( adv_w_ux, adv_w_uy, adv_w_uz, adv_w_vx, adv_w_vy, adv_w_vz )

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!  S  T  R  E  T  C  H  I  N  G       T  E  R  M
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		CALL compute_strain_tensor
    ! REF-> <<< system_advfunctions >>>

		! First getting the real vorticity
		CALL fft_c2r_vec( w_vx_pred, w_vy_pred, w_vz_pred, w_ux, w_uy, w_uz )

		sth_u_x = str_xx * w_ux + str_xy * w_uy + str_zx * w_uz
		sth_u_y = str_xy * w_ux + str_yy * w_uy + str_yz * w_uz
		sth_u_z = str_zx * w_ux + str_yz * w_uy + str_zz * w_uz

		! Calculate the stretching term in spectral space by doing iFFT
		CALL fft_r2c_vec( sth_u_x, sth_u_y, sth_u_z, sth_v_x, sth_v_y, sth_v_z )

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla w' and the stretching term 'w\cdot \Nabla u ' in spectral space
		w_vx_dot_m3 = truncator * ( sth_v_x - adv_w_vx )
		w_vy_dot_m3 = truncator * ( sth_v_y - adv_w_vy )
		w_vz_dot_m3 = truncator * ( sth_v_z - adv_w_vz )

	END
	! </f>

	SUBROUTINE deallocate_solver_AB4
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE( w_vx_pred, w_vy_pred, w_vz_pred )
		DEALLOCATE( w_vx_dot_m1, w_vx_dot_m2, w_vx_dot_m3 )
		DEALLOCATE( w_vy_dot_m1, w_vy_dot_m2, w_vy_dot_m3 )
		DEALLOCATE( w_vz_dot_m1, w_vz_dot_m2, w_vz_dot_m3 )

	END
	! </f>
! --------------------------------------------------------------
	SUBROUTINE deallocate_solver_main
	! <f
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE( adv_w_vx, adv_w_vy, adv_w_vz )
		DEALLOCATE( adv_w_ux, adv_w_uy, adv_w_uz )
		DEALLOCATE( sth_v_x, sth_v_y, sth_v_z )
		DEALLOCATE( sth_u_x, sth_u_y, sth_u_z )
		DEALLOCATE( grad_x, grad_y, grad_z )

	END
	! </f>

 END MODULE system_solver
