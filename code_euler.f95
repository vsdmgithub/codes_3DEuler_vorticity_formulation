module parameters
! All the parameters for the simulation space are declared here
	implicit none
	double precision::lth,vol,vcy,dx,dvol
	double precision::dt,time_total,time_save
	integer (kind=4)::ln2_N,N,N_z
	integer (kind=4)::t_step_total,t_step_save
	integer (kind=4)::i_x,i_y,i_z
	double complex(parameter)::iota
	double precision(parameter)::two_pi
	
	contains 
	subroutine init_space_parameters
		implicit none
		two_pi=atan(1.0)*8.0
		iota=COMPLEX(0.0,1.0)
		! Defining the value of 2\pi and iota
		lth=two_pi
		vol=(two_pi)**(3.0)
		! Length and Volume of periodic box (cube)
		ln2_N=4
		N=2**(ln2_N)
		! No of collocation points in physical space.
		N_z=N/2+1
		! No of collocation points in Fourier space
		dx=two_pi/DFLOAT(N)
		dvol=dx**(3.0)
		! Grid distance
	end subroutine init_space_parameters
    subroutine init_time_parameters
		implicit none
		dt=0.0001
		t_step_total=10**6
		t_step_save=100
		step_to_time_convert(t_step_total,time_total)
		step_to_time_convert(t_step_save,time_save)
		! Time steps, No of time steps,Total time to simulate,time to save
	end subroutine init_time_parameters
	subroutine init_fluid_parameters
		implicit none
		vcy=0.0001
		! Viscosity of Fluid
	end subroutine init_fluid_parameters
	subroutine step_to_time_convert(step,time)
	! Converts time step into actual time
		implicit none
		integer (kind=4),intent(in)::step
		double precision,intent(out)::time
		time=DFLOAT(step)*dt
	end subroutine step_to_time_convert
end module parameters

module initial_condition
! Declares the initial condition for the simulation in fourier space
	use parameters
	implicit none
	double complex,dimension(:,:,:),allocatable::v_x,v_y,v_z
	double precision,dimension(3)::temp
	double precision::phi,theta
	double precision::pwr_k,k0,pre_factor,k_mod
	contains
    subroutine init_initial_velocity
    ! Exponent in the initial condition k^(pwr_k)
	pwr_k=2.0
	! Shell radius,till the energy is concentrated initially
	k0=SQRT(3.0)*DFLOAT(N_half)/(10.0)
	do i_x=1,N
	do i_y=1,N
	do i_z=1,N_z
		call initial_velocity_function
		v_x(i_x,i_y,i_z)=temp(1)
		v_y(i_x,i_y,i_z)=temp(2)
		v_z(i_x,i_y,i_z)=temp(3)
	end do
	end do
	end do
	end subroutine init_initial_velocity
	subroutine initial_velocity_function
	call random_number(phi)
	call random_number(theta)
	phi=two_pi*phi
	theta=two_pi*0.5*theta
	call l2_norm
    pre_factor=(k_mod**(pwr_k))*(exp(-(k_mod**(2.0))/(2.0*k0**(2.0)))
	temp(1)=pre_factor*SIN(theta)*COS(phi)
	temp(2)=pre_factor*SIN(theta)*SIN(phi)
	temp(3)=pre_factor*COS(theta)
	end subroutine initial_velocity_function
	subroutine l2_norm
	implicit none
	k_mod=SQRT(DFLOAT(i_x**2+i_y**2+i_z**2))
	end subroutine l2_norm
end module initial_condition

module solver
! takes the initial condition and solves till time 't' 
	use initial_condition
	use fft
	implicit none
	double complex,dimension(:,:,:),allocatable::u_x,u_y,u_z
	integer(kind=4)::t_step
	
	contains
	subroutine time_forward
	implicit none
	do t_step=1,t_step_total
		call evolve_algorithm
		! Saving output		
		if (MOD(t_step,t_step_save)==0) then
			call save_file
		end if
	end do
	end subroutine time_forward
	subroutine evolve_algorithm
	end subroutine evolve_algorithm
	subroutine evolve_step
	implicit none
	call convection_term_dft_r2c
	call dissipation_term
	v_x_dot=v_k_x_diss+u_grad_u_k_x
	v_y_dot=v_k_x_diss+u_grad_u_k_x
	v_z_dot=v_k_x_diss+u_grad_u_k_x
	end subroutine evolve_step
	subroutine convection_term
	end subroutine convection_term
	subroutine dissipation_term
	end subroutine dissipation_term
end module solver

module fft
	implicit none
	use 'fftw.h'
	contains
	subroutine fft
	end subroutine fft
end module

module output

end module output
program euler
	use solver
	implicit none
	! Initialize all the parameters for the simulation
	call init_space_parameters
	call init_time_parameters
	call init_fluid_parameters
	allocate(u_x(N,N,N),u_y(N,N,N),u_z(N,N,N),v_x(N,N,N_z),v_y(N,N,N_z),v_z(N,N,N_z))
	! Call and initialize initial velocity in fourier space
	call init_initial_velocity
	! Call solver module
	call solver
	
end program euler
