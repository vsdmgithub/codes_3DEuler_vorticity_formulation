module parameters
! All the parameters for the simulation space are declared here
	implicit none
	double precision::lth,vol,vcy,dx,dvol
	double precision::dt,time_total,time_save
	integer (kind=4)::ln2_N,N,N_half
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
		N_half=N/2
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
! Declares the initial condition for the simulation
	use parameters
	implicit none
	double precision,dimension(:,:,:),allocatable::u_r
	double precision,dimension(:,:,:),allocatable::u_k_x,u_k_y,u_k_z
	double precision::phi,theta
	double precision::pwr_k,k0,pre_factor
	contains
	subroutine init_initial_velocity
	pwr_k=2.0
	! Exponent in the initial condition k^(pwr_k)
	k0=SQRT(3.0)*DFLOAT(N_half)/(10.0)
	! Shell radius,till the energy is concentrated initially
	do i_x=1,N_half
		do i_y=1,N_half
			do i_z=1,N_half
				call random_number(phi)
				call random_number(theta)
				phi=two_pi*phi
				theta=two_pi*theta
				call distance(i_x,i_y,i_z,k_mod)
				pre_factor=(k_mod**(pwr_k))*(exp(-(k_mod**(2.0))/(2.0*k0**(2.0)))
				u_k_x(i_x,i_y,i_z)=pre_factor*SIN(theta)*COS(phi)
				u_k_y(i_x,i_y,i_z)=pre_factor*SIN(theta)*SIN(phi)
				u_k_z(i_x,i_y,i_z)=pre_factor*COS(theta))
			end do
		end do
	end do
	end subroutine init_initial_velocity
end module initial_condition

program euler
	use initial_condition
	implicit none
	call init_space_parameters
	call init_time_parameters
	call init_fluid_parameters
	allocate(u_k_x(N_half,N_half,N_half),u_k_y(N_half,N_half,N_half),u_k_z(N_half,N_half,N_half))
	call init_initial_velocity
end program euler
