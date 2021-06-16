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
! MODULE: system_analysis
! LAST MODIFIED: 3 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! RUN AND OUTPUT MODULE TO RUN THE TIME EVOLUTION FOR 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_analysis
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major system_analysis for the Code to run.
! The standard analysis are done in system_functions module
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_variables
  USE system_fftw
  USE system_output
  USE system_functions

  IMPLICIT NONE
  ! _________________________
  ! MOMENTS VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  :: m_ind, m_ind_max, m_ind_min
  DOUBLE PRECISION :: circulation
  DOUBLE PRECISION :: loc_stretching, vx_stretching_max, loc_vx_stretching_max
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::m_arr,m2_inv_arr,alpha_arr,alpha_by_2m_arr
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vx_D_moment,vx_O_moment
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vx_dot_moment
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::loc_vx_dot_moment
  ! DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::bck_vx_dot_moment
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::str_xx,str_yy,str_zz,str_xy,str_yz,str_zx
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_str_xx,bck_str_yy,bck_str_zz,bck_str_xy,bck_str_yz,bck_str_zx
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_str_opr,w_mod_2
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::bck_vx_stretching,vx_stretching

  ! StraIN tensor IN real space

  CONTAINS

  SUBROUTINE allocate_vorticity_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the moments array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    m_ind_min = 1
    m_ind_max = 9

    CALL compute_vorticity

    CALL compute_enstrophy

    circulation = DSQRT( two * enstrophy )
    ! Calculates a time scale fron initial condition

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( m_arr( m_ind_min            : m_ind_max ) )
    ALLOCATE( m2_inv_arr( m_ind_min       : m_ind_max ) )
    ALLOCATE( alpha_arr( m_ind_min        : m_ind_max ) )
    ALLOCATE( alpha_by_2m_arr( m_ind_min  : m_ind_max ) )
    ALLOCATE( vx_O_moment( m_ind_min : m_ind_max ) )
    ALLOCATE( vx_D_moment( m_ind_min : m_ind_max ) )
    ALLOCATE( w_mod_2( 0 : N - 1, 0 : N - 1, 0 : N - 1 ) )

    m_arr =  (/1.0D0, 1.25D0, 1.5D0, 5.0D0 / 3.0D0, 2.0D0, 2.5D0, 3.0D0, 4.0D0, 6.0D0 /)
    file_name = TRIM( ADJUSTL( file_address ) ) // 'moments_exponent_list.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(unit = 8880, file = file_name )

    DO m_ind = m_ind_min, m_ind_max

      ! m_arr( m_ind )           = DBLE( m_ind )
      m2_inv_arr( m_ind )      = one / ( two * m_arr( m_ind ) )
      alpha_arr( m_ind )       = two * m_arr( m_ind ) / ( two * two * m_arr( m_ind ) - thr )
      alpha_by_2m_arr( m_ind ) = one / ( two * two * m_arr( m_ind ) - thr )

      WRITE(8880,f_d8p4,ADVANCE ='yes') m_arr( m_ind )
      ! SAVING THE MOMENT EXPONENTS

    END DO

    CLOSE(8880)

    file_name = TRIM( ADJUSTL( file_address ) ) // 'vorticity_moments_O.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(unit = 8888, file = file_name )

    file_name = TRIM( ADJUSTL( file_address ) ) // 'vorticity_moments_D.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(unit = 8889, file = file_name )

  END

  SUBROUTINE compute_vorticity_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the L-p moments and Rescaled Dimensionless moments of
  ! the vorticity field
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   V  O  R  T  I  C  I  T  Y      M  O  M  E  N  T  S
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    w_mod_2       = w_ux ** two + w_uy ** two + w_uz ** two

    DO m_ind = m_ind_min, m_ind_max

      vx_O_moment( m_ind ) = SUM( w_mod_2 ** m_arr( m_ind ) )
      vx_O_moment( m_ind ) = vx_O_moment( m_ind ) / N3
      vx_D_moment( m_ind ) = vx_O_moment( m_ind ) ** alpha_by_2m_arr( m_ind )
      vx_D_moment( m_ind ) = vx_D_moment( m_ind ) / ( circulation ** alpha_arr( m_ind ) )
      vx_O_moment( m_ind ) = vx_O_moment( m_ind ) ** m2_inv_arr( m_ind )

    END DO

    CALL write_vorticity_moments

  END

  SUBROUTINE write_vorticity_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to write the moments
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

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

  SUBROUTINE allocate_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the strain tensor array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(str_xx(0:N-1,0:N-1,0:N-1),str_yy(0:N-1,0:N-1,0:N-1),str_zz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(str_zx(0:N-1,0:N-1,0:N-1),str_xy(0:N-1,0:N-1,0:N-1),str_yz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(vx_stretching(0:N-1,0:N-1,0:N-1))
    ALLOCATE(vx_dot_moment( m_ind_min : m_ind_max ) )

    sub_dir  =   'section_vx_stretching/'
    ! Sub directory name to store spectral data

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( file_address ) ) // TRIM( ADJUSTL( sub_dir ) ) )

    file_name = TRIM( ADJUSTL( file_address ) ) // 'vorticity_stretching_moments.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(unit = 7777, file = file_name )

  END

  SUBROUTINE allocate_bck_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the bckground strain operator array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::arg, r_local

    r_local = 3.5D0 * dx
    ! Size of the region to integrate and consider as local

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(bck_str_opr(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(bck_str_xx(0:N-1,0:N-1,0:N-1),bck_str_yy(0:N-1,0:N-1,0:N-1),bck_str_zz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(bck_str_zx(0:N-1,0:N-1,0:N-1),bck_str_xy(0:N-1,0:N-1,0:N-1),bck_str_yz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(bck_vx_stretching(0:N-1,0:N-1,0:N-1))
    ALLOCATE(loc_vx_dot_moment( m_ind_min : m_ind_max ) )
    ! ALLOCATE(bck_vx_dot_moment( m_ind_min : m_ind_max ) )

    DO i_x = 0, Nh
    DO i_y = -Nh, Nh - 1
    DO i_z = -Nh, Nh - 1

      arg                 = r_local * DSQRT( k_2( i_x, i_y, i_z ) )

      bck_str_opr( i_x, i_y, i_z ) = thr * ( DSIN( arg ) - arg * DCOS( arg ) ) /  (arg ** thr)

    END DO
    END DO
    END DO

    bck_str_opr( 0, 0, 0 ) = one

    file_name = TRIM( ADJUSTL( file_address ) ) // 'vorticity_stretching_loc_moments.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! OPEN(unit = 7778, file = file_name )

    ! file_name = TRIM( ADJUSTL( file_address ) ) // 'vorticity_stretching_bck_moments.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(unit = 7779, file = file_name )

  END

  SUBROUTINE compute_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the strain tensor array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  T  R  A  I  N        T  E  N  S  O  R        C  A  L  C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL fft_c2r(i*k_x*v_x,hf*i*(k_y*v_x+k_x*v_y),i*k_z*v_z,N,Nh,str_xx,str_xy,str_zz)
    CALL fft_c2r(i*k_y*v_y,hf*i*(k_y*v_z+k_z*v_y),hf*i*(k_x*v_z+k_z*v_x),N,Nh,str_yy,str_yz,str_zx)
    ! All six components of strain tensor.

  END

  SUBROUTINE compute_bck_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the background strain tensor array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  T  R  A  I  N        T  E  N  S  O  R        C  A  L  C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL fft_c2r(i*bck_str_opr*k_x*v_x,hf*i*bck_str_opr*(k_y*v_x+k_x*v_y),i*bck_str_opr*k_z*v_z,N,Nh,&
                  bck_str_xx,bck_str_xy,bck_str_zz)
    CALL fft_c2r(i*bck_str_opr*k_y*v_y,hf*i*bck_str_opr*(k_y*v_z+k_z*v_y),hf*i*bck_str_opr*(k_x*v_z+k_z*v_x),N,Nh,&
                  bck_str_yy,bck_str_yz,bck_str_zx)
    ! All six components of strain tensor.

  END

  SUBROUTINE compute_vortex_stretching
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the scalar for vortex stretching
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   V  O  R  T  E  X        S  T  R  E  T  C  H  I  N  G.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    vx_stretching = ( str_xx * w_ux * w_ux + str_yy * w_uy * w_uy + str_zz * w_uz * w_uz + two * &
                    ( str_xy * w_ux * w_uy + str_yz * w_uy * w_uz + str_zx * w_ux * w_uz ) ) / w_mod_2


  END

  SUBROUTINE compute_bck_vortex_stretching
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the scalar for background vortex stretching
  ! local vortex stretcing can be obtained by subracting this from vortex stretching
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! ( BCK ) V  O  R  T  E  X        S  T  R  E  T  C  H  I  N  G.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    bck_vx_stretching = ( bck_str_xx * w_ux * w_ux + bck_str_yy * w_uy * w_uy + bck_str_zz * w_uz * w_uz + two * &
                        ( bck_str_xy * w_ux * w_uy + bck_str_yz * w_uy * w_uz + bck_str_zx * w_ux * w_uz ) ) / w_mod_2


  END

  SUBROUTINE compute_vorticity_dot_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the L-p moments and Rescaled Dimensionless moments of
  ! the derivatives of moments in terms of the stretching of vorticity
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   V  O  R  T  I  C  I  T  Y      M  O  M  E  N  T  S
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO m_ind = m_ind_min, m_ind_max

      vx_dot_moment( m_ind ) = SUM( ( w_mod_2 ** m_arr( m_ind ) ) * vx_stretching )
      vx_dot_moment( m_ind ) = vx_dot_moment( m_ind ) / ( N3 * ( vx_O_moment( m_ind ) ** ( two * m_arr( m_ind ) ) ) )

    END DO

    CALL write_vorticity_dot_moments

    CALL write_vx_stretching_section
    ! writes the section file

  END

  SUBROUTINE write_vorticity_dot_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to write the moments for the vorticity_dot
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(7777,f_d8p4,ADVANCE   ='no')  time_now

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

    SUBROUTINE compute_vorticity_dot_loc_moments
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL this to get the L-p moments and Rescaled Dimensionless moments of
    ! the derivatives of moments in terms of the local stretching of vorticity
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IMPLICIT NONE

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !   V  O  R  T  I  C  I  T  Y      M  O  M  E  N  T  S
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO m_ind = m_ind_min, m_ind_max

        ! bck_vx_dot_moment( m_ind ) = SUM( ( w_mod_2 ** m_arr( m_ind ) ) * bck_vx_stretching )
        ! bck_vx_dot_moment( m_ind ) = bck_vx_dot_moment( m_ind ) / ( N3 * ( vx_O_moment( m_ind ) ** ( two * m_arr( m_ind ) ) ) )
        loc_vx_dot_moment( m_ind ) = SUM( ( w_mod_2 ** m_arr( m_ind ) ) * ( vx_stretching - bck_vx_stretching ) )
        loc_vx_dot_moment( m_ind ) = loc_vx_dot_moment( m_ind ) / ( N3 * ( vx_O_moment( m_ind ) ** ( two * m_arr( m_ind ) ) ) )

      END DO

      CALL write_vorticity_dot_loc_moments

      CALL write_loc_vx_stretching_section
      ! writes the section file

    END

  SUBROUTINE write_vorticity_dot_loc_moments
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to write the moments for the vorticity_dot
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

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

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir ) ) &
                // 'vx_stretching_section_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 788, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    vx_stretching_max = MAXVAL( DABS( vx_stretching ) )

       i_y = Nh - 1
    DO i_x = 0, N - 1
    DO i_z = 0, N - 1
         ! WRITE(788,'(F32.17)',ADVANCE='no') vx_stretching(i_x,i_y,i_z) / vx_dot_moment( 1 )
         ! Normalized by the vortex stretching moment weighted with enstrophy (m=1 moment)
         WRITE(788,'(F32.17)',ADVANCE='yes') vx_stretching(i_x,i_y,i_z) / vx_stretching_max
         ! Normalized by the infinite norm of the vx stretching factor.
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

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir ) ) &
                // 'loc_vx_stretching_section_t_'//TRIM( ADJUSTL( file_time ) ) // '.dat'

    OPEN( UNIT = 789, FILE = file_name)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE-section
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    loc_vx_stretching_max = MAXVAL( DABS( vx_stretching - bck_vx_stretching ) )

       i_y = Nh - 1
    DO i_x = 0, N - 1
    DO i_z = 0, N - 1
         loc_stretching = vx_stretching(i_x,i_y,i_z) - bck_vx_stretching(i_x,i_y,i_z)
         ! WRITE(789,'(F32.17)',ADVANCE='no')  loc_stretching / loc_vx_dot_moment( 1 )
         ! WRITE(789,'(F32.17)',ADVANCE='no')  loc_stretching / vx_dot_moment( 1 )
         ! WRITE(789,'(F32.17)',ADVANCE='no')  loc_stretching / vx_stretching_max
         WRITE(789,'(F32.17)',ADVANCE='yes') loc_stretching / loc_vx_stretching_max
    END DO
    END DO

    CLOSE(789)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

	SUBROUTINE deallocate_vorticity_moments
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(w_mod_2)

	END

	SUBROUTINE deallocate_strain_tensor
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(str_xx,str_yy,str_zz)
		DEALLOCATE(str_xy,str_yz,str_zx)
		DEALLOCATE(vx_stretching)

	END

	SUBROUTINE deallocate_bck_strain_tensor
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE(bck_str_xx,bck_str_yy,bck_str_zz)
    DEALLOCATE(bck_str_xy,bck_str_yz,bck_str_zx)
		DEALLOCATE(bck_vx_stretching)
		DEALLOCATE(bck_str_opr)

	END

END MODULE system_analysis
