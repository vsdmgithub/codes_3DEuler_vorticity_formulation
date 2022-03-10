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
! MODULE: system_advfunctions
! LAST MODIFIED: 22 JAN 2022
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED FUNCTIONS MODULE FOR 3D EULER ANALYSIS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_advfunctions
! <f Info
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major advanced functions involving analysis.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! </f>

! <f Glob Dec
  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicfunctions
  USE system_advvariables
  USE system_advoutput

  IMPLICIT NONE

! </f>
  CONTAINS

  SUBROUTINE compute_strain_tensor
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the strain tensor array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  T  R  A  I  N        T  E  N  S  O  R        C  A  L  C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL fft_c2r_vec(i*k_x*v_x,hf*i*(k_y*v_x+k_x*v_y),i*k_z*v_z,str_xx,str_xy,str_zz)
    CALL fft_c2r_vec(i*k_y*v_y,hf*i*(k_y*v_z+k_z*v_y),hf*i*(k_x*v_z+k_z*v_x),str_yy,str_yz,str_zx)
    ! All six components of strain tensor.

    w_mod_2 = w_ux ** two + w_uy ** two + w_uz ** two

  END
! </f>

  SUBROUTINE compute_bck_strain_tensor
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the background strain tensor array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  T  R  A  I  N        T  E  N  S  O  R        C  A  L  C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL fft_c2r_vec(i*bck_str_opr*k_x*v_x,hf*i*bck_str_opr*(k_y*v_x+k_x*v_y),i*bck_str_opr*k_z*v_z,&
                  bck_str_xx,bck_str_xy,bck_str_zz)
    CALL fft_c2r_vec(i*bck_str_opr*k_y*v_y,hf*i*bck_str_opr*(k_y*v_z+k_z*v_y),hf*i*bck_str_opr*(k_x*v_z+k_z*v_x),&
                  bck_str_yy,bck_str_yz,bck_str_zx)
    ! All six components of bckground strain tensor.

  END
! </f>

  SUBROUTINE compute_filtered_strain_tensor
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the filtered strain tensor array
  ! based on the decomposition of strain and connecting the local and nonlocal strain
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  F  I  L  T  E  R  I  N  G      S  T  R  A  I  N        C  A  L  C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  1. First compute the actual strain tensor
    !  2. Then, decompose it
    !  3. Compute the thermalisation coefficient
    !  4. Combine the local and non-local strain using the thermalisation coefficient

    CALL compute_strain_tensor

    ! CALL write_linedata("lin_unf_xy",str_xy(:,0,0))
    ! CALL write_linedata("lin_unf_xx",str_xx(:,0,0))

    CALL compute_bck_strain_tensor

    CALL write_linedata("BS_xy",bck_str_xy(:,0,0))
    CALL write_linedata("LS_xy",str_xy(:,0,0)-bck_str_xy(:,0,0))
    CALL write_linedata("S_xy",str_xy(:,0,0))

    th_coeff =        str_xx ** two + str_yy ** two + str_zz ** two + &
              two * ( str_xy ** two + str_yz ** two + str_zx ** two )

    ! th_coeff = th_coeff ** ( th_coeff_exponent / two )

    th_norm  = SUM( th_coeff ) / N3
    ! th_norm  = th_norm ** ( two / th_coeff_exponent )

    th_coeff = hf + ( th_coeff / th_norm ) ** th_coeff_exponent
    th_coeff = DERF( th_coeff )
    ! th_coeff = 0.5D0

    CALL write_linedata("f",th_coeff(:,0,0))

    str_xx   = bck_str_xx + th_coeff * ( str_xx - bck_str_xx )
    str_yy   = bck_str_yy + th_coeff * ( str_yy - bck_str_yy )
    str_zz   = bck_str_zz + th_coeff * ( str_zz - bck_str_zz )
    str_xy   = bck_str_xy + th_coeff * ( str_xy - bck_str_xy )
    str_yz   = bck_str_yz + th_coeff * ( str_yz - bck_str_yz )
    str_zx   = bck_str_zx + th_coeff * ( str_zx - bck_str_zx )

    CALL write_linedata("Sf_xy",str_xy(:,0,0))

    ! CALL write_linedata("lin_fil_xy",str_xy(:,0,0))
    ! CALL write_linedata("lin_w_z",w_uz(:,0,0))
    ! CALL write_linedata("lin_fil_xx",str_xx(:,0,0))

  END
! </f>

  SUBROUTINE purge_vorticity
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the filtered strain tensor and generate velocity
  ! and vorticity from it.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    ! DOUBLE PRECISION::kx,ky,kz
    ! DOUBLE COMPLEX  ::sxx,syy,szz,sxy,syz,szx
    ! DOUBLE COMPLEX  ::wx,wy,wz
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   ::str_k_xx,str_k_yy,str_k_zz
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   ::str_k_xy,str_k_yz,str_k_zx
    ALLOCATE( str_k_xx( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( str_k_yy( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( str_k_zz( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( str_k_xy( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( str_k_yz( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    ALLOCATE( str_k_zx( kMin_x : kMax_x, kMin_y : kMax_y, kMin_z : kMax_z ) )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  U  R  G  I  N  G
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! CALL write_linedata("unf_vor_z",w_uz(:,0,0))
    ! CALL write_linedata("unf_vel_z",u_y(:,0,0))

    CALL compute_filtered_strain_tensor

    CALL fft_r2c_vec(str_xx,str_yy,str_zz,str_k_xx,str_k_yy,str_k_zz)
    CALL fft_r2c_vec(str_xy,str_yz,str_zx,str_k_xy,str_k_yz,str_k_zx)
    ! Computing the fourier transform of strain

    ! v_x = -i * two * ( k_x * str_k_xx + k_y * str_k_xy + k_z * str_k_zx ) / k_2
    ! v_y = -i * two * ( k_x * str_k_xy + k_y * str_k_yy + k_z * str_k_yz ) / k_2
    ! v_z = -i * two * ( k_x * str_k_zx + k_y * str_k_yz + k_z * str_k_zz ) / k_2

    w_vx = two * ( k_x*k_y*str_k_zx+k_y*k_y*str_k_yz+k_y*k_z*str_k_zz - &
                  k_x*k_z*str_k_xy+k_y*k_z*str_k_yy+k_z*k_z*str_k_yz ) / k_2
    w_vy = two * ( k_x*k_z*str_k_xx+k_y*k_z*str_k_xy+k_z*k_z*str_k_zx - &
                  k_x*k_x*str_k_zx+k_x*k_y*str_k_yz+k_z*k_x*str_k_zz ) / k_2
    w_vz = two * ( k_x*k_x*str_k_xy+k_x*k_y*str_k_yy+k_x*k_z*str_k_yz - &
                  k_x*k_y*str_k_xx+k_y*k_y*str_k_xy+k_y*k_z*str_k_zx ) / k_2

    DEALLOCATE(str_k_xx,str_k_yy,str_k_zz,str_k_xy,str_k_yz,str_k_zx)

    w_vx( 0, 0, 0 ) = c_zero
    w_vy( 0, 0, 0 ) = c_zero
    w_vz( 0, 0, 0 ) = c_zero

    ! CALL compute_vorticity
    ! ! REF-> <<< system_basicfunctions >>>
    !
    CALL write_linedata("W1_z",w_uz(:,0,0))

    CALL compute_velocity
    ! ! REF-> <<< system_basicfunctions >>>

    ! CALL compute_projected_velocity
    ! REF-> <<< system_basicfunctions >>>

    CALL normalize_velocity
    ! REF-> <<< system_basicfunctions >>>

    CALL compute_vorticity
    ! REF-> <<< system_basicfunctions >>>

    ! CALL write_linedata("fil_vel_z",u_y(:,0,0))

    CALL write_linedata("W2_z",w_uz(:,0,0))

  END
! </f>

  SUBROUTINE compute_vortex_stretching
! <f
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
    ! Rate of vortex stretcing scalar

    ! CALL write_vx_stretching_section
    ! REF-> <<< system_advoutput >>>

  END
! </f>

  SUBROUTINE compute_bck_vortex_stretching
! <f
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
    ! Rate of vortex stretching scalar due to background strain,

    ! CALL write_loc_vx_stretching_section
    ! REF-> <<< system_advoutput >>>

  END
! </f>

  SUBROUTINE compute_energy_filter
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute energy in the truncation wavenumber filter modes
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! E N E R G Y     I N     F I L T E R
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    energy_filter = SUM( ( CDABS( v_x ) ** two + CDABS( v_y ) ** two + CDABS( v_z ) ** two ) * tr_wave_filter )

    energy_filter_spectral = SUM( spectral_energy( k_P : ) )

    CALL write_energy_filter
    ! REF-> <<< system_advoutput >>>

  END
! </f>

END MODULE system_advfunctions
