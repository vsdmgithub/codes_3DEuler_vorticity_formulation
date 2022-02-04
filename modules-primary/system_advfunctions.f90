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

    CALL write_linedata("lin_unf_xy",str_xy(:,0,0))
    CALL write_linedata("lin_unf_xx",str_xx(:,0,0))

    CALL compute_bck_strain_tensor

    var_str  = SUM( str_xx ** ( two * th_coeff_exponent ) ) / N3
    th_coeff = DERF( ( str_xx ** ( two * th_coeff_exponent ) / var_str ) )
    str_xx   = bck_str_xx + th_coeff * ( str_xx - bck_str_xx )

    var_str  = SUM( str_yy ** ( two * th_coeff_exponent ) ) / N3
    th_coeff = DERF( ( str_yy ** ( two * th_coeff_exponent ) / var_str ) )
    str_yy   = bck_str_yy + th_coeff * ( str_yy - bck_str_yy )

    var_str  = SUM( str_zz ** ( two * th_coeff_exponent ) ) / N3
    th_coeff = DERF( ( str_zz ** ( two * th_coeff_exponent ) / var_str ) )
    str_zz   = bck_str_zz + th_coeff * ( str_zz - bck_str_zz )

    var_str  = SUM( str_xy ** ( two * th_coeff_exponent ) ) / N3
    th_coeff = DERF( ( str_xy ** ( two * th_coeff_exponent ) / var_str ) )
    str_xy   = bck_str_xy + th_coeff * ( str_xy - bck_str_xy )

    var_str  = SUM( str_yz ** ( two * th_coeff_exponent ) ) / N3
    th_coeff = DERF( ( str_yz ** ( two * th_coeff_exponent ) / var_str ) )
    str_yz   = bck_str_yz + th_coeff * ( str_yz - bck_str_yz )

    var_str  = SUM( str_zx ** ( two * th_coeff_exponent ) ) / N3
    th_coeff = DERF( ( str_zx ** ( two * th_coeff_exponent ) / var_str ) )
    str_zx   = bck_str_zx + th_coeff * ( str_zx - bck_str_zx )

    CALL write_linedata("lin_fil_xy",str_xy(:,0,0))
    CALL write_linedata("lin_fil_xx",str_xx(:,0,0))

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
