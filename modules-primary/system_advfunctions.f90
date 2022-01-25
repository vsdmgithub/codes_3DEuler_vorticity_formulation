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

MODULE system_advfunctions
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major advanced functions involving analysis.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicfunctions
  USE system_advvariables
  USE system_advoutput

  IMPLICIT NONE

  CONTAINS

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
    CALL fft_c2r_vec(i*k_x*v_x,hf*i*(k_y*v_x+k_x*v_y),i*k_z*v_z,str_xx,str_xy,str_zz)
    CALL fft_c2r_vec(i*k_y*v_y,hf*i*(k_y*v_z+k_z*v_y),hf*i*(k_x*v_z+k_z*v_x),str_yy,str_yz,str_zx)
    ! All six components of strain tensor.

    w_mod_2 = w_ux ** two + w_uy ** two + w_uz ** two

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
    CALL fft_c2r_vec(i*bck_str_opr*k_x*v_x,hf*i*bck_str_opr*(k_y*v_x+k_x*v_y),i*bck_str_opr*k_z*v_z,&
                  bck_str_xx,bck_str_xy,bck_str_zz)
    CALL fft_c2r_vec(i*bck_str_opr*k_y*v_y,hf*i*bck_str_opr*(k_y*v_z+k_z*v_y),hf*i*bck_str_opr*(k_x*v_z+k_z*v_x),&
                  bck_str_yy,bck_str_yz,bck_str_zx)
    ! All six components of bckground strain tensor.

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
    ! Rate of vortex stretcing scalar

    CALL write_vx_stretching_section
    ! REF-> <<< system_advoutput >>>

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
    ! Rate of vortex stretching scalar due to background strain,

    CALL write_loc_vx_stretching_section
    ! REF-> <<< system_advoutput >>>

  END

  SUBROUTINE compute_energy_filter
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

END MODULE system_advfunctions
