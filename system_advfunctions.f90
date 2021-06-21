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
! LAST MODIFIED: 21 JUNE 2021
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
    ! REF-> <<< system_advoutput >>>

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
    ! REF-> <<< system_advoutput >>>

    CALL write_vx_stretching_section
    ! REF-> <<< system_advoutput >>>

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

      loc_vx_dot_moment( m_ind ) = SUM( ( w_mod_2 ** m_arr( m_ind ) ) * ( vx_stretching - bck_vx_stretching ) )
      loc_vx_dot_moment( m_ind ) = loc_vx_dot_moment( m_ind ) / ( N3 * ( vx_O_moment( m_ind ) ** ( two * m_arr( m_ind ) ) ) )

    END DO

    CALL write_vorticity_dot_loc_moments
    ! REF-> <<< system_advoutput >>>

    CALL write_loc_vx_stretching_section
    ! REF-> <<< system_advoutput >>>

  END

END MODULE system_advfunctions
