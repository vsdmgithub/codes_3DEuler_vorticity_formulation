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

  SUBROUTINE compute_thermalised_vorticity_field
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the background strain tensor array
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  L O C A L    V O R T I C I T Y    C A L C
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL fft_c2r_vec(loc_opr*w_vx,loc_opr*w_vy,loc_opr*w_vz,w_ux_loc,w_uy_loc,w_uz_loc)
    ! Local vorticity field

    w_mod_2      = w_ux ** two + w_uy ** two + w_uz ** two
    CALL write_linedata("lin_es",w_mod_2(:,0,0))
    CALL write_section_xy("sec_es",w_mod_2(:,:,0))

    w_k_mod_2    = CDABS( w_vx ) ** two + CDABS( w_vy ) ** two + CDABS( w_vz ) ** two
    w_k_mod_2    = w_k_mod_2 * th_filter
    th_threshold = DBLE( N3 ) * MAXVAL( w_k_mod_2 )
    ! Finding a threshold level of the oscillations
    print*,th_threshold,MAXVAL( w_k_mod_2 )

    w_loc_mod_2  = w_ux_loc ** two + w_uy_loc ** two + w_uz_loc ** two
    CALL write_linedata("lin_es_loc",w_loc_mod_2(:,0,0))
    CALL write_section_xy("sec_es_loc",w_loc_mod_2(:,:,0))

    w_fil_mod_2  = w_loc_mod_2 * hf * ( one + DERF( th_threshold - w_k_mod_2 ) )
    CALL write_linedata("lin_es_fil",w_fil_mod_2(:,0,0))
    CALL write_section_xy("sec_es_fil",w_fil_mod_2(:,:,0))

  END
! </f>

END MODULE system_advfunctions
