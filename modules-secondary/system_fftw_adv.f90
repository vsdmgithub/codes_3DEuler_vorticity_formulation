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

! #########################
! MODULE: system_fftw_adv
! LAST MODIFIED: 5 JULY 2021.
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! FFTW MODULE ADVANCED
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_fftw_adv
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This module is to compute the forward and backward DFT for vectors
! and scalars
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! ----------------------
  ! HEADER FILES/MODULES INCLUSION
  ! ----------------------
  USE,INTRINSIC::ISO_C_BINDING
  ! STANDARD MODULE WHICH DEFINES THE EQUIVALENT OF C TYPES IN FORTRAN

  IMPLICIT NONE
  INCLUDE 'fftw3.f03'
  ! FORTRAN INTERFACE FILES FOR ALL OF THE C ROUTINES FOR FFTW OPERATION

  ! _________________________
  !  VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! -----------------------
  ! C VARIABLES DECLARATION
  ! -----------------------
  INTEGER(C_INT)::Z_x,Z_y,Z_z
  INTEGER(C_INT)::q_x,q_y,q_z
  INTEGER(C_INT)::qMax_x,qMax_y,qMax_z
  INTEGER(C_INT)::qMin_x,qMin_y,qMin_z
  INTEGER(C_INT)::real_pts,fourier_pts
  DOUBLE PRECISION::Z_factor

  TYPE(C_PTR)::plan_r2c,plan_c2r
  TYPE(C_PTR)::cdata_r2c_in,cdata_r2c_out
  TYPE(C_PTR)::cdata_c2r_in,cdata_c2r_out
  ! ALL FFTW PLANS ARE OF THIS DATATYPE IN FORTRAN

  COMPLEX(C_DOUBLE_COMPLEX),POINTER::data_r2c_out(:,:,:)
  COMPLEX(C_DOUBLE_COMPLEX),POINTER::data_c2r_in(:,:,:)
  ! FOURIER SPACE DATA

  REAL(C_DOUBLE),POINTER::data_c2r_out(:,:,:)
  REAL(C_DOUBLE),POINTER::data_r2c_in(:,:,:)
  ! REAL SPACE DATA

  ! -----------------------------------------------
  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  ! INFO: WE HAVE FOUR SUBROUTINES
  ! -----------------------------------------------
  ! 1. fft_r2c_vec - Vector DFT Forward
  ! 2. fft_cr2_vec - Vector DFT Backward
  ! 3. fft_r2c - Scalar DFT Forward
  ! 4. fft_cr2 - Scalar DFT Backward
  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

  CONTAINS

  SUBROUTINE init_fft_size( siz_x, siz_y, siz_z )
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS AT THE BEGINNNING OF PROG TO
  ! INITIALIZE THE SIZES FOR FFT COMPUTATION
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4),INTENT(IN)::siz_x,siz_y,siz_z

    ! ---------------------------
    ! ALLOTING THE DOMAIN SIZES
    ! ---------------------------
    Z_x         = siz_x
    Z_y         = siz_y
    Z_z         = siz_z

    ! ---------------------------
    ! MAXIMUM WAVENUMBERS ALLOWED
    ! ---------------------------
    qMin_x      = + 0
    qMax_x      = + INT( Z_x / 2 )
    qMin_y      = - INT( Z_y / 2 )
    qMax_y      = + INT( Z_y / 2 ) - 1
    qMin_z      = - INT( Z_z / 2 )
    qMax_z      = + INT( Z_z / 2 ) - 1

    real_pts    = Z_x * Z_y * Z_z
    ! NO OF COLLACATION PTS IN REAL SPACE

    fourier_pts = ( ( Z_x / 2 ) + 1 ) * Z_y * Z_z
    ! NO OF COLLACATION PTS IN FOURIER SPACE

    Z_factor = DBLE( real_pts )
    ! NORMALIZATION FACTOR IN THE FFT

END

SUBROUTINE fft_r2c_vec( in_x, in_y, in_z, out_x, out_y, out_z )
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! CALL THIS WITH REAL INPUT ARRAY 'IN'
! AND GET SPECTRAL OUTPUT ARRAY 'OUT'
! ------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,DIMENSION( 0 : Z_x - 1, 0 : Z_y - 1, 0 : Z_z - 1 ),INTENT(IN)::in_x,in_y,in_z
    DOUBLE COMPLEX,DIMENSION( qMin_x : qMax_x, qMin_y : qMax_y, qMin_z : qMax_z ),INTENT(OUT)::out_x,out_y,out_z

    ! ---------------------------------------
    ! ALLOCATE ARRAYS - DYNAMIC
    ! NOTE:- The array dimensions are in reverse order for FORTRAN
    ! ---------------------------------------
    cdata_r2c_in  = FFTW_ALLOC_REAL( INT( real_pts, C_SIZE_T ) )
    CALL C_F_POINTER( cdata_r2c_in, data_r2c_in,[ Z_x, Z_y, Z_z ] )

    cdata_r2c_out = FFTW_ALLOC_COMPLEX( INT( fourier_pts, C_SIZE_T ) )
    CALL C_F_POINTER( cdata_r2c_out, data_r2c_out,[ ( Z_x / 2 + 1 ), Z_y, Z_z ] )

    ! -----------------------------------
    ! PLAN FOR OUT-PLACE FORWARD DFT R2C
    ! -----------------------------------
    plan_r2c      = FFTW_PLAN_DFT_R2C_3D( Z_z, Z_y, Z_x, data_r2c_in, data_r2c_out, FFTW_ESTIMATE)
    ! Note the *** REVERSE *** order of the array (unlike in C)

    !_____________________________________________________________
    ! 'x' comp
    ! ============================================================

    ! ---------------------
    ! INITIALIZE INPUT DATA
    ! ---------------------
    data_r2c_in = in_x

    ! -----------
    ! EXECUTE DFT
    ! -----------
    CALL FFTW_EXECUTE_DFT_R2C( plan_r2c, data_r2c_in, data_r2c_out )

    ! -------------------------------------------------------
    ! WRITE OUTPUT (in format of first Brillouin zone format)
    ! -------------------------------------------------------
    DO q_x = qMin_x, qMax_x ! Along 'q_x' direction

    DO q_y = qMin_y, - 1
    DO q_z = qMin_z, - 1

            out_x( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + Z_y + 1, q_z + Z_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = 0, qMax_z

            out_x( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + 1, q_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = qMin_z, - 1

            out_x( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + 1, q_z + Z_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = qMin_y, - 1
    DO q_z = 0, qMax_z

            out_x( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + Z_y + 1, q_z + 1 ) / Z_factor

    END DO
    END DO

    END DO

    !_____________________________________________________________
    ! 'y' comp
    ! ============================================================

    ! ---------------------
    ! INITIALIZE INPUT DATA
    ! ---------------------
    data_r2c_in = in_y

    ! -----------
    ! EXECUTE DFT
    ! -----------
    CALL FFTW_EXECUTE_DFT_R2C( plan_r2c, data_r2c_in, data_r2c_out )

    ! -------------------------------------------------------
    ! WRITE OUTPUT (in format of first Brillouin zone format)
    ! -------------------------------------------------------
    DO q_x = qMin_x, qMax_x ! Along 'q_x' direction

    DO q_y = qMin_y, - 1
    DO q_z = qMin_z, - 1

            out_y( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + Z_y + 1, q_z + Z_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = 0, qMax_z

            out_y( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + 1, q_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = qMin_z, - 1

            out_y( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + 1, q_z + Z_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = qMin_y, - 1
    DO q_z = 0, qMax_z

            out_y( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + Z_y + 1, q_z + 1 ) / Z_factor

    END DO
    END DO

    END DO

    !_____________________________________________________________
    ! 'z' comp
    ! ============================================================

    ! ---------------------
    ! INITIALIZE INPUT DATA
    ! ---------------------
    data_r2c_in = in_z

    ! -----------
    ! EXECUTE DFT
    ! -----------
    CALL FFTW_EXECUTE_DFT_R2C( plan_r2c, data_r2c_in, data_r2c_out )

    ! -------------------------------------------------------
    ! WRITE OUTPUT (in format of first Brillouin zone format)
    ! -------------------------------------------------------
    DO q_x = qMin_x, qMax_x ! Along 'q_x' direction

    DO q_y = qMin_y, - 1
    DO q_z = qMin_z, - 1

            out_z( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + Z_y + 1, q_z + Z_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = 0, qMax_z

            out_z( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + 1, q_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = qMin_z, - 1

            out_z( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + 1, q_z + Z_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = qMin_y, - 1
    DO q_z = 0, qMax_z

            out_z( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + Z_y + 1, q_z + 1 ) / Z_factor

    END DO
    END DO

    END DO

    ! -------------
    ! DESTROY PLANS
    ! -------------
    CALL FFTW_DESTROY_PLAN( plan_r2c )
    CALL FFTW_FREE( cdata_r2c_in )
    CALL FFTW_FREE( cdata_r2c_out )

  END

  SUBROUTINE fft_c2r_vec( in_x, in_y, in_z, out_x, out_y, out_z )
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS WITH COMPLEX INPUT ARRAY 'IN'
  ! AND GET REAL OUTPUT ARRAY 'OUT'
  ! ------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX,DIMENSION( qMin_x : qMax_x, qMin_y : qMax_y, qMin_z : qMax_z ),INTENT(IN)::in_x,in_y,in_z
    DOUBLE PRECISION,DIMENSION( 0 : Z_x - 1, 0 : Z_y - 1, 0 : Z_z - 1 ),INTENT(OUT)::out_x,out_y,out_z

    ! ---------------------------------------
    ! ALLOCATE ARRAYS - DYNAMIC
    ! NOTE:- The array dimensions are in reverse order for FORTRAN
    ! ---------------------------------------
    cdata_c2r_in  = FFTW_ALLOC_COMPLEX( INT( fourier_pts, C_SIZE_T ) )
    CALL C_F_POINTER( cdata_c2r_in, data_c2r_in,[ ( Z_x / 2 + 1 ), Z_y, Z_z ] )

    cdata_c2r_out = FFTW_ALLOC_REAL( INT( real_pts, C_SIZE_T ) )
    CALL C_F_POINTER( cdata_c2r_out, data_c2r_out,[ Z_x, Z_y, Z_z ] )

    ! -----------------------------------
    ! PLAN FOR OUT-PLACE FORWARD DFT R2C
    ! -----------------------------------
    plan_c2r      = FFTW_PLAN_DFT_C2R_3D( Z_z, Z_y, Z_x, data_c2r_in, data_c2r_out, FFTW_ESTIMATE)
    ! Note the *** REVERSE *** order of the array (unlike in C)

    !_____________________________________________________________
    ! 'x' comp
    ! ============================================================

    ! ---------------------
    ! INITIALIZE INPUT DATA  (in format of FFTW from first Brillouin zone format)
    ! ---------------------
    DO q_x = qMin_x, qMax_x ! Along 'q_x' direction

    DO q_y = qMin_y, - 1
    DO q_z = qMin_z, - 1

            data_c2r_in( q_x + 1, q_y + Z_y + 1, q_z + Z_z + 1 ) = in_x( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = 0, qMax_z

            data_c2r_in( q_x + 1, q_y + 1, q_z + 1 ) = in_x( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = qMin_z, - 1

            data_c2r_in( q_x + 1, q_y + 1, q_z + Z_z + 1 ) = in_x( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = qMin_y, - 1
    DO q_z = 0, qMax_z

            data_c2r_in( q_x + 1, q_y + Z_y + 1, q_z + 1 ) = in_x( q_x, q_y, q_z )

    END DO
    END DO

    END DO

    ! -----------
    ! EXECUTE DFT
    ! -----------
    CALL FFTW_EXECUTE_DFT_C2R( plan_c2r, data_c2r_in, data_c2r_out )

    ! ------------
    ! WRITE OUTPUT
    ! ------------
    out_x = data_c2r_out

    !_____________________________________________________________
    ! 'y' comp
    ! ============================================================

    ! ---------------------
    ! INITIALIZE INPUT DATA  (in format of FFTW from first Brillouin zone format)
    ! ---------------------
    DO q_x = qMin_x, qMax_x ! Along 'q_x' direction

    DO q_y = qMin_y, - 1
    DO q_z = qMin_z, - 1

            data_c2r_in( q_x + 1, q_y + Z_y + 1, q_z + Z_z + 1 ) = in_y( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = 0, qMax_z

            data_c2r_in( q_x + 1, q_y + 1, q_z + 1 ) = in_y( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = qMin_z, - 1

            data_c2r_in( q_x + 1, q_y + 1, q_z + Z_z + 1 ) = in_y( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = qMin_y, - 1
    DO q_z = 0, qMax_z

            data_c2r_in( q_x + 1, q_y + Z_y + 1, q_z + 1 ) = in_y( q_x, q_y, q_z )

    END DO
    END DO


    END DO

    ! -----------
    ! EXECUTE DFT
    ! -----------
    CALL FFTW_EXECUTE_DFT_C2R( plan_c2r, data_c2r_in, data_c2r_out )

    ! ------------
    ! WRITE OUTPUT
    ! ------------
    out_y = data_c2r_out

    !_____________________________________________________________
    ! 'z' comp
    ! ============================================================

    ! ---------------------
    ! INITIALIZE INPUT DATA  (in format of FFTW from first Brillouin zone format)
    ! ---------------------
    DO q_x = qMin_x, qMax_x ! Along 'q_x' direction

    DO q_y = qMin_y, - 1
    DO q_z = qMin_z, - 1

            data_c2r_in( q_x + 1, q_y + Z_y + 1, q_z + Z_z + 1 ) = in_z( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = 0, qMax_z

            data_c2r_in( q_x + 1, q_y + 1, q_z + 1 ) = in_z( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = qMin_z, - 1

            data_c2r_in( q_x + 1, q_y + 1, q_z + Z_z + 1 ) = in_z( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = qMin_y, - 1
    DO q_z = 0, qMax_z

            data_c2r_in( q_x + 1, q_y + Z_y + 1, q_z + 1 ) = in_z( q_x, q_y, q_z )

    END DO
    END DO

    END DO

    ! -----------
    ! EXECUTE DFT
    ! -----------
    CALL FFTW_EXECUTE_DFT_C2R( plan_c2r, data_c2r_in, data_c2r_out )

    ! ------------
    ! WRITE OUTPUT
    ! ------------
    out_z = data_c2r_out

    ! -------------
    ! DESTROY PLANS
    ! -------------
    CALL FFTW_DESTROY_PLAN( plan_c2r )
    CALL FFTW_FREE( cdata_c2r_in )
    CALL FFTW_FREE( cdata_c2r_out )

  END

  SUBROUTINE fft_r2c( in, out )
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS WITH REAL INPUT ARRAY 'IN'
  ! AND GET SPECTRAL OUTPUT ARRAY 'OUT'
  ! ------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,DIMENSION( 0 : Z_x - 1, 0 : Z_y - 1, 0 : Z_z - 1 ),INTENT(IN)::in
    DOUBLE COMPLEX,DIMENSION( qMin_x : qMax_x, qMin_y : qMax_y, qMin_z : qMax_z ),INTENT(OUT)::out

    ! ---------------------------------------
    ! ALLOCATE ARRAYS - DYNAMIC
    ! NOTE:- The array dimensions are in reverse order for FORTRAN
    ! ---------------------------------------
    cdata_r2c_in  = FFTW_ALLOC_REAL( INT( real_pts, C_SIZE_T ) )
    CALL C_F_POINTER( cdata_r2c_in, data_r2c_in,[ Z_x, Z_y, Z_z ] )

    cdata_r2c_out = FFTW_ALLOC_COMPLEX( INT( fourier_pts, C_SIZE_T ) )
    CALL C_F_POINTER( cdata_r2c_out, data_r2c_out,[ ( Z_x / 2 + 1 ), Z_y, Z_z ] )

    ! -----------------------------------
    ! PLAN FOR OUT-PLACE FORWARD DFT R2C
    ! -----------------------------------
    plan_r2c      = FFTW_PLAN_DFT_R2C_3D( Z_z, Z_y, Z_x, data_r2c_in, data_r2c_out, FFTW_ESTIMATE)
    ! Note the *** REVERSE *** order of the array (unlike in C)

    ! ---------------------
    ! INITIALIZE INPUT DATA
    ! ---------------------
    data_r2c_in = in

    ! -----------
    ! EXECUTE DFT
    ! -----------
    CALL FFTW_EXECUTE_DFT_R2C( plan_r2c, data_r2c_in, data_r2c_out )

    ! -------------------------------------------------------
    ! WRITE OUTPUT (in format of first Brillouin zone format)
    ! -------------------------------------------------------
    DO q_x = qMin_x, qMax_x ! Along 'q_x' direction

    DO q_y = qMin_y, - 1
    DO q_z = qMin_z, - 1

            out( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + Z_y + 1, q_z + Z_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = 0, qMax_z

            out( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + 1, q_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = qMin_z, - 1

            out( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + 1, q_z + Z_z + 1 ) / Z_factor

    END DO
    END DO

    DO q_y = qMin_y, - 1
    DO q_z = 0, qMax_z

            out( q_x, q_y, q_z ) = data_r2c_out( q_x + 1, q_y + Z_y + 1, q_z + 1 ) / Z_factor

    END DO
    END DO

    END DO

    ! -------------
    ! DESTROY PLANS
    ! -------------
    CALL FFTW_DESTROY_PLAN( plan_r2c )
    CALL FFTW_FREE( cdata_r2c_in )
    CALL FFTW_FREE( cdata_r2c_out )

  END

  SUBROUTINE fft_c2r( in, out )
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS WITH COMPLEX INPUT ARRAY 'IN'
  ! AND GET REAL OUTPUT ARRAY 'OUT'
  ! ------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX,DIMENSION( qMin_x : qMax_x, qMin_y : qMax_y, qMin_z : qMax_z ),INTENT(IN)::in
    DOUBLE PRECISION,DIMENSION( 0 : Z_x - 1, 0 : Z_y - 1, 0 : Z_z - 1 ),INTENT(OUT)::out

    ! ---------------------------------------
    ! ALLOCATE ARRAYS - DYNAMIC
    ! NOTE:- The array dimensions are in reverse order for FORTRAN
    ! ---------------------------------------
    cdata_c2r_in  = FFTW_ALLOC_COMPLEX( INT( fourier_pts, C_SIZE_T ) )
    CALL C_F_POINTER( cdata_c2r_in, data_c2r_in,[ ( Z_x / 2 + 1 ), Z_y, Z_z ] )

    cdata_c2r_out = FFTW_ALLOC_REAL( INT( real_pts, C_SIZE_T ) )
    CALL C_F_POINTER( cdata_c2r_out, data_c2r_out,[ Z_x, Z_y, Z_z ] )

    ! -----------------------------------
    ! PLAN FOR OUT-PLACE FORWARD DFT R2C
    ! -----------------------------------
    plan_c2r      = FFTW_PLAN_DFT_C2R_3D( Z_z, Z_y, Z_x, data_c2r_in, data_c2r_out, FFTW_ESTIMATE)
    ! Note the *** REVERSE *** order of the array (unlike in C)

    !_____________________________________________________________
    ! 'x' comp
    ! ============================================================

    ! ---------------------
    ! INITIALIZE INPUT DATA  (in format of FFTW from first Brillouin zone format)
    ! ---------------------
    DO q_x = qMin_x, qMax_x ! Along 'q_x' direction

    DO q_y = qMin_y, - 1
    DO q_z = qMin_z, - 1

            data_c2r_in( q_x + 1, q_y + Z_y + 1, q_z + Z_z + 1 ) = in( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = 0, qMax_z

            data_c2r_in( q_x + 1, q_y + 1, q_z + 1 ) = in( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = 0, qMax_y
    DO q_z = qMin_z, - 1

            data_c2r_in( q_x + 1, q_y + 1, q_z + Z_z + 1 ) = in( q_x, q_y, q_z )

    END DO
    END DO

    DO q_y = qMin_y, - 1
    DO q_z = 0, qMax_z

            data_c2r_in( q_x + 1, q_y + Z_y + 1, q_z + 1 ) = in( q_x, q_y, q_z )

    END DO
    END DO

    END DO

    ! -----------
    ! EXECUTE DFT
    ! -----------
    CALL FFTW_EXECUTE_DFT_C2R( plan_c2r, data_c2r_in, data_c2r_out )

    ! ------------
    ! WRITE OUTPUT
    ! ------------
    out = data_c2r_out

    ! -------------
    ! DESTROY PLANS
    ! -------------
    CALL FFTW_DESTROY_PLAN( plan_c2r )
    CALL FFTW_FREE( cdata_c2r_in )
    CALL FFTW_FREE( cdata_c2r_out )

  END

END MODULE system_fftw_adv
