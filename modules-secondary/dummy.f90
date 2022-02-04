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

! #########################
! MODULE:
! LAST MODIFIED:
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
!
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! </f>

MODULE system_dummy
! <f Info
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
!
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! </f>

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
! <f Glob Dec
  USE system_auxilaries
	! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

	IMPLICIT  NONE
  ! _________________________
  ! REAL SPACE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN = 10) 							::
	  INTEGER  (KIND= 4 ) 		 					::
		DOUBLE PRECISION									::
		DOUBLE PRECISION,INTENT(IN)			  ::
		DOUBLE PRECISION,INTENT(OUT)			::
		DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE   ::
	  DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE		::
	  INTEGER(KIND=4) ,DIMENSION(:),    ALLOCATABLE 	::
! </f>

  CONTAINS

	SUBROUTINE name
! <f
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN = 10) 							::
	  INTEGER  (KIND= 4 ) 		 					::
		DOUBLE PRECISION									::
		DOUBLE PRECISION,INTENT(IN)			  ::
		DOUBLE PRECISION,INTENT(OUT)			::
		DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE   ::
	  DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE		::
	  INTEGER(KIND=4) ,DIMENSION(:),    ALLOCATABLE 	::

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( ax( 0 : N ) )

    OPEN( UNIT = 1, FILE = TRIM( ADJUSTL( file_name ) ) )

	    READ( 1, f_i4,  ADVANCE ='yes')  N
	    READ( 1, f_d8p4, ADVANCE ='yes') 	x

    CLOSE(1)

    WRITE( char, f_i4 ) N

    DO WHILE(  )

    END DO

    IF(  ) THEN

		ELSE

    END IF

		DEALLOCATE( ax )

	END
! </f>

END MODULE system_dummy
