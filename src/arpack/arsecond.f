      SUBROUTINE ARSECOND( T )
*
      REAL       T
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     July 26, 1991
*
*  Purpose
*  =======
*
*  ARSECOND returns the user time for a process in seconds.
*  This version gets the time from the system function ETIME.
*
*  DRR - ETIME is not officially part of the Fortran language
*        and so fails with certain compilers e.g. Cray. Use
*        CPU_TIME instead.
*
*     .. Local Scalars ..
*     REAL               T1
*     ..
*     .. Local Arrays ..
*     REAL               TARRAY( 2 )
*     ..
*     .. External Functions ..
      REAL               ETIME
*      EXTERNAL           ETIME
*     ..
*     .. Executable Statements ..
*

      CALL CPU_TIME( T )
*     T1 = ETIME( TARRAY )
*     T  = TARRAY( 1 )

      RETURN
*
*     End of ARSECOND
*
      END
