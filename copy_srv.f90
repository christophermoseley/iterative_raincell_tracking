! This example program just copies an input file in SRV format to an output file
! to demonstrate how to handle SRV files in Fortran. Note that the REAL values
! are 32 bit floating point value by default. Some CDO versions use 64 bit,
! which will cause the Fortran programs to crash. In this case, use the option
! cdo -b F32 
! in the CDO commands.
PROGRAM copy_srv

IMPLICIT NONE

INTEGER                 :: xsize, ysize
REAL, ALLOCATABLE       :: field(:,:)
INTEGER                 :: header(8)

OPEN (10,FILE='irt_objects_input_00.srv', FORM='unformatted', STATUS='old', ACTION='read')
OPEN (20,FILE='irt_objects_input_00_copy.srv',FORM='unformatted', ACTION='write')

   ! read in first record header, to get the domain size: 
   READ(10) header
   ! header(1) : Code ID
   ! header(3) : Date in the form: YYYYMMDDhhmmss
   ! header(4) : Time increment
   ! header(5) : domain size in x direction
   ! header(6) : domain size in y direction
   ! header(7) = 0
   ! header(8) = 0

   xsize = header(5)
   ysize = header(6)
   ALLOCATE (field(xsize,ysize))

DO
   ! read in SRV record
   READ(10) field

   ! write out header+record into output file
   WRITE(20) header
   WRITE(20) field

   ! read in header of next record
   READ(10,END=300) header
ENDDO

300 CONTINUE ! input file has reached end

CLOSE(10)
CLOSE(20)

END

