! needs the sorted output and the objects mask file as input
! generate sorted output with:
! > sort -n -k2 irt_tracks_nohead_output.txt > irt_tracks_sorted.txt
! writes out a mask file with track id's
! Compile: ifort -no-wrap-margin -o irt_trackmask_release.x irt_trackmask_release.f90

PROGRAM irt_trackmask

USE irt_parameters, ONLY: domainsize_x, domainsize_y, max_no_of_tracks

IMPLICIT NONE
INTEGER              :: ii, ij,n,ncell
INTEGER              :: time_step, time_step_event
INTEGER              :: track_id(max_no_of_tracks)
INTEGER              :: cell_id(max_no_of_tracks)
INTEGER              :: status_beginning(max_no_of_tracks)
INTEGER              :: status_end(max_no_of_tracks)
REAL                 :: in_field(domainsize_x,domainsize_y)
REAL                 :: out_field(domainsize_x,domainsize_y)

INTEGER              :: srv_header(8)

CHARACTER (len=90), PARAMETER   :: input_filename1 = "irt_objects_mask.srv"
CHARACTER (len=90), PARAMETER   :: input_filename2 = "irt_tracks_sorted.txt"
CHARACTER (len=90), PARAMETER   :: output_filename = "irt_tracks_mask.srv"

OPEN(10,FILE=trim(input_filename1),FORM='unformatted', ACTION='read')
OPEN(20,FILE=trim(input_filename2),FORM='formatted', ACTION='read')
OPEN(30,FILE=trim(output_filename),FORM='unformatted', ACTION='write')

READ (20,*) track_id(1),time_step_event,cell_id(1),status_beginning(1),status_end(1)
!WRITE (*,*) time_step_event,track_id

time_step=-1

! beginning of main loop
DO

WRITE (*,*) 'timestep:',time_step

READ (10,END=200) srv_header
READ (10) in_field
time_step = time_step+1

IF (time_step<time_step_event) CYCLE

ncell = 1
DO WHILE (time_step==time_step_event)
  ncell = ncell+1
  READ (20,*,END=200) track_id(ncell),time_step_event,cell_id(ncell),status_beginning(ncell),status_end(ncell)
ENDDO

!WRITE (*,*) cell_id(1:ncell-1)

out_field = 0

DO ii=1,domainsize_x
  DO ij=1,domainsize_y
    IF (NINT(in_field(ii,ij)) .GT. 0) THEN
      out_field(ii,ij) = -1   ! default for cells which are not assigned to a track
      DO n=1,ncell-1
        IF (NINT(in_field(ii,ij)) .EQ. cell_id(n)) THEN
	  out_field(ii,ij) = track_id(n)
	  
	  !This code segement can be modified if e.g. tracks should be
	  !marked according to their status at the track beginning and end
	  !Here is an example:
	  !IF (status_beginning(n) .EQ. 0 .AND. status_end(n) .EQ. 0) THEN
	  !  out_field(ii,ij) = 1  ! solitary
	  !ELSEIF (status_beginning(n) .EQ. -1 .OR. status_end(n) .EQ. -1) THEN
	  !  out_field(ii,ij) = 6  ! contaminated by missing values
	  !ELSEIF (status_beginning(n) .EQ. 2 .AND. status_end(n) .EQ. 2) THEN
	  !  out_field(ii,ij) = 2  ! merging and then splitting up
	  !ELSEIF (status_beginning(n) .EQ. 1) THEN
	  !  out_field(ii,ij) = 3  ! other mergers
	  !ELSEIF (status_beginning(n) .EQ. 2) THEN
	  !  out_field(ii,ij) = 4  ! fragments
	  !ELSE
	  !  out_field(ii,ij) = 5
	  !ENDIF
	  EXIT
	ENDIF
      ENDDO
    ELSEIF (NINT(in_field(ii,ij)) .LT. 0) THEN
      out_field(ii,ij) = -1
    ENDIF
  ENDDO
ENDDO

track_id(1) = track_id(ncell)
cell_id(1) = cell_id(ncell)
status_beginning(1) = status_beginning(ncell)
status_end(1) = status_end(ncell)

WRITE (30) srv_header
WRITE (30) REAL(out_field)

! end main loop
ENDDO

200 CONTINUE

CLOSE(30)
CLOSE(20)
CLOSE(10)

END PROGRAM irt_trackmask
