! needs the output of "irt objects" as input file
! writes out single cell tracks, ignoring joining and splitting
! Compile: ifort -no-wrap-margin -o irt_tracks_release.x irt_tracks_release.f90

PROGRAM irt_tracks

USE irt_parameters, ONLY: max_no_of_tracks, max_length_of_track, &
                          termination_sensitivity, n_fields

IMPLICIT NONE

! cell data storage
INTEGER              :: cell_timestep(max_length_of_track,max_no_of_tracks)
INTEGER              :: cell_number(max_length_of_track,max_no_of_tracks)
INTEGER              :: cell_id(max_length_of_track,max_no_of_tracks)   ! cell_id begins with 1 at every time step
REAL                 :: totarea(max_length_of_track,max_no_of_tracks)
REAL                 :: field_mean(max_length_of_track,max_no_of_tracks,n_fields+1)
REAL                 :: field_min(max_length_of_track,max_no_of_tracks,n_fields+1)
REAL                 :: field_max(max_length_of_track,max_no_of_tracks,n_fields+1)
REAL                 :: center_of_mass_x(max_length_of_track,max_no_of_tracks)
REAL                 :: center_of_mass_y(max_length_of_track,max_no_of_tracks)
REAL                 :: velocity_x(max_length_of_track,max_no_of_tracks)
REAL                 :: velocity_y(max_length_of_track,max_no_of_tracks)
INTEGER              :: actual_cell_number(max_no_of_tracks)
INTEGER              :: next_cell_number(max_no_of_tracks)
INTEGER              :: track_length(max_no_of_tracks)
INTEGER              :: status_beginning(max_no_of_tracks),status_end(max_no_of_tracks)

! status_beginning=0 : Track emerges by itself
! status_beginning=1 : Track is fragment of splitting event
! status_beginning=2 : Track is result of merging event
! status_beginning=-1: Track started by contact with missing values

! status_end=0 : Track dissipates
! status_end=1 : Track terminates by merging event
! status_end=2 : Track terminates by splitting event
! status_end=-1: Track terminates by contact with missing values

! dummies for reading
INTEGER              :: cell_timestep0,cell_number0,cell_id0
REAL                 :: field_mean0(n_fields+1),field_min0(n_fields+1),field_max0(n_fields+1)
REAL                 :: totarea0
INTEGER              :: xfirst0,xlast0,yfirst0,ylast0
REAL                 :: center_of_mass_x0,center_of_mass_y0
REAL                 :: velocity_x0,velocity_y0
CHARACTER (len=1)    :: sternchen

! link numbers and sizes
INTEGER              :: largest_forward_link,second_largest_forward_link
INTEGER              :: largest_backward_link,second_largest_backward_link
REAL                 :: largest_forward_link_size,second_largest_forward_link_size
REAL                 :: largest_backward_link_size,second_largest_backward_link_size

INTEGER              :: status_forward,status_backward  ! these ones will not be read in, but diagnosed

! indices
INTEGER              :: i,j,counter,actual_time
INTEGER              :: no_of_tracks, longest_track
!INTEGER              :: no_of_links_backward,no_of_links_forward

CHARACTER (len=90), PARAMETER   :: input_filename  = "irt_objects_output.txt"
CHARACTER (len=90), PARAMETER   :: output_filename = "irt_tracks_output.txt"
CHARACTER (len=90), PARAMETER   :: output2_filename = "irt_tracks_nohead_output.txt"

OPEN(10,FILE=trim(input_filename),FORM='formatted', ACTION='read')
OPEN(20,FILE=trim(output_filename),FORM='formatted', ACTION='write')
OPEN(30,FILE=trim(output2_filename),FORM='formatted', ACTION='write')

track_length=0 ! set all tracks to the beginning
counter=0
no_of_tracks=0
actual_time=0
longest_track=0

! beginning of main loop
DO

READ(10,*,END=200) cell_timestep0,cell_id0,cell_number0,totarea0, &
                   field_mean0(:),field_min0(:),field_max0(:), &
                   xfirst0,xlast0,yfirst0,ylast0,center_of_mass_x0,center_of_mass_y0,velocity_x0,velocity_y0, &
		   largest_forward_link,largest_forward_link_size, &
		   second_largest_forward_link,second_largest_forward_link_size, &
		   largest_backward_link,largest_backward_link_size, &
		   second_largest_backward_link,second_largest_backward_link_size

IF (cell_timestep0 .NE. actual_time) THEN
  actual_time = cell_timestep0
  WRITE (*,*) actual_time,cell_number0,no_of_tracks, longest_track
ENDIF

! diagnose ratio_forward
IF (largest_forward_link .EQ. -1) THEN  ! interrupted by missing values --> terminating the track
  status_forward = -1
ELSEIF (largest_forward_link .EQ. 0) THEN ! no forward link
  status_forward = 0
ELSEIF (second_largest_forward_link_size/largest_forward_link_size .LE. termination_sensitivity) THEN 
  ! track is continued
  status_forward = 1
ELSE   ! major splitting event, terminating the track
  status_forward = 2
ENDIF

! diagnose ratio_backward
IF (largest_backward_link .EQ. -1) THEN  ! interrupted by missing values --> initiating new track
  status_backward = -1
ELSEIF (largest_backward_link .EQ. 0) THEN ! no forward link
  status_backward = 0
ELSEIF (second_largest_backward_link_size/largest_backward_link_size .LE. termination_sensitivity) THEN 
  ! track is continued
  status_backward = 1
ELSE   ! major merging event, initiating new track
  status_backward = 2
ENDIF

! If cell_number0 terminates tracks because of a merging event...
IF (status_backward .EQ. 2) THEN
  ! Find the tracks which have cell_number0 as successors
  DO i=1,max_no_of_tracks
    IF (next_cell_number(i) .EQ. cell_number0) THEN
      ! terminate the track i
      IF (track_length(i) .GT. 1) THEN
        ! Mark this case with a value of 1 for status_end(i): Merging event
	status_end(i) = 1
        CALL write_output(max_length_of_track,longest_track,counter,n_fields, &
                    track_length(i),cell_timestep(:,i),cell_number(:,i),cell_id(:,i),&
		    status_beginning(i),status_end(i),totarea(:,i), &
		    field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
		    center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i))
      ENDIF
      track_length(i) = 0 ! give the track space free for a new one
      no_of_tracks = no_of_tracks-1
    ENDIF
  ENDDO
ENDIF

! If cell_number0 continuates a track, it could terminate smaller ones, dependent on termination sensitivity
IF (status_backward .EQ. 1) THEN
  ! find out if there is a predecessor
  DO i=1,max_no_of_tracks
    IF (next_cell_number(i) .EQ. cell_number0 .AND. actual_cell_number(i) .NE. largest_backward_link) THEN
      ! terminate the cell i
      IF (track_length(i) .GT. 1) THEN
        ! Mark this case with a value of 1 for status_end(i): Merging event
	status_end(i) = 1
        CALL write_output(max_length_of_track,longest_track,counter,n_fields, &
                    track_length(i),cell_timestep(:,i),cell_number(:,i),cell_id(:,i),&
		    status_beginning(i),status_end(i),totarea(:,i),&
		    field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
		    center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i))
      ENDIF
      track_length(i) = 0 ! give the track space free for a new one
      no_of_tracks = no_of_tracks-1
    ENDIF
  ENDDO
ENDIF

! Continuation of a track?
IF (status_backward .EQ. 1 .AND. status_forward .EQ. 1) THEN
  ! find the track to which the cell belongs to
  DO i=1,max_no_of_tracks
    IF (next_cell_number(i) .EQ. cell_number0 .AND. &
        actual_cell_number(i) .EQ. largest_backward_link) EXIT
  ENDDO
  IF (i .GT. max_no_of_tracks) THEN
    ! No track found. --> This cell initiates a new track which is a fragment of a splitting event
    ! Mark this case with a value of -2
    status_backward=-2
    GOTO 100
  ENDIF
  track_length(i)=track_length(i)+1
  j=track_length(i)
  IF (j .GE. max_length_of_track) THEN
    WRITE (*,*) 'ERROR: track too long, increase max_length_of_track!'
    STOP
  ENDIF
  actual_cell_number(i)=cell_number0
  next_cell_number(i)=largest_forward_link
  cell_timestep(j,i)=cell_timestep0
  cell_number(j,i)=cell_number0
  cell_id(j,i)=cell_id0
  totarea(j,i)=totarea0
  field_mean(j,i,:)=field_mean0(:)
  field_min(j,i,:)=field_min0(:)
  field_max(j,i,:)=field_max0(:)
  center_of_mass_x(j,i)=center_of_mass_x0
  center_of_mass_y(j,i)=center_of_mass_y0
  velocity_x(j,i)=velocity_x0
  velocity_y(j,i)=velocity_y0
ENDIF

100 CONTINUE

! beginning of a new single track?
IF (status_backward .NE. 1 .AND. status_forward .EQ. 1) THEN
  ! find first free track space
  DO i=1,max_no_of_tracks
    IF (track_length(i) .EQ. 0) EXIT
  ENDDO
  IF (i .GT. max_no_of_tracks) THEN
    WRITE (*,*) 'ERROR: no free track space, increase max_no_of_tracks!'
    STOP
  ENDIF
  track_length(i)=1
  actual_cell_number(i)=cell_number0
  next_cell_number(i)=largest_forward_link
  IF (status_backward .EQ. -2) THEN
    status_beginning(i)=1  ! splitting event
  ELSE
    status_beginning(i)=status_backward
  ENDIF
  cell_timestep(1,i)=cell_timestep0
  cell_number(1,i)=cell_number0
  cell_id(1,i)=cell_id0
  totarea(1,i)=totarea0
  field_mean(1,i,:)=field_mean0(:)
  field_min(1,i,:)=field_min0(:)
  field_max(1,i,:)=field_max0(:)
  center_of_mass_x(1,i)=center_of_mass_x0
  center_of_mass_y(1,i)=center_of_mass_y0
  velocity_x(1,i)=velocity_x0
  velocity_y(1,i)=velocity_y0
  no_of_tracks=no_of_tracks+1
ENDIF

! end of track?
IF (status_backward .EQ. 1 .AND. status_forward .NE. 1) THEN
  ! find the track which the cell belongs to
  DO i=1,max_no_of_tracks
    IF (next_cell_number(i) .EQ. cell_number0 .AND. &
        actual_cell_number(i) .EQ. largest_backward_link) EXIT
  ENDDO
  IF (i .GT. max_no_of_tracks) THEN
    ! In this case, ignore the cell since it would lead to
    ! a track length of 1.
    CYCLE
  ENDIF
  track_length(i)=track_length(i)+1
  j=track_length(i)
  IF (j .GT. max_length_of_track) THEN
    WRITE (*,*) 'ERROR: track too long, increase max_length_of_track!'
    STOP
  ENDIF
  actual_cell_number(i)=cell_number0
  status_end(i)=status_forward
  cell_timestep(j,i)=cell_timestep0
  cell_number(j,i)=cell_number0
  cell_id(j,i)=cell_id0
  totarea(j,i)=totarea0
  field_mean(j,i,:)=field_mean0(:)
  field_min(j,i,:)=field_min0(:)
  field_max(j,i,:)=field_max0(:)
  center_of_mass_x(j,i)=center_of_mass_x0
  center_of_mass_y(j,i)=center_of_mass_y0
  velocity_x(j,i)=velocity_x0
  velocity_y(j,i)=velocity_y0
  ! write output
  CALL write_output(max_length_of_track,longest_track,counter,n_fields, &
  	      track_length(i),cell_timestep(:,i),cell_number(:,i),cell_id(:,i),&
	      status_beginning(i),status_end(i),totarea(:,i), &
	      field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
  	      center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i))
  track_length(i)=0 ! give the track space free for a new one
  no_of_tracks=no_of_tracks-1
ENDIF

! end main loop
ENDDO

200 CONTINUE

! write out all remaining tracks
DO i=1,max_no_of_tracks
  IF (track_length(i) .GT. 1) THEN
    status_end(i) = -1    ! irregular termination
    CALL write_output(max_length_of_track,longest_track,counter,n_fields, &
        	track_length(i),cell_timestep(:,i),cell_number(:,i),cell_id(:,i),&
		status_beginning(i),status_end(i),totarea(:,i), &
		field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
    		center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i))
  ENDIF
ENDDO

WRITE (*,*) actual_time,cell_number0,no_of_tracks, longest_track

CLOSE(30)
CLOSE(20)
CLOSE(10)

END PROGRAM irt_tracks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE write_output(max_length_of_track,longest_track,counter,n_fields,track_length,cell_timestep, &
                        cell_number,cell_id,status_beginning,status_end,totarea, &
			field_mean,field_min,field_max, &
			center_of_mass_x,center_of_mass_y,velocity_x,velocity_y)

INTEGER, PARAMETER    :: area_bins=1

INTEGER, INTENT(IN)    :: max_length_of_track
INTEGER, INTENT(INOUT) :: counter, longest_track
INTEGER, INTENT(IN)    :: n_fields
INTEGER, INTENT(IN)    :: track_length
INTEGER, INTENT(IN)    :: status_beginning
INTEGER, INTENT(IN)    :: status_end
INTEGER, INTENT(IN)    :: cell_timestep(max_length_of_track)
INTEGER, INTENT(IN)    :: cell_number(max_length_of_track)
INTEGER, INTENT(IN)    :: cell_id(max_length_of_track)
REAL, INTENT(IN)       :: totarea(max_length_of_track)
REAL, INTENT(IN)       :: field_mean(max_length_of_track,n_fields+1)
REAL, INTENT(IN)       :: field_min(max_length_of_track,n_fields+1)
REAL, INTENT(IN)       :: field_max(max_length_of_track,n_fields+1)
REAL, INTENT(IN)       :: center_of_mass_x(max_length_of_track)
REAL, INTENT(IN)       :: center_of_mass_y(max_length_of_track)
REAL, INTENT(IN)       :: velocity_x(max_length_of_track)
REAL, INTENT(IN)       :: velocity_y(max_length_of_track)

counter=counter+1

WRITE (20,*) '*'
WRITE (20,*) counter,cell_timestep(1),track_length,status_beginning,status_end

DO j=1,track_length
  WRITE(20,*) counter,cell_timestep(j),cell_id(j),cell_number(j),totarea(j), &
  field_mean(j,:),field_min(j,:),field_max(j,:), &
  center_of_mass_x(j),center_of_mass_y(j),velocity_x(j),velocity_y(j)
  WRITE(30,*) counter,cell_timestep(j),cell_id(j),status_beginning,status_end, &
              center_of_mass_x(j),center_of_mass_y(j), totarea(j), field_mean(j,1)
ENDDO

IF (track_length .GT. longest_track) longest_track = track_length

END SUBROUTINE
