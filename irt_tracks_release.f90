! needs the output of "irt objects" as input file
! writes out single cell tracks, ignoring joining and splitting
! Compile: gfortran -O0 -o irt_tracks_release.x irt_tracks_release.f90 irt_parameters.f90

PROGRAM irt_tracks

USE irt_parameters, ONLY: max_no_of_tracks, max_length_of_track, &
                          termination_sensitivity, n_fields

IMPLICIT NONE

! cell data storage
INTEGER              :: cell_timestep(max_length_of_track,max_no_of_tracks)
INTEGER              :: cell_id(max_length_of_track,max_no_of_tracks)   ! cell_id begins with 1 at every time step
INTEGER              :: cell_number(max_length_of_track,max_no_of_tracks)
INTEGER              :: cell_age1(max_length_of_track,max_no_of_tracks)
INTEGER              :: cell_age2(max_length_of_track,max_no_of_tracks)
REAL                 :: totarea(max_length_of_track,max_no_of_tracks)
REAL                 :: field_mean(max_length_of_track,max_no_of_tracks,n_fields+1)
REAL                 :: field_min(max_length_of_track,max_no_of_tracks,n_fields+1)
REAL                 :: field_max(max_length_of_track,max_no_of_tracks,n_fields+1)
INTEGER              :: xfirst(max_length_of_track,max_no_of_tracks)
INTEGER              :: xlast(max_length_of_track,max_no_of_tracks)
INTEGER              :: yfirst(max_length_of_track,max_no_of_tracks)
INTEGER              :: ylast(max_length_of_track,max_no_of_tracks)
INTEGER              :: lfl(max_length_of_track,max_no_of_tracks) ! largest forward link
INTEGER              :: slfl(max_length_of_track,max_no_of_tracks) ! second largest forward link
INTEGER              :: lbl(max_length_of_track,max_no_of_tracks) ! largest backward link
INTEGER              :: slbl(max_length_of_track,max_no_of_tracks) ! second largest backward link
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
INTEGER              :: cell_timestep0,cell_number0,cell_id0,cell_age10,cell_age20
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

!for debugging
LOGICAL              :: cell_used

CHARACTER (len=90), PARAMETER   :: input_filename  = "irt_objects_output.txt"
CHARACTER (len=90), PARAMETER   :: output_filename = "irt_tracks_output.txt"
CHARACTER (len=90), PARAMETER   :: output2_filename = "irt_tracks_nohead_output.txt"
CHARACTER (len=90), PARAMETER   :: output3_filename = "irt_tracks_nohead_begin_output.txt" ! temporary
!CHARACTER (len=90), PARAMETER   :: output4_filename = "irt_tracks_nohead_end_output.txt"
!CHARACTER (len=90), PARAMETER   :: output5_filename = "irt_tracks_nohead_maxarea_output.txt"

OPEN(10,FILE=trim(input_filename),FORM='formatted', ACTION='read')
OPEN(20,FILE=trim(output_filename),FORM='formatted', ACTION='write')
OPEN(30,FILE=trim(output2_filename),FORM='formatted', ACTION='write')
OPEN(40,FILE=trim(output3_filename),FORM='formatted', ACTION='write') !temporary
!OPEN(50,FILE=trim(output4_filename),FORM='formatted', ACTION='write')
!OPEN(60,FILE=trim(output5_filename),FORM='formatted', ACTION='write')

track_length=0 ! set all tracks to the beginning
counter=0
no_of_tracks=0
actual_time=0
longest_track=0

! beginning of main loop
DO

READ(10,*,END=200) cell_timestep0,cell_id0,cell_number0,cell_age10,cell_age20,totarea0, &
                   field_mean0(:),field_min0(:),field_max0(:), &
                   xfirst0,xlast0,yfirst0,ylast0,center_of_mass_x0,center_of_mass_y0,velocity_x0,velocity_y0, &
		   largest_forward_link,largest_forward_link_size, &
		   second_largest_forward_link,second_largest_forward_link_size, &
		   largest_backward_link,largest_backward_link_size, &
		   second_largest_backward_link,second_largest_backward_link_size

IF (cell_id0 .EQ. 0) CYCLE

cell_used = .FALSE. ! debug

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
      IF (track_length(i) .GT. 0) THEN
        ! Mark this case with a value of 1 for status_end(i): Merging event
	status_end(i) = 1
        !CALL write_output(max_length_of_track,longest_track,counter,n_fields, &
        !            track_length(i),cell_timestep(:,i),cell_age(:,i),cell_number(:,i),cell_id(:,i),&
	!	    status_beginning(i),status_end(i),totarea(:,i), &
	!	    field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
	!	    center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i))
        CALL write_output(max_length_of_track,longest_track,counter,n_fields, track_length(i), &
                    cell_timestep(:,i),cell_id(:,i),cell_number(:,i), &
                    cell_age1(:,i),cell_age2(:,i),totarea(:,i), &
		    field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
                    xfirst(:,i),xlast(:,i),yfirst(:,i),ylast(:,i), &
		    center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i), &
                    lfl(:,i),slfl(:,i),lbl(:,i),slbl(:,i), &
		    status_beginning(i),status_end(i))
        !ENDIF
        !IF (track_length(i) .LE. 1) WRITE(*,*) "merging2",track_length(i) ! debug
        track_length(i) = 0 ! give the track space free for a new one
        no_of_tracks = no_of_tracks-1
      ENDIF
    ENDIF
  ENDDO
ENDIF

! If cell_number0 continues a track, it could terminate smaller ones, dependent on termination sensitivity
IF (status_backward .EQ. 1) THEN
  ! find out if there is a predecessor
  DO i=1,max_no_of_tracks
    IF (next_cell_number(i) .EQ. cell_number0 .AND. actual_cell_number(i) .NE. largest_backward_link) THEN
      ! terminate the cell i
      IF (track_length(i) .GT. 0) THEN
        ! Mark this case with a value of 1 for status_end(i): Merging event
	status_end(i) = 1
        !CALL write_output(max_length_of_track,longest_track,counter,n_fields, &
        !            track_length(i),cell_timestep(:,i),cell_age(:,i),cell_number(:,i),cell_id(:,i),&
	!	    status_beginning(i),status_end(i),totarea(:,i),&
	!	    field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
	!	    center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i))
        CALL write_output(max_length_of_track,longest_track,counter,n_fields, track_length(i), &
                    cell_timestep(:,i),cell_id(:,i),cell_number(:,i), &
                    cell_age1(:,i),cell_age2(:,i),totarea(:,i), &
		    field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
                    xfirst(:,i),xlast(:,i),yfirst(:,i),ylast(:,i), &
		    center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i), &
                    lfl(:,i),slfl(:,i),lbl(:,i),slbl(:,i), &
		    status_beginning(i),status_end(i))
        !ENDIF
        !IF (track_length(i) .LE. 1) WRITE(*,*) "merging1",track_length(i) ! debug
        track_length(i) = 0 ! give the track space free for a new one
        no_of_tracks = no_of_tracks-1
      ENDIF
    ENDIF
  ENDDO
ENDIF

! Continuation of a track?
!IF (status_backward .EQ. 1 .AND. status_forward .EQ. 1) THEN
IF (status_backward .EQ. 1) THEN
  ! find the track to which the cell belongs to
  DO i=1,max_no_of_tracks
    IF (next_cell_number(i) .EQ. cell_number0 .AND. &
        actual_cell_number(i) .EQ. largest_backward_link .AND. track_length(i) .GT. 0) EXIT
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
  cell_id(j,i)=cell_id0
  cell_number(j,i)=cell_number0
  cell_age1(j,i)=cell_age10
  cell_age2(j,i)=cell_age20
  totarea(j,i)=totarea0
  field_mean(j,i,:)=field_mean0(:)
  field_min(j,i,:)=field_min0(:)
  field_max(j,i,:)=field_max0(:)
  xfirst(j,i)=xfirst0
  xlast(j,i)=xlast0
  yfirst(j,i)=yfirst0
  ylast(j,i)=ylast0
  center_of_mass_x(j,i)=center_of_mass_x0
  center_of_mass_y(j,i)=center_of_mass_y0
  velocity_x(j,i)=velocity_x0
  velocity_y(j,i)=velocity_y0
  lfl(j,i)=largest_forward_link
  slfl(j,i)=second_largest_forward_link
  lbl(j,i)=largest_backward_link
  slbl(j,i)=second_largest_backward_link
  cell_used = .TRUE. ! debug
  IF (status_forward .NE. 1) THEN ! track terminates
    status_end(i)=status_forward
    ! write output
    !CALL write_output(max_length_of_track,longest_track,counter,n_fields, &
    !	      track_length(i),cell_timestep(:,i),cell_age(:,i),cell_number(:,i),cell_id(:,i),&
    !	      status_beginning(i),status_end(i),totarea(:,i), &
    !         field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
    !         center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i))
    CALL write_output(max_length_of_track,longest_track,counter,n_fields, track_length(i), &
              cell_timestep(:,i),cell_id(:,i),cell_number(:,i), &
              cell_age1(:,i),cell_age2(:,i),totarea(:,i), &
	      field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
              xfirst(:,i),xlast(:,i),yfirst(:,i),ylast(:,i), &
	      center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i), &
              lfl(:,i),slfl(:,i),lbl(:,i),slbl(:,i), &
	      status_beginning(i),status_end(i))
    track_length(i)=0 ! give the track space free again
    no_of_tracks=no_of_tracks-1
  ENDIF
ENDIF

100 CONTINUE

! beginning of a new single track?
IF (status_backward .NE. 1) THEN
!IF (status_backward .NE. 1 .AND. status_forward .EQ. 1) THEN
  ! find first free track space
  DO i=1,max_no_of_tracks
    IF (track_length(i) .EQ. 0) EXIT
  ENDDO
  IF (i .GT. max_no_of_tracks) THEN
    WRITE (*,*) 'ERROR: no free track space, increase max_no_of_tracks!'
    STOP
  ENDIF
  actual_cell_number(i)=cell_number0
  next_cell_number(i)=largest_forward_link
  IF (status_backward .EQ. -2) THEN
    status_beginning(i)=1  ! splitting event
  ELSE
    status_beginning(i)=status_backward
  ENDIF
  !IF (status_backward .EQ. -1) WRITE(*,*) "status_backward=-1"
  cell_timestep(1,i)=cell_timestep0
  cell_id(1,i)=cell_id0
  cell_number(1,i)=cell_number0
  cell_age1(1,i)=cell_age10
  cell_age2(1,i)=cell_age20
  totarea(1,i)=totarea0
  field_mean(1,i,:)=field_mean0(:)
  field_min(1,i,:)=field_min0(:)
  field_max(1,i,:)=field_max0(:)
  xfirst(1,i)=xfirst0
  xlast(1,i)=xlast0
  yfirst(1,i)=yfirst0
  ylast(1,i)=ylast0
  center_of_mass_x(1,i)=center_of_mass_x0
  center_of_mass_y(1,i)=center_of_mass_y0
  velocity_x(1,i)=velocity_x0
  velocity_y(1,i)=velocity_y0
  lfl(1,i)=largest_forward_link
  slfl(1,i)=second_largest_forward_link
  lbl(1,i)=largest_backward_link
  slbl(1,i)=second_largest_backward_link
  cell_used = .TRUE. ! debug
  IF (status_forward .EQ. 1) THEN
    ! occupy slot for a new track in memory
    track_length(i)=1
    no_of_tracks=no_of_tracks+1
  ELSE ! track is only one timestep long, write immediately
    !CALL write_output(max_length_of_track,longest_track,counter,n_fields, &
    !	        1,cell_timestep(:,i),cell_age(:,i),cell_number(:,i),cell_id(:,i),&
    !	        status_beginning(i),status_end(i),totarea(:,i), &
    !        field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
    !        center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i))
    !cell_timestep(1,i)=999 ! for debugging
    status_end(i)=status_forward
    CALL write_output(max_length_of_track,longest_track,counter,n_fields, 1, &
              cell_timestep(:,i),cell_id(:,i),cell_number(:,i), & 
              cell_age1(:,i),cell_age2(:,i),totarea(:,i), &
	      field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
              xfirst(:,i),xlast(:,i),yfirst(:,i),ylast(:,i), &
	      center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i), &
              lfl(:,i),slfl(:,i),lbl(:,i),slbl(:,i), &
	      status_beginning(i),status_end(i))   ! write track with length 1
  ENDIF
ENDIF

! debug
IF (.NOT. cell_used) THEN ! should never happen
  WRITE(*,*) "Cell not assigned to any track!"
  WRITE(*,*) cell_timestep0,cell_id0,cell_number0,cell_age10,cell_age20,totarea0, &
                   field_mean0(:),field_min0(:),field_max0(:), &
                   xfirst0,xlast0,yfirst0,ylast0,center_of_mass_x0, &
                   center_of_mass_y0,velocity_x0,velocity_y0, &
		   largest_forward_link,largest_forward_link_size, &
		   second_largest_forward_link,second_largest_forward_link_size, &
		   largest_backward_link,largest_backward_link_size, &
		   second_largest_backward_link,second_largest_backward_link_size
  STOP
ENDIF

! end main loop
ENDDO

200 CONTINUE

! write out all remaining tracks
DO i=1,max_no_of_tracks
  IF (track_length(i) .GT. 0) THEN
    status_end(i) = -1    ! irregular termination
    !CALL write_output(max_length_of_track,longest_track,counter,n_fields, &
    !    	track_length(i),cell_timestep(:,i),cell_age(:,i),cell_number(:,i),cell_id(:,i),&
    !		status_beginning(i),status_end(i),totarea(:,i), &
    !		field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
    !		center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i))
    CALL write_output(max_length_of_track,longest_track,counter,n_fields, track_length(i), &
              cell_timestep(:,i),cell_id(:,i),cell_number(:,i), &
              cell_age1(:,i),cell_age2(:,i),totarea(:,i), &
	      field_mean(:,i,:),field_min(:,i,:),field_max(:,i,:), &
              xfirst(:,i),xlast(:,i),yfirst(:,i),ylast(:,i), &
	      center_of_mass_x(:,i),center_of_mass_y(:,i),velocity_x(:,i),velocity_y(:,i), &
              lfl(:,i),slfl(:,i),lbl(:,i),slbl(:,i), &
	      status_beginning(i),status_end(i))

  ENDIF
ENDDO

WRITE (*,*) actual_time,cell_number0,no_of_tracks, longest_track

 !CLOSE(60)
 !CLOSE(50)
 CLOSE(40)  ! temporary
 CLOSE(30)
 CLOSE(20)
 CLOSE(10)

END PROGRAM irt_tracks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!SUBROUTINE write_output(max_length_of_track,longest_track,counter,n_fields,track_length,cell_timestep, &
!                        cell_age,cell_number,cell_id,status_beginning,status_end,totarea, &
!			field_mean,field_min,field_max, &
!			center_of_mass_x,center_of_mass_y,velocity_x,velocity_y)

SUBROUTINE write_output(max_length_of_track, longest_track, counter, n_fields, track_length, &
                        cell_timestep, cell_id, cell_number,cell_age1,cell_age2,totarea, &
                        field_mean, field_min, field_max, &
                        xfirst, xlast, yfirst, ylast, &
                        center_of_mass_x, center_of_mass_y, velocity_x, velocity_y, &
                        lfl, slfl, lbl, slbl, &
                        status_beginning, status_end)

INTEGER, PARAMETER    :: area_bins=1

INTEGER, INTENT(IN)    :: max_length_of_track
INTEGER, INTENT(INOUT) :: counter, longest_track
INTEGER, INTENT(IN)    :: n_fields
INTEGER, INTENT(IN)    :: track_length
INTEGER, INTENT(IN)    :: status_beginning
INTEGER, INTENT(IN)    :: status_end
INTEGER, INTENT(IN)    :: cell_timestep(max_length_of_track)
INTEGER, INTENT(IN)    :: cell_age1(max_length_of_track)
INTEGER, INTENT(IN)    :: cell_age2(max_length_of_track)
INTEGER, INTENT(IN)    :: cell_number(max_length_of_track)
INTEGER, INTENT(IN)    :: cell_id(max_length_of_track)
REAL, INTENT(IN)       :: totarea(max_length_of_track)
REAL, INTENT(IN)       :: field_mean(max_length_of_track,n_fields+1)
REAL, INTENT(IN)       :: field_min(max_length_of_track,n_fields+1)
REAL, INTENT(IN)       :: field_max(max_length_of_track,n_fields+1)
INTEGER, INTENT(IN)    :: xfirst(max_length_of_track)
INTEGER, INTENT(IN)    :: xlast(max_length_of_track)
INTEGER, INTENT(IN)    :: yfirst(max_length_of_track)
INTEGER, INTENT(IN)    :: ylast(max_length_of_track)
REAL, INTENT(IN)       :: center_of_mass_x(max_length_of_track)
REAL, INTENT(IN)       :: center_of_mass_y(max_length_of_track)
REAL, INTENT(IN)       :: velocity_x(max_length_of_track)
REAL, INTENT(IN)       :: velocity_y(max_length_of_track)
INTEGER, INTENT(IN)    :: lfl(max_length_of_track)
INTEGER, INTENT(IN)    :: slfl(max_length_of_track)
INTEGER, INTENT(IN)    :: lbl(max_length_of_track)
INTEGER, INTENT(IN)    :: slbl(max_length_of_track)

REAL     :: max_area
INTEGER  :: max_j

max_area = 0.
counter=counter+1

!IF ((status_beginning.EQ.0 .AND. lbl(1).NE.0) .OR. &
!    (status_beginning.EQ.1 .AND. lbl(1).EQ.0) .OR. &
!     status_beginning.EQ.2 .AND. slbl(1).EQ.0) THEN
!  WRITE(*,*) counter,cell_timestep(1),track_length,status_beginning,status_end,lbl(1),slbl(1),lfl(1),slfl(1)
!ENDIF

!IF (status_end.EQ.2 .AND. slfl(track_length).EQ.0) THEN
!  WRITE(*,*) counter,cell_timestep(1),track_length,status_beginning,status_end,lbl(1),slbl(1),lfl(1),slfl(1)
!ENDIF

WRITE (20,*) '*'
WRITE (20,*) counter,cell_timestep(1),track_length,status_beginning,status_end

DO j=1,track_length
  IF (totarea(j) .GT. max_area) THEN
    max_area = totarea(j)
    max_j = j
  ENDIF
  !WRITE(20,*) counter,cell_timestep(j),cell_age(j),cell_id(j),cell_number(j),totarea(j), &
  !field_mean(j,:),field_min(j,:),field_max(j,:), &
  !center_of_mass_x(j),center_of_mass_y(j),velocity_x(j),velocity_y(j)
  WRITE(20,*) counter,cell_timestep(j),cell_id(j),cell_number(j), &
              cell_age1(j),cell_age2(j),totarea(j), &
	      field_mean(j,:),field_min(j,:),field_max(j,:), &
              xfirst(j),xlast(j),yfirst(j),ylast(j), &
	      center_of_mass_x(j),center_of_mass_y(j),velocity_x(j),velocity_y(j), &
              lfl(j),slfl(j),lbl(j),slbl(j)
  WRITE(30,*) counter,cell_timestep(j),cell_age1(j),cell_age2(j),cell_id(j), &
              status_beginning,status_end, track_length, j, &
              center_of_mass_x(j),center_of_mass_y(j), totarea(j), field_mean(j,1)
ENDDO

! only track beginnings
IF (status_beginning .EQ. 0 .AND. status_end .EQ. 0 .AND. track_length .GE. 6) THEN
  WRITE(40,*) counter,cell_timestep(1),center_of_mass_x(1),center_of_mass_y(1)
ENDIF

! only track ends
!WRITE(50,*) counter,cell_timestep(track_length),cell_id(track_length),status_beginning,status_end, &
!	    center_of_mass_x(track_length),center_of_mass_y(track_length), totarea(max_j), field_mean(track_length,1)

! only maximum area of track
!WRITE(60,*) counter,cell_timestep(max_j),cell_id(max_j),status_beginning,status_end, &
!	    center_of_mass_x(max_j),center_of_mass_y(max_j), totarea(max_j), field_mean(max_j,1)

IF (track_length .GT. longest_track) longest_track = track_length

END SUBROUTINE
