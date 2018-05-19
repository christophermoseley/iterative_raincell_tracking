! needs the output of irt_objects as input.
! writes out coarse advection field
! Compile: ifort -no-wrap-margin -o irt_advection_field_release.x irt_advection_field_release.f90

PROGRAM irt_advection_field

USE irt_parameters, ONLY: domainsize_x, domainsize_y, lperiodic, &
    n_fields, time_steps, nt_bins, nx_bins, ny_bins, max_velocity, &
    min_cells,miss

IMPLICIT NONE

CHARACTER (len=90), PARAMETER   :: input_filename  = "irt_objects_output.txt"
CHARACTER (len=90), PARAMETER   :: output_filename = "irt_advection_field.srv"

! cell data storage
INTEGER              :: cell_timestep0
REAL                 :: totarea0
REAL                 :: center_of_mass_x0,center_of_mass_y0
INTEGER              :: cell_number0,cell_id0
REAL                 :: field_mean0(n_fields+1),field_min0(n_fields+1),field_max0(n_fields+1)
INTEGER              :: xfirst0,xlast0,yfirst0,ylast0
INTEGER              :: largest_forward_link0,second_largest_forward_link0
INTEGER              :: largest_backward_link0,second_largest_backward_link0
REAL                 :: largest_forward_link_size0,second_largest_forward_link_size0
REAL                 :: largest_backward_link_size0,second_largest_backward_link_size0
REAL                 :: velocity_x0, velocity_y0

! variables on the coarse velocity field
REAL                 :: coarse_vel_x(nt_bins,nx_bins,ny_bins)
REAL                 :: coarse_vel_y(nt_bins,nx_bins,ny_bins)
REAL                 :: area_weight(nt_bins,nx_bins,ny_bins)
REAL                 :: no_of_data_points(nt_bins,nx_bins,ny_bins)

! indices
INTEGER              :: nt,nx,ny

INTEGER              :: srv_header(8)

srv_header(1) = 1	    ! Code
srv_header(2) = 1	    ! Level
srv_header(3) = 20000101    ! Datum
srv_header(4) = 0000	    ! Zeitinkrement
srv_header(5) = nx_bins
srv_header(6) = ny_bins
srv_header(7) = 0
srv_header(8) = 0

OPEN(10,FILE=trim(input_filename),FORM='formatted', ACTION='read')
OPEN(30,FILE=trim(output_filename),FORM='unformatted', ACTION='write')

! initialize coarse field
coarse_vel_x=0
coarse_vel_y=0
area_weight=0
no_of_data_points=0

! beginning of main loop
DO

   READ(10,*,END=200) cell_timestep0,cell_id0,cell_number0,totarea0, &
                      field_mean0(:),field_min0(:),field_max0(:), &
                      xfirst0,xlast0,yfirst0,ylast0, &
		      center_of_mass_x0,center_of_mass_y0, &
		      velocity_x0,velocity_y0, & 
		      largest_forward_link0,largest_forward_link_size0, &
		      second_largest_forward_link0,second_largest_forward_link_size0, &
		      largest_backward_link0,largest_backward_link_size0, &
		      second_largest_backward_link0,second_largest_backward_link_size0

   ! filter out missing values and velocities that exceed max_velocity
   IF (velocity_x0 .LT. miss+0.1 .OR. velocity_y0 .LT. miss+0.1 .OR. &
       velocity_x0*velocity_x0+velocity_y0*velocity_y0 .GT. max_velocity*max_velocity) THEN
     CYCLE
   ENDIF
   
   IF (cell_id0 .EQ. 0) THEN
        velocity_x0 = 0.
        velocity_y0 = 0.
   ENDIF

   nt = MIN(FLOOR((REAL(cell_timestep0-1)/REAL(time_steps))*nt_bins)+1,nt_bins)
   nx = MIN(FLOOR((REAL(center_of_mass_x0)/REAL(domainsize_x))*nx_bins)+1,nx_bins)
   ny = MIN(FLOOR((REAL(center_of_mass_y0)/REAL(domainsize_y))*ny_bins)+1,ny_bins)

   coarse_vel_x(nt,nx,ny) = coarse_vel_x(nt,nx,ny)+velocity_x0*totarea0
   coarse_vel_y(nt,nx,ny) = coarse_vel_y(nt,nx,ny)+velocity_y0*totarea0
   area_weight(nt,nx,ny)  = area_weight(nt,nx,ny)+totarea0
   no_of_data_points(nt,nx,ny) = no_of_data_points(nt,nx,ny)+1

! end main loop
ENDDO

200 CONTINUE

! divide by weight
DO nt=1,nt_bins

   DO nx=1,nx_bins
     DO ny=1,ny_bins
       IF (no_of_data_points(nt,nx,ny) .GE. min_cells) THEN
         coarse_vel_x(nt,nx,ny) = coarse_vel_x(nt,nx,ny)/area_weight(nt,nx,ny)
         coarse_vel_y(nt,nx,ny) = coarse_vel_y(nt,nx,ny)/area_weight(nt,nx,ny)
       ELSE
         coarse_vel_x(nt,nx,ny) = miss
         coarse_vel_y(nt,nx,ny) = miss
       ENDIF
     ENDDO
   ENDDO
  
   CALL fill_missing_values(lperiodic,miss,nx_bins,ny_bins,coarse_vel_x(nt,:,:))
   CALL fill_missing_values(lperiodic,miss,nx_bins,ny_bins,coarse_vel_y(nt,:,:))

   ! write SRV format
   srv_header(4)=nt

   srv_header(1)=1  ! vx
   WRITE(30) srv_header
   WRITE(30) coarse_vel_x(nt,:,:)

   srv_header(1)=2  ! vy
   WRITE(30) srv_header
   WRITE(30) coarse_vel_y(nt,:,:)

   srv_header(1)=3  ! number of cells
   WRITE(30) srv_header
   WRITE(30) no_of_data_points(nt,:,:)

ENDDO

CLOSE(30)
CLOSE(10)

END PROGRAM irt_advection_field


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! fill up missing values by interpolating between neighbours


SUBROUTINE fill_missing_values(lperiodic,miss,nx_bins,ny_bins,coarse_field)

IMPLICIT NONE

INTEGER, INTENT(IN)    :: nx_bins, ny_bins
REAL , INTENT(IN)      :: miss
LOGICAL, INTENT(IN)    :: lperiodic        ! periodic boundaries?

! the indices for the track data.
INTEGER               :: nx,ny,iteration

INTEGER               :: no_of_neighbours

REAL, INTENT(INOUT)   :: coarse_field(nx_bins,ny_bins)

REAL                  :: coarse_field_bck(nx_bins,ny_bins)

REAL                  :: coarse_field_w,coarse_field_e,coarse_field_n,coarse_field_s

LOGICAL               :: iterate  ! are there still missing values? --> repeat iteration
LOGICAL               :: only_missing_values ! if all values are missing, set them to zero

! iteration
iterate = .TRUE.
DO
  IF (.NOT. iterate) EXIT
  iterate = .FALSE.
  coarse_field_bck = coarse_field
  
  only_missing_values = .TRUE.

  DO nx=1,nx_bins
    DO ny=1,ny_bins

      ! is the value already defined?
      IF (coarse_field_bck(nx,ny) .GT. miss+1) THEN
	only_missing_values = .FALSE.
        CYCLE
      ENDIF
      
      IF (lperiodic) THEN
        IF (nx==1) THEN ! west
	  coarse_field_w = coarse_field_bck(nx_bins,ny)
	ELSE
	  coarse_field_w = coarse_field_bck(nx-1,ny)
	ENDIF
        IF (nx==nx_bins) THEN ! east
	  coarse_field_e = coarse_field_bck(1,ny)
	ELSE
	  coarse_field_e = coarse_field_bck(nx+1,ny)
	ENDIF
        IF (ny==1) THEN ! south
	  coarse_field_s = coarse_field_bck(nx,ny_bins)
	ELSE
	  coarse_field_s = coarse_field_bck(nx,ny-1)
	ENDIF
        IF (ny==ny_bins) THEN ! north
	  coarse_field_n = coarse_field_bck(nx,1)
	ELSE
	  coarse_field_n = coarse_field_bck(nx,ny+1)
	ENDIF
      ELSE
        IF (nx==1) THEN ! west
	  coarse_field_w = miss
	ELSE
	  coarse_field_w = coarse_field_bck(nx-1,ny)
	ENDIF
        IF (nx==nx_bins) THEN ! east
	  coarse_field_e = miss
	ELSE
	  coarse_field_e = coarse_field_bck(nx+1,ny)
	ENDIF
        IF (ny==1) THEN ! south
	  coarse_field_s = miss
	ELSE
	  coarse_field_s = coarse_field_bck(nx,ny-1)
	ENDIF
        IF (ny==ny_bins) THEN ! north
	  coarse_field_n = miss
	ELSE
	  coarse_field_n = coarse_field_bck(nx,ny+1)
	ENDIF
      ENDIF
    
      coarse_field(nx,ny)=0
      no_of_neighbours=0

      ! west neighbour defined?
      IF (coarse_field_w .GT. miss+1) THEN
        coarse_field(nx,ny) = coarse_field(nx,ny)+coarse_field_w
        no_of_neighbours=no_of_neighbours+1
      ENDIF

      ! right neighbour defined?
      IF (coarse_field_e .GT. miss+1) THEN
        coarse_field(nx,ny) = coarse_field(nx,ny)+coarse_field_e
        no_of_neighbours=no_of_neighbours+1
      ENDIF

      ! south neighbour defined?
      IF (coarse_field_s .GT. miss+1) THEN
        coarse_field(nx,ny) = coarse_field(nx,ny)+coarse_field_s
        no_of_neighbours=no_of_neighbours+1
      ENDIF

      ! north neighbour defined?
      IF (coarse_field_n .GT. miss+1) THEN
        coarse_field(nx,ny) = coarse_field(nx,ny)+coarse_field_n
        no_of_neighbours=no_of_neighbours+1
      ENDIF

      ! is at least one neighbour defined?
      IF (no_of_neighbours .GT. 0) THEN
	coarse_field(nx,ny) = coarse_field(nx,ny)/no_of_neighbours
      ELSE
	coarse_field(nx,ny) = miss
	iterate = .TRUE.
      ENDIF
      
    ENDDO
  ENDDO
  
  IF (only_missing_values) THEN
    coarse_field = 0.
    RETURN
  ENDIF
  
ENDDO ! iteration

RETURN

END SUBROUTINE
