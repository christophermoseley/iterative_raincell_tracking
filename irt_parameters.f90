MODULE irt_parameters

INTEGER, PARAMETER    :: domainsize_x = 960
INTEGER, PARAMETER    :: domainsize_y = 960

LOGICAL, PARAMETER    :: lperiodic = .TRUE.

INTEGER, PARAMETER    :: n_fields = 1   ! number of additional averaging fields

! bins of coarse velocity field
INTEGER, PARAMETER    :: time_steps = 120    ! total number of timesteps
INTEGER, PARAMETER    :: nt_bins = 5         ! 2 hourly
INTEGER, PARAMETER    :: nx_bins = 2
INTEGER, PARAMETER    :: ny_bins = 2

REAL, PARAMETER       :: threshold = 1.0            ! for rain intensity
REAL, PARAMETER       :: minimum_size = 4.       ! events smaller than that will be sorted out

REAL, PARAMETER       :: termination_sensitivity=0.5      ! Choose value between 0.0 and 1.0

REAL, PARAMETER       :: max_velocity = 10.   ! adjust acordingly
                                              ! velocities>max_velocity will be ignored to remove outliers
! define a minimal number of cells required for a coarse grained coordinate to be evaluated 
! if there are less, missing value will be assigned to that coarse cell
INTEGER, PARAMETER    :: min_cells = 10

INTEGER, PARAMETER    :: max_no_of_cells=5000  ! buffer size, increase if necessary
INTEGER, PARAMETER    :: max_no_of_tracks=5000    ! buffer size, increase if necessary
INTEGER, PARAMETER    :: max_length_of_track=1000  ! buffer size, increase if necessary

REAL, PARAMETER       :: miss=-9999.           ! value<miss ==> missing_value

END MODULE irt_parameters
