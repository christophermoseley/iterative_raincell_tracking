#!/bin/bash
set -ex

#cdo -f srv copy $DATAFILE irt_objects_input_00.srv

#cp irt_parameters.f90_gpm irt_parameters.f90
./compile.sh
./compile.sh

# iterate object files
./irt_objects_release.x 1
./irt_advection_field_release.x
cp irt_advection_field.srv irt_advection_field_it1.srv
for ITERATION in 2 3; do
echo iteration $ITERATION
./irt_objects_release.x 2
./irt_advection_field_release.x
cp irt_advection_field.srv irt_advection_field_it${ITERATION}.srv
done

cp irt_objects_output.txt irt_objects_output_it3.txt

# build tracks
./irt_tracks_release.x

# generate field of track IDs
sort -n -k2 irt_tracks_nohead_output.txt > irt_tracks_sorted.txt
./irt_trackmask_release.x

# Save output files.
# Otherwise they will be overwritten when the next iteration is performed
cp irt_tracks_output.txt irt_tracks_output_it3.txt 
cdo -f nc copy irt_tracks_mask.srv irt_tracks_mask_it3.nc

exit

