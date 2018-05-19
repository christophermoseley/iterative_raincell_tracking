#!/bin/ksh
set -ex

typeset -Z2 iteration

# Set in here: Which iteration should be done?
iteration=1

# identify objects
if [ $iteration -eq 1 ] ; then
   ./irt_objects_v1.x 1
else
   ./irt_objects_v1.x 2
fi

# generate coarse velocity field
./irt_advection_field_v1.x

# build tracks
./irt_tracks_v1.x

# generate field of track IDs
sort -n -k2 irt_tracks_nohead_output.txt > irt_tracks_sorted.txt
./irt_trackmask_v1.x

# Save output files.
# Otherwise they will be overwritten when the next iteration is performed
cp irt_objects_output.txt irt_objects_output_it_${iteration}.txt
#cp irt_objects_mask.srv irt_objects_mask_it_${iteration}.srv #comment in if needed
cp irt_advection_field.srv irt_advection_field_it_${iteration}.srv
cp irt_tracks_output.txt irt_tracks_output_it_${iteration}.txt 
cp irt_tracks_mask.srv irt_tracks_mask_it_${iteration}.srv

exit
