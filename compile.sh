#!/bin/bash
set -ex

#COMPILE_COMMAND='ifort -no-wrap-margin'
COMPILE_COMMAND='gfortran -O0'

# clean:
set +e
rm irt_parameters.mod
rm irt_objects_release.x
rm irt_advection_field_release.x
rm irt_tracks_release.x
rm irt_trackmask_release.x
rm irt_agemask_release.x
rm irt_tracklinks_release.x
set -ex

# compile:
${COMPILE_COMMAND} -o irt_objects_release.x irt_parameters.f90 irt_objects_release.f90 
${COMPILE_COMMAND} -o irt_advection_field_release.x irt_advection_field_release.f90 irt_parameters.f90
${COMPILE_COMMAND} -o irt_tracks_release.x irt_tracks_release.f90 irt_parameters.f90
${COMPILE_COMMAND} -o irt_trackmask_release.x irt_trackmask_release.f90 irt_parameters.f90
${COMPILE_COMMAND} -o irt_agemask_release.x irt_agemask_release.f90 irt_parameters.f90
${COMPILE_COMMAND} -o irt_tracklinks_release.x irt_tracklinks_release.f90 irt_parameters.f90

exit
