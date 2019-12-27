#!/bin/bash
set -ex

#COMPILE_COMMAND='ifort -no-wrap-margin'
COMPILE_COMMAND='gfortran -O0'

${COMPILE_COMMAND} -o irt_objects_release.x irt_objects_release.f90 irt_parameters.f90
${COMPILE_COMMAND} -o irt_advection_field_release.x irt_advection_field_release.f90 irt_parameters.f90
${COMPILE_COMMAND} -o irt_tracks_release.x irt_tracks_release.f90 irt_parameters.f90
${COMPILE_COMMAND} -o irt_trackmask_release.x irt_trackmask_release.f90 irt_parameters.f90

exit
