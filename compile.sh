#!/bin/ksh
set -ex

COMPILE_COMMAND='ifort -no-wrap-margin'

${COMPILE_COMMAND} -o irt_objects_v1.x irt_objects_v1.f90 irt_parameters.f90
${COMPILE_COMMAND} -o irt_advection_field_v1.x irt_advection_field_v1.f90 irt_parameters.f90
${COMPILE_COMMAND} -o irt_tracks_v1.x irt_tracks_v1.f90 irt_parameters.f90
${COMPILE_COMMAND} -o irt_trackmask_v1.x irt_trackmask_v1.f90 irt_parameters.f90

exit
