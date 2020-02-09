#!/bin/bash
#SBATCH --job-name=0529          # Specify job name
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1             # Specify max. number of tasks to be invoked
#SBATCH --cpus-per-task=1      # Specify number of CPUs per task
#SBATCH --time=08:00:00        # Set a limit on the total run time
#SBATCH --account=bm0982       # Charge resources on this project account
#SBATCH --output=track0529.o%j    # File name for standard output
#SBATCH --error=track0529.o%j     # File name for standard error output
#shopt -s expand_aliases
#source ~/.bashrc

set -ex

#typeset -Z2 iteration

#caselist=(tpe20060508 tpe20060721 tpe20070830 tpe20090707 tpe20090827 tpe20100629 tpe20100630 tpe20100803 tpe20100912)
#caselist=(tpe20060508)
total=${#caselist[*]}

#for ((i=0;i<=$(($total-1));i++)); do
#for ((i=0;i<=0;i++)); do
#echo $i
#echo ${caselist[$i]}

   ### Definition
   SUFFIX=tpe20060508
   PFAD=`pwd`
#   EXPT=cln 
   EXPT=cln
 
   cd $PFAD

   ### File Preparation ###
#   cd /home/C.brobbtyh/tracking_RR/${SUFFIX}"cln"
   #cd /home/C.brobbtyh/tracking_RR/${SUFFIX}"dty"
   #cp ${PFAD}/compile.sh .
   #cp ${PFAD}/irt_advection_field_release.f90 .
   #cp ${PFAD}/irt_objects_release.f90 .
   #cp ${PFAD}/irt_parameters.f90 .
   #cp ${PFAD}/irt_parameters.f90_yuhung .
   #cp ${PFAD}/irt_parameters.mod .
   #cp ${PFAD}/IRT_tracking_documentation.pdf .
   #cp ${PFAD}/irt_trackmask_release.f90 .
   #cp ${PFAD}/irt_tracks_release.f90 .
   cp irt_parameters.f90_yuhung irt_parameters.f90

   ### Compile ###
   ./compile.sh

   #cd $PFAD

   #cd /home/C.brobbtyh/tracking_RR/${SUFFIX}"cln"
   cd /data/moseley/yuhung/

   ### Iterate Object Files ###
   #PFAD=`pwd`
   ${PFAD}/irt_objects_release.x 1
   ${PFAD}/irt_advection_field_release.x
   cp irt_advection_field.srv irt_advection_field_${SUFFIX}_it1.srv
   for ITERATION in 2 3; do
        echo iteration $ITERATION
        ${PFAD}/irt_objects_release.x 2
        ${PFAD}/irt_advection_field_release.x
        cp irt_advection_field.srv irt_advection_field_${SUFFIX}_it${ITERATION}.srv
   done
   cp irt_objects_output.txt irt_objects_output_${SUFFIX}_it3.txt

   ### Build Tracks ###
   ${PFAD}/irt_tracks_release.x

   ### Generate Field of Track IDs ###
   sort -n -k2 irt_tracks_nohead_output.txt > irt_tracks_sorted.txt
   ${PFAD}/irt_trackmask_release.x

   ### add links to preceding and successing tracks to track headers
   ${PFAD}/irt_tracklinks_release.x

   ### Save Output Files ###
   # Otherwise they will be overwritten when the next iteration is performed
   cp irt_tracks_output.txt irt_tracks_output_${SUFFIX}_it3.txt 
   cdo -f nc copy irt_tracks_mask.srv irt_tracks_mask_${SUFFIX}_it3.nc
   cdo -f nc copy irt_objects_mask.srv irt_objects_mask_${SUFFIX}_${EXPT}.nc
   cdo -f nc copy irt_tracks_mask.srv irt_tracks_mask_${SUFFIX}_${EXPT}.nc

   ### Change File Name
   cp irt_tracklinks_output.txt irt_tracklinks_${SUFFIX}_${EXPT}.txt 
   cp irt_objects_output.txt irt_objects_${SUFFIX}_${EXPT}.txt
   
#done 

exit
