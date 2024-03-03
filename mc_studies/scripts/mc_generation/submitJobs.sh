#!/bin/bash

FILES=/home/carpinteiroinacio/RAT/miniProd_RAT-7-0-14_ASCI_RATHS/condor/jobSimulation_runMacroBiPo212*.sh

for f in $FILES
do
  echo "Processing $f file..."
  # take action on each file. $f store current file name
  #  qsub -q 'lipq' $f

  condor_submit -batch-name miniProd $f

  sleep 1
done