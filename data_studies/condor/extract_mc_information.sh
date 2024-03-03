#!/bin/bash

# source /home/hunt-stokes/rat4_env.sh
source /home/carpinteiroinacio/RAT/env_rat-7-0-14-recoord.sh

cd /home/hunt-stokes/multisite_analysis

python3 extract_mc_information.py ${RUN_NUMBER} ${ISOTOPE}