#!/bin/bash

source /home/hunt-stokes/rat_fresh.sh

cd /data/snoplus3/hunt-stokes/multisite_clean/data_studies/scripts/background_model

python3 extract_mc_information.py ${RUN_NUMBER} ${ISOTOPE}