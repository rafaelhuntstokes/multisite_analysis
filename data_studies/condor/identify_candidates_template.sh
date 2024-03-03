#!/bin/bash
source /home/hunt-stokes/rat_fresh.sh
cd /data/snoplus3/hunt-stokes/multisite_clean/data_studies/scripts

./identify_solar_candidates ${RUN_NUMBER} ${FV_CUT} ${Z_CUT} ${INPUT} ${OUTPUT}
