#!/bin/bash
source /home/hunt-stokes/rattus_novus.sh
cd /data/snoplus3/hunt-stokes/multisite_clean/mc_studies/scripts

./create_pdfs ${ISOTOPE} ${RUN_NUMBER} ${FV_CUT} ${Z_CUT} ${INPUT} ${OUTPUT}