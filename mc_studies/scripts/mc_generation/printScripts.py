import numpy as np
import os

type = ['BiPo212']
runs = np.loadtxt("/home/hunt-stokes/multisite_analysis/mc_list.txt", dtype = int)

for j in range(len(type)):
  for i in range(len(runs)):

    outName = '/home/carpinteiroinacio/RAT/miniProd_RAT-7-0-14_ASCI_RATHS/sh/runMacro'+type[j]+'_%d.sh' % ( i )

    # Create output file
    outputfile = open(outName, 'w')
    os.chmod(outName, 0o0777)
    print('Creating macro ', outName)

    outputfile.write('#!/bin/bash')
    outputfile.write('\n')
    outputfile.write('source /home/carpinteiroinacio/RAT/env_rat-7-0-14-recoord.sh\n')
    outputfile.write('cd /home/carpinteiroinacio/RAT/miniProd_RAT-7-0-14_ASCI_RATHS/mac/')
    outputfile.write('\n')
    outputfile.write('rat -n '+str(runs[i])+' -N 1000 '+type[j]+'.mac -o /data/snoplus3/SNOplusData/production/miniProd_RAT-7-0-14_ASCI_RATHS_newRecoordination/simulation'+type[j]+'_'+str(runs[i])+'\n')

    outputfile.write('\n')
    outputfile.write('exit 0')