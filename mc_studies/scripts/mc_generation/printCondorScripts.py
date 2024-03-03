import glob

dataSet = []

for file in glob.glob("/home/carpinteiroinacio/RAT/miniProd_RAT-7-0-14_ASCI_RATHS/sh/runMacro*.sh"):
    dataSet.append(file)
print(dataSet)
print(len(dataSet))
name = "/home/carpinteiroinacio/RAT/miniProd_RAT-7-0-14_ASCI_RATHS/sh/"
print(dataSet[0][len(name):])
for t in range(len(dataSet)):

    outName = '/home/carpinteiroinacio/RAT/miniProd_RAT-7-0-14_ASCI_RATHS/condor/jobSimulation_'+dataSet[t][len(name):-3]+'.sh'

    # Create output file
    outputfile = open(outName, 'w')
    print('Creating script ', outName)

    outputfile.write('executable = /home/carpinteiroinacio/RAT/miniProd_RAT-7-0-14_ASCI_RATHS/sh/'+dataSet[t][len(name):]+'\n')
    outputfile.write('getenv = false')
    outputfile.write('\n')
    outputfile.write('output = /home/carpinteiroinacio/RAT/miniProd_RAT-7-0-14_ASCI_RATHS/logs/results.output.$(ClusterId)\n')
    outputfile.write('error = /home/carpinteiroinacio/RAT/miniProd_RAT-7-0-14_ASCI_RATHS/logs/results.error.$(ClusterId)\n')
    outputfile.write('log = /home/carpinteiroinacio/RAT/miniProd_RAT-7-0-14_ASCI_RATHS/logs/results.log.$(ClusterId)\n')
    outputfile.write('queue 1\n')