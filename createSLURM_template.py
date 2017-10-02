#!/usr/bin/python

# -----------------------------------------------------------------------------------
# Create PBS template
# Author:    Kelly Quek
# Date started:  30-Jul-2015
# Last updated:  01-Aug-2015
# Objective:  This is a script to generate PBS submission file when given a list of patient id[1] and a template[2]
# -----------------------------------------------------------------------------------

import string
import sys        

inFile = open(sys.argv[1], 'r')

key = ['patientID']                             # Name of the key; this is subject to changes depending on how you name your variable in the template file

d = {}                                          # Create empty dictionary

fh = open (sys.argv[2], 'r')
template = fh.read() 
tpl = string.Template(template)

for id in inFile :              
    patients = id.split()                       # Iterate each patient ID
    
    for v, k in zip(patients, key) :            # Assign key to each of the patient ID
        d[k] = v                                # example - {'patientID':'ADI_0050'}    
        result = tpl.safe_substitute(d)         # Substitute $patientID (notation for key) in template with value
    
    for line in patients :                      # Write output using patient id[1]
        with open(line.rstrip() + '_STAR.SLURM', 'w') as outFile :
            outFile.write(result)
            outFile.close()
            
