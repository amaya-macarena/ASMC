#!/bin/bash
source /home/mamaya/load_python.scr

# To submit to the AMD nodes (nodes 01-10), use OpenLava.
# Submit the job writing: bsub < submit_python.sh 

# Uncomment this if you want to use the high-memory node (node15)
# #BSUB -q highmem 

# Uncomment this if you want to use nodes 11 to 14
 # #BSUB -q intel 

# Job name
#BSUB -J 20steps_test1B
# Redirect screen output to output.txt
#BSUB -o /home/mamaya/test_AIS/Uncertainty_quatinfication/CM1/20steps_B/test1/20steps_test1B.out

# Error File
#BSUB -e /home/mamaya/test_AIS/Uncertainty_quatinfication/CM1/20steps_B/test2/20steps_test1B.err

# Mail notification
#BSUB -u macarena.amaya@unil.ch

# Number of processes. <- number of tasks=cpu_cores to reserve for the job
#BSUB -n 5,20

#BSUB -m "node02"

# Add the line for pytorch 1.3

export LD_PRELOAD=/soft/glibc/glibc-2.14/lib/libc-2.14.so

# Launch the Parrot batch job
python run_mcmc.py
#python Plotgraphs_DREAM_Output.py
#python Generate_lastmodelchain.py
