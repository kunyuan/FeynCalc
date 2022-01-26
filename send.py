#!/usr/bin/python
import random
import os
import sys

##### Modify parameters here  ###############
# Cluster="PBS"
Cluster = "local"
# Cluster = "condor"

############################################

rootdir = os.getcwd()
execute = "feyncalc.exe"

assert len(sys.argv) == 2, "Number of jobs is needed."

Number = int(sys.argv[1])
print "Creating {0} jobs...".format(Number)

PIDList = range(Number)

# if int(para[-2])==0:
#     title="freq"
# elif int(para[-2])==1:
#     title="eqTime"
# else:
#     print "Not yet implemented!"
#     break

homedir = os.getcwd() + "/data"
if(os.path.exists(homedir) != True):
    os.system("mkdir "+homedir)

os.system("cp -r groups* "+homedir)
os.system("cp {0} {1}".format(execute, homedir))
os.system("cp *.data "+homedir)
# os.system("cp reweight.data "+homedir)
os.system("cp parameter "+homedir)

outfilepath = homedir+"/outfile"
if(os.path.exists(outfilepath) != True):
    os.system("mkdir "+outfilepath)
jobfilepath = homedir+"/jobfile"
if(os.path.exists(jobfilepath) != True):
    os.system("mkdir "+jobfilepath)

for pid in PIDList:

    ### terminal output goes here #############
    outfile = "_out"+str(pid)
    ### job file to submit to cluster  ########
    jobfile = "_job"+str(pid)+".sh"

    if Cluster == "local":
        os.chdir(homedir)
        os.system("./"+execute+" > "+outfilepath+"/"+outfile+" &")
        os.chdir("..")
    elif Cluster == "condor":
        with open(jobfilepath+"/"+jobfile, "w") as fjob:
            fjob.write("executable = {0}\n".format(execute))
            fjob.write("output ={0}/{1}\n".format(outfilepath, outfile))
            fjob.write("initialdir ={0}\n".format(homedir))
            fjob.write("queue")

        os.chdir(homedir)
        os.system("condor_submit {0}/{1}".format(jobfilepath, jobfile))
        os.system("rm "+jobfilepath + "/"+jobfile)
        os.chdir("..")
    elif Cluster == "PBS":
        with open(jobfilepath+"/"+jobfile, "w") as fjob:
            fjob.write("#!/bin/sh\n"+"#PBS -N "+jobfile+"\n")
            fjob.write("#PBS -o "+homedir+"/Output\n")
            fjob.write("#PBS -e "+homedir+"/Error\n")
            fjob.write("#PBS -l walltime=2000:00:00\n")
            fjob.write("echo $PBS_JOBID >>"+homedir+"/id_job.log\n")
            fjob.write("cd "+homedir+"\n")
            fjob.write("./"+execute+" > "+outfilepath+"/"+outfile)

        os.chdir(homedir)
        os.system("qsub "+jobfilepath + "/"+jobfile)
        os.system("rm "+jobfilepath + "/"+jobfile)
        os.chdir("..")
    else:
        print("I don't know what is {0}".format(Cluster))
        break

print("Jobs manage daemon is ended")
sys.exit(0)
