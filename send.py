#!/usr/bin/python
import random
import os
import sys

##### Modify parameters here  ###############
# Cluster="PBS"
Cluster = "local"
# Cluster="condor"

############################################
inlist = open("./inlist","r")
rootdir = os.getcwd()
# execute = "feyncalc_eqTime.exe"
# assert len(sys.argv)==2, "Number of jobs is needed."
# Number=int(sys.argv[1])

suffix = ["_eqTime", "_freq", "_freqTau", "_Ek"]

for index, eachline in enumerate(inlist):
    para = eachline.split()
    if len(para)==0:
        print ("All submitted!")
        break
    
    order= int(para[0])
    beta = float(para[1])
    rs   = float(para[2])
    lam  = float(para[4])
    kmin  = float(para[5])

    print ("Creating {0} jobs...".format(int(para[-2])))

    PIDList=range(int(para[-2]))
    title = suffix[int(para[-1])]
    execute = "feyncalc"+title+".exe"
    fname = "beta{0}_rs{1}_lam{2}".format(beta,rs,lam)
    homedir = os.getcwd() +"/"+fname+"_o{0}".format(order)+title
#    if kmin == 0.02:
#        homedir = os.getcwd() +"/"+fname+"_o{0}".format(order)+title+"_1"
#    else:
#        homedir = os.getcwd() +"/"+fname+"_o{0}".format(order)+title+"_2"

    if(os.path.exists(homedir) != True):
        os.system("mkdir "+homedir)
    else:
        print (homedir+" alreadly exists!")
        # homedir = homedir + "_6kf"
        # os.system("mkdir "+homedir)
        break

    with open("./parameter", "w") as file:
        parameters = ' '.join(para[:-2])
        file.write(parameters+"\n\n")
        file.write("#Order, Beta, rs, Mass2, Lambda, MinExtMom(*kF), MaxExtMom, TotalStep(*1e6)")

    os.system("cp -r groups* "+homedir)
    os.system("cp {0} {1}".format(execute, homedir))
    os.system("cp reweight.data "+homedir)
    os.system("cp parameter "+homedir)
    os.system("cp ph_10000_1e-8.dlr "+homedir)
    os.system("cp ./sigma/sigma3D_"+fname+".txt "+homedir+"/sigma3D.txt")

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
