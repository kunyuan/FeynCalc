#!/bin/sh
julia fock_static.jl
mkdir data
cp -r ./groups* ./data
cp -r ./*.exe ./data
cp parameter ./data
cp sigma.* ./data
cd data

sbatch -N1 -p ccq -t 7-00:00:00 -J 3dueg disBatch.py ../task
#sbatch -N1 -p ccq -t 7-00:00:00 -J 3dueg disBatch.py ../task
#sbatch -N1 -p ccq -t 7-00:00:00 -J rs1 disBatch.py ../task
#sbatch -N1 -p ccq -t 7-00:00:00 -J rs1 disBatch.py ../task
#sbatch -n 28 -p gen -t 7-00:00:00 --ntasks-per-node 28 --exclusive --wrap="./feyncalc.exe"
#sbatch -n 28 -p gen -t 7-00:00:00 --ntasks-per-node 28 --exclusive --wrap="disBatch.py ../task"
#sbatch -a 0-28 -p ccq -t 7-00:00:00 -N1 --wrap="./feyncalc.exe"
#sbatch -N1 -p genx -t 7-00:00:00 disBatch.py ../task


#sbatch -N 1 -a 0-28 -p ccq -t 7-00:00:00 -n 1 --ntasks-per-node=28 --wrap="./feyncalc.exe"
#sbatch -a 0-64 -p ccq -t 7-00:00:00 -n 1 --wrap="./feyncalc.exe"
