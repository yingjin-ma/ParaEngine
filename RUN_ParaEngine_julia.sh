#!/bin/bash

#SBATCH -J Julia_TEST
#SBATCH -N 1
#SBATCH --ntasks-per-node=5
#SBATCH --gres=dcu:4
#SBATCH --mem=100G
#SBATCH --gres-flags=disable-binding
#SBATCH --exclusive
#SBATCH -p debug
#SBATCH -o out-julia
#SBATCH -e err-julia

module  use /public/home/mayj/Quantum_modulefiles
module load apps/julia/1.7.2

module load IntelMKL/2021.1.1
module load apps/nwchem-DCU/6.8.1/hpcx-2.7.4-intel2017

date1=$(date +%Y%m%d%H%M%S)
#echo $date1

scontrol show hostname > nodelist

julia --machine-file=nodelist ParaEngine.jl > out-$date1

