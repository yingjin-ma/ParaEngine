# Self-written ParaFrags (Julia) version
# 2021.3.14  Initial REM-MPI version refering Xpol, but it was totally rewritten (due to X-pol is not easy to debug for us)
#            REM-MPI was also used as ParaFrags as the computating driver for GridMol2.0
# 2021.11.28 Testing version use the Fortran + Shell style
# 2022.2.28  Initial pass for Julia version, using Julia to refactering the ParaFrags 
#            A lot of works need to do ...  

using Distributed
using BenchmarkTools

include("PrintParaFrags.jl")
include("GlobalParameters.jl")
include("ReadInputs.jl")
include("QCTasks.jl")

# launch worker processes
addprocs(0)

startinfo()
slurminfo()

readinp("input")

gentask()
distributingtask()

try
    global NN = parse(Int, ENV["SLURM_NNODES"])
    global IFSLURM = true 
catch err
    println("ENV(SLURM_NNODES) is not found, use only 1 local process")
    global NN = 1 
    global IFSLURM = false
end

IFDYNA= false
println("NN : ",NN, " IFSLURM : ",IFSLURM," IFDYNAMIC : ",IFDYNA)

runtask(NN,IFSLURM,IFDYNA)

finishinfo()

