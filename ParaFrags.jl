# Self-written ParaFrags (Julia) version
# 2021.3.14  Initial REM-MPI version refering Xpol, but it was totally rewritten (due to X-pol is not easy to debug for us)
#            REM-MPI was also used as ParaFrags as the computating driver for GridMol2.0
# 2021.11.28 Testing version use the Fortran + Shell style
# 2022.2.28  Initial pass for Julia version, using Julia to refactering the ParaFrags 
#            A lot of works need to do ...  
# 2022.3.14  Both dynamic and static load-balancing ways can work
# 2022.3.17  Multi-nodes case can work


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

try
    global NN = parse(Int, ENV["SLURM_NNODES"])
    global IFSLURM = true 
catch err
    println("ENV(SLURM_NNODES) is not found, use only 1 local process")
    global NN = 1 
    global IFSLURM = false
end
flush(stdout)

gentask()
distributingtask(NN)

ISPAWN = 1
println("NN : ",NN, " IFSLURM : ",IFSLURM," ISPAWN : ",ISPAWN, " (1:@spawn 2:@spawnat)")
flush(stdout)

runtask(NN,IFSLURM,ISPAWN)

finishinfo()

