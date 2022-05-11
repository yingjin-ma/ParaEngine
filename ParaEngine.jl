# Self-written ParaEngine/ParaFrags (Julia) version
# 2021.3.14  Initial REM-MPI version refering Xpol, but it was totally rewritten (due to X-pol is not easy to debug for us)
#            REM-MPI was also used as ParaFrags as the computating driver for GridMol2.0
# 2021.11.28 Testing version use the Fortran + Shell style
# 2022.2.28  Initial pass for Julia version, using Julia to refactering the ParaFrags 
#            A lot of works need to do ...  
# 2022.3.14  Both dynamic and static load-balancing ways can work
# 2022.3.17  Multi-nodes case can work
# 2022.4.12  Molecular suits can work
# 2022.4.14  Rename as ParaEngine, because not only fragments can be calculated


using Distributed
using BenchmarkTools

include("PrintParaEngine.jl")
include("GlobalParameters.jl")
include("ReadInputs.jl")
include("QCTasks.jl")
include("Monitor.jl")

# launch worker processes
addprocs(0)

startinfo()
slurminfo()

readinp(ARGS[1])

try
    global NN = parse(Int, ENV["SLURM_NNODES"])
    global IFSLURM = true 
catch err
    println("ENV(SLURM_NNODES) is not found, use only 1 local process")
    global NN = 10 
    global IFSLURM = false
end
flush(stdout)

gentask()
distributingtask(NN)

ISPAWN = 1
println("NN : ",NN, " IFSLURM : ",IFSLURM," ISPAWN : ",ISPAWN, " (1:@spawn 2:@spawnat)")
flush(stdout)

println(" workdir : ", workdir)
if isfile(workdir)
    run(`rm $(workdir)"/hostlock"`)
end 

# Monitor the states in each nodes
police = @spawnat 2 nodes_monitor(workdir,"hostlock")
@async fetch(police)

runtask(NN,IFSLURM,ISPAWN)

println("After nodes_monitor")

finishinfo()

