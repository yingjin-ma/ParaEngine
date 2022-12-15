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
# 2022.7.14  Both Slurm (DongSheng-No.1) and LSF (ERA) HPC clusters are supported
# 2022.7.19  Start to support Gaussian and LSF
# 2022.8.10  More for Gaussian 

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

# Slurm case 
try
    global NN = parse(Int, ENV["SLURM_NNODES"])
    global IFSLURM = true 
catch err
    println("ENV(SLURM_NNODES) is not found ")
    global NN = 1 
    global IFSLURM = false
end

println(" == ",NN)

# LSF case 
try
    global LSB_MCPU_HOSTS = ENV["LSB_MCPU_HOSTS"]
    global IFLSF = true
catch err
    println("ENV(LSB_MCPU_HOSTS) is not found ")
    global NN = 1
    global IFLSF = false
end

if IFLSF
    Ntmp = length(split(LSB_MCPU_HOSTS))
    #println(" =Ntmp= ",Ntmp)
    global NN = Int(Ntmp/2)
end 
#println(" == ",NN)
flush(stdout)

gentask()
distributingtask(NN)

ISPAWN = 1
println("NN : ",NN, " IFSLURM : ",IFSLURM, "IFLSF", IFLSF,  " ISPAWN : ",ISPAWN, " (1:@spawn 2:@spawnat)")
flush(stdout)

println(" workdir : ", workdir)
if isfile(workdir)
    run(`rm $(workdir)"/hostlock"`)
end 

# Monitor the states in each nodes
police = @spawnat 2 nodes_monitor(workdir,"hostlock")
@async fetch(police)

if IFSLURM
    runtask(NN,IFSLURM,ISPAWN)
end 
if IFLSF
    LSB_JOBID=ENV["LSB_JOBID"]    
    cp(string("hosts.",LSB_JOBID),"nodelist",force=true)
    runtask(NN,IFLSF,ISPAWN)
end 


println("After nodes_monitor")

finishinfo()

