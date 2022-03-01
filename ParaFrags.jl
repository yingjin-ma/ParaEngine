# Self-written ParaFrags (Julia) version
# 2021.3.14  Initial REM-MPI version refering Xpol, but it was totally rewritten (due to X-pol is not easy to debug for us)
#            REM-MPI was also used as ParaFrags as the computating driver for GridMol2.0
# 2021.11.28 Testing version use the Fortran + Shell style
# 2022.2.28  Initial pass for Julia version, using Julia to refactering the ParaFrags 
#            A lot of works need to do ...  

using Distributed
using BenchmarkTools

include("PrintParaFrags.jl")
include("ReadInputs.jl")

# launch worker processes
addprocs(0)

startinfo()
slurminfo()

readinp("./input")

finishinfo()

