
using Dates

function slurminfo()
    println(" ")
    println("-----------------------------------------------------------")
    println("-- LSF/Slurm info : ---------------------------------------")
    println("-----------------------------------------------------------")
    println("Number of processes : ", nprocs())
    println("Number of workers   : ", nworkers())
    println("-----------------------------------------------------------")
    println(" ")
end

function startinfo()
    println("")
    println("Start : ",Dates.Date(Dates.now())," ",Dates.Time(Dates.now()))
    println("===========================================================")
    println("-----------------------------------------------------------")
    println("  The ParaEngine (Julia) code   ")
    println("-----------------------------------------------------------")
    println("===========================================================")
    println("")
end

function finishinfo()
    println("")
    println("===========================================================")
    println("-----------------------------------------------------------")
    println("  Done the ParaEngine calculations   ")
    println("-----------------------------------------------------------")
    println("===========================================================")
    println("Finish : ",Dates.Date(Dates.now())," ",Dates.Time(Dates.now()))
    println("")
end

