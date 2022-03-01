
using Dates

function startinfo()
    println("")
    println("Start : ",Dates.Date(Dates.now())," ",Dates.Time(Dates.now()))
    println("===========================================================")
    println("-----------------------------------------------------------")
    println("  The ParaFrags (Julia) code   ")
    println("-----------------------------------------------------------")
    println("===========================================================")
    println("")
end

function finishinfo()
    println("")
    println("===========================================================")
    println("-----------------------------------------------------------")
    println("  Done the ParaFrags calculations   ")
    println("-----------------------------------------------------------")
    println("===========================================================")
    println("Finish : ",Dates.Date(Dates.now())," ",Dates.Time(Dates.now()))
    println("")
end

