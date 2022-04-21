
target=ARGS[1]
outlog=ARGS[2]

if !isdir(target)
    println("The target suit $(target) is not exist")
    exit("Stopped. Reason: $(target) is not exist.")
else
    println("Monitor folder : ", target )
end

filelist = readdir(target)
println("")
println(filelist)
println("")

open(outlog,"w") do wlog
    i = 0
    for ifile in filelist
        i      = i + 1  
        ifile  = string(target,"/",ifile)
        println("ifile : ",ifile)
        nlines = countlines(ifile)
        println(wlog,i," ",nlines)
    end
end 

