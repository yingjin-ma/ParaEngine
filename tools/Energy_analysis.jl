
target=ARGS[1]
outlog=ARGS[2]
nnodes=parse(Int64,ARGS[3])

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

times=[]
percentage=0.0
open(outlog,"w") do wlog
    neff = 0
    i = 0
    for ifile in filelist
        i      = i + 1 
        ifile  = string(target,"/",ifile)
        println("ifile : ",ifile)
        nlines = countlines(ifile)
        neff   = neff + nlines 
        println(wlog,i," ",nlines)
        push!(times,nlines)
    end
    promote(i,nnodes,neff)
    nidl = i * nnodes
    percentage = neff / nidl
    println(wlog,"")
    println(wlog,"The effective   mins (active) are       : ", neff)
    println(wlog,"The ideal nodes mins (active) should be : ", nidl)
    println(wlog,"The utilization percentage for the nodes in life cycle : ", percentage )
    println(wlog,"")
end 

println("times :", times)
nt=length(times)
for i in 1:length(times)-1
   if times[nt-i] < times[nt-i+1]
      times[nt-i] = times[nt-i+1]
   end  
end 

percentage=0.0
outlog=string(outlog,"_correct")
open(outlog,"w") do wlog
    neff = 0
    for i in 1:length(times)
        nlines=times[i]
        neff   = neff + nlines
        println(wlog,i," ",nlines)
        push!(times,nlines)
    end
    promote(nt,nnodes,neff)
    nidl = nt * nnodes
    percentage = neff / nidl
    println(wlog,"")
    println(wlog,"The effective   mins (active) are       : ", neff)
    println(wlog,"The ideal nodes mins (active) should be : ", nidl)
    println(wlog,"The utilization percentage for the nodes in life cycle : ", percentage )
    println(wlog,"")
end


