using Printf

mutable struct FRAGS
      idx::Int
     name::String
     node::String
 nodeidx0::Int
 nodeidx1::Int
 nodeidx2::Int
  timeCPU::Float64
  timeELP::Float64
end

# e.g.
#       julia  Timing_analysis_sequeue.jl  ../ParaEngine-git-Gaussian/P38-M062x-ALL   MLSLB-ONLY

folder=ARGS[1]
target=ARGS[2]

if !isdir(folder)
    println("The folder $(folder) is not exist")
    exit("Stopped. Reason: $(folder) is not exist.")
else
    println("Monitor folder : ", folder )
end

if !isdir(target)
    println("The target suit $(target) is not exist")
    exit("Stopped. Reason: $(target) is not exist.")
else
    println("Monitor folder : ", target )
end

filelist = readdir(folder)
println("")
#println(filelist)
println("")

ifrag = 0
names = []
fraglist = [] 
for ifile in filelist
    push!(names,split(ifile,".")[1])
    id1 = string(folder,"/",ifile)
    global ifrag = ifrag + 1
    push!(fraglist,FRAGS(ifrag,names[ifrag],"",0,0,0,0.0,0.0))
end

#println(fraglist)

ftime1=string(target,"/Timing1.txt")
ftime2=string(target,"/Timing2.txt")
ftime3=string(target,"/Timing3.txt")

# Timing2.txt
open(ftime2,"r") do read1
    global nodetask1=[]
    global nodetask2=[]
    while !eof(read1)
        line=readline(read1)
        sline=split(line)
        ntmp=length(sline)
        if ntmp > 1
           TL1=[] 
           for i in 2:length(sline)
               push!(TL1,sline[i])
           end  
           push!(nodetask1,sline[1]) 
           push!(nodetask2,TL1) 
        end 
    end
    # println(nodetask1) 
end

# Timing1.txt
open(ftime1,"r") do read1
    global nodetask3=[]
    while !eof(read1)
        line=readline(read1)
        sline=split(line)
        ntmp=length(sline)
        if ntmp > 1
           TL1=[]
           for i in 2:length(sline)
               dv = parse(Float64,sline[i])
               push!(TL1,dv)
           end
           push!(nodetask3,TL1)
        end
    end
    println(nodetask3)
end 

# Timing3.txt
open(ftime3,"r") do read1
    global nodetask4=[]
    while !eof(read1)
        line=readline(read1)
        sline=split(line)
        ntmp=length(sline)
        if ntmp > 1
           TL1=[]
           for i in 2:length(sline)
               dv = parse(Float64,sline[i])
               push!(TL1,dv)
           end
           push!(nodetask4,TL1)
        end
    end
    println(nodetask4)
end

if length(nodetask4) == length(nodetask3) == length(nodetask2) == length(nodetask1)
     println(" ==>  DIM match ")
     println(" ==>  DIM match ")
     println(" ==>  DIM match ")
end 

for i in 1:length(fraglist)
    idx1name = string(fraglist[i].name,".log")
    idx0 = 0 
    idx1 = 0 
    for j in 1:length(nodetask2)
        #println(idx1name)
        #println(nodetask2[j])
        idx1pos = findfirst(==(idx1name),nodetask2[j])
        if idx1pos != nothing  
           idx0 = j
           idx1 = idx1pos
           #println("break ",j, "  ", idx1pos)
           break
        end 
    end
    # println("idx0, idx1 : ", idx0, " ", idx1)
    fraglist[i].timeCPU = nodetask3[idx0][idx1]
    fraglist[i].timeELP = nodetask4[idx0][idx1]
    fraglist[i].node    = nodetask1[idx0]

    fraglist[i].nodeidx0    =  idx0
    fraglist[i].nodeidx1    =  idx1

    #println(fraglist[i])
end 

nlbdef = trunc(Int,length(fraglist)/length(nodetask1)+1)

println("nlbdef : ",nlbdef )

MNelems = Array{FRAGS,2}(undef,length(nodetask1),10*(nlbdef)) 

println(" ============== auto check when no LB file shoud be used ==============  ")
icx = []
for i in 1:length(fraglist)
    if fraglist[i].nodeidx1 == 1
        push!(icx,i)
    end 
end 

for i in 1:length(icx)
    println(icx[i]," ",fraglist[icx[i]])
end

idxnode = Dict(fraglist[icx[i]].node => fraglist[icx[i]].nodeidx0 for i in 1:length(icx))

println(idxnode)

for i in 1:length(fraglist)
    fraglist[i].nodeidx2 = idxnode[fraglist[i].node]
    println(fraglist[i])
end 

println(length(nodetask1), " - ",(nlbdef))
for i in 1:length(fraglist)
    idx0=fraglist[i].nodeidx0
    idx1=fraglist[i].nodeidx1
    idx2=fraglist[i].nodeidx2
    #println("idx1 ,idx2 : ", idx1, " ",idx2 )

    if idx0 == idx2                           # No DLB was involved 
        MNelems[idx2,idx1]=fraglist[i]
    else                                      # DLB was involved
        offset = (nlbdef)
        println(offset)
        itmp = -1
        while true
            try 
                tmp = MNelems[idx2,idx1+offset]
                offset = offset + 1
            catch err 
                itmp = 1
            end 
            if itmp == 1
                break
            end
        end 
        MNelems[idx2,idx1+offset]=fraglist[i]
        # println(offset)
    end 
end 

open("Opted_Tasks.txt","w") do write1
open("Opted_Tasks_GNUPLOT-CPU.dat","w") do write2
open("Opted_Tasks_GNUPLOT-ELP.dat","w") do write3
lnmax=0
for i in 1:length(nodetask1)
    iln=0 
    for j in 1:10*(nlbdef)
        try 
            if MNelems[i,j].timeCPU > 0 
                iln=iln+1
            end 
        catch err
        end
    end
    if iln > lnmax
        lnmax = iln
    end
end
println(" lnmax ",lnmax) 
print(write2,"node ")
print(write3,"node ")
for i in 1:lnmax
    print(write2," task$(i)")
    print(write3," task$(i)")
end 
println(write2,"")
println(write3,"")
for i in 1:length(nodetask1)
    print(write2,MNelems[i,1].node) 
    print(write3,MNelems[i,1].node)
    icount = 0 
    for j in 1:10*nlbdef
        try 
          println(write1,i," ", j, " - ",MNelems[i,j]) 
          print(write2," ",MNelems[i,j].timeCPU)
          print(write3," ",MNelems[i,j].timeELP)
          #@printf(write2," %.1f",MNelems[i,j].timeCPU)
          #@printf(write3," %.1f",MNelems[i,j].timeELP)
          icount = icount +1
        catch err
          println(i," ",j," ") 
          # print(write2," "," - ")
        end
    end
    for j in icount+1:lnmax
        print(write2," 0.0")
        print(write3," 0.0")
    end 
    println(write2,"") 
    println(write3,"") 
end
end
end
end




exit(0)

