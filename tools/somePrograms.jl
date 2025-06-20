#----对比两个文件，将相同的复制出来-----

"""global opt = "C:\\Users\\22103\\Desktop\\input\\opt"
global freq = "C:\\Users\\22103\\Desktop\\input\\freq"

filelist1 = readdir(opt)
filelist2 = readdir(freq)

opt1 = string(opt,"\\opt1")
mkdir(opt1)

for i in filelist1
    for j in filelist2
        if i == j
            ifiles = string(opt,"\\",i)
            ifiles1 = string(opt1,"\\",i)
            cp(ifiles,ifiles1)
            break
        end
    end
end
"""





#----按照顺序文件修改名字-----

#=
global opt1 = "C:\\Users\\22103\\Desktop\\input\\opt\\opt1"
global sFile = "C:\\Users\\22103\\Desktop\\input\\sorted_result.txt"
global freq = "C:\\Users\\22103\\Desktop\\input\\freq"

freqs = string(freq,"\\freqs")
mkdir(freqs)

mutable struct SORTS
      idx::Int
    names::String
end
sortlist = []
#global i = 0
open(sFile,"r") do molread
    global i = 0
    while !eof(molread)
        i = i + 1
        line = readline(molread)
        sp = split(line,",")
        spl = split(sp[2],"'")
        push!(sortlist,SORTS(i,spl[2]))
    end
end

count = length(sortlist)
"""for j in 1:count
    println(sortlist[j].idx,"  ",sortlist[j].names)
end"""

filelist3 = readdir(opt1)
for ifiel in filelist3
    ifiel1 = string(opt1,"\\",ifiel)
    oldRoot = string(freq,"\\",ifiel)
    #println("ifiel1: ",ifiel1)
    #println("oldRoot: ",oldRoot)
    open(ifiel1,"r") do molread
        #println("open")
        while !eof(molread)
            #println("while")
            #readline(molread)
            #readline(molread)
            readline(molread)
            readline(molread)
            line = readline(molread)
            sp1 = split(line,"/")
            if length(sp1) > 2
                #println("sp1 ")
                for j in 1:count
                    if occursin(sp1[3],sortlist[j].names)
                        sline = split(sortlist[j].names,".")
                        newName = string(sortlist[j].idx,"-",split(ifiel,".")[1],"-",sline[1],".gjf")
                        newRoot = string(freqs,"\\",newName)
                        cp(oldRoot,newRoot)
                        break
                    end
                end
            end
        end
    end
end

=#





#----将gjf转换成xyz格式------

"""global freq = "C:\\Users\\22103\\Desktop\\input\\freq\\freqs"

mutable struct ATOM
      elem::String
    coordx::Float64
    coordy::Float64
    coordz::Float64
end

global xyz = string(freq,"\\xyz")
mkdir(xyz)

filelist = readdir(freq)
global atomlist = []
println(filelist)

for ifile in filelist
    ifiles = string(freq,"\\",ifile)
    if isfile(ifiles)
        global name = split(ifile,".")[1]       
        global names =string(xyz,"\\",name,".xyz")
        #println("names: ",names)
        open(ifiles,"r") do molread
            println("0000")
            readline(molread)
            readline(molread)
            readline(molread)
            readline(molread)
            readline(molread)
            while !eof(molread)
                #println("1111")
                line = readline(molread)
                sp1 = split(line)
                #println("sp1: ",sp1)
                if length(sp1) == 4
                    #println("if")
                    dx = parse(Float64,sp1[2])
                    dy = parse(Float64,sp1[3])
                    dz = parse(Float64,sp1[4])
                    #println("push")
                    push!(atomlist,ATOM(sp1[1],dx,dy,dz))
                end
            end
        end
        counts = length(atomlist)
        open(names,"a") do steams
            #println("22222")
            println(steams,counts)
            println(steams,name)
        end
        for k in 1:counts
            open(names,"a") do steams
                #println("22222")
                println(steams,atomlist[k].elem," ",atomlist[k].coordx," ",atomlist[k].coordy," ",atomlist[k].coordz)
            end
        end
        for k in 1:counts
            pop!(atomlist)
        end
    end  
end
"""




#----提取红外和拉曼坐标----

#=
mutable struct COORD
    coordX::String 
    coordY::String
end

global freq = "C:\\Users\\22103\\Desktop\\output\\freq"

filelist = readdir(freq)

IR = string(freq,"\\IR")
Raman = string(freq,"\\Raman")
mkdir(IR)
mkdir(Raman)

global IRList = []
global RamanList = []
for ifile in filelist
    ifiles = string(freq,"\\",ifile)
    global name = split(ifile,".")[1]       
    global IRNames =string(IR,"\\IR_",name,".txt")   
    global RaNames =string(Raman,"\\Raman_",name,".txt")
    open(ifiles,"r") do molread
        while !eof(molread)
            line = readline(molread)
            if occursin("Simulated IR Spectrum",line)
                readline(molread)
                frags = true
                while(frags)
                    ats = readline(molread)                    
                    if length(ats) == 0 || occursin("IR Spectrum",ats)
                        #println("end ")
                        break
                    end
                    sp = split(ats)
                    if length(sp) > 1
                        push!(IRList,COORD(sp[1],sp[2]))
                    end                                        
                end
            end
        end
    end
    open(ifiles,"r") do molread
        while !eof(molread)
            line = readline(molread)
            if occursin("Simulated Raman Spectrum",line)
                readline(molread)
                frags = true
                while(frags)
                    ats = readline(molread)                    
                    if length(ats) == 0 || occursin("Raman Spectrum",ats) || occursin("Axes restored to original set",ats)
                        #println("end")
                        break
                    end
                    sp = split(ats)
                    if length(sp) > 1
                        #println("1111")
                        push!(RamanList,COORD(sp[1],sp[2]))
                        #println("2222")
                    end                                        
                end
            end
        end
    end
    counts1 = length(IRList)
    counts2 = length(RamanList)
    #println("c1: ",counts1,"c2: ",counts2)
    for k in 1:counts1
        open(IRNames,"a") do steams
            println(steams,IRList[k].coordX," ",IRList[k].coordY)
        end
    end
    for k in 1:counts1
        pop!(IRList )
    end
    for k in 1:counts2
        open(RaNames,"a") do steams
            println(steams,RamanList[k].coordX," ",RamanList[k].coordY)
        end
    end
    for k in 1:counts2
        pop!(RamanList)
    end
end
=#





#----提取排序文件前100个pdb结构小分子----

"""mutable struct SORTS
    idx::Int
  names::String
end

global pdb_all = "C:\\Users\\22103\\Desktop\\pdbs_withH"
global softFile = "C:\\Users\\22103\\Desktop\\pdbs_withH\\sorted_result.txt"

filelist = readdir(pdb_all)
pdb_100 = string(pdb_all,"\\pdb_100")
mkdir(pdb_100)

sortlist = []
open(softFile,"r") do molread
    global i = 0
    while !eof(molread)
        i = i + 1
        line = readline(molread)
        sp = split(line,",")
        spl = split(sp[2],"'")
        push!(sortlist,SORTS(i,spl[2]))
    end
end

for ifile in filelist
    pdbName = split(ifile,".")[1] 
    ifile1 = string(pdb_all,"\\",ifile)    
    for j in 1:100
        if occursin(pdbName,sortlist[j].names)
            newName = string(sortlist[j].idx,"-",ifile)
            newRoot = string(pdb_100,"\\",newName)
            cp(ifile1,newRoot)
            break
        end
    end            
end"""





#----pdb转gjf----

global pdb_400 = "/work1/scquant/ParaEngine-git-Drug_small_molecule/pdb_400"
filelist = readdir(pdb_400)
for ifile in filelist
    gjfName = split(ifile,".")[1]
    newName = string(gjfName,".gjf")
    run(`obabel $(pdb_400)/$(ifile) -ogjf -O $(pdb_400)/$(newName)`)
end
