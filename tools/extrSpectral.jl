# Yingjin, pickup and rewritten from Ling Li's someProgram.jl code

mutable struct COORD
    coordX::String 
    coordY::String
end

global freqout  = "/work1/scquant/ParaEngine-git-CASinfo/data-CASinfo/top100drugs-withH-freq-results"
global freqout2 = string(freqout,"/data")

if !isdir(freqout2)
    mkdir(freqout2)
end 

filelist = readdir(freqout)

IR    = string(freqout2,"/IR")
Raman = string(freqout2,"/Raman")
if !isdir(IR)
    mkdir(IR)
end 
if !isdir(Raman)
    mkdir(Raman)
end

global IRList = []
global RamanList = []
for ifile in filelist
    taskdir = string(freqout,"/",ifile)
    println("taskdir : ", taskdir)
    filelist1 = readdir(taskdir)
    # println("filelist1 : ", filelist1)
    #exit(1) 
    for ifile1 in filelist1
        taskdir1 = string(taskdir,"/",ifile1)
        #println("taskdir1 :",taskdir1)
        if isdir(taskdir1)
            filelist2 = readdir(taskdir1)
            #println("filelist2 :",filelist2)
            for ifile2 in filelist2 
                if uppercase(split(ifile2,".")[2]) == "GJF"
                    global name = split(ifile2,".")[1]       
                    global IRNames =string(IR,"/",ifile,"/IR_",name,".txt")   
                    global RaNames =string(Raman,"/",ifile,"/Raman_",name,".txt")                 
                    if !isdir(string(IR,"/",ifile))
                        mkdir(string(IR,"/",ifile))
                    end
                    if !isdir(string(Raman,"/",ifile))
                        mkdir(string(Raman,"/",ifile))
                    end 
                    open(IRNames,"w") do blank
                        print(blank) 
                    end  
                    open(RaNames,"w") do blank
                        print(blank) 
                    end  
                    ifiles = string(taskdir1,"/",name,".log")
                    println("ifiles : ", ifiles)
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
            end 
        end 
    end
end


