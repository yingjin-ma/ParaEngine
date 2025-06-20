# e.g.
#       julia  Timing_analysis.jl  ../../ParaEngine-git-Gaussian/Test_P38-M062x-ALL-MLSLB-DLB

target=ARGS[1]
cpgjfdir=ARGS[2]
cplogdir=ARGS[3]

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

global fraglist=[]
#icount = 0
open("fcst_raw_wild00.dat","w") do record1
for ifile in filelist
    print("ifile : ",ifile)
    # print(record1,ifile)
    ifolder  = string(target,"/",ifile)
    if isdir(ifolder)
        dirlist = readdir(ifolder)
        #println("dirlist : ",dirlist)
        for id in dirlist
            sline=split(id,".")
            #println("id(00) : ", id)
            #print(id," ")
            Eiters=[]
            if uppercase(sline[2])=="LOG"
                println(" id(11) :  ",id)
                sline=split(split(id,".")[1],"-")
                #global icount = icount + 1
                energy = 0
                outidx = string(ifolder,"/",id)
                gjfidx = string(ifolder,"/",split(id,".")[1],".gjf")

                Month1 = ""
                Month2 = ""
                 iday1 = 0
                 iday2 = 0
                timeS1 = ""
                timeS2 = ""
                timeCPU = 0.0
                timeELP = 0.0
                 nbasis = 0
                  basis = ""
                 energy = 0.0
                 errgau = -1
                println("outidx : ",outidx)

                open(outidx,"r") do readlog
                    while !eof(readlog)
                        line=readline(readlog)
                        sline=split(line)
                        ntmp=length(sline)
                        if ntmp > 8
                            if sline[1] == "Leave" && sline[2] == "Link" && sline[3] == "1"
                                Month1 = sline[6]
                                 iday1 = parse(Int32,sline[7])
                                timeS1 = sline[8]
                            end
                            if sline[1] == "Normal" && sline[2] == "termination" && sline[4] == "Gaussian"
                                Month2 = sline[8]
                                 iday2 = parse(Int32,sline[9])
                                timeS2 = sline[10]
                            end
                            if sline[1] == "Job" && sline[2] == "cpu" && sline[3] == "time:"
                                 idv1 = parse(Int32,sline[4])
                                 idv2 = parse(Int32,sline[6])
                                 idv3 = parse(Int32,sline[8])
                                 idv4 = parse(Float64,sline[10])
                                 timeCPU = idv1*86400 + idv2*3600 + idv3*60 +idv4
                            end
                            if occursin("Delta-E=",line)
                                dv1 = parse(Float64,sline[4])
                                push!(Eiters,dv1)
                            end
                            if occursin("basis functions,",line)
                                nbasis = parse(Int32,sline[1])
                            end
                            if occursin("Standard basis:",line)
                                basis = sline[3]
                            end
                            if occursin("SCF Done:",line)
                                energy = parse(Float64,sline[5])
                            end
                            if occursin("Error termination",line)
                                errgau = 1
                            end
                        end
                    end
                end

                if errgau == -1

                    run(`cp $(gjfidx) $(cpgjfdir)`)
                    run(`cp $(outidx) $(cplogdir)`)

                    #println("Time1 : ",Month1," ",iday1," ",timeS1)
                    #println("Time2 : ",Month2," ",iday2," ",timeS2)

                    timeSS1=split(timeS1,":")
                    timeSS2=split(timeS2,":")

                    iidv1 = parse(Int32,timeSS1[1])
                    iidv2 = parse(Int32,timeSS1[2])
                    iidv3 = parse(Int32,timeSS1[3])

                    jjdv1 = parse(Int32,timeSS2[1])
                    jjdv2 = parse(Int32,timeSS2[2])
                    jjdv3 = parse(Int32,timeSS2[3])

                    if Month1 == Month2
                        timeELP = (iday2-iday1)*86400+(jjdv1-iidv1)*3600+(jjdv2-iidv2)*60+(jjdv3-iidv3)
                    else
                        println("Not implemented")
                    end

                    niter = length(Eiters) + 1

                    # println("timeCPU  ", timeCPU, "  ||  timeELP  ", timeELP)
                    println(record1, nbasis, " 0 ", timeCPU, " ", timeELP, " ", split(id,".")[1], " ", niter, " ", timeCPU/niter, " 999999999 ", basis)
                end

            end
        end
    end
    #println(record1)
end
end

