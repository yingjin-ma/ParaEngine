using Distributed
using BenchmarkTools

@everywhere function read_last(file)
    r=""
    open(file) do io
        seekend(io)
        seek(io, position(io) - 2)
        while Char(peek(io)) != '\n'
            seek(io, position(io) - 1)
        end
        Base.read(io, Char)
        r=Base.read(io, String)
        return r
    end
end

@everywhere function check_last_NWChem(out)
    while true
      sleep(5)
      tailout=read_last(out)
      println(tailout)
      if occursin("Total times  cpu:",tailout)
        break
      end
    end
end

function nwchemtask(task,frag,atoms,par)

    # Generating NWChem headers
   
    fragidx=lpad(frag.idx,8,"0")
    #run(`cp Template.nwchem $(fragidx).nw`) 
 
    # 
    open("Frag-$(fragidx).nw","w") do nwinp
        txt = read("Template.nwchem",String)
        txt = replace(txt,"FRAG-i" => "Frag-$(fragidx)","X-library-basis" => "* library 6-31g","dftxc" => "xc m06-2x")
        print(nwinp,txt)
    end     

    open("XYZ.tmp","w") do xyztmp
        for iatom in atoms
            println(xyztmp,iatom.elem," ",iatom.coord[1]," ",iatom.coord[2]," ",iatom.coord[3])
        end
    end 

    run(`sed -i "4r XYZ.tmp" Frag-$(fragidx).nw`)

    task.infile  = "Frag-$(fragidx).nw"
    task.outfile = "Frag-$(fragidx).out"

end 

function gentask()

    @everywhere global tasklist=[]

    println("")
    if isdir(workdir)
        println("WorkDir : ",workdir)
    else 
        mkdir(workdir)
        println("WorkDir : ",workdir)
    end

    if qcdriver == "NWCHEM"
        println("Generating the NWChem tasks ... ")
        println("")        
        for i in 1:total_frags
             istart = fraglist[i].iatom
            ifinish = fraglist[i].iatom+fraglist[i].natoms-1
            print(" frags-",i," : atomlist ",istart,"-",ifinish)           
            push!(tasklist,TASKS(i,fraglist[i].idx,workdir,"","","NWCHEM")) 
            nwchemtask(tasklist[i],fraglist[i],atomlist[istart:ifinish],"DEFAULT")
            println(" ... done ")
        end  
        println("")
    end
 
    println(tasklist)

end 

function distributingtask()

    nnodes = -1
    iflag  = -1
    icount =  0 

    println("")
    println("Distributing the QC tasks ... ")
    if LBfile != ""        
        print("Assigning tasks with --> ")
        print(LBfile)
        println(" <-- ")
        # Read in the loadbalance file 
        open(LBfile,"r") do LBread
            while !eof(LBread)
                line=readline(LBread)
                sline=split(line)
                ntmp=length(sline)
                if ntmp > 0                   
                    if uppercase(sline[1]) == "NODES"
                        nnodes = parse(Int32,sline[2])
                        LBmat=[] 
                        for i in 1:nnodes
                            line=readline(LBread)
                            sline=split(line)
                            ntmp=length(sline)
                            LBvec=[]
                            for i in 1:ntmp
                                push!(LBvec,parse(Int32,sline[i]))
                            end 
                            push!(LBmat,LBvec)
                        end  
                        push!(LBcube,LBmat)
                    end 
                end
            end
        end
    else
        println("Assigning tasks with --> Default <-- ")
    end  
    println("")

end


@everywhere function NWChemRUN(tvec,tlist,id)
 
#    println("taskvec ",tvec) 
#    println("tasklist",tlist) 

    for itask in tvec
        if itask != 0
            print("id",id,gethostname()," itask ",itask," ",(tlist[itask].infile)," ",(tlist[itask].outfile))
            RUNXX=`time ../../NWChemRUN $(tlist[itask].infile) $(tlist[itask].outfile) `
            tlist[itask].folder=string(tlist[itask].folder,"/",gethostname())
            println(" ",tlist[itask].folder) 
            if !isdir(tlist[itask].folder)
                mkpath(tlist[itask].folder)
            end 
            run(`mv  $(tlist[itask].infile)  $(tlist[itask].folder)`)
            run(Cmd(RUNXX,dir=tlist[itask].folder,detach=true))
            ccheck=string(tlist[itask].folder,"/",tlist[itask].outfile)
            rdlast=@spawnat id check_last_NWChem(ccheck)
            fetch(rdlast)
        end 
    end 

    #dir0="/work1/mayj/Test_CODES/Test_Julia/Frag-000000000"
    #for i in i1:i2
    #    AXX=`time ../NWChemRUN Frag-000000000$(i).nw out$(i).txt `
    #    cdir=string(dir0,string(i))
    #    run(Cmd(AXX,dir=cdir,detach=true))
    #    #b=@spawnat id run(Cmd(AXX,dir=cdir,detach=true))
    #    ccheck=string("./Frag-000000000",string(i),"/out",string(i),".txt")
    #    rdlast=@spawnat id check_last_NWChem(ccheck)
    #    fetch(rdlast)
    #end
    println("Done in NWChemRUN with inode-",id)
end



function runtask(NNSLURM, IFSLURM)

    println("")
    for i in 1:length(LBcube) 
        if i == NNSLURM
            println("Matched the number of NODES between SLURM and Task-Pre-Assignment ") 
            println("    NSLURM_NODES : ", NNSLURM, "        NPRE_NODES : ", i)  
            println("                  ------------           ")  
            for j in 1:length(LBcube[i])
                println("Tasks on node-",j,"  :  ", LBcube[i][j])
            end  
            println("                  ------------           ")  
        end  
    end
    println("Running the QC tasks ... ")
    println("    Number of   cores : ", nprocs())
    println("    Number of workers : ", nworkers())
    if IFSLURM 
        open("nodelist","r") do stream
            for line in eachline(stream)
                println(line)
            end
        end
    end 

    for i in workers()
        id, pid, host = fetch(@spawnat i (myid(), getpid(), gethostname()))
        println( i, " ", id, " ", pid, " ", host, " done " )
    end
 
    @everywhere println("    process: $(myid()) on host $(gethostname())")
    println("")
    println("")

    run=[]
    for i in 1:NNSLURM
        if IFSLURM
          inode = i+1
        else
          inode = i
        end  
        push!(run,@spawnat inode NWChemRUN(LBcube[NNSLURM][i],tasklist,inode))
    end 

    for i in 1:NNSLURM
        fetch(run[i])
    end  

end

 
