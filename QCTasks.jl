function nwchemtask(frag,atoms,par)

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

end 

function gentask()

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
            nwchemtask(fraglist[i],atomlist[istart:ifinish],"DEFAULT")
            println(" ... done ")
        end  
        println("")
    end

end 

function distributingtask(NNSLURM)

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
            for line in eachline(LBread)
                sline=split(line)
                ntmp=length(sline)
                if iflag == 1
                    icount = icount +1 
                    LBvec=[]
                    for i in 1:ntmp 
                        push!(LBvec,parse(Int32,sline[i]))
                    end   
                    push!(LBmat,LBvec)
                end
                if icount == nnodes
                    break
                end
                if ntmp > 0                    
                    if uppercase(sline[1]) == "NODES"
                        nnodes = parse(Int32,sline[2])
                        if nnodes == NNSLURM
                            println("Matched the number of NODES between SLURM and Task-Pre-Assignment ") 
                            println("NSLURM_NODES : ", NNSLURM, " NPRE_NODES : ", nnodes)  
                            iflag=1
                        end
                    end 
                end

            end
        end
        println(LBmat)
    else
        println("Assigning tasks with --> Default <-- ")
    end  
    println("")

end


function runtask()

    println("")
    println("Running the QC tasks ... ")
    println("")
    println("")

end 
