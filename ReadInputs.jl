function getmult(frag,atoms,ifcs=false)

    elemdict=Dict( "H" => 1,
               "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8, "F" => 9,
               "P" => 15, "S" => 16, "Cl" => 17 )

    istart  = frag.iatom 
    ifinish = frag.iatom + frag.natoms - 1

    #println(" istart, ifinish : ", istart, "  " , ifinish)

    nele = 0
    imod = 0
    for i in istart:ifinish
        #println(" ELEM : ",atoms[i].elem)
        ielem = elemdict[strip(atoms[i].elem)]
        nele = nele + ielem
    end 

    #println("00 : ifcs ",ifcs, " | nele, Multi, charge", nele, " ",frag.multiple, " ",frag.icharge)

    nele = nele - frag.icharge
    imod = mod(nele,2)
    if ifcs  # assign closed Shell
        if imod != 0
            frag.icharge = frag.icharge - 1
            nele = nele + 1      # "-" electron
            imod = mod(nele,2)
        end 
    end 

    frag.multiple = 2 * imod +1
    #println("11 : nele, Multi, charge", nele, " ",frag.multiple," ", frag.icharge)
 
    return frag.icharge

end

function readpdbfile(inpdb)

    # Tested, not used
    pdb=open(inpdb,"r")
    nlines=countlines(pdb)
    # println(nlines)    

    #
    global atomlist=[]  
    natoms=0 
    open(inpdb,"r") do stream
        for line in eachline(stream)
            sline=split(line)
            ntmp=length(sline)
            if ntmp > 1
                if uppercase(line[1:6]) == "HETATM" || uppercase(line[1:4]) == "ATOM"
                    icharge= 0
                    natoms = natoms+1    
                    #println(line[23:26],line[31:38],line[39:46],line[47:54]) 
                    #println(line)
                    ifrag = parse(Int32,line[23:26])
                       dx = parse(Float64,line[31:38])              
                       dy = parse(Float64,line[39:46])              
                       dz = parse(Float64,line[47:54])              
                    if line[79:79] != " "
                        # println("line : ",line[79:79])  
                        icharge = parse(Int32,line[79:79])                        
                    end
                    if length(line) > 80
                        icharge = parse(Int32,line[81:82])
                        push!(atomlist,ATOMS(natoms,ifrag,icharge,line[13:16],(dx,dy,dz),0.0))
                    else                             
                        if line[80:80] == "-"
                            icharge = -1 * icharge   
                        end
                        push!(atomlist,ATOMS(natoms,ifrag,icharge,line[13:16],(dx,dy,dz),0.0))
                    end   
                end 
            end 
        end 
    end 
    global total_atoms=natoms 

    ifrag=1
    for i in 1:total_atoms-1
        idiff = -1
        if atomlist[i].ifrag != atomlist[i+1].ifrag
            idiff = 1
        end         
        atomlist[i].ifrag = ifrag
        if idiff == 1
            ifrag = ifrag + 1
        end
    end 
    atomlist[total_atoms].ifrag = ifrag

    for i in 1:100
        if i > total_atoms  
            break
        else
            println(atomlist[i])
        end  
    end
    if total_atoms > 100
        println("... (more) ")
    end 

    ifrag=0
    global fraglist=[]
    icharge    = 0 
    multiple   = 1
    natoms_frg = 0
    for i in 1:total_atoms
        #println("atomlist[i][1] : ",atomlist[i].ifrag)
        if atomlist[i].ifrag != ifrag
            if ifrag != 0
                istart=atomlist[i].idx-natoms_frg
                push!(fraglist,FRAGS(ifrag,istart,natoms_frg,icharge,1,"unnamed",0.0))
            end 
            ifrag = ifrag+1
            icharge    = 0 
            natoms_frg = 0
        end 
        if atomlist[i].ifrag == ifrag
            natoms_frg=natoms_frg+1
        end 
        if atomlist[i].icharge != 0
            #println("atomlist[i][3] : ",atomlist[i].icharge)
            icharge=icharge+atomlist[i].icharge
        end
        if i == total_atoms
            istart=atomlist[i].idx-natoms_frg+1 
            push!(fraglist,FRAGS(ifrag,istart,natoms_frg,icharge,1,"unnamed",0.0))
        end 
    end
    global total_frags=ifrag

    println("total_frags : ",total_frags) 
    for i in 1:total_frags
        fraglist[i].icharge = getmult(fraglist[i],atomlist,closedshell)
        print(" ==  closedshell($(closedshell)) ")
        println(fraglist[i]) 
    end

    #exit(0)

end

function readxyz2gauss(suit)

    if !isdir(suit)
        println("The target suit $(suit) is not exist")
        exit("Stopped. Reason: $(suit) is not exist.")
    end

    println("suit : ", suit)
    filelist = readdir(suit)
    println(filelist)
    println("")

    global atomlist=[]
    natoms= 0
    ifrag = 0
    names = []
    for ifile in filelist
        if uppercase(split(ifile,".")[2]) == "XYZ" && uppercase(split(ifile,".")[1]) != "LIGAND0"
            push!(names,split(ifile,".")[1])
            ifile = string(suit,"/",ifile)
            ifrag = ifrag + 1
            println("ifile : ", ifile, " || XYZ  ==>  GJF " )
            i = 0
            open(ifile,"r") do xyzread
                readline(xyzread)
                cline = readline(xyzread)
                icharge = parse(Int32,cline)
                readline(xyzread)
                while !eof(xyzread)
                    cline  = readline(xyzread)
                    sline  = split(cline)
                    ntmp   = length(sline)
                    #println("ntmp : ", ntmp)
                    if ntmp > 3
                        i = i + 1
                        dx = parse(Float64,sline[2])
                        dy = parse(Float64,sline[3])
                        dz = parse(Float64,sline[4])
                        push!(atomlist,ATOMS(i+natoms,ifrag,icharge,sline[1],(dx,dy,dz),0.0))
                    end                    
                end 
            end 
            natoms = natoms + i
            println(natoms," of  natoms")
        end
    end 
    global total_atoms=natoms

    for i in 1:100
        if i > total_atoms
            break
        else
            println(atomlist[i])
        end
    end
    if total_atoms > 100
        println("... (more) ")
    end

    ifrag=0
    global fraglist=[]
    icharge    = 0
    multiple   = 1
    natoms_frg = 0
    for i in 1:total_atoms
        #println("atomlist[i][1] : ",atomlist[i].ifrag)
        if atomlist[i].ifrag != ifrag
            if ifrag != 0
                icharge = icharge/natoms_frg
                istart=atomlist[i].idx-natoms_frg
                push!(fraglist,FRAGS(ifrag,istart,natoms_frg,icharge,1,names[ifrag],0.0))
            end
            ifrag = ifrag+1
            icharge    = 0
            natoms_frg = 0
        end
        if atomlist[i].ifrag == ifrag
            natoms_frg=natoms_frg+1
        end
        if atomlist[i].icharge != 0
            #println("atomlist[i][3] : ",atomlist[i].icharge)
            icharge=icharge+atomlist[i].icharge
        end
        if i == total_atoms
            icharge = icharge/natoms_frg
            istart=atomlist[i].idx-natoms_frg+1
            push!(fraglist,FRAGS(ifrag,istart,natoms_frg,icharge,1,names[ifrag],0.0))
        end
    end
    global total_frags=ifrag

    println("total molecules in suit : ",total_frags)
    for i in 1:total_frags
        if closedshell        
            fraglist[i].icharge = getmult(fraglist[i],atomlist,closedshell)
            print(" ==  closedshell($(closedshell))  ")
        end  
        println(fraglist[i])
    end

#    exit(0)
end

function readgauss(suit)

    if !isdir(suit)
        println("The target suit $(suit) is not exist")
        exit("Stopped. Reason: $(suit) is not exist.")
    end

    println("suit : ", suit)
    filelist = readdir(suit)
    println(filelist)
    println("")

    global atomlist=[]
    natoms= 0
    ifrag = 0
    names = []
    for ifile in filelist
        push!(names,split(ifile,".")[1])
        ifile = string(suit,"/",ifile)
        ifrag = ifrag + 1
        println("ifile : ", ifile, " || change proc to 24 cores in ERA"  )

        txt = read(ifile,String)
        txt = replace(txt,"%mem=2GB" => "%mem=108GB", "%nproc=20" => "%nproc=24", "%nproc=10" => "%nproc=24", "MP2/" => "MP2=FullDirect " )
        open(ifile,"w") do nwinp
            print(nwinp,txt)
        end
 
        i = 0 
        open(ifile,"r") do molread
            # skip the 7 lines that generated by gridmol
            readline(molread)
            readline(molread)
            readline(molread)
            readline(molread)
            readline(molread)
            readline(molread)
            while !eof(molread)
                readline(molread)
                cline  = readline(molread)
                sline  = split(cline)
                ntmp   = length(sline)
                #println("ntmp : ", ntmp)
                icharge=0
                if ntmp > 3
                    i = i + 1
                    dx = parse(Float64,sline[2])
                    dy = parse(Float64,sline[3])
                    dz = parse(Float64,sline[4])
                    push!(atomlist,ATOMS(i+natoms,ifrag,icharge,sline[1],(dx,dy,dz),0.0))
                end  
            end
            natoms = natoms + i 
        end 
    end 
    global total_atoms=natoms 

    for i in 1:100
        if i > total_atoms
            break
        else
            println(atomlist[i])
        end
    end
    if total_atoms > 100
        println("... (more) ")
    end

    ifrag=0
    global fraglist=[]
    icharge    = 0
    multiple   = 1
    natoms_frg = 0
    for i in 1:total_atoms
        #println("atomlist[i][1] : ",atomlist[i].ifrag)
        if atomlist[i].ifrag != ifrag
            if ifrag != 0
                istart=atomlist[i].idx-natoms_frg
                push!(fraglist,FRAGS(ifrag,istart,natoms_frg,icharge,1,names[ifrag],0.0))
            end
            ifrag = ifrag+1
            icharge    = 0
            natoms_frg = 0
        end
        if atomlist[i].ifrag == ifrag
            natoms_frg=natoms_frg+1
        end
        if atomlist[i].icharge != 0
            #println("atomlist[i][3] : ",atomlist[i].icharge)
            icharge=icharge+atomlist[i].icharge
        end
        if i == total_atoms
            istart=atomlist[i].idx-natoms_frg+1
            push!(fraglist,FRAGS(ifrag,istart,natoms_frg,icharge,1,names[ifrag],0.0))
        end
    end
    global total_frags=ifrag

    println("total molecules in suit : ",total_frags)
    for i in 1:total_frags
        fraglist[i].icharge = getmult(fraglist[i],atomlist,closedshell)
        print(" ==  closedshell($(closedshell))  ")
        println(fraglist[i])
    end

end



function readsuits(suit)
 
    if !isdir(suit)
        println("The target suit $(suit) is not exist") 
        exit("Stopped. Reason: $(suit) is not exist.")
    end

    println("suit : ", suit)
    filelist = readdir(suit)
    println(filelist)
    println("")

    global atomlist=[]
    natoms= 0
    ifrag = 0 
    names = []
    for ifile in filelist
        ifile = string(suit,"/",ifile) 
        ifrag = ifrag + 1
        println("ifile : ", ifile)
        open(ifile,"r") do molread
            name  = readline(molread)
            #if name == ""
            sline  = split(ifile,"/")
            ntmp   = length(sline)
            ssline = split(sline[ntmp],".")
            name   = ssline[1] 
            #end
            push!(names,name)
            readline(molread)          
            readline(molread)          
            line  = readline(molread)          
            #sline = split(line)
            #println("line :",line)
            natom = parse(Int32,strip(line[1:3])) 
            nconn = parse(Int32,strip(line[4:6]))
            for i in 1:natom
                line=readline(molread)
                sline=split(line)
                icharge=0 
                dx = parse(Float64,sline[1])
                dy = parse(Float64,sline[2])
                dz = parse(Float64,sline[3])
                push!(atomlist,ATOMS(i+natoms,ifrag,icharge,sline[4],(dx,dy,dz),0.0))
            end
            for i in 1:nconn
                readline(molread)
            end 
            line  = readline(molread)
            sline = split(line)
            if sline[2] == "CHG"
                nelec = parse(Int32,sline[3])
                for i in 0:nelec-1
                    iatom   = parse(Int32,sline[2*i+4])
                    icharge = parse(Int32,sline[2*i+5])
                    atomlist[iatom+natoms].icharge=icharge
                end 
            end
            natoms=natoms+natom
        end  
    end 
    global total_atoms=natoms
    
    for i in 1:100
        if i > total_atoms
            break
        else
            println(atomlist[i])
        end
    end
    if total_atoms > 100
        println("... (more) ")
    end

    ifrag=0
    global fraglist=[]
    icharge    = 0
    multiple   = 1
    natoms_frg = 0
    for i in 1:total_atoms
        #println("atomlist[i][1] : ",atomlist[i].ifrag)
        if atomlist[i].ifrag != ifrag
            if ifrag != 0
                istart=atomlist[i].idx-natoms_frg
                push!(fraglist,FRAGS(ifrag,istart,natoms_frg,icharge,1,names[ifrag],0.0))
            end
            ifrag = ifrag+1
            icharge    = 0
            natoms_frg = 0
        end
        if atomlist[i].ifrag == ifrag
            natoms_frg=natoms_frg+1
        end
        if atomlist[i].icharge != 0
            #println("atomlist[i][3] : ",atomlist[i].icharge)
            icharge=icharge+atomlist[i].icharge
        end
        if i == total_atoms
            istart=atomlist[i].idx-natoms_frg+1
            push!(fraglist,FRAGS(ifrag,istart,natoms_frg,icharge,1,names[ifrag],0.0))
        end
    end
    global total_frags=ifrag

    println("total molecules in suit : ",total_frags)
    for i in 1:total_frags
        fraglist[i].icharge = getmult(fraglist[i],atomlist,closedshell)
        print(" ==  closedshell($(closedshell))  ")
        println(fraglist[i])
    end

    #exit()

end 

function readinp(infile)

    global iffifo = true
    global closedshell = false

    open(infile,"r") do stream
        while !eof(stream) 
            line=readline(stream)
            sline=split(line)
            ntmp=length(sline) 
            if ntmp > 0 

#                println(" sline : ", sline, "")

                if uppercase(sline[1]) == "TASK"
#                   println(uppercase(sline[1])) 
                    if uppercase(sline[2]) == "DFT"
                        global runtype = "DFT"
                    end  
                    if uppercase(sline[2]) == "GAUSSIAN"
                        global runtype = "GAU"
                        if uppercase(sline[3]) == "XYZ"
                            global runtype2 = "XYZGAU"
                        end  
                    end  
                    if uppercase(sline[3]) == "FRAG"
                        global runtype2 = "FRAG"
                    end  
                    if uppercase(sline[3]) == "GJF" || uppercase(sline[3]) == "COM"
                        global runtype2 = "G09"
                    end  
                    #println(" sline[1] [2] [3] : ", sline[1], sline[2], sline[3])
                    #println(" uppercase[2] [3] : ", uppercase(sline[2]), uppercase(sline[3]))
                end

                #println(" runtype ( 0 ) : ", runtype, " runtype2 ( 0 ) : ", runtype2 )

                # QC package as the engine
                if uppercase(sline[1]) == "ENGINE"
                    if uppercase(sline[2]) == "NWCHEM"
                        global qcdriver = "NWCHEM"
                    elseif uppercase(sline[2]) == "GAMESS"
                        global qcdriver = "GAMESS"
                    elseif uppercase(sline[2]) == "G09"
                        global qcdriver = "GAUSSIAN"
                    else
                        println("EXIT()")
                        exit()
                    end 
                    if length(sline) > 2 
                        if uppercase(sline[3]) == "CLOSEDSHELL"  # force closed-shell
                            global closedshell = true 
                        end 
                    end 
                end
                if uppercase(sline[1]) == "COORD"
                    global pdbfile = sline[2]
                    readpdbfile(pdbfile)
                end
                if uppercase(sline[1]) == "LOADBALANCE"
                    global LBfile  = sline[2] 
                end
                if uppercase(sline[1]) == "OVERLOAD"
                    global overload = parse(Float64,sline[2])
                    if length(sline) >= 3 
                        global nrepeat = parse(Int32,sline[3])
                    else 
                        global nrepeat = -1
                    end 
                    if length(sline) >= 4 
                        fifo = uppercase(sline[4])
                        if fifo == "FIFO"
                           global iffifo = true
                        else
                           global iffifo = false
                        end 
                    end
                end
                if uppercase(sline[1]) == "WORKDIR"
                    global workdir = sline[2] 
                end
                if uppercase(sline[1]) == "SUIT"
                    # println(" runtype2 : ",runtype2)
                    if runtype2 == "G09" || runtype2 == "G16"
                        global targetsuit = sline[2] 
                        readgauss(targetsuit)
                    elseif runtype2 == "XYZGAU"
                        global targetsuit = sline[2] 
                        readxyz2gauss(targetsuit)                      
                    else
                        global targetsuit = sline[2] 
                        readsuits(targetsuit)
                    end 
                end

                if sline[1] == "!>"
                    if uppercase(sline[2]) == "NWCHEM"
                        ctmp=""
                        while ctmp != "!>"
                            line=readline(stream)
                            sline=split(line)
                            ntmp=length(sline) 
                            if ntmp > 0
                                push!(QCpara,line)
                                ctmp=sline[1]
                            else
                                push!(QCpara," ")
                            end  
                        end
                        QCpara[length(QCpara)]=" "
                    end 
                    if uppercase(sline[2]) == "GAUSSIAN"
                        ctmp=""
                        while ctmp != "!>"
                            line=readline(stream)
                            sline=split(line)
                            ntmp=length(sline)
                            if ntmp > 0
                                push!(QCpara,line)
                                ctmp=sline[1]
                            end
                        end
                    end
                end
            end 
        end
    end


    println("runtype,runtype2,qcdriver,pdbfile,LBfile")
    println(runtype,runtype2,qcdriver,pdbfile,LBfile)
    #return runtype,runtype2,qcdriver,pdbfile,LBfile

    #println(QCpara)
    println("closedshell  ",closedshell)
    #exit() 

end 



function readPDB(flag)

    println(" ")
    println("-----------------------------------------------------------")
    if uppercase(flag) == "GRIDMOL"
        println("PDBs from GridMol package ") 
    else
        println("PDBs (standard) ")
    end 
    println("-----------------------------------------------------------")
    println(" ")

end 


