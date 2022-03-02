function readpdbfile(inpdb)

    # Tested, not used
    pdb=open(inpdb,"r")
    nlines=countlines(pdb)
    println(nlines)    

    #
    atomlist=[]  
    natoms=0 
    open(inpdb,"r") do stream
        for line in eachline(stream)
            sline=split(line)
            ntmp=length(sline)
            if ntmp > 0
                if uppercase(sline[1]) == "HETATM" || uppercase(sline[1]) == "ATOM"
                   natoms = natoms+1    
                   # println(line[23:26],line[31:38],line[39:46],line[47:54]) 
                    ifrag = parse(Int32,line[23:26])
                  icharge = 0 
                       dx = parse(Float64,line[31:38])              
                       dy = parse(Float64,line[39:46])              
                       dz = parse(Float64,line[47:54])              
                    if length(line) > 80
                        icharge = parse(Int32,line[81:82])
                        push!(atomlist,ATOMS(natoms,ifrag,icharge,line[13:16],(dx,dy,dz),0.0))
                    else                              
                        push!(atomlist,ATOMS(natoms,ifrag,icharge,line[13:16],(dx,dy,dz),0.0))
                    end   
                end 
            end 
        end 
    end 
    total_atoms=natoms  

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


end 

function readinp(infile)
    open(infile,"r") do stream
        for line in eachline(stream)
            sline=split(line)
            ntmp=length(sline) 
            if ntmp > 0 
                if uppercase(sline[1]) == "TASK"
                    if uppercase(sline[2]) == "DFT"
                        runtype = "DFT"
                    end  
                    if uppercase(sline[3]) == "FRAG"
                        runtype2 = "FRAG"
                    end  
                end
                if uppercase(sline[1]) == "ENGINE"
                    if uppercase(sline[2]) == "NWCHEM"
                        qcdriver = "NWCHEM"
                    end  
                    if uppercase(sline[3]) == "NWCHEM"
                        qcdriver = "NWCHEM"
                    end  
                end
                if uppercase(sline[1]) == "COORD"
                    pdbfile = sline[2]
                    readpdbfile(pdbfile)
                end
                if uppercase(sline[1]) == "LOADBALANCE"
                    LBfile  = sline[2] 
                end
            end 
        end
    end
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


