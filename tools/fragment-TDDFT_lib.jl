using Printf

mutable struct ATOMS
      idx::Int
    ifrag::Int
  icharge::Int
     elem::String
    coord::Tuple{Float64,Float64,Float64}
       ZZ::Float64
  residue::String
end

mutable struct FRAGS
      idx::Int
    iatom::Int
   natoms::Int
  icharge::Int
 multiple::Int
     name::String
   energy::Float64
end

function getmult(frag,atoms,ipdb)

    # ipdb :  1   PDB style
    # ipdb : -1   XYZ style

    elemdict=Dict( "H" => 1,
               "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8, "F" => 9,
               "P" => 15, "S" => 16, "Cl" => 17 )

    istart  = frag.iatom
    ifinish = frag.iatom + frag.natoms - 1

    #println(" istart, ifinish : ", istart, "  " , ifinish)


    nele=0
    for i in istart:ifinish
        #println(" ELEM : ",atoms[i].elem)
        if ipdb == 1
            elem = strip(atoms[i].elem[2:2])
        elseif ipdb == -1
            elem = strip(atoms[i].elem)
        end
        ielem = elemdict[elem]
        nele = nele + ielem
    end

    nele = nele - frag.icharge
    imod = mod(nele,2)

    frag.multiple = 2 * imod +1

    #println("nele, Multi", nele, frag.multiple)

end

function residueADDH()

    println("Add the H to the residues")

    i  = 0
    ii = 0
    iatom = 0
    names1 = []
    icharges1 = []
    for id in 1:length(ResID2)
        i = i + 1
        if ResID2[id] in ResID
                  ii = ii + 1
             istart1 = fraglist[i-1].iatom
            ifinish1 = fraglist[i-1].iatom+fraglist[i-1].natoms-1
             istartx = fraglist[i].iatom
            ifinishx = fraglist[i].iatom+fraglist[i].natoms-1
             istart3 = fraglist[i+1].iatom
            ifinish3 = fraglist[i+1].iatom+fraglist[i+1].natoms-1
              atoms1 = atomlist[istart1:ifinish1]
              atomsx = atomlist[istartx:ifinishx]
              atoms3 = atomlist[istart3:ifinish3]
             id1,id3,v1,v3 = minIDFRAG(atomsx,atoms1,atoms3)
            #println("addH1[",i-1,"] : ", id1, "   addH2[",i+1,"] : ", id3)

            name = fraglist[i].name
            push!(names1,name)
            push!(icharges1,fraglist[i].icharge)

            # push!(atomlist1,ATOMS(1+iatom,ii,0," H s",atoms1[id1].coord,0.0,atomsx[1].residue))
            push!(atomlist1,ATOMS(1+iatom,ii,0," H s",v1,0.0,atomsx[1].residue))
            natomsx = length(atomsx)
            println("natomsx", natomsx , " ResID2[id]",ResID2[id] )
            for j in 1:natomsx
                push!(atomlist1,ATOMS(1+j+iatom,ii,atomsx[j].icharge,atomsx[j].elem,atomsx[j].coord,0.0,atomsx[j].residue))
            end
            # push!(atomlist1,ATOMS(2+natomsx+iatom,ii,0," H f",atoms3[id3].coord,0.0,atomsx[1].residue))
            push!(atomlist1,ATOMS(2+natomsx+iatom,ii,0," H f",v3,0.0,atomsx[1].residue))
            iatom = iatom + 2 + natomsx

        end
    end
    global total_atoms1=iatom
    for i in 1:20
        if i > total_atoms1
            break
        else
            println(atomlist1[i])
        end
    end
    if total_atoms1 > 20
        println("... (more) for CAPPED residues ")
    end

    ifrag=0
    global fraglist1=[]
    icharge    = 0
    multiple   = 1
    natoms_frg = 0
    for i in 1:total_atoms1
        #println("atomlist[i][1] : ",atomlist1[i].ifrag)
        if atomlist1[i].ifrag != ifrag
            if ifrag != 0
                istart=atomlist1[i].idx-natoms_frg
                icharge=icharges1[ifrag]
                idx = parse(Int32,split(names1[ifrag],"-")[3])
                push!(fraglist1,FRAGS(idx,istart,natoms_frg,icharge,1,names1[ifrag],0.0))
            end
            ifrag = ifrag+1
            natoms_frg = 0
        end
        if atomlist1[i].ifrag == ifrag
            natoms_frg=natoms_frg+1
        end
        if atomlist1[i].icharge != 0
            #println("atomlist1[i][3] : ",atomlist1[i].icharge)
            icharge=icharge+atomlist1[i].icharge
        end
        if i == total_atoms1
            istart=atomlist1[i].idx-natoms_frg+1
            idx = parse(Int32,split(names1[ifrag],"-")[3])
            push!(fraglist1,FRAGS(idx,istart,natoms_frg,icharge,1,names1[ifrag],0.0))
        end
    end
    global total_frags1=ifrag

    println("total molecules in capped residue suit : ",total_frags1)
    for i in 1:total_frags1
        println(" fraglist1[",i,"]",fraglist1[i])
        getmult(fraglist1[i],atomlist1,1)
        #println(" ==== ")
        #println(fraglist1[i])
        #println(" ==== 2 ==== ")
    end
end


function readpdbfile(inpdb,flag)

    global atomlist1=[]
    global fraglist1=[]

    if uppercase(flag) == "ADDH"
        ResTMP=[]
        for id in 1:length(ResID)
            push!(ResTMP,ResID[id]-1)
            push!(ResTMP,ResID[id])
            push!(ResTMP,ResID[id]+1)
        end
        global ResID2 = unique(ResTMP)
    else
        ResID2 = ResID
    end

    println("ResID extended : ", ResID2)

    global modelpdb=[]
    global model_natoms=[]
    imodel=0 
    open(inpdb,"r") do pdbstream
        for line in eachline(pdbstream)
            sline=split(line)
            ntmp=length(sline)
            if ntmp > 0
                # One time one residual/part/fragment
                if uppercase(line[1:5])==flagMODEL1 || flagMODEL1 ==""
                    global natoms=0
                    global atomlist=[]
                end                  
                if ntmp > 10 
                        #println(line)
                        icharge= 0
                        natoms = natoms + 1
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
                            #push!(atomlist,ATOMS(natoms,ifrag,icharge,line[13:16],(dx,dy,dz),0.0,line[18:26]))
                            push!(atomlist,ATOMS(natoms,ifrag,icharge,line[77:80],(dx,dy,dz),0.0,line[18:26]))
                        else
                            if line[80:80] == "-"
                                icharge = -1 * icharge
                            end
                            #push!(atomlist,ATOMS(natoms,ifrag,icharge,line[13:16],(dx,dy,dz),0.0,line[18:26]))
                            push!(atomlist,ATOMS(natoms,ifrag,icharge,line[77:80],(dx,dy,dz),0.0,line[18:26]))
                        end
                end
                if uppercase(line[1:5])==flagMODEL2 || flagMODEL2==""
                    imodel = imodel + 1 
                    push!(modelpdb,atomlist) 
                    push!(model_natoms,natoms)
                end                  
            end
        end
    end
    #println(model)
    #println(model_natoms)
    #exit(1)

    global model_nfrags=[]
    global modelfrag = []
    for imodel in 1:length(modelpdb)
        ifrag=1
        total_atoms = model_natoms[imodel]
        for i in 1:total_atoms-1
            idiff = -1
            if modelpdb[imodel][i].ifrag != modelpdb[imodel][i+1].ifrag
                idiff = 1
            end
            modelpdb[imodel][i].ifrag = ifrag
            if idiff == 1
                ifrag = ifrag + 1
            end
        end
        modelpdb[imodel][total_atoms].ifrag = ifrag

        if imodel == length(modelpdb)
            println(" === Print the final models for checking === ")
            for i in 1:20
                if i > total_atoms
                    break
                else
                    println(modelpdb[imodel][i])
                end
            end
            if total_atoms > 20
                println("... (more) ")
            end
        end 

        ifrag=0
        fraglist=[]
        icharge    = 0
        multiple   = 1
        natoms_frg = 0
        for i in 1:total_atoms
            #println("model[imodel][i].ifrag : ",modelpdb[imodel][i].ifrag)
            if modelpdb[imodel][i].ifrag != ifrag
                if ifrag != 0
                    istart=modelpdb[imodel][i].idx-natoms_frg            
                    #println(modelpdb[imodel][i-1])        
                    idx = parse(Int32,split(modelpdb[imodel][i-1].residue)[2])
                    push!(fraglist,FRAGS(idx,istart,natoms_frg,icharge,1,replace(modelpdb[imodel][i-1].residue," "=>"-"),0.0))
                end
                ifrag = ifrag+1
                icharge    = 0
                natoms_frg = 0
            end
            if modelpdb[imodel][i].ifrag == ifrag
                natoms_frg=natoms_frg+1
            end
            if modelpdb[imodel][i].icharge != 0
                #println("atomlist[i][3] : ",modelpdb[imodel][i].icharge)
                icharge=icharge+modelpdb[imodel][i].icharge
            end
            if i == total_atoms
                istart=modelpdb[imodel][i].idx-natoms_frg+1
                idx = parse(Int32,split(modelpdb[imodel][i].residue)[2])
                push!(fraglist,FRAGS(idx,istart,natoms_frg,icharge,1,replace(modelpdb[imodel][i].residue," "=>"-"),0.0))
            end
        end
        global total_frags=ifrag

        if imodel == length(modelpdb)
            println("total_frags : ",total_frags," of the ",imodel, "-th model")
        end 
        for i in 1:total_frags
            getmult(fraglist[i],modelpdb[imodel],1)
            if imodel == length(modelpdb)
                # println(fraglist[i])
                # println(" ==== ")
            end 
        end
        push!(model_nfrags,total_frags)
        push!(modelfrag,fraglist)
    end   
end


function distanceMIN(atoms1,atoms2)

    mindist=10000.0
    for i in 1:length(atoms1)
        for j in 1:length(atoms2)
            dvx = atoms1[i].coord[1] - atoms2[j].coord[1]
            dvy = atoms1[i].coord[2] - atoms2[j].coord[2]
            dvz = atoms1[i].coord[3] - atoms2[j].coord[3]
            dv  = dvx^2 + dvy^2 + dvz^2
            dv  = sqrt(dv)
            if dv < mindist
                mindist = dv
            end
        end
    end
    return mindist
end


function combineFRAG(natoms,nij,frag1,atoms1,frag2,atoms2)

    name = string(frag1.name,"_",frag2.name)
    for i in 1:length(atoms1)
        push!(atomlist12,ATOMS(i+natoms,nij,atoms1[i].icharge,atoms1[i].elem[2:2],atoms1[i].coord,0.0,name))
    end
    for i in 1:length(atoms2)
        push!(atomlist12,ATOMS(i+natoms,nij,atoms2[i].icharge,atoms2[i].elem,atoms2[i].coord,0.0,name))
    end
    push!(fraglist12,FRAGS(nij,natoms+1,frag1.natoms+frag2.natoms,frag1.icharge+frag2.icharge,1,name,0.0))
end


function MFCC_TDDFT_PRE(fraglist,atomlist,exID,threshold1,threshold2,fraglist12,atomlist12)

    #println(" fraglist ", fraglist)
    #println(" atomlist ", atomlist)

    ijfrag = 0
    #global atomlist12=[]
    #global fraglist12=[]
    for i in 1:length(fraglist)-1
         istart = fraglist[i].iatom
        ifinish = fraglist[i].iatom+fraglist[i].natoms-1
        for j in i+1:length(fraglist)
             jstart = fraglist[j].iatom
            jfinish = fraglist[j].iatom+fraglist[j].natoms-1
            distance1 = distanceMIN(atomlist[istart:ifinish],atomlist[jstart:jfinish])
            if i in exID || j in exID
                threshold = threshold1
            else
                threshold = threshold2
            end 
            if distance1 < threshold
                println("distance[",i,"][",j,"] : ", distance1)              
                # println("fraglist[i].natoms, fraglist[j].natoms : ", fraglist[i].natoms , "  ",  fraglist[j].natoms)              
                ijfrag = ijfrag + 1
                combineFRAG(natoms,ijfrag,fraglist[i],atomlist[istart:ifinish],fraglist[j],atomlist[jstart:jfinish])
                natoms12 = natoms12 + fraglist[i].natoms + fraglist[j].natoms
            end
        end
    end
    #println("natoms12 : ", natoms12 )             
    return natoms12  

end 

function wffragPDB(ibase1,ibase2,frag,atoms,pdbfile,pdbinfo,stamp,ipdb)

    # ipdb :  1   PDB style
    # ipdb : -1   XYZ style

    println("ibase1 : ", ibase1)
    println("ibase2 : ", ibase2)
    println("frag   : ", frag)
    println("atoms  : ", atoms)

    open(pdbfile,"a+") do wfpdb
        for i in 1:length(frag)
            ibase1 = ibase1 + 1
            i1 = frag[i].iatom
            i2 = frag[i].iatom + frag[i].natoms -1
            for j in i1:i2
                ibase2 = ibase2 + 1
                println("atoms[",j,"].elem : ",atoms[j].elem)
                if ipdb == 1
                    elem = atoms[j].elem[2:2]
                else
                    elem = atoms[j].elem
                end
                idx1 = lpad("$ibase1",4)
                idx2 = lpad("$ibase2",7)
                # println("idx1 : ",idx1, "   idx2 : ", idx2)
                dx = atoms[j].coord[1]
                dy = atoms[j].coord[2]
                dz = atoms[j].coord[3]

                if atoms[j].icharge > 0
                    sign = "+"
                end
                if atoms[j].icharge < 0
                    sign = "-"
                end

                ielec = abs(atoms[j].icharge)
                if ielec == 0
                    elec = " "
                    sign = " "
                else
                    elec = string(ielec)
                end

                print(wfpdb,"ATOM",idx2,"  ",elem," ","  ",stamp,idx1,"    ")
                @printf(wfpdb,"%8.3f%8.3f%8.3f",dx,dy,dz)
                print(wfpdb,"      ","      ","           ",elem,elec,sign)

                println(wfpdb)
            end
        end
    end

    println("before pdbinfo in wffragPDB")

    open(pdbinfo,"a+") do infopdb
        for i in 1:length(frag)
            println(infopdb,stamp," frag[",i,"] : ",frag[i])
        end
    end

    return ibase1,ibase2

end


function MFCC_TDDFT_GEN(fmt,pdbfile,fraglist,atomlist,fraglist12,atomlist12)

    if uppercase(fmt)=="PDB"

        pdbinfo = replace(pdbfile,".pdb"=>".info")

        open(pdbfile,"w") do wfpdb
            println(wfpdb,"HEADER")
            println(wfpdb,"TITLE     Built with ParaEngine")
            println(wfpdb,"REMARK   ParaEngine generated pdb file")
            println(wfpdb,"REMARK   GITHUB")
            println(wfpdb,"REMARK   ")
        end

        # 
        # fraglist    -->    fraglist12
        # fraglist        : Monomers
        # fraglist12      : Dimers
        # 

        open(pdbinfo,"w") do infopdb
            println(infopdb," Fragments TDDFT with N-body interactions ")
        end

        ibase1 = 0
        ibase2 = 0
        ibase1,ibase2 = wffragPDB(ibase1,ibase2,fraglist,atomlist,pdbfile,pdbinfo,"MONOM",1)
        ibase1,ibase2 = wffragPDB(ibase1,ibase2,fraglist12,atomlist12,pdbfile,pdbinfo,"DIMER",1)

        open(pdbfile,"a+") do wfpdb
            println(wfpdb,"END")
        end
    end

end

