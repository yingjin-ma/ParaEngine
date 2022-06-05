# using Distributed
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

function minIDFRAG(atoms1,atoms2,atoms3)

    mindist=10000.0
    pos1   =1

    minID2 = -1
    minID3 = -1
    stdV2  = (0.0,0.0,0.0)
    stdV3  = (0.0,0.0,0.0)

    for i in 1:length(atoms1)
        for j in 1:length(atoms2)
            dvx = atoms1[i].coord[1] - atoms2[j].coord[1]
            dvy = atoms1[i].coord[2] - atoms2[j].coord[2]
            dvz = atoms1[i].coord[3] - atoms2[j].coord[3]
            dv  = dvx^2 + dvy^2 + dvz^2
            dv  = sqrt(dv)
            dvx = atoms1[i].coord[1]-dvx/dv
            dvy = atoms1[i].coord[2]-dvy/dv
            dvz = atoms1[i].coord[3]-dvz/dv
            if dv < mindist
                mindist = dv 
                pos1    = i
                minID2  = j
                stdV2   = (dvx,dvy,dvz)
            end
        end
    end

    mindist=10000.0
    for i in 1:length(atoms1)
        if i != pos1
            for j in 1:length(atoms3)
                dvx = atoms1[i].coord[1] - atoms3[j].coord[1]
                dvy = atoms1[i].coord[2] - atoms3[j].coord[2]
                dvz = atoms1[i].coord[3] - atoms3[j].coord[3]
                dv  = dvx^2 + dvy^2 + dvz^2
                dv  = sqrt(dv)
                dvx = atoms1[i].coord[1]-dvx/dv
                dvy = atoms1[i].coord[2]-dvy/dv
                dvz = atoms1[i].coord[3]-dvz/dv
                if dv < mindist
                    mindist = dv
                    minID3  = j
                    stdV3   = (dvx,dvy,dvz) 
                end
            end
        end
    end

    return minID2,minID3,stdV2,stdV3

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

    global atomlist=[]
    natoms=0
    open(inpdb,"r") do pdbstream
        for line in eachline(pdbstream)
            sline=split(line)
            ntmp=length(sline)
            if ntmp > 0
                if uppercase(line[22:22])==flagSpike

                    # One time one residual
                    idres = parse(Int32,line[23:26])
                    for id in 1:length(ResID2)
                        if ResID2[id] == idres
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

    for i in 1:20
        if i > total_atoms
            break
        else
            println(atomlist[i])
        end
    end
    if total_atoms > 20
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
                idx = parse(Int32,split(atomlist[i-1].residue)[3])
                push!(fraglist,FRAGS(idx,istart,natoms_frg,icharge,1,replace(atomlist[i-1].residue," "=>"-"),0.0))
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
            idx = parse(Int32,split(atomlist[i].residue)[3])
            push!(fraglist,FRAGS(idx,istart,natoms_frg,icharge,1,replace(atomlist[i].residue," "=>"-"),0.0))
        end
    end
    global total_frags=ifrag

    println("total_frags : ",total_frags)
    for i in 1:total_frags
        getmult(fraglist[i],atomlist,1)
        #println(" ==== ")
        println(fraglist[i])
    end

end

function readXYZfile(suit,flag)

    filelist = readdir(suit)
    #println(filelist)
    println("")

    global atomlist2=[]
    natoms= 0
    ifrag = 0
    names2 = []
    icharges2 = []
    for ifile in filelist

        sline = split(ifile,".")
        ntmp  = length(sline)
        ifxyz = sline[ntmp]

        if uppercase(ifxyz) == "XYZ" && uppercase(sline[1][1:6])==string("CHAIN",flag) 

            name  = sline[1]
            push!(names2,name)
            ifile = string(suit,"/",ifile)
            ifrag = ifrag + 1
            # println("ifile : ", ifile)      
            nlines = countlines(ifile)
            natom = 0

            open(ifile,"r") do molread
                readline(molread)
                icharge=parse(Int32,readline(molread))
                push!(icharges2,icharge) 
                while !eof(molread)
                    line=readline(molread)
                    sline=split(line)
                    ntmp=length(sline)
                    if ntmp > 3
                        natom = natom + 1
                        dx = parse(Float64,sline[2])
                        dy = parse(Float64,sline[3])
                        dz = parse(Float64,sline[4])
                        push!(atomlist2,ATOMS(natom+natoms,ifrag,0,sline[1],(dx,dy,dz),0.0,name))
                    end 
                end  
                #atomlist2[1+natoms].icharge=icharge  # keep the charge in the first atom
            end
            natoms=natoms+natom
        end
    end
    global total_atoms2=natoms

    global atomlist3=[]
    natoms= 0
    ifrag = 0
    names3 = []
    icharges3 = []
    for ifile in filelist

        sline = split(ifile,".")
        ntmp  = length(sline)
        ifxyz = sline[ntmp]

        if uppercase(ifxyz) == "XYZ" && uppercase(sline[1][1:4])==string("CAP",flag)

            name  = sline[1]
            push!(names3,name)
            ifile = string(suit,"/",ifile)
            ifrag = ifrag + 1
            # println("ifile : ", ifile)
            nlines = countlines(ifile)
            natom = 0

            open(ifile,"r") do molread
                readline(molread)
                icharge=parse(Int32,readline(molread))
                push!(icharges3,icharge)
                while !eof(molread)
                    line=readline(molread)
                    sline=split(line)
                    ntmp=length(sline)
                    if ntmp > 3
                        natom = natom + 1
                        dx = parse(Float64,sline[2])
                        dy = parse(Float64,sline[3])
                        dz = parse(Float64,sline[4])
                        push!(atomlist3,ATOMS(natom+natoms,ifrag,0,sline[1],(dx,dy,dz),0.0,name))
                    end 
                end
                #atomlist2[1+natoms].icharge=icharge  # keep the charge in the first atom
            end
            natoms=natoms+natom
        end
    end
    global total_atoms3=natoms


    for i in 1:20
        if i > total_atoms2
            break
        else
            println(atomlist2[i])
        end
    end
    if total_atoms2 > 20
        println("... (more) for Chain ")
    end

    ifrag=0
    global fraglist2=[]
    icharge    = 0
    multiple   = 1
    natoms_frg = 0
    for i in 1:total_atoms2
        #println("atomlist[i][1] : ",atomlist2[i].ifrag)
        if atomlist2[i].ifrag != ifrag
            if ifrag != 0
                istart=atomlist2[i].idx-natoms_frg
                icharge=icharges2[ifrag]
                idx = parse(Int32,names2[ifrag][7:10])
                push!(fraglist2,FRAGS(idx,istart,natoms_frg,icharge,1,names2[ifrag],0.0))
            end
            ifrag = ifrag+1
            natoms_frg = 0
        end
        if atomlist2[i].ifrag == ifrag
            natoms_frg=natoms_frg+1
        end
        if atomlist2[i].icharge != 0
            #println("atomlist2[i][3] : ",atomlist2[i].icharge)
            icharge=icharge+atomlist2[i].icharge
        end
        if i == total_atoms2
            istart=atomlist2[i].idx-natoms_frg+1
            idx = parse(Int32,names2[ifrag][7:10])
            push!(fraglist2,FRAGS(idx,istart,natoms_frg,icharge,1,names2[ifrag],0.0))
        end
    end
    global total_frags2=ifrag

    println("total molecules in Chain suit : ",total_frags2)
    for i in 1:total_frags2
        #println(" fraglist2[",i,"]",fraglist2[i])
        getmult(fraglist2[i],atomlist2,-1)
        #println(" ==== ")
        #println(fraglist2[i])
        #println(" ==== 2 ==== ")
    end   



    for i in 1:20
        if i > total_atoms3
            break
        else
            println(atomlist3[i])
        end
    end
    if total_atoms3 > 20
        println("... (more) for CAP ")
    end

    ifrag=0
    global fraglist3=[]
    icharge    = 0
    multiple   = 1
    natoms_frg = 0
    for i in 1:total_atoms3
        #println("atomlist[i][1] : ",atomlist3[i].ifrag)
        if atomlist3[i].ifrag != ifrag
            if ifrag != 0
                istart=atomlist3[i].idx-natoms_frg
                icharge=icharges3[ifrag]
                idx = parse(Int32,names3[ifrag][5:8])
                push!(fraglist3,FRAGS(idx,istart,natoms_frg,icharge,1,names3[ifrag],0.0))
            end
            ifrag = ifrag+1
            natoms_frg = 0
        end
        if atomlist3[i].ifrag == ifrag
            natoms_frg=natoms_frg+1
        end
        if atomlist3[i].icharge != 0
            #println("atomlist3[i][3] : ",atomlist3[i].icharge)
            icharge=icharge+atomlist3[i].icharge
        end
        if i == total_atoms3
            istart=atomlist3[i].idx-natoms_frg+1
            idx = parse(Int32,names3[ifrag][5:8])
            push!(fraglist3,FRAGS(idx,istart,natoms_frg,icharge,1,names3[ifrag],0.0))
        end
    end
    global total_frags3=ifrag

    println("total molecules in CAP suit : ",total_frags3)
    for i in 1:total_frags3
        #println(" fraglist3[",i,"]",fraglist3[i])
        getmult(fraglist3[i],atomlist3,-1)
        #println(" ==== ")
        #println(fraglist3[i])
        #println(" ==== 2 ==== ")
    end



end 


function distanceFRAG(atoms1,atoms2)
# @everywhere function distanceFRAG(atoms1,atoms2)

    maxdist=10000.0
    for i in 1:length(atoms1)
        for j in 1:length(atoms2)
            dvx = atoms1[i].coord[1] - atoms2[j].coord[1] 
            dvy = atoms1[i].coord[2] - atoms2[j].coord[2] 
            dvz = atoms1[i].coord[3] - atoms2[j].coord[3] 
            dv  = dvx^2 + dvy^2 + dvz^2 
            dv  = sqrt(dv)
            if dv < maxdist
                maxdist = dv
            end 
        end 
    end
    return maxdist 

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

function generateMFCC_PRE(threshold)

    #println(atomlist)

    #println(fraglist)
    #println(fraglist2)

    global atomlist12=[] 
    global fraglist12=[]

    ijfrag = 0 
    natoms = 0
    for i in 1:length(fraglist)
         istart = fraglist[i].iatom
        ifinish = fraglist[i].iatom+fraglist[i].natoms-1        
        for j in 1:length(fraglist2)
             jstart = fraglist2[j].iatom
            jfinish = fraglist2[j].iatom+fraglist2[j].natoms-1                     
            distance1 = distanceFRAG(atomlist[istart:ifinish],atomlist2[jstart:jfinish])
            # println("distance[",i,"][",j,"] :", distance1 )
            if distance1 < threshold
                ijfrag = ijfrag + 1
                combineFRAG(natoms,ijfrag,fraglist[i],atomlist[istart:ifinish],fraglist2[j],atomlist2[jstart:jfinish]) 
                natoms = natoms + fraglist[i].natoms + fraglist2[j].natoms
            end 
        end
    end 
    global total_atoms12 = natoms
    global total_frags12 = ijfrag

    println("total_frags of Res-Chain : ",total_frags12)
    for i in 1:total_frags12
        getmult(fraglist12[i],atomlist12,-1)
        println("fraglist12[",i,"] : ",fraglist12[i])
        #println(" ==== ")
        #println(fraglist[i])
    end

#    for i in 1:length(atomlist12)
#        println("atomlist12[",i,"] : ",atomlist12[i])
#    end

    global res3in12=[]
    global res2in12=[]

    for i in 1:length(fraglist12)
        # println("fraglist12[",i,"] : ",fraglist12[i])
        sline = split(fraglist12[i].name,"_") 
        idx = parse(Int32,sline[2][7:10])
        push!(res2in12,idx)
    end 

    res2in12 = sort(unique(res2in12))
    if last(res2in12) > total_frags3 
        res3in12 = res2in12[1:length(res2in12)-1]
    else 
        res3in12 = res2in12
    end 
    println("res2in12(Chain ) : ",res2in12)
    println("res3in12(Capped) : ",res3in12)

    # For matching 
    IDfrags1 = Array{String}(undef, length(fraglist1))
    IDfrags2 = Array{String}(undef, length(fraglist2))
    IDfrags3 = Array{String}(undef, length(fraglist3))
    for i in 1:length(fraglist1)
        IDfrags1[i] = split(fraglist1[i].name,"-")[3]
    end 
    for i in 1:length(fraglist2)
        IDfrags2[i] = fraglist2[i].name[7:10]
    end 
    for i in 1:length(fraglist3)
        IDfrags3[i] = fraglist3[i].name[5:8]
    end 

#    println("IDfrags1 : ",IDfrags1)
#    println("IDfrags3 : ",IDfrags3)

    global atomlist13=[]
    global fraglist13=[]

    ifrag = 0 
    natoms = 0 
    for i in 1:length(fraglist12)
        #println("fraglist12[",i,"] : ",fraglist12[i])
        sline = split(fraglist12[i].name,"_") 
        idx1 = split(sline[1],"-")[3]
        idx3 = sline[2][7:10]

        idx1pos = findfirst(==(idx1),IDfrags1)
        idx3pos = findfirst(==(idx3),IDfrags3)

        # println(idx1,"'s pos : ",idx1pos,"     ", idx3,"'s pos : ",idx3pos)
     

        ifrag = ifrag + 1
        natoms13  = fraglist1[idx1pos].natoms+fraglist3[idx3pos].natoms
        icharge13 = fraglist1[idx1pos].icharge+fraglist3[idx3pos].icharge
        name13    = string(fraglist1[idx1pos].name,"_",fraglist3[idx3pos].name)

        jstart1   = fraglist1[idx1pos].iatom
        jfinish1  = fraglist1[idx1pos].iatom+fraglist1[idx1pos].natoms-1
        atoms1    = atomlist1[jstart1:jfinish1]

        jstart3   = fraglist3[idx3pos].iatom
        jfinish3  = fraglist3[idx3pos].iatom+fraglist3[idx3pos].natoms-1
        atoms3    = atomlist3[jstart3:jfinish3]

        push!(fraglist13,FRAGS(ifrag,natoms+1,natoms13,icharge13,1,name13,0.0)) 
        for i in 1:length(atoms1)        
            push!(atomlist13,ATOMS(i+natoms,ifrag,atoms1[i].icharge,atoms1[i].elem[2:2],atoms1[i].coord,0.0,name13))
        end
        natoms = natoms + fraglist1[idx1pos].natoms
        for i in 1:length(atoms3)        
            push!(atomlist13,ATOMS(i+natoms,ifrag,atoms3[i].icharge,atoms3[i].elem,atoms3[i].coord,0.0,name13))
        end
        natoms = natoms + fraglist3[idx3pos].natoms

    end
    global total_atoms13 = natoms
    global total_frags13 = ijfrag

    for i in 1:20
        if i > total_atoms13
            break
        else
            println(atomlist13[i])
        end
    end
    if total_atoms13 > 20
        println("... (more) for Res-CAP ")
    end

    println("total_frags of Res-Cap : ",total_frags13)
    for i in 1:total_frags13
        getmult(fraglist13[i],atomlist13,-1)
        println("fraglist13[",i,"] : ",fraglist13[i])
        #println(" ==== ")
        #println(fraglist[i])
    end


    global atomlist2in12=[]
    global fraglist2in12=[]

    ifrag = 0
    natoms = 0
    for i in 1:length(fraglist12)
        sline = split(fraglist12[i].name,"_")
        idx = sline[2][7:10]
        idxpos = findfirst(==(idx),IDfrags2)

        ifrag = ifrag + 1
        natoms2x  = fraglist2[idxpos].natoms
        icharge2x = fraglist2[idxpos].icharge
        name2x    = fraglist2[idxpos].name

        jstart   = fraglist2[idxpos].iatom
        jfinish  = fraglist2[idxpos].iatom+fraglist2[idxpos].natoms-1
        atoms2x  = atomlist2[jstart:jfinish]

        push!(fraglist2in12,FRAGS(ifrag,natoms+1,natoms2x,icharge2x,1,name2x,0.0))
        for i in 1:length(atoms2x)
            push!(atomlist2in12,ATOMS(i+natoms,ifrag,atoms2x[i].icharge,atoms2x[i].elem,atoms2x[i].coord,0.0,name2x))
        end
        natoms = natoms + fraglist2[idxpos].natoms

    end


    global atomlist3in13=[]
    global fraglist3in13=[]

    ifrag = 0
    natoms = 0
    for i in 1:length(fraglist13)
        sline = split(fraglist13[i].name,"_")
        idx = sline[2][5:8]
        idxpos = findfirst(==(idx),IDfrags3)

        ifrag = ifrag + 1
        natoms3x  = fraglist3[idxpos].natoms
        icharge3x = fraglist3[idxpos].icharge
        name3x    = fraglist3[idxpos].name

        jstart   = fraglist3[idxpos].iatom
        jfinish  = fraglist3[idxpos].iatom+fraglist3[idxpos].natoms-1
        atoms3x  = atomlist3[jstart:jfinish]

        push!(fraglist3in13,FRAGS(ifrag,natoms+1,natoms3x,icharge3x,1,name3x,0.0))
        for i in 1:length(atoms3x)
            push!(atomlist3in13,ATOMS(i+natoms,ifrag,atoms3x[i].icharge,atoms3x[i].elem,atoms3x[i].coord,0.0,name3x))
        end
        natoms = natoms + fraglist3[idxpos].natoms

    end


end 

function wffragPDB(ibase1,ibase2,frag,atoms,pdbfile,pdbinfo,stamp,ipdb) 

    # ipdb :  1   PDB style
    # ipdb : -1   XYZ style

    open(pdbfile,"a+") do wfpdb
        for i in 1:length(frag)
            ibase1 = ibase1 + 1
            i1 = frag[i].iatom
            i2 = frag[i].iatom + frag[i].natoms -1
            for j in i1:i2
                ibase2 = ibase2 + 1
                #println("atoms[",j,"].elem : ",atoms[j].elem)
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

    open(pdbinfo,"a+") do infopdb
        for i in 1:length(frag)
            println(infopdb,stamp," frag[",i,"] : ",frag[i])
        end 
    end 

    return ibase1,ibase2

end

function generateMFCC(fmt,pdbfile)

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
        # fraglist1 --> fraglist(2in12) --> fraglist(3in13) --> fraglist12 --> fraglist13
        # fraglist1       : Selected spike residues
        # fraglist(2in12) : So-called Chain  in hACE
        # fraglist(3in13) : So-called Capped in hACE
        # fraglist12      : Spike's residue - Chain  --<  dimmers
        # fraglist13      : Spike's residue - Capped --<  dimmers
        # 

        open(pdbinfo,"w") do infopdb
            println(infopdb," Spike - hACE  ")        
        end 

        ibase1 = 0 
        ibase2 = 0 
        ibase1,ibase2 = wffragPDB(ibase1,ibase2,fraglist1,atomlist1,pdbfile,pdbinfo,"SPK A",1)
        ibase1,ibase2 = wffragPDB(ibase1,ibase2,fraglist2in12,atomlist2in12,pdbfile,pdbinfo,"ACE A",-1)
        ibase1,ibase2 = wffragPDB(ibase1,ibase2,fraglist3in13,atomlist3in13,pdbfile,pdbinfo,"ACE B",-1)
        ibase1,ibase2 = wffragPDB(ibase1,ibase2,fraglist12,atomlist12,pdbfile,pdbinfo,"DIM A",-1)
        ibase1,ibase2 = wffragPDB(ibase1,ibase2,fraglist13,atomlist13,pdbfile,pdbinfo,"DIM B",-1)
        
        open(pdbfile,"a+") do wfpdb
            println(wfpdb,"END")
        end   
    end 

end 



