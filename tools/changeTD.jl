global suits = "C:\\Users\\22103\\Desktop\\protein\\6wzo\\opt\\output"

mutable struct ATOM
       idx::Int
      elem::String
    coordx::Float64
    coordy::Float64
    coordz::Float64
     flags::Int
end


function TDchange(suit)

    if !isdir(suit)
        println("The target suit $(suit) is not exist")
        exit("Stopped. Reason: $(suit) is not exist.")
    end
    
    println("suit : ", suit)
    filelist = readdir(suit)
    println(filelist)
    println("")
    
    global names = []
    global atomlist = []
    global atomic = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe"]
    ifrag = 0
    Charges = "0"
    Mul = "0"
     

    td = string(suit,"\\","input-pdb50-freq")
    mkdir(td)

    for ifile in filelist
        name = split(ifile,".")[1]    
        push!(names,name)
        ifile = string(suit,"\\",ifile)
        ifrag = ifrag + 1
        println("ifile : ", ifile)
        tdfile = string(td,"\\",name,".gjf")
        println(tdfile)
        i = 0
        frags = true
        open(ifile,"r") do molread
            while !eof(molread)
                line = readline(molread)
                #println("line")
                if occursin("%mem",line)
                    open(tdfile,"a") do steams
                        write(steams,line,"\n")
                    end
                end
                if occursin("%nproc",line)
                    open(tdfile,"a") do steams
                        write(steams,line,"\n")
                    end
                end
                if occursin("Multiplicity",line)
                    Charges = split(line)[3]
                    Mul = split(line)[6]
                end



                if occursin("Standard orientation:",line)
                    i = i + 1
                    j = 0
                    println("=============",i)
                    readline(molread)
                    readline(molread)
                    readline(molread)
                    readline(molread)                    
                    while(frags)
                        ats = readline(molread)
                        if occursin("------",ats)
                            #frags = false                    
                            break
                        end
                        #println(ats)
                        sp = split(ats)
                        j = j +1
                        id1 = parse(Int,sp[1])
                        ele = parse(Int,sp[2])
                        dx = parse(Float64,sp[4])
                        dy = parse(Float64,sp[5])
                        dz = parse(Float64,sp[6])
                        #println(id1,atomic[ele],dx,dy,dz,i)
                        push!(atomlist,ATOM(id1,atomic[ele],dx,dy,dz,i))
                        #atomlist[j] = ATOMS(sp[1],atomic[ele],dx,dy,dz)
                        
                        
                    end
                    ifrag = j
                    #println("j: ",j)    
                end                
            end
        end
        counts = length(atomlist)
        for k in 1:counts
            if atomlist[k].flags == i
                #for c in 1:ifrag
                    println(atomlist[k])
                #end
            end
        end
        ncount = length(names)  
        for k in 1:ncount
            println(names[k])
        end

        println("Charges: ",Charges,"Mul: ",Mul)
        open(tdfile,"a") do steams
            write(steams,"#P\n\n",name,"\n\n",Charges," ",Mul,"\n")
            #write(steams,charges," ",mul,"\n")
        end
        for k in 1:counts
            if atomlist[k].flags == i
                open(tdfile,"a") do steams
                    println(steams,atomlist[k].elem," ",atomlist[k].coordx," ",atomlist[k].coordy," ",atomlist[k].coordz)
                    #write(steams,charges," ",mul,"\n")
                end
            end
        end
        open(tdfile,"a") do steams
            write(steams,"\n")
        end
        for k in 1:counts
            pop!(atomlist)
        end


    end



end   

TDchange(suits)