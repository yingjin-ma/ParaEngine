mutable struct RESI
      idx::Int
     Eint::Float64
end

mutable struct FRAGS
      idx1::Int
      idx2::Int
   residue
end

# e.g.
#       julia Energy_analysis.jl  ??? 

target=ARGS[1]

data=["ACE2_Ab-Wild-type", "ACE2_Ab-Omicron_Q493R"]
nresi = 15
threshold = 0.2

if !isdir(target)
    println("The target suit $(target) is not exist")
    exit("Stopped. Reason: $(target) is not exist.")
else
    println(" folder : ", target )
end

filelist = readdir(target)
println("")
#println(filelist)
println("")

VVV=[]
for i in 1:length(data)
    println(i," ", data[i])
    VV=[]
    for ifile in filelist
        if occursin(data[i],ifile)
            id = string(target,"/",ifile)
            println("id : ",id)
            sline=split(id,"_")
            ntmp=length(sline)
            ssline=sline[ntmp-1]
            idx1 = parse(Int32,sline[ntmp-2])
            idx2 = parse(Int32,split(ssline,"-")[1])
            println("id : ", id, " ", idx1, " ", idx2)
 
            V1=[] 
            open(id,"r") do readlog
                for i in 1:nresi
                    line  = readline(readlog)
                    #println(line)
                    sline = split(line)
                    ntmp  = length(sline)
                    dv    = parse(Float64,sline[3])
                    push!(V1,dv)
                end 
            end 
            push!(VV,V1)
            println(i,"  ",ifile,"  ",V1)
        end  
    end
    push!(VVV,VV)    
end 


Vres =[] 
VVres=[]
for i in 1:length(data)
    #println("  ",i)
    Vresi = Array{Float64}(undef,nresi)
    VresI = Array{Int32}(undef,nresi)
    v1=[] 
    for k in 1:nresi
        Vresi[k] = 0.0
        VresI[k] = 0
        push!(v1,[])
    end 
    for j in 1:length(VVV[i])
        # println("  ",VVV[i])
        for k in 1:nresi
            if abs(VVV[i][j][k]) <= threshold
                Vresi[k] = Vresi[k] + VVV[i][j][k]
                VresI[k] = VresI[k] + 1               
                push!(v1[k],VVV[i][j][k]*627.5094) 
            else
                println("remove VVV[",i,"][",j,"][",k,"] : ", VVV[i][j][k]) 
            end  
        end 
    end 
    push!(VVres,v1)
    #println(v1)
    for k in 1:nresi
        Vresi[k]=Vresi[k]/VresI[k]
    end
    println("count : ",VresI) 
    push!(Vres,Vresi)  
end 

for i in 1:length(data)
    #println(" length(data[i]) : ",length(data[i]))
    for j in 1:15
        vvv = sort(VVres[i][j])  
        lenvvv = length(vvv) 
        lenmid1 = round(Int,lenvvv/2)  
        lenmid2 = round(Int,lenvvv/2) + 1
        #println("lenmid1 : ", lenmid1," lenmid2 : ",lenmid2) 
        if lenvvv%2 == 0
            dvmid2 = (vvv[lenmid1]+vvv[lenmid2])/2
            vvv1 = vvv[1:lenmid1]
            vvv2 = vvv[lenmid2:lenvvv]
        else
            dvmid2 = vvv[lenmid2]
            vvv1 = vvv[1:lenmid1]
            vvv2 = vvv[lenmid2+1:lenvvv]
        end  
 
        lenvvv1 = length(vvv1)
        lenmid1 = round(Int,lenvvv1/2)
        lenmid2 = round(Int,lenvvv1/2) + 1
        if lenvvv1%2 == 0
            dvmid1 = (vvv1[lenmid1]+vvv1[lenmid2])/2
        else
            dvmid1 = vvv1[lenmid2]
        end 

        lenvvv2 = length(vvv2)
        lenmid1 = round(Int,lenvvv2/2)
        lenmid2 = round(Int,lenvvv2/2) + 1
        if lenvvv2%2 == 0
            dvmid3 = (vvv2[lenmid1]+vvv2[lenmid2])/2
        else
            dvmid3 = vvv2[lenmid2]
        end              
        if i == 1
            println(j-0.15,"  ", vvv[1], " ", dvmid1, " ",dvmid2 ," ", dvmid3," ", vvv[lenvvv]," ",Vres[i][j]*627.5094)         
        end 
#        if i == 2 
#            println(j+0.00,"  ", vvv[1], " ", dvmid1, " ",dvmid2 ," ", dvmid3," ", vvv[lenvvv]," ",Vres[i][j]*627.5094)         
#        end 
        if i == 2
            println(j+0.15,"  ", vvv[1], " ", dvmid1, " ",dvmid2 ," ", dvmid3," ", vvv[lenvvv]," ",Vres[i][j]*627.5094)         
        end 
    end
end 

for i in 1:length(data)
    println(" length(data[i]  ",length(data[i]))
    for j in 1:15
        println(j,"  ",Vres[i][j]*627.5094,)
    end
end


 
