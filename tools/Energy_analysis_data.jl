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
#       julia Energy_analysis.jl ../ACE2_Ab-Omicron_Q493R_0_0-8.0A-MFCC.info ../COVID-10nodes-TEST-DFT-TRY4-S  

target=ARGS[1]

data=["ACE2_Ab-Wild-type", "ACE2_Ab-Omicron_Q493K", "ACE2_Ab-Omicron_Q493R"]
nresi = 15
threshold = 1.0

if !isdir(target)
    println("The target suit $(target) is not exist")
    exit("Stopped. Reason: $(target) is not exist.")
else
    println(" folder : ", target )
end

filelist = readdir(target)
println("")
println(filelist)
println("")

VVV=[]
for i in 1:length(data)
    println(i,data[i])
    VV=[]
    for ifile in filelist
        if occursin(data[i],ifile)
            id = string(target,"/",ifile)
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


Vres=[] 
for i in 1:length(data)
    #println("  ",i)
    Vresi = Array{Float64}(undef,nresi)
    VresI = Array{Int32}(undef,nresi)
    for k in 1:nresi
        Vresi[k] = 0.0
        VresI[k] = 0
    end 
    for j in 1:length(VVV[i])
        # println("  ",VVV[i])
        for k in 1:nresi
            if abs(VVV[i][j][k]) <= 5.0        
                Vresi[k] = Vresi[k] + VVV[i][j][k]
                VresI[k] = VresI[k] + 1               
            else
                println("remove VVV[",i,"][",j,"][",k,"] : ", VVV[i][j][k]) 
            end  
        end 
    end
    for k in 1:nresi
        Vresi[k]=Vresi[k]/VresI[k]
    end
    println("count : ",VresI) 
    push!(Vres,Vresi)  
end 

for i in 1:length(data)
    println(Vres[i]*627.5094)
end 
