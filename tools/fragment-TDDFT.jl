# This code belongs to the ParaEngine/Fcst_sys package
#    - Application of fragment TDDFT calculation (1 excition)
#      with N-body interactions
# 

include("fragment-TDDFT_lib.jl")

# e.g. usage
# julia   fragment-TDDFT.jl  Target-PDB[XXX.pdb]   N-nody interaction[1,2,3]    

target  =             ARGS[1]
nbody   = parse(Int32,ARGS[2])

ResID      = [1]                   # The excited frags
PreNAME    = "GFPA-snapshot-"      # folder's prename
flagMODEL1 = "MODEL"
flagMODEL2 = "ENDMD"
flagADDH   = ""

threshold1 = 5.0  # Considered interaction range for excited frags
threshold2 = 3.0  # Considered interaction range for ground state frags

if !isfile(target)
    println("The target file $(target) is not exist")
    exit("Stopped. Reason: $(target) is not exist.")
else
    println("The input PDB      : ", target )
end

if nbody != 2
    println("Only support 2-body interactions")
    exit("Stopped. Reason: 2 is valid.")
else
    println("N-body interaction : ", nbody )
end

readpdbfile(target,flagADDH)
if flagADDH == "ADDH"
    residueADDH()
    atomlist=[]
    fraglist=[]
    atomlist=atomlist1
    fraglist=fraglist1
end

for i in 1:length(modelfrag)
    idx  = lpad(i,8,"0")
    ifilepdb = string(PreNAME,idx,".pdb")
    #nfrags = model_nfrags[i]
    global natoms12 = model_natoms[i] 
    global fraglist12 = [] 
    global atomlist12 = [] 
    natoms12 = MFCC_TDDFT_PRE(modelfrag[i],modelpdb[i],ResID,threshold1,threshold2,fraglist12,atomlist12)
    println("natoms12 : ",natoms12)
    #println("fraglist12 : ",fraglist12)
    #println("atomlist12 : ",atomlist12)
    MFCC_TDDFT_GEN("PDB",ifilepdb,modelfrag[i],modelpdb[i],fraglist12,atomlist12) 
end 




