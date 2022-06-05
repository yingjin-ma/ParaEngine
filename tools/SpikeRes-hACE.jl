#  
# This code belongs to the ParaEngine package
#    - Application of spikeRes-hACE interactions
# 

include("SpikeRes-hACE_lib.jl")

# e.g. usage
# julia   SpikeRes-hACE.jl   ACE2_Ab-Omicron_Q493R_0_0.pdb[Spike]   ACE2_Ab-Omicron_Q493R_0_0[MPCC_FOLDER]

target    = ARGS[1]
hACE_subs = ARGS[2]

flagSpike = "E"
flaghACE  = "A"
flagADDH  = "ADDH"

ResID=[339, 371, 373, 375, 417, 440, 446, 477, 478, 484, 493, 496, 498, 501, 505]

threshold = 8.0  # Considered interaction range 

if !isfile(target)
    println("The target file $(target) is not exist")
    exit("Stopped. Reason: $(target) is not exist.")
else
    println("input PDB : ", target )
end

if !isdir(hACE_subs)
    println("The hACE_subs dir $(hACE_subs) is not exist")
    exit("Stopped. Reason: $(hACE_subs) is not exist.")
else
    println("input hACE_subs : ", hACE_subs )
end

readpdbfile(target,flagADDH)
if flagADDH == "ADDH"
    residueADDH()
    atomlist=[]
    fraglist=[]
    atomlist=atomlist1
    fraglist=fraglist1
end 
readXYZfile(hACE_subs,flaghACE) # flag of hACE : only pick up the hACE residues


sline0=split(target,"/")
ntmp0=length(sline0)
ssline0=split(sline0[ntmp0],".")
pdbfile0=string(ssline0[1],"-",string(threshold),"A-MFCC.pdb")
println("MFCC pdbfile generated : ", pdbfile0)

generateMFCC_PRE(threshold)
generateMFCC("PDB",pdbfile0)




