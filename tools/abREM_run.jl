#  
# This code belongs to the ParaEngine package
#    - ab inito REM calculation
#      
# 

include("abREM_lib.jl")
include("abREM_comp.jl")

# e.g. usage
# julia   abREM_run.jl   INFORMATION[.info]   RESULTS-FOLDER[FOLDER]

info  = ARGS[1]
rdir  = ARGS[2]

if_REM = true

if !isfile(info)
    println("The target file $(info) is not exist")
    exit("Stopped. Reason: $(info) is not exist.")
else
    println("The input PDB info. file : ", info )
end

if !isdir(rdir)
    println("The result folder $(rdir) is not exist")
    exit("Stopped. Reason: $(rdir) is not exist.")
else
    println("Result folder : ", rdir)
end

# =========== INFORMATION  =============

readINFO(info)

# ============   RESULTS   =============

EXdata = readELST(rdir,"FULL","Overlap",if_REM)

# ==========  ab initio REM  ===========
# note : updated fraglist + EXdata

REM_comp(fraglist,EXdata,if_REM)


