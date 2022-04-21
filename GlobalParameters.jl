runtype  = ""
runtype2 = ""
qcdriver = ""
pdbfile  = ""
LBfile   = ""
LBcube   = []
MNcube   = []
QCpara   = []
@everywhere workdir  = "./workdir"
@everywhere IOrecord = ""
@everywhere IFDONE   = false


mutable struct ATOMS
      idx::Int
    ifrag::Int
  icharge::Int 
     elem::String
    coord::Tuple{Float64,Float64,Float64}
       ZZ::Float64
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

@everywhere mutable struct TASKS
      idx::Int
   fragid::Int
   folder::String 
   infile::String
  outfile::String
   qcsoft::String
   nnodes::Int
end

