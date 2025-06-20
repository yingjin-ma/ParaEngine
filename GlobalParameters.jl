runtype  = ""
runtype2 = ""
qcdriver = ""
pdbfile  = ""
LBfile   = ""
global LBcube   = []
MNcube   = []
QCpara   = []
closedshell = false

@everywhere pepath   = "./"
@everywhere workdir  = "./workdir"
@everywhere IOrecord = ""
@everywhere IFDONE   = false
@everywhere overload = 0.0
@everywhere nrepeat  = -1
@everywhere iffifo = true

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

CDcube = []
CDtask = []
mutable struct CODED
    inode::Int
     ipos::Int
    itask::Int
    icopy::Int
end 


