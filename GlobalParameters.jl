runtype  = ""
runtype2 = ""
qcdriver = ""
pdbfile  = ""
LBfile   = ""
LBmat    = [[]]
workdir  = "./workdir"

struct ATOMS
      idx::Int
    ifrag::Int
  icharge::Int 
     elem::String
    coord::Tuple{Float64,Float64,Float64}
       ZZ::Float64
end 

struct FRAGS
      idx::Int
    iatom::Int
   natoms::Int
  icharge::Int
 multiple::Int
   energy::Float64
end 


