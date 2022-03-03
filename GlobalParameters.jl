global runtype,runtype2
global qcdriver
global pdbfile
global LBfile

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
   natoms::Int
  icharge::Int
 multiple::Int
   energy::Float64
end 
