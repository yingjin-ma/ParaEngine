using LinearAlgebra

mutable struct EXREM2
    idx1::Int
    idx2::Int
    eigs::Any
    vecs::Any
end



function schmidt_orthogonalization(vectors)
    n = length(vectors)
    ortho_vectors = Vector{Vector{Float64}}(undef, n)
    ortho_vectors[1]=vectors[1]
    for i in 2:n
        v = vectors[i]
        for j in 1:i-1
            temp=dot(v, ortho_vectors[j])/dot(ortho_vectors[j],ortho_vectors[j])
            v-=temp*ortho_vectors[j]
        end
        ortho_vectors[i] = v
    end
    for i in 1:n
        ortho_vectors[i] = normalize(ortho_vectors[i])  # 
    end
    return ortho_vectors
end

function REM_proj1(EX1,EX2)
# --> project EX1 into Ex2

    Cmat1 = EX1.MO[1].C
    Cmat2 = EX2.MO[1].C 
    Smat1 = EX1.MO[1].S 
    Smat2 = EX2.MO[1].S 
 
    Smo = similar(Cmat1)
 
    # println("Cmat1 : ", Cmat1) 
    # println("transpose(Cmat1) : ", transpose(Cmat1)) 
 
    Tmat1 = transpose(Cmat1) * Smat1
    Smo   = Tmat1 * Cmat2 

    # println("Smo : ", Smo )

    Vocc1 = []
    for i in 1:EX1.MO[1].nele/2
        push!(Vocc1,Int(i))  
    end      
    Vocc2 = []
    for i in 1:EX2.MO[1].nele/2
        push!(Vocc2,Int(i))
    end  

#    println("Vocc1 : ", Vocc1)
#    println("Vocc2 : ", Vocc2)

    CCproj = zeros(Float64,length(EX1.Ex),length(EX2.Ex))
    for i in 1:length(EX1.Ex)

        v1=zeros(Int,length(Vocc1))
        for ii in 1:length(Vocc1)
            v1[ii] = Vocc1[ii]
        end 
        #println("EX1.Ex[", i, "].exst : ", EX1.Ex[i].exst)

        Cprojs = []
        for i1 in 1:length(EX1.Ex[i].exst)

            vv1=zeros(Int,length(v1))
            for ii in 1:length(v1)
                vv1[ii] = v1[ii]
            end

            imo1 = EX1.Ex[i].exst[i1].exmo1  
            imo2 = EX1.Ex[i].exst[i1].exmo2  
            cci1 = EX1.Ex[i].exst[i1].exci
            #println(" vv1 : ", vv1) 
            for i2 in 1:length(vv1)
                if vv1[i2] == imo1
                    vv1[i2] = imo2
                end 
            end
            #println("imo1 = ",imo1, "imo2 = ", imo2, " ; vv1 : ", vv1) 

 
            Cproj = []
            for j in 1:length(EX2.Ex)
                v2=zeros(Int,length(Vocc2))
                for jj in 1:length(v1)
                    v2[jj] = Vocc2[jj]
                end
                #println("EX2.Ex[", j, "].exst : ", EX2.Ex[j].exst)

                dv2 = 0.0
                for j1 in 1:length(EX2.Ex[j].exst)

                    vv2=zeros(Int,length(v2))
                    for jj in 1:length(v2)
                        vv2[jj] = v2[jj]
                    end 
                    jmo1 = EX2.Ex[j].exst[j1].exmo1
                    jmo2 = EX2.Ex[j].exst[j1].exmo2
                    ccj1 = EX2.Ex[j].exst[j1].exci
                    for j2 in 1:length(vv2)
                        if vv2[j2] == jmo1
                            vv2[j2] = jmo2
                        end
                    end
                    #println("jmo1 = ",jmo1, "jmo2 = ", jmo2, " ;      vv2 : ", vv2) 

                    dvv2 = 0.0
                    SmoREM = zeros(Float64,length(vv1),length(vv2))
                    for ii in 1:length(vv1)      
                        for jj in 1:length(vv2)      
                            SmoREM[ii,jj] = Smo[vv1[ii],vv2[jj]]
                        end 
                    end 
                    dvv2 = det(SmoREM)
                    dv2 = dv2 + ccj1 * dvv2

                end 
                push!(Cproj,dv2)
            end 
       
            Vdvv1 = []
            for j in 1:length(EX2.Ex)
                dvv1 = cci1 * Cproj[j]  
                push!(Vdvv1, dvv1)
            end 
            push!(Cprojs,Vdvv1)
        end
 
        for j in 1:length(EX2.Ex)
            dtmp = 0.0
            for i1 in 1:length(EX1.Ex[i].exst)
                dtmp = dtmp + Cprojs[i1][j]
            end 
            CCproj[i,j] = dtmp            
        end  

    end   

    CCvec = Vector{Vector{Float64}}(undef,length(EX2.Ex))
    for i in 1:length(EX1.Ex)
        CCvec[i] = CCproj[i,:]
#        println("CCvec[i]  :  ", CCvec[i]) 
    end 

#   println("CCproj         :  ", CCproj)
    ortho_CCvec = schmidt_orthogonalization(CCvec)
#    println("CCvec (orth)  :  ", ortho_CCvec)

    return ortho_CCvec

end 

function REM_proj2(EX1,EX2,EX3)
# --> project EX1 into Ex2(*)Ex3

    Cmat1 = EX1.MO[1].C
    Cmat2 = EX2.MO[1].C
    Cmat3 = EX3.MO[1].C
    Smat  = EX1.MO[1].S

    m1,n1 = size(Cmat1)
    m2,n2 = size(Cmat2)
    m3,n3 = size(Cmat3)

    #println("m1,n1,m2,n2,m3,n3 ", m1,n1,m2,n2,m3,n3)

    Cmat23 = similar(Cmat1) 
    zero23 = zeros(Float64,m2,n3)
    zero32 = zeros(Float64,m3,n2)
    Cmat23[1:m2,1:n2] = Cmat2
    Cmat23[m2+1:m2+m3,1:n2] = zero32
    Cmat23[1:m2,n2+1:n2+n3] = zero23
    Cmat23[m2+1:m2+m3,n2+1:n2+n3] = Cmat3

    println(" Q inner-product : TBM --> Smat: ", Smat)
    Tmat1 = transpose(Cmat1) * Smat
    println(" Q inner-product : TBM --> Tmat1 : ", Tmat1)
    Smo   = Tmat1 * Cmat23
    println(" Q inner-product : check --> Smo : ", Smo)

    # Smo = similar(Cmat1)     
    #println("Cmat23", size(Cmat23))

    Vocc1 = []
    for i in 1:EX1.MO[1].nele/2
        push!(Vocc1,Int(i))
    end
    Vocc2 = []
    for i in 1:EX2.MO[1].nele/2
        push!(Vocc2,Int(i))
    end 
    Vocc3 = []
    for i in 1:EX3.MO[1].nele/2
	push!(Vocc3,Int(i))
    end

#    println("Vocc1 : ", Vocc1)
#    println("Vocc2 : ", Vocc2)
#    println("Vocc3 : ", Vocc3)

    println("length(EX1.Ex) ", length(EX1.Ex))
    println("length(EX2.Ex) ", length(EX2.Ex))
    println("length(EX3.Ex) ", length(EX3.Ex))

    Cproj = zeros(Float64,length(EX1.Ex),length(EX2.Ex),length(EX3.Ex))
    for i in 1:length(EX1.Ex)
        for j in 1:length(EX2.Ex)
            for k in 1:length(EX3.Ex)
                dijk = 0                 
                for i1 in 1:length(EX1.Ex[i].exst)
                    vv1 = zeros(Int,length(Vocc1))
                    for ii in 1:length(Vocc1)
                        vv1[ii] = Vocc1[ii]
                    end
                    imo1 = EX1.Ex[i].exst[i1].exmo1
                    imo2 = EX1.Ex[i].exst[i1].exmo2
                    cci1 = EX1.Ex[i].exst[i1].exci
                    for i2 in 1:length(vv1)
                        if vv1[i2] == imo1
                            vv1[i2] = imo2
                        end
                    end
 
                    for j1 in 1:length(EX2.Ex[j].exst)
                        vv2=zeros(Int,length(Vocc2))
                        for jj in 1:length(Vocc2)
                            vv2[jj] = Vocc2[jj]
                        end                        
                        jmo1 = EX2.Ex[j].exst[j1].exmo1
                        jmo2 = EX2.Ex[j].exst[j1].exmo2
                        ccj1 = EX2.Ex[j].exst[j1].exci
                        for j2 in 1:length(vv2)
                            if vv2[j2] == jmo1
                                vv2[j2] = jmo2
                            end
                        end 

                        for k1 in 1:length(EX3.Ex[k].exst) 
                            vv3=zeros(Int,length(Vocc3))
                            for kk in 1:length(Vocc3)
                                vv3[kk] = Vocc3[kk] 
                            end
                            kmo1 = EX3.Ex[k].exst[k1].exmo1
                            kmo2 = EX3.Ex[k].exst[k1].exmo2
                            cck1 = EX3.Ex[k].exst[k1].exci                        
                            for k2 in 1:length(vv3)
                                if vv3[k2] == kmo1
                                    vv3[k2] = kmo2   # plus the offset
                                end
                                vv3[k2] = vv3[k2] + m2
                            end
                         
                            #println("vv2 : ", vv2)
                            #println("vv3 : ", vv3)
                            vv23 = vcat(vv2,vv3)
                            #println("vv23 : ", vv23)
                            SmoREM = zeros(Float64,length(vv1),length(vv23))
                            for ii in 1:length(vv1)
                                for jj in 1:length(vv23)
                                    SmoREM[ii,jj] = Smo[vv1[ii],vv23[jj]]
                                end
                            end
                            #println("det(SmoREM) : ", det(SmoREM))
                            dv = det(SmoREM) * cci1 * ccj1 * cck1
                            dijk = dijk + dv
                        end
                    end
                end
                Cproj[i,j,k] = dijk
            end 
        end
    end 

    #println("Cproj ", Cproj)
    #exit(1) 

    CCproj = zeros(Float64,length(EX1.Ex),length(EX2.Ex)*length(EX3.Ex))
    for i in 1:length(EX1.Ex)
        jk=0
        for j in 1:length(EX2.Ex)
            for k in 1:length(EX3.Ex)
                jk = jk + 1
                CCproj[i,jk] = Cproj[i,j,k]
            end 
        end 
    end

    CCvec = Vector{Vector{Float64}}(undef,length(EX1.Ex)-1)
    for i in 2:length(EX1.Ex)
        CCvec[i-1] = CCproj[i,2:3]
#        println("CCvec[i]  :  ", CCvec[i])
    end

#   println("CCproj         :  ", CCproj)
    ortho_CCvec = schmidt_orthogonalization(CCvec)
#    println("CCvec (orth)  :  ", ortho_CCvec)

    return ortho_CCvec

end


function REM_comp(fraglist,EXdata,if_REM=false)

    nmonomers = 0 
    idx_dimer=[]
    # idex0 = zeros(Int,length(fraglist))  # record the index
    idex0 = []
    for i in 1:length(fraglist) 
        println("fraglist[",i,"] in REM_comp : ", fraglist[i])
        if occursin("_",fraglist[i].name) 
            push!(idx_dimer,i)
        else 
            nmonomers = nmonomers + 1
        end 
        push!(idex0,EXdata[i].info.idx)
    end
    Vperm = sortperm(idex0)
    println("idex0 : ", idex0)
    println("Vperm : ", Vperm)

    ndimers = length(fraglist) - nmonomers
    println("idx_dimer : ", idx_dimer)
   
    idx_mon1=[]
    idx_mon2=[]
    for i in 1:length(idx_dimer)
        idx = idx_dimer[i]
        sline=split(fraglist[idx].name,"_")
        for j in 1:length(fraglist) 
            if occursin(sline[1],fraglist[j].name)
                if (j in idx_dimer) == false # not in
                    push!(idx_mon1,j)
                end
            end 
            if occursin(sline[2],fraglist[j].name)
                if (j in idx_dimer) == false # not in
                    push!(idx_mon2,j)
                end 
            end 
        end 
    end

    println("idx_mon1 : ", idx_mon1)
    println("idx_mon2 : ", idx_mon2)

    REMlist = []
    for i in 1:length(idx_dimer)
        idx1 = 0 
        idx2 = 0 
        idx3 = 0 
        for j in 1:length(EXdata)
            if idx_dimer[i] == EXdata[j].info.idx
                # println(i," ",j)
                idx1 = j
            end 
            if idx_mon1[i]  == EXdata[j].info.idx
                # println(i," ",j)
                idx2 = j
            end 
            if idx_mon2[i]  == EXdata[j].info.idx
                # println(i," ",j)
                idx3 = j
            end 
        end
             # REM_proj1(EXdata[idx3],EXdata[idx3])                 # Test the code
        Pvec = REM_proj2(EXdata[idx1],EXdata[idx2],EXdata[idx3])    # Dimer --> Monomers

        push!(REMlist,EXREM2(idx_mon1[i],idx_mon2[i],[EXdata[idx1].Ex[2].ex_ene EXdata[idx1].Ex[3].ex_ene],Pvec))

    end

#    println("REMlist : ", REMlist)   
#    for i in 1:length(EXdata)
#        println(EXdata[i].info)
#    end 
  
    if if_REM

        HREM = zeros(Float64,nmonomers,nmonomers) 
        for i in 1:nmonomers
            idx = Vperm[nmonomers - i + 1]
            HREM[i,i] = EXdata[idx].Ex[2].ex_ene 
            println(i," ",HREM[i,i])
        end 
        
        for i in 1:length(idx_dimer)
            #println("REMlist[",i,"] : ", REMlist[i])   
            # Get the index
            idx1 = REMlist[i].idx1
            idx2 = REMlist[i].idx2

            Htmp = zeros(2,2)
            Ctmp = zeros(2,2)
            Etmp = zeros(2,2)

            # H = C*e*inv(C)
            Etmp[1,1] = REMlist[i].eigs[1] 
            Etmp[2,2] = REMlist[i].eigs[2] 
            for j in 1:2
                for k in 1:2
                    Ctmp[k,j] = REMlist[i].vecs[j][k]
                end 
            end
            Htmp = Ctmp * Etmp * inv(Ctmp)
            #println("Htmp : ", Htmp) 

            # H^eff(ij) - H^eff(i) - H^eff(j)
            Htmp[1,1] = Htmp[1,1] - EXdata[idx1].Ex[2].ex_ene
            Htmp[2,2] = Htmp[2,2] - EXdata[idx2].Ex[2].ex_ene
            for j in 1:2
                jdx = nmonomers - idx1 + 1
                for k in 1:2
                    kdx = nmonomers - idx2 + 1
                    HREM[jdx,kdx] = HREM[jdx,kdx] + Htmp[j,k]
                end 
            end
 
        end 

        EREM, VREM = eigen(HREM)

        for i in 1:5
            println(" EREM[",i,"], ", EREM[i], "  ",length(VREM[:,i]) )
            for j in 1:length(VREM[:,i])
                if abs(VREM[j,i]) >  0.01
                    println(" VREM[",length(VREM[:,i]) - j + 1,"] : ", VREM[j,i] )
                end 
            end 
        end 

    else
        println("Multiple excitations are not implemented yet")
        exit(1)
    end

end 


