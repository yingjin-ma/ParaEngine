# e.g.
#       julia Gau2xyz.jl   GAU_Folders   XYZ_Folder

gaufolder=ARGS[1]
xyzfolder=ARGS[2]
sdffolder=ARGS[3]

if !isdir(gaufolder)
    println("The folder $(gaufolder) is not exist")
    exit("Stopped. Reason: $(gaufolder) is not exist.")
else
    println("Gaussian folder : ", gaufolder )
end

if !isdir(xyzfolder)
    println("The folder $(xyzfolder) is not exist")
    exit("Stopped. Reason: $(xyzfolder) is not exist.")
else
    println("XYZfolder folder : ", xyzfolder )
end

filelist = readdir(gaufolder)
println("")
println(filelist)
println("")

gaulist=[]
#icount = 0
for ifile in filelist
    sline=split(ifile,".")	 
    if length(sline) > 1 
        if uppercase(sline[2])=="GJF" || uppercase(sline[2])=="COM"
            name=sline[1]
	    jfile=string(name, ".xyz") 
	    kfile=string(name, ".smi") 
            id1 = string(gaufolder,"/",ifile)
            id2 = string(xyzfolder,"/",jfile)
            id3 = string(sdffolder,"/",kfile)

            idx = 0
	    open(id1,"r") do readgau
		while !eof(readgau)     
                    line=readline(readgau)
		    sline=split(line)
		    ntmp=length(sline)
		    if occursin("#",line)
                         readline(readgau)
                         readline(readgau)
                         readline(readgau)
                         while !eof(readgau)
                             line=readline(readgau)
  		             if length(split(line))==4
			        idx = idx +1			    
	     		     end
			 end 

	            end     
		end
            end

            open(id1,"r") do readgau
                while !eof(readgau)
                    line=readline(readgau)
                    sline=split(line)
                    ntmp=length(sline)
                    if occursin("#",line)
                         readline(readgau)
                         readline(readgau)
                         readline(readgau)
                         open(id2,"w") do wxyz
			     println(wxyz,idx)
			     print(wxyz," ",name," charge-and-mult ")
                             while !eof(readgau)
                                 line=readline(readgau)
                                 println(wxyz,line)
                             end
                         end

			 run(`obabel -i xyz $(id2) -o smi -O $(id3)`)
                    end
                end
            end
	    

 	end
    end
end    

# All convert
#obabel -i xyz test1/*.xyz -o smi -O iii
#
