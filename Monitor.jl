

@everywhere function nodes_monitor(folder, lkfile)

    inum     =  0
    lockfile = string(folder,"/",lkfile)
    println(" lockfile : ",lockfile) 

    monidir  = string(folder,"/","monitor")
    
    if !isdir(monidir)
        mkdir(monidir)
    end   

    while true
        # Sample interval is 60s 
        tint=60

        if isfile(lockfile)

            println("lockfile was locked!")

            idx=lpad(inum,6,"0") 
            LKfile = string("hostlock.",idx)
            LKfile = string(monidir,"/",LKfile)

            println("LKfile ", LKfile)
            run(`cp $(lockfile) $(LKfile)`) 

            nlines = countlines(LKfile)

            println("Checkpoint : ",inum," th, active nodes : ",nlines)  

            sleep(tint)
            println("  ")

            if nlines == 0 
                break
            end 
 
            inum = inum + 1

        else

            println("lockfile was not yet generated!")
            sleep(tint)

        end
    end

end 

