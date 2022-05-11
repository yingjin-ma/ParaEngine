

@everywhere function nodes_monitor(folder, lkfile)

    inum     =  0
    lockfile = string(folder,"/",lkfile)
    println(" lockfile : ",lockfile) 

    monidir  = string(folder,"/","monitor")
    
    if !isdir(monidir)
        mkdir(monidir)
    end   

    idone = 0
    while true
        # Sample interval is 60s 
        tint=60

        if isfile(lockfile)

            println("lockfile was locked!")

            idx=lpad(inum,6,"0") 
            LKfile = string("hostlock.",idx)
            LKfile = string(monidir,"/",LKfile)

            println("LKfile ", LKfile)

            global ifcp = -1
            while ifcp !=1
                try 
                    run(`cp $(lockfile) $(LKfile)`) 
                    ifcp = 1 
                    println("success of sync of LKfile") 
                catch err
                    sleep(0.1) 
                    println("extra 0.1s for waiting sync of LKfile") 
                end  
            end  

            nlines = countlines(LKfile)

            println("Checkpoint : ",inum," th, active nodes : ",nlines)  

            sleep(tint)
            println("  ")

            if nlines == 0 
                idone = idone + 1  
                if idone == 5
                    break
                end 
            end 
 
            inum = inum + 1

        else

            println("lockfile was not yet generated!")
            sleep(tint)

        end
    end

end 

