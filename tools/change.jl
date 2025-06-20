println("000")
global inputdir = "./input-test10"
global pdbfile = "./original_document/test10.sdf"
function changesdf()
    println("111")	
    run(`obabel $(pdbfile) -ogjf -O $(inputdir)/test.gjf -m`)
    println("222")
end
println("333")
changesdf()
println("444")
