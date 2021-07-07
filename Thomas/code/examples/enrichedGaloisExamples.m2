--file with examples of enriched Galois groups.
restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2")

--Esterov's example
A = {matrix{{0,0,1,2,2},{0,2,1,0,2}},matrix{{0,0,1,2,2},{0,2,1,0,2}}}
M = sparseMonodromy(A,Solver=>M2)
monodromyLoop(M,100)
printCycleTally M












