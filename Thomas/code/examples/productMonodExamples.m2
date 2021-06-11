--small file illustrating when Galois groups are products.
restart
loadPackage("Monodromy",Reload=>true)
loadPackage("cycleCombinatorics",Reload=>true)

A = {matrix{{0,0,1,1},{0,1,0,1},{0,0,0,0}},matrix{{0,0,1,1},{0,1,0,1},{0,0,0,0}},matrix{{0,0,0},{0,0,0},{0,1,2}}}
M = sparseMonodromy A
monodromyLoop(M,30)
printCycleTally M
printCycleTypes M
productTypes(2,2)

A = {matrix{{0,0,0,1,1,2},{0,1,2,0,1,0},{0,0,0,0,0,0}},matrix{{0,0,0,1,1,2},{0,1,2,0,1,0},{0,0,0,0,0,0}},matrix{{0,0,0},{0,0,0},{0,1,2}}}
M = sparseMonodromy A
monodromyLoop(M,30)
printCycleTally M
printCycleTypes M
productTypes(4,2)

A = {matrix{{0,0,0,1,1,2},{0,1,2,0,1,0},{0,0,0,0,0,0}},matrix{{0,0,0,1,1,2},{0,1,2,0,1,0},{0,0,0,0,0,0}},matrix{{0,2,4},{0,1,2},{0,1,2}}}
M = sparseMonodromy A
monodromyLoop(M,30)
printCycleTally M
printCycleTypes M
productTypes(4,2)

