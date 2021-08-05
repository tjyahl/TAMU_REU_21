--Esterov's Example and variants
restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,1,2,2},{0,2,1,0,2}},matrix{{0,0,1,2,2},{0,2,1,0,2}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,2},{0,2,0}},matrix{{0,0,2,2},{0,2,0,2}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,1,2},{0,2,1,0}},matrix{{1,0,2,2},{1,2,0,2}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


--scaling everything by 2 seems to always work
restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {2*matrix{{-1,0,2,3},{0,3,-1,2}},2*matrix{{1,0},{0,1}}*matrix{{-1,0,2,3},{0,3,-1,2}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,2,4},{0,2,4,0}},matrix{{1,0},{0,1}}*matrix{{0,0,2,4},{0,2,4,0}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M

--some triangular examples
restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,2},{0,2,0}},matrix{{0,1,2},{0,2,4}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,2,2},{0,2,0,2}},matrix{{0,1,2},{0,1,2}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,2},{0,2,0}},matrix{{0,1,2},{0,2,4}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,2},{0,2,0}},matrix{{0,random(-6,6)},{0,random(-6,6)}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{-1,-1,-1,-1,1,1,1,1},{-1,-1,1,1,-1,-1,1,1},{-1,1,-1,1,-1,1,-1,1}},matrix{{0,1,2,3,4},{-1,1,3,5,7},{0,0,0,0,0}},matrix{{0,3,4},{0,3,4},{0,6,8}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{-1,-1,1,1},{-1,1,-1,1},{0,0,0,0}},matrix{{-1,-1,1,1},{-1,1,-1,1},{0,0,0,0}},matrix{{0,0,0},{0,0,0},{0,1,2}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


--I like this example
restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,2},{0,2,0},{0,0,0}},matrix{{0,2,2},{2,2,0},{0,0,0}},matrix{{0,0,0,0},{0,0,0,0},{0,1,2,3}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


--binomial example, but still interesting. monodromy is isomorphic to cokernel of exponent matrix.
--equivalent to A = {matrix{{0,3}}}
restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,1},{0,2}},matrix{{0,2},{0,1}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,1,1},{0,0,1},{0,0,0}},matrix{{0,2},{0,0},{0,0}},matrix{{0,0,0},{0,0,0},{0,1,2}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,2,1},{0,2,0,1},{0,0,0,0}},matrix{{0,0,2,1},{0,2,0,1},{0,0,0,0}},matrix{{0,0,0},{0,0,0},{0,1,2}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,2},{0,2,0}},matrix{{0,1,2},{0,2,4}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


--odd examples
restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{2,0,0,0},{0,2,0,0},{0,0,2,0}},matrix{{2,0,0,2},{0,2,0,2},{0,0,2,2}},matrix{{0,1,2,3},{0,2,1,3},{0,1,1,2}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M


restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)
A = {matrix{{0,0,1,2},{0,2,1,0}},matrix{{0,1,2},{0,1,2}}}
M = sparseMonodromy(A)
monodromyLoop(M,30)
printCycleTally M
testAlternatingMonodromy M

--use to check lattice index
first smithNormalForm fold(A,(M,N)->M|N)

--vertices of minkowski sum
vertices sum apply(A,convexHull)

--
loadPackage("Polyhedra",Reload=>true)
P = sum (A/convexHull)
vertices P
M = first facets P


