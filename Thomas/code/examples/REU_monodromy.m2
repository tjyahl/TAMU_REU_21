--Small file to illustrate numerical homotopy methods and Monodromy.m2

--M2 doesn't have built in methods for path tracking/solving.
----We use the package NumericalAlgebraicGeometry.m2 (NAG)
restart
loadPackage("NumericalAlgebraicGeometry")



--always make your ring
R = CC[x,y]

--start system and solutions
----"solveSystem" is a method from NAG to well.. solve systems.
F = { x^3-1 , y^2-2 }
startSolns = solveSystem F

--target system we want solutions to
G = { x^2-y^2+2*x-1 , x*y+2*y^2+y+1 }



--track solutions from start system to target system
----"track" (from NAG) tracks solutions!
----the option "gamma=>random(CC)" isn't necessary.. but it kind of is.
targetSolns = track(F,G,startSolns,gamma=>random(CC))

--notice that targetSolns is a list of "Point"s. 
----these points have different attributes like "status" and "coordinates"
targetSolns#0
status targetSolns#0
coordinates targetSolns#0

targetSolns#2
status targetSolns#2
coordinates targetSolns#2

apply(targetSolns,status)
apply(targetSolns,coordinates)



--we can use path-tracking as above to compute monodromy permutations.
----notice "startSolns" and "trackedSolns" are permutations of one another.
F = { x^3-1 , y^2-2 }
startSolns = solveSystem F

G1 = { x^3+1-ii , y^2-3+ii }
G2 = { x^3+1+ii , y^2-1+ii }

trackedSolns = track(F,G1,startSolns)
trackedSolns = track(G1,G2,trackedSolns)
trackedSolns = track(G2,F,trackedSolns)

startSolns
trackedSolns



--luckily, smart people make packages that does the above for us.
----this loads the package.
restart
loadPackage("Monodromy")

--recall we represent sparse systems by supports, points in the plane determining exponents.
----here the supports for each equation are written as a matrix.
A = {matrix{{0,0,1,2,2},{0,2,1,0,2}},matrix{{0,0,1,2,2},{0,2,1,0,2}}}



--"sparseMonodromy" creates a "Monodromy" object.
----this contains information needed to compute monodromy as well as count our permutations.
----("Solver=>M2" is needed unless you have phc)
M = sparseMonodromy(A,Solver=>M2)



--this is the information a "Monodromy" object keeps track of.
----"M#family" is the corresponding family of polynomial systems.
----"M#basePoint" and "M#baseSolutions" are a start system and solutions.
----"M#group" is empty since we haven't told it to compute any permutations!
M#family
M#basePoint
M#baseSolutions
M#group



--compute monodromy loops with.. "monodromyLoop". 
----Can give it an explicit loop or tell it to do a number of random loops.
monodromyLoop(M,10)

--to see permutations, use this.
printCycleTally M


--OTHER STUFF--
--can save "Monodromy" objects to a file with "saveMonodromy"
--can also load "Monodromy" objects from a file with "loadMonodromy"
--"testAlternatingMonodromy" and "printForGAP" will likelyl be useful.
--we can add other functionality as needed!
