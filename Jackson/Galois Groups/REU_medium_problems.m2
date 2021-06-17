--An ongoing list of medium problems to examine for understanding.

-----------------
----Problem 1----
-----------------
--Consider the system 'F' below.
----a) Determine the number of solutions to 'F' using 'solveSystem'.
----
----b) Draw the support for each polynomial in 'F'.
----
----c) Show that the Galois group corresponding to these supports is the symmetric group S_3.
restart
loadPackage("Monodromy")
R = CC[x,y]
F = {1+2*x+3*y+5*x*y+7*y^2, 1-3*x+7*y-13*x*y}
solveSystem F -- there are 3 solutions
A = {matrix{{0,1,0,1,0},{0,0,1,1,2}}, matrix{{0,1,0,0},{0,0,1,1}}}--write the supports here as matrices separated by a comma
M = sparseMonodromy(A,Solver=>M2)
monodromyLoop(M,100)
printCycleTally M
-- I believe the Galois group is S_2 not S_3. This is because the solutions contain two complex solutions
-- (one original and a complex conjugate of both the x and y values) and a single real solution. You can permute the complex conjugate
-- with the original but the real solution cannot be permuted.


-----------------
----Problem 2----
-----------------
--Consider the new polynomial system 'G' obtained by replacing 'x' by 'x^2' in 'F'.
----a) Show that 'G' has twice as many solutions as 'F'. (Note: the solutions to 'G' can be computed
------by computing solutions to 'F' and taking the square root of the x-coordinate)
----
----b) Draw the supports for 'G', how do they differ from the supports of 'F'?
----
----c) Show that the Galois group is the wreath product (S_3 wr S_2) (3 blocks with 2 elements).
----
----Note: Polynomial systems like 'G' that are really polynomials in different variables (here 'x^2' and 'y') are
------called "lacunary". The Galois group of lacunary systems is always imprimitive (there are blocks).
------
------The blocks correspond to solutions to 'F' (giving 3 blocks). The elements of the blocks correspond to the 
------2 solutions of 'G' obtained by taking the square root of the x-coordinate of a solution of 'F'.
G = apply(F,f->sub(f,x=>x^2))

B = {}--enter new supports here
N = sparseMonodromy(B,Solver=>M2)
monodromyLoop(N,30)
printCycleTally N



-----------------
----Problem 3----
-----------------
--Consider the polynomial 'F'. Show there are 3 solutions and the Galois group corresponding
----to the supports is S_3.
restart
loadPackage("Monodromy")
R = CC[x]
F = {1+3*x+5*x^2+7*x^3}



-----------------
----Problem 4----
-----------------
--Consider the polynomial system 'G' obtained by adding another polynomial to the system 'F' above.
----a) Show there are 6 solutions to 'G'. (Note: the solutions to 'G' can be computed by computing
------solutions to the first polynomial and substituting those x-coordinates into the second polynomial)
----
----b) Draw the supports of the polynomials in 'G'.
----
----c) Show the Galois group corresponding to the supports is equal to (S_3 wr S_2) (again 3 blocks with 2 elements).
----
----Note: Systems like 'G' that have a subsystem that can be solved are called "triangular". The Galois
------group of triangular systems is always imprimitive (there are blocks).
------
------The blocks correspond to the solutions of 'F' (giving 3 blocks) and the elements of the blocks 
------correspond to the 2 solutions of the second polynomial after substituting a solution from the first.
restart
loadPackage("Monodromy")
R = CC[x,y]
F = {1+3*x+5*x^2+7*x^3, 2-5*x*y+4*y^2}


