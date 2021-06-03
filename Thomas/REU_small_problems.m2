--An ongoing list of small problems to examine for understanding.

--------------------------------------------------------------------------
--Problem 1
----Show that the wreath product (S_2 wr S_2) has 8 elements and determine
----the different cycle types of elements.
----
----(This is one to do by hand)



--------------------------------------------------------------------------
--Problem 2
----The matrices 'A' below represent supports for a polynomial system.
----The columns of the first matrix are the exponents for monomials in
----the first polynomial and the columns of the second matrix are the 
----exponents for monomials in the second polynomial.
----
----Draw the supports to see what they look like. Create a random system 
----with these supports and determine the number of solutions using 'solveSystem'.
loadPackage("NumericalAlgebraicGeometry",Reload=>true)
A = {matrix{{0,1,2,0},{0,0,0,2}},matrix{{0,1,2,0},{0,0,0,2}}}

R = CC[x,y]
F = {}--type your system here
solveSystem F



--------------------------------------------------------------------------
--Problem 3
----The Galois group corresponding to the supports 'A' above is
----the wreath product (S_2 wr S_2). (Note: this is because they 
----are really polynomials in the variables x and y^2, or x_1 and x_2^2)
----
----Using 'A' from above, use 'monodromyLoop' to generate elements
----of the Galois group (may need more than 10). Check that you
----can generate 8 distinct elements and that their cycle types are
----the cycle types of elements in the wreath product (S_2 wr S_2).
----
----You might try to see if you can determine which sets of numbers
----form the blocks of the wreath product.
----(remember in (S_2 wr S_2) there are 2 blocks of size 2)
loadPackage("Monodromy",Reload=>true)
M = sparseMonodromy A
M#family
monodromyLoop(M,10)
printCycleTally M
printCycleTypes M



--------------------------------------------------------------------------
--Problem 4
----The supports 'A' below are the supports used above but with one point removed.
----Use 'monodromyLoop' as above to show you can only generate 4 elements
----of the Galois group. 
----
----It can be shown that this Galois group is the product (S_2 x S_2), which
----is smaller than the wreath product (S_2 wr S_2). Can you determine why the
----Galois group is smaller in this case?
A = {matrix{{0,2,0},{0,0,2}},matrix{{0,2,0},{0,0,2}}}
M = sparseMonodromy A
monodromyLoop(M,10)
printCycleTally M
printCycleTypes M
