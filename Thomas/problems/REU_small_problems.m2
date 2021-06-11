--An ongoing list of small problems to examine for understanding.

--------------------------------------------------------------------------
--Problem 1
----Show that the wreath product (S_2 wr S_2) has 8 elements and determine
----the different cycle types of elements.
----
----(This is one to do by hand)

----------------
----SOLUTION----
----------------

--There are 8 elements of the group since there are 2 blocks with 2 elements..
----In each of the 2 blocks, there are 2! possibilities, switch the elements or dont.
----Then you can decide to switch the blocks or not.
----Together, there are (2!)^2*2=8 possibilities

----The cycle types are (1,1,1,1),(1,1,2),(2,2),(4)

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
F = {2+3*y^2-x*y+7*x^2+11*x^2*y^2,-1+3*y^2-6*x*y+11*x^2+13*x^2*y^2}--type your system here
solveSystem F

----------------
----SOLUTION----
----------------

--Notice there are 8 solutions! (You might not ALWAYS get 8 solutions, but almost always).
----The supports look like triangles with length and width = 2. (Notice also both triangles
----have all 3 points on the base filled in, but not on along the vertical edge)




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
loadPackage("Monodromy",FileName=>"../code/Monodromy.m2",Reload=>true)
M = sparseMonodromy(A,Solver=>M2)
M#family
monodromyLoop(M,10)
printCycleTally M
printCycleTypes M


----------------
----SOLUTION----
----------------

--Really just check that the above code agrees with the wreath product S_2 wr S_2.
----The blocks will vary since the computer might not list the solutions the same every time.
----You can determine which elements make a block since a simple transposition like (1,2) can only occur
----in a block! I see transpositions (1,2) and (3,4), so my blocks are the first and second solutions and the
----third and fourth solutions.



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


----------------
----SOLUTION----
----------------

--You can generate as many monodromy loops as you want, you'll see you'll only ever get 4 permutations!
----The reason for this is a bit difficult to see at first. The supports are again triangles as before,
----but now the middle point of the base is removed. This shows that there's more going on than just
----the shape of the supports. 

----Here, the polynomials are polynomials in the variables x^2 and y^2 (or 
----x_1^2 and x_2^2). This is important, we'll talk about this sometime.


