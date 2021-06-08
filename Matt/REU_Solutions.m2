--Problems for REU students

--I would recommend at least doing 1,2, 3, (bonus problems only if you want to)
-- feel free to look at 4 and 5, do 6 if you really want to.
--7 if you want to start looking at the actual research problem


--I have not done 4 yet, but am interested in it.
--Also feel free to check other properties of these examples or look at specialized versions.

--Note:
--Do not worry about definitions like "coverings", "topological crystals", or "quotient graphs", if you are
--interested in this project we will go over all that later if necessary.
--These examples are just to get you comfortable working with these objects.

--Finally: If stuck, as in no idea what to do next, for more than a few minutes, please feel free to email me. 
--I probably neglected to give some vital information.

--refer to my slides, Frank's slides, or Do Kuchment Sottile paper for sources 



--1
--Problem:  Write out the matrix representation for the laplace beltrami operator over the graphene by hand.
--Check to make sure this is similar to the code that was provided (it's okay if x and x^{-1} are swapped)
--Problem: Draw the newton polytope and calculate the volume by hand.
--(I have already given the polytope, I just want you to verify you understand why that is the correct newton polytope)\

------------------------------------------------------------------------------------

--Solution: See symbolic code and images



--Base of polytope has area 3, height is 2. Then the volume of the polytope is 2. (b*h/dimension)




------------------------------------------------------------------------------------




--2 Run the adjacentdensematrix command.
--AdjacentDensePeriodicMatrix(2,3); -- this a graph as before but now with 3 vertices in the fundamental domain and 2 actions
--Problem: Can you draw the periodic graph this is an operator over? 
--hint: add another vertex to the 2,2 dense periodic graph fundamental domain, leave actions the same.


--Problem: What does the newton polytope of the laplace beltrami operator look like? What is its volume?

-- Bonus Problem: Can you come up with a general volume for the newton polytope of the operator given by AdjacentDensePeriodicMatrix(2,n)? 

------------------------------------------------------------------------------------

--Solution: See crude drawing adjecentdense23



--Solutions: Is a scaled version of the adjacentdenseperiodic 2,2 polytope. Volume: base has area 18, height 3, so 18*3/3 = 18. (same formular as in 1).


--Bonus problem: As we add more vertices, we have the height = #vertices and the base = (#vertices)^2 * 2, we are always in 3 dimensions. Thus the area is given by 2*(n)^3/3

------------------------------------------------------------------------------------



--3.
--Problem: Find the matrix represenation of the laplace beltrami operator of the diamond crystal. 
--Note: (this is a 3D crystal so should have 3 actions)
-- see images the pink vertices are the vertices of the fundamental domain, so the matrix should be 2 by 2
--feel free to search the full image online search something like "diamond crystal periodic graph"


--Side note and questions:
--This might remind you of the graphene. These two operators are very closely related.
--All periodic graphs are also called topological crystals
--Topological crystals have quotient graphs, and they are called covering graphs of these quotient graphs
--the diamond crystal has the same quotient graph as graphene, but with one additional edge

--Bonus Problems
--If diamond and graphene were part of a family, where graphene is the first member and diamond is the second
--Bonus problems 1:  Can you come up with the closed form of the matrix representation of a Laplace Beltrami operator over this family?
--Bonus problems 2:  Can you come up with a closed form for the volume of the newton polytope of laplace Beltrami operators over this family?

--Hint: You may want to use software to guess a closed form

------------------------------------------------------------------------------------

--Solution: 
R = QQ[x_1,x_2,x_3,y_1,y_2,y_3,z,a,b,c,d]
diamondcrystal = matrix{{a+b+c+d,-(d + a*y_1 + b*y_2 + c*y_3)},{-(d + a*x_1 + b*x_2 + c*x_3), a+b+c+d}}

diamondcrystalcharpoly = matrix{{a+b+c+d-z,-(d + a*y_1 + b*y_2 + c*y_3)},{-(d + a*x_1 + b*x_2 + c*x_3), a+b+c+d-z}}


-- Bonus Solution 1: 
--Let the 2nd member be the graphene, 3rd member be diamond, then the nth member is given by:
-- graphenediamondn = matrix{{a_1+a_2+...+ a_n + a_(n+1),-(a_1*y_1+a_2*y_2+...+ a_n*y_n + a_(n+1)},{-(a_1*x_1+a_2*x_2+...+ a_n*x_n + a_(n+1)), a_1+a_2+...+ a_n + a_(n+1)}}
-- a_1, \dots, a_n , a_(n+1) are edge weights

--Bonus Solution  2: 

--We know the volume of the graphene is 2.

-- Lets calculate diamond
--We will use modulated version, that is we multiply everything by x_1*x_2*x_3
-- you could have also done this in sage or anything else with many other methods
-- We will just set all edges to weight 1
load "functions.m2"
R = QQ[x_1,x_2,z]
graphenecharpoly = matrix{{x_1*x_2*(3-z),-(x_1*x_2 + x_2 + x_1)},{-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(3-z)}}    
F = sub(det(graphenecharpoly,Strategy => Cofactor),R); --get char polynomial
FsN = convexHull (transpose matrix ( exponents(F)));
3*2*volume FsN --normalized volume

R = QQ[x_1,x_2,x_3,z]
diamondcrystalcharpoly = matrix{{x_1*x_2*x_3*(4-z),-(x_1*x_2*x_3*1 + x_2*x_3 + x_1*x_3 + x_1*x_2)},{-x_1*x_2*x_3*(1 + 1*x_1 + 1*x_2 + 1*x_3), x_1*x_2*x_3*(4-z)}}    
F = sub(det(diamondcrystalcharpoly,Strategy => Cofactor),R); --get char polynomial
FsN = convexHull (transpose matrix ( exponents(F)));
4*3*2*volume FsN

R = QQ[x_1,x_2,x_3,x_4,z]
graphenediamondcharpoly4 = matrix{{x_1*x_2*x_3*x_4*(5-z),-(x_1*x_2*x_3*x_4*1 + x_2*x_3*x_4 + x_1*x_3*x_4 + x_1*x_2*x_4 + x_1*x_2*x_3)},{-x_1*x_2*x_3*x_4*(1 + 1*x_1 + 1*x_2 + 1*x_3 + x_4), x_1*x_2*x_3*x_4*(5-z)}}    
F = sub(det(graphenediamondcharpoly4,Strategy => Cofactor),R); --get char polynomial
FsN = convexHull (transpose matrix ( exponents(F)));
5*4*3*2*volume FsN

--we get 12, 40, 140, ...
--plugging the values into oeis we see this is probably twice the central binomial coefficent.
--solving this in closed form we get 2*(2(n) choose n) for the volume of the newton polytope of the characteristic polynomial of the nth member
------------------------------------------------------------------------------------




--If you want more:


--4.
--Background: If you did 3 (CORRECTION from 2), the family you were looking at were "maximal abelian coverings" of graphs with 2 vertices
--and n edges between them (starting at n=3, this is graphene).
--Graphene is part of another family however, we will call this the graphene-dice family as the first two memeber
--are the graphene and the dice lattice, these also are related by their topological crystals

--See the image of the dice lattice, remark dice lattice has 3 vertex fundamental domain
--so the matrix representation should be 3 by 3.

--Problem: Calculate the volume of the newton polytope of a laplace Beltrami operator over the dice lattice?
--Hint: How does the Newton polytope of the dice crystal relate to that of the graphene?
--Problem: What might the volume of the nth memeber of the larger family mentioned above be?




--------------------------------------------------------------------------------------

--Solution:
load "functions.m2"
R = QQ[x_1,x_2,z]
dicecharpoly = matrix{{x_1*x_2*(3-z),-(x_1*x_2 + x_2 + x_1),0},{-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(6-z),-(x_1*x_2 + x_2 + x_1)}, {0, -x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(3-z)}}    
F = sub(det(dicecharpoly,Strategy => Cofactor),R); --get char polynomial
FsN = convexHull (transpose matrix ( exponents(F)));
3*2*volume FsN --normalized volume


load "functions.m2"
R = QQ[x_1,x_2,z]
dicecharpoly3 = matrix{{x_1*x_2*(3-z),-(x_1*x_2 + x_2 + x_1),0,0},{-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(6-z),-(x_1*x_2 + x_2 + x_1),0}, {0, -x_1*x_2*(1 + 1*x_1 + 1*x_2),x_1*x_2*(6-z),-(x_1*x_2 + x_2 + x_1)},{0,0,-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(3-z) }}    
F = sub(det(dicecharpoly3,Strategy => Cofactor),R); --get char polynomial
F
FsN = convexHull (transpose matrix ( exponents(F)));
3*2*volume FsN --normalized volume


load "functions.m2"
R = QQ[x_1,x_2,z]
dicecharpoly4 = matrix{{x_1*x_2*(3-z),-(x_1*x_2 + x_2 + x_1),0,0,0},{-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(6-z),-(x_1*x_2 + x_2 + x_1),0,0}, {0, -x_1*x_2*(1 + 1*x_1 + 1*x_2),x_1*x_2*(6-z),-(x_1*x_2 + x_2 + x_1),0},{0,0,-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(6-z),-(x_1*x_2 + x_2 + x_1)},{0,0,0,-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(3-z)}}    
F = sub(det(dicecharpoly4,Strategy => Cofactor),R); --get char polynomial
FsN = convexHull (transpose matrix ( exponents(F)));
3*2*volume FsN --normalized volume

load "functions.m2"
R = QQ[x_1,x_2,z]
dicecharpoly5 = matrix{{x_1*x_2*(3-z),-(x_1*x_2 + x_2 + x_1),0,0,0,0},{-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(6-z),-(x_1*x_2 + x_2 + x_1),0,0,0}, {0, -x_1*x_2*(1 + 1*x_1 + 1*x_2),x_1*x_2*(6-z),-(x_1*x_2 + x_2 + x_1),0,0},{0,0,-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(6-z),-(x_1*x_2 + x_2 + x_1),0},{0,0,0,-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(6-z),-(x_1*x_2 + x_2 + x_1)},{0,0,0,0,-x_1*x_2*(1 + 1*x_1 + 1*x_2), x_1*x_2*(3-z)}}    
F = sub(det(dicecharpoly5,Strategy => Cofactor),R); --get char polynomial
FsN = convexHull (transpose matrix ( exponents(F)));
3*2*volume FsN --normalized volume

-- This is not a known sequence.
-- my bad
-- 12, 30, 96, 168, 324

--Can we figure out what is going on though?
-- clearly, dice lattice is a graphene polytope with its base now having height 1, thus having volume 3+2
-- in the next shape, we have a graphene with all the base points scaled by 2, thus we multiply the base by 4, our height is now 4 as well so multiply by 2, so 8*2 = 16.

--Now we have a general formula in mind:
-- odd members: the (2n-1)th element has same base as graphene scaled by n, it has a height 2n, thus the volume is 3n^2 * (2n/3) = 2n^3
-- Even members: the (2n)th element is the (2n-1)th element with a base of height 1, thus it is 2n^3 + 3n^2.

--See image of the dice lattice polytope for reference and what I mean when I say a base of height 1.
----------------------------------------------------------------------------------------------





--5.
-- Investigate the K4 crystal. I haven't actually done this yet, and constructing it is a little hard 
--with not much resources out there
--So I will go ahead and give you the matrix representation.
--Next time we meet I will try to have written some code to generate what the fundamental domain looks like
-- (its in 3d space with 4 vertices)
--you might want to eliminate the y_i variables, see how I did this in graphenesym.m2

R = QQ[x_1,x_2,x_3,y_1,y_2,y_3,z,e_1,e_2,e_3,e_4,e_5,e_6]
k4cry = matrix{{e_1+e_2+e_3,-e_1,-e_2,-e_3},{-e_1,e_1+e_4+e_5,-e_4*x_1,-e_5*x_2},{-e_2,-e_4*y_1, e_2+e_4+e_6, -e_6*y_3},{-e_3,-e_5*y_2,-e_6*x_3 ,e_3+e_4+e_6}}
--do not forget we are looking at characteristic polynomials, so you will want to subtract I*z from k4cry
--Also do not forget to set x_i and y_i as inverses
inverseideal = ideal(x_1*y_1 -1, x_2*y_2-1,x_3*y_3 -1)


--Calculate the volume of the newton polytope of the characteristic polynomial:
-- First specialize the edges to some nonzero values and then take the newton polytope of this and get the volume
-- This polytope is 4 dimensional, so you will definitely want to use software for this

--only way to visualize this would be to look at the 3 dimensional base of the polytope,
-- that is the newton polytope of the characteristic polynomial (when z= 0)

------------------------------------------------------------------------------------------------------------------

--Solution:

load "functions.m2"
R = QQ[x_1,x_2,x_3,z]
k4charpoly = matrix{{3*x_1*x_2*x_3,-x_1*x_2*x_3,-x_1*x_2*x_3,-x_1*x_2*x_3},{-x_1*x_2*x_3,3*x_1*x_2*x_3,-x_1*x_2*x_3*x_1,-x_1*x_2*x_3*x_2},{-x_1*x_2*x_3,-x_2*x_3, 3*x_1*x_2*x_3, -x_1*x_2},{-x_1*x_2*x_3,-x_1*x_3,-x_1*x_2*x_3*x_3 ,3*x_1*x_2*x_3}}    
F = sub(det(k4charpoly,Strategy => Cofactor),R); --get char polynomial
FsN = convexHull (transpose matrix ( exponents(F)));
4*3*2*volume FsN --normalized volume


--If anyone actually checked solutions for some generic edge weights let me know. Otherwise I will do this later.

------------------------------------------------------------------------------------------------------------------




--6.
--Problem : Calculate the volume of the newton polytope of a laplace beltrami operator on the kagome lattice.
--see images, I put a fundamental domain image in the directory, you can also just search for a full image
--Hint: should only be 2 actions (as kagome lattice is embedded in K2)


------------------------------------------------------------------------------------------------------------------
--Solution: 


load "functions.m2"
R = QQ[x_1,x_2,z]
dicecharpoly = matrix{{x_1*x_2*(4-z),-(x_1*x_2 +x_1),-(x_1*x_2 +x_2)},{-x_1*x_2*(1 + 1*x_2), x_1*x_2*(4-z),-(x_1*x_2 + x_2*x_2)}, {-(x_1*x_2+ x_1*x_1*x_2), -(x_1*x_2 + 1*x_1*x_1), x_1*x_2*(4-z)}}    
F = sub(det(dicecharpoly,Strategy => Cofactor),R); --get char polynomial
FsN = convexHull (transpose matrix ( exponents(F)));
3*2*volume FsN --normalized volume

--Kagome lattice has the same polytope as the dice lattice.


------------------------------------------------------------------------------------------------------------------



--Problem if you want start thinking about the actual research question:
--7.
--Can you give any conjectures on what might be a minimal subgraph of the periodic graph that
--AdjacentDensePeriodicMatrix(2,3) is an operator over?
--what about AdjacentDensePeriodicMatrix(2,n)?

--Recall by minimal subgraph, we mean one that has as many critical points as the operator over the dense periodic graph
--but that doesn't contain any smaller subgraph with the same number of critical points.


-------------------------------------------------------------------------------------------------------------------

--Solution:

--I do not know, let me know if one of you have done this. I have some conjectures but that is it.

--I would guess removing only edges u to x_i u for each i and each u will result in one of the minimal subgraphs.






-----------------------------------------------------------------------------------------------------------------------