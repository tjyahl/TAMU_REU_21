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





--2 Run the adjacentdensematrix command.
--AdjacentDensePeriodicMatrix(2,3); -- this a graph as before but now with 3 vertices in the fundamental domain and 2 actions
--Problem: Can you draw the periodic graph this is an operator over? 
--hint: add another vertex to the 2,2 dense periodic graph fundamental domain, leave actions the same.
load "../functions.m2"
M=AdjacentDensePeriodicMatrix(2,2)
M_0
M_1
F=M_0 - id_(M_1^2)*x_1*x_2*z --produce the matrix and subtract z from the diagonal
f= det F; -- get the determinant
--for i from 1 to 21 do f = sub(f,e_i=>random(ZZ))
--f
M
L = exponents f;
p = length L - 1
L 
l = new List
for i from 0 to p do l = append(l,{L_i_0-L_i_2,L_i_1-L_i_3,L_i_4})
l
l = unique l;
length l
volume convexHull transpose matrix l
-- by hand these points give 18 but this gives 27
--Problem: What does the newton polytope of the laplace beltrami operator look like? What is its volume?

-- Bonus Problem: Can you come up with a general volume for the newton polytope of the operator given by AdjacentDensePeriodicMatrix(2,n)? 




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

--Problem 4 code:
load("functions.m2")
R = QQ[x_1,x_2,y_1,y_2,a,b,c,d,e,f,z]
I = ideal(x_1*y_1 - 1,x_2*y_2 - 1)
operator = matrix{{e+f+a-z,-e*y_2-f*y_1-a,0},{-a-f*x_1-e*x_2,a+b+e+c+d+f-z,-b-d*y_2+c*y_1},{0,-b-c*x_1-d*x_2,b+c+d-z}}
g = det operator
h = exponents g
l = new List
for i from 0 to length h - 1 do l = append(l,{h_i_0-h_i_2,h_i_1-h_i_3,h_i_10})
length l
l = unique l
length l
volume convexHull transpose matrix l
-- this calculation gives 14/3 when it should give 5 (reason is that a single point is missing from the 
-- supports but I am not sure why it is missing)
--If you want more:


--4.
--Background: If you did 2, the family you were looking at were "maximal abelian coverings" of graphs with 2 vertices
--and n edges between them (starting at n=3, this is graphene).
--Graphene is part of another family however, we will call this the graphene-dice family as the first two memeber
--are the graphene and the dice lattice, these also are related by their topological crystals

--See the image of the dice lattice, remark dice lattice has 3 vertex fundamental domain
--so the matrix representation should be 3 by 3.

--Problem: Calculate the volume of the newton polytope of a laplace Beltrami operator over the dice lattice?
--Hint: How does the Newton polytope of the dice crystal relate to that of the graphene?
--Problem: What might the volume of the nth memeber of the larger family mentioned above be?

-- the vloume of the newton polytope for the dice lattice is 5

-- I am not sure what the volume of the nth member of the graphene-dice family should have. It is hard
-- to make any conclusions with only two data points. Although, it may be the sum of 2 to n.
-- graphene is 2 and the dice lattice is 2 + 3.



--5.
-- Investigate the K4 crystal. I haven't actually done this yet, and constructing it is a little hard 
--with not much resources out there
--So I will go ahead and give you the matrix representation.
--Next time we meet I will try to have written some code to generate what the fundamental domain looks like
-- (its in 3d space with 4 vertices)
--you might want to eliminate the y_i variables, see how I did this in graphenesym.m2
load("functions.m2")
k4cry = AdjacentDensePeriodicMatrix(3,2)
R = QQ[x_1,x_2,x_3,y_1,y_2,y_3,z,e_1,e_2,e_3,e_4,e_5,e_6]
k4cry = matrix{{e_1+e_2+e_3,-e_1,-e_2,-e_3},{-e_1,e_1+e_4+e_5,-e_4*x_1,-e_5*x_2},{-e_2,-e_4*y_1, e_2+e_4+e_6, -e_6*y_3},{-e_3,-e_5*y_2,-e_6*x_3 ,e_3+e_5+e_6}}
--do not forget we are looking at characteristic polynomials, so you will want to subtract I*z from k4cry
--Also do not forget to set x_i and y_i as inverses
inverseideal = ideal(x_1*y_1 -1, x_2*y_2-1,x_3*y_3 -1)
--M = k4cry*x_1*x_2*x_3
--Calculate the volume of the newton polytope of the characteristic polynomial:
-- First specialize the edges to some nonzero values and then take the newton polytope of this and get the volume
-- This polytope is 4 dimensional, so you will definitely want to use software for this

M = k4cry - matrix{{z,0,0,0},{0,z,0,0},{0,0,z,0},{0,0,0,z}}
f = det M
--sub(f,z=>0)
--for i from 1 to 6 do f = sub(f,e_i=>random(ZZ))
E = exponents f
L =  new List
length E
for i from 0 to length E - 1 do L = append(L,{E_i_0-E_i_3,E_i_1-E_i_4,E_i_2-E_i_5,E_i_6})
L = unique L
K = transpose matrix L
CH = convexHull K
volume CH -- outputs 6 for e_i = 1 and 16/3, 5, or 6 for random edge weights. Matt has 4 as answer.
-- interesting to note the ratio of this to the given answer is the same as in problem 2 (27 but should be 18)
--only way to visualize this would be to look at the 3 dimensional base of the polytope,
-- that is the newton polytope of the characteristic polynomial (when z= 0)



--6.
--Problem : Calculate the volume of the newton polytope of a laplace beltrami operator on the kagome lattice.
--see images, I put a fundamental domain image in the directory, you can also just search for a full image
--Hint: should only be 2 actions (as kagome lattice is embedded in K2)







--Problem if you want start thinking about the actual research question:
--7.
--Can you give any conjectures on what might be a minimal subgraph of the periodic graph that
--AdjacentDensePeriodicMatrix(2,3) is an operator over?
--what about AdjacentDensePeriodicMatrix(2,n)?

--Recall by minimal subgraph, we mean one that has as many critical points as the operator over the dense periodic graph
--but that doesn't contain any smaller subgraph with the same number of critical points.
