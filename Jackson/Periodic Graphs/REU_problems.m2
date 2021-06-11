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
load "functions.m2"
M=AdjacentDensePeriodicMatrix(2,3)
F=M_0 - matrix{{z,0,0},{0,z,0},{0,0,z}} --produce the matrix and subtract z from the diagonal
f= det F -- get the determinant
--for i from 1 to 21 do f = sub(f,e_i=>random(ZZ))
--f
L = exponents f
p = length L - 1
L = new MutableList from L
for i from 0 to p do L#i = {L#i_0,L#i_1,L#i_2,L#i_3,L#i_4}
L = new List from L
length L
--length unique L
--for i from 0 to length L - 1 if L_i_2 != 0 print i
l = new List
for i from 0 to p do l = append(l,{L_i_0-L_i_2,L_i_1-L_i_3,L_i_4})
l
l = unique l
length l
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
