--See graphene picture for labeling to see what is going on
--Here we look at a laplace-beltrami opertor over the graphene
--Laplace-beltrami operator assigns edges weights
--In periodic graph, after a floquet transform, this operators are 
--determined fully by the vertices of a chosen fundamental domain, u and v in the image

--The graphene is embedded into Z^2 so we have two free abelian actions
--Our free actions have inverses so we have Laurent Polynomials
--That K[x,x^{-1}] would be one variable Laurent Polynomials over K
--We are interested in the spectrum, so we have an additional variable


clearAll
needsPackage "NumericalAlgebraicGeometry"


R = CC[x_1,x_2,z,y_1,y_2,a,b,c] 
--We have 3 variables defined x_1, x_2 , and z, the y_i will be the inverse
---of x_i
--a,b,c are edge weights
--To establish this inverse behavior we need to define a quotient ring
I = ideal(x_1*y_1 - 1, x_2*y_2 - 1)

Q = R/I

operator = matrix{{a+b+c-z, -a-b*x_1-c*x_2},{-a-b*y_1-c*y_2, a+b+c-z}}

--build F, system of critical point equations    
F = {};    
F = append(F, sub(det(operator,Strategy => Cofactor),Q));
for i from 1 to 2 do (
    tempI = ideal(x_i^2 * F_0);
    F = append(F, y_i*(diff(x_i, tempI_0) - 2*x_i * F_0));
);
edges = {0};
--play with edge specializations
for i from 1 to 3 do( --add random aspect maybe
    edges = append(edges, random (1, 1321));--no particular reason its 11 just feel thats probably enough
 );

--Rs ring where we drop inverses 
Rs = CC[t_1,t_2,l,e_1,e_2,e_3]
--RsEd ring where we can fix the edges 
RsEd =  CC[s_1,s_2,w]

--specialized system, will send the system F to one with fixed edges
specMap = map(Rs,Q,{t_1,t_2,l,t_1,t_2,e_1,e_2,e_3})
specMapEd = map(RsEd,Rs,{s_1,s_2,w,edges_1,edges_2,edges_3})

--Fp will just allow us to ignore the inverse variables
Fp = {F_0*x_1*x_2,F_1*x_1*x_2,F_2*x_1*x_2}
--we removes y_i variables of F that value can be viewed in CC[x_1,x_2,z,a,b,c]
Fs = apply(Fp, f -> specMap f)
FsEd = apply(Fs, f -> specMapEd(f))

--numerical solve with fixed edges
ss = solveSystem FsEd
#ss
--want to remove zero solutions
--first remove solutions 0 in the first
sift = {}
for i from 0 to #(ss)-1 do(
    curp = coordinates ss_i;
    if (round(realPart(curp_0*conjugate(curp_0))*10000000000) >=  1) then (
	sift = append(sift, curp)
	)
    )
--then remove solutions 0 in the second
sift
#sift

--for generic choices, the solve is missing some solutions

--try to solve symbolically

--When the discriminant is 0 then singular points exist. We homogenize to look at the discriminant

--rebuild as ring of coefficients extended, just going to resolve solve directly

RH = QQ[r_1,r_2,r_3][h_1,h_2,p,o] 
fv = (r_1*o+r_2*o+r_3*o-p)*(r_1*o+r_2*o+r_3*o-p)*h_1*h_1*h_2*h_2 - ((-r_1*o*(h_1*h_2)-r_2*h_1*(h_1*h_2)-r_3*h_2*(h_1*h_2))*(-r_1*o*h_1*h_2-r_2*h_2*o^2-r_3*h_1*o^2))
--fv = homogenize(fv,o,{1,1,1,1})-- doesnt want to work? will just homogenize by hand

fv = (ideal fv)_0


discriminant fv

--gives 0... makes sense because all choices of coefficients give a degenerate system


--... okay, again, use cleaning variables
clearAll
RH = ZZ/8191[a,b,c][x,y,z] 
fv = (a+b+c-z)*(a+b+c-z)*x^2*y^2 - ((-a*(x*y)-b*x*(x*y)-c*y*(x*y))*(-a*x*y-b*y-c*x))

ifv = ideal(fv, (x*diff(x, fv) - 2*fv), y*diff(y,fv) - 2*fv)

--If a,b,c non zero, then we have for the derivatives that - a*x^2*y - c*x^2 + c*y^2 + a*y = 0 and - a*x*y^2 + b*x^2 - b*y^2 + a*x = 0 

--can redefine ideal

ifv =  ideal(fv, - a*x^2*y - c*x^2 + c*y^2 + a*y , -a*x*y^2 + b*x^2 - b*y^2 + a*x)


-- we can fully solve this... but I am unsure how to do this symbolically in macaulay2 without it being more of a hassal than writing it out
