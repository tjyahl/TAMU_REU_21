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


R = QQ[x_1,x_2,z,y_1,y_2,a,b,c] 
--We have 3 variables defined x_1, x_2 , and z, the y_i will be the inverse
---of x_i
--a,b,c are edge weights
--To establish this inverse behavior we need to define a quotient ring
I = ideal(x_1*y_1 - 1, x_2*y_2 - 1)

--Q = R/I

operator = matrix{{(a+b+c-z)*x_1*x_2, (-a-b*x_1-c*x_2)*x_1*x_2},{-a*x_1*x_2-b*x_2-c*x_1, (a+b+c-z)*x_1*x_2}}

--build F, system of critical point equations    
F = {};    
F = append(F, sub(det(operator,Strategy => Cofactor),R));
for i from 1 to 2 do (
    tempI = ideal(x_i * F_0);
    F = append(F, y_1*(diff(x_i, tempI_0) - F_0));
);

J = ideal (F) + I
dim J
degree J 
`
