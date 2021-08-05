--Code to take a small loop around discriminant.
restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2",Reload=>true)

A = {matrix{{0,0,1,2,2},{0,2,1,0,2}},matrix{{0,0,1,2,2},{0,2,1,0,2}}}
M = sparseMonodromy A
monodromyLoop(M,20,Verbosity=>2)
numSolns = #(M#baseSolutions)

F = transpose (M#family).PolyMap
R = ring F;
C = coefficientRing R;

n = #(gens R)
m = #(gens C)

--Create a system with a singular solution.
singSoln = apply(n,i->random(CC));
phi = map(C,R,singSoln);

S = CC[t];
A = transpose (sub(gens kernel transpose last coefficients phi(F),S)*random(CC^(m-n),CC^2));
B = (t*A^{0} + (1-t)*A^{1});
rho = map(S,C,B);

--compute singular value of t by Newton's method.
JacF = rho phi jacobian F;
a = sub(JacF,t=>0)
b = sub(JacF-a,t=>1)
singT = -1/(1+trace(inverse(a)*b))
j = 0;
while (j < 40) and (abs det(a+singT*b) > 1e-8) do (
    j = j+1;
    singT = singT - 1/(1+trace(inverse(a+singT*b)*b))
    );
det(a+singT*b)


--check sing system is sing.
singBase = flatten entries ((map(CC,S,{singT}))(B));
f = map(CC,R,singSoln|singBase);
norm f(F)
det f(jacobian F)

--linear space containing singBase and M#basePoint
B = matrix {(1-t)*(coordinates M#basePoint) + t*singBase}

--might need to play with eps. 
eps = .001
t0 = 1 + eps;
t1 = 1 + eps*exp(2*pi*ii/3);
t2 = 1 + eps*exp(4*pi*ii/3);

newBase = point {flatten entries ((map(CC,S,{t0}))(B))};
pt1 = point {flatten entries ((map(CC,S,{t1}))(B))};
pt2 = point {flatten entries ((map(CC,S,{t2}))(B))};

N = changeBasePoint(M,newBase)
isWellDefined N

N = monodromyLoop(N,{pt1,pt2},Verbosity=>2)







