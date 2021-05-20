--The code below generates a polynomial system of the given supports,
----finds a univariate polynomial with the same Galois group, and
----uses the Frobenius method to compute cycle types of elements
----of the Galois group.
restart
loadPackage("Polyhedra",Reload=>true)

R = ZZ[x,y] --make a polynomial ring in 2 variables

supports = {matrix{{0,0,1,2},{0,1,0,1}},matrix{{0,0,0,1},{0,1,2,0}}} --make a list of supports for each polynomial (written as a matrix for each)
P1 = convexHull(supports#0)
P2 = convexHull(supports#1)
numSols = volume(P1+P2) - volume(P1) - volume(P2) --mixed volume = number of solutions

F = for A in supports list (
    	terms := for j from 0 to numcols(A)-1 list (
	    random(-30,30)*x^(A_(0,j))*y^(A_(1,j))
	    );
	sum(terms)
	) --makes our polynomials
    
f = (gens eliminate(ideal F,{y}))_(0,0) --univariate polynomial with same Galois group

p = 31
S = (ZZ/p)[t]
phi = map(S,R,{t,0})

fReduced = phi(f) --reducing mod p
factors = apply(toList factor fReduced,toList) --list of factors and their respective powers

if (first degree fReduced == numSols) and all(factors,p-> last p === 1) then (
    delete(0,sort factors/first@@degree@@first)
    ) --outputs the cycle type for this prime









--This code applies the Frobeius method (modified slightly from above) to various systems 
----of a given support for different primes p.
restart
n = 2 --number of variables
N = 5 --number of primes to use
q = 10 --starting prime (smallest prime larger than q)
k = 25 --number of systems for each prime

supports = n:(matrix toList(n:{0})|matrix table(n,n+3,(i,j)->random(0,3))) --same (random) support for each polynomial

loadPackage("Polyhedra",Reload=>true)
P = convexHull(supports#0)
numSols = n!*volume(P) --mixed volume (has an nice formula when supports are all the same)

primesList = new MutableList from N:nextPrime(q)
scan(N-1,i->primesList#(i+1) = nextPrime(primesList#i+1))
primesList = toList primesList --list of primes to use

L = tally flatten for p in primesList list (
    for l from 1 to k list (
	R := ZZ/p[x_0..x_(n-1)];
    	
    	F := apply(supports,A->sum(numcols(A),j->random(ZZ/p)*product(n,i->x_i^(A_(i,j)))));
    	
	I := eliminate(ideal F,toList(x_1..x_(n-1)));
    	f := if (numgens I > 0) then (gens I)_(0,0) else 1_R;
    	
    	factors := (toList factor f)/toList;
    	
    	if (first degree f == numSols) and all(factors,p-> last p === 1) then (
    	    delete(0,sort factors/first@@degree@@first)
    	    ) else (
	    {}
	    )
    	)
    ) --list of cycle types that showed up using the Frobenius method

T = apply(sort keys L,l->(l,L#l,sub(N*k/(product(l)*product(values tally l,i->i!)),RR))) --list of cycle types, the number of them that appeared, the number of them that should of appeared.


