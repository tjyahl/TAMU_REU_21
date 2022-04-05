
-- installPackage "Polyhedra";
needsPackage "Polyhedra"
needsPackage "NumericalAlgebraicGeometry"

--literally just map
Specialize = (P,Q,M) -> map(P,Q,M)

--Define pruneList
--We want to reduce the dimension to get rid of a and b for constructing our system of newton polytope
--first list will be a given list of unpruned vertices of equal dimension
--second list is a given set of what indices to prune, the first corresponding to 0 and so on
--basically just a more customizable list drop
--could have just applied drop
--must give p list in reverse order (I should rewrite this to reorder list in descending order regardless of input)
pruneList = method(TypicalValue => List)
pruneList (List,List) := List => (initList, pList) -> (
    local prunedList; local listLength; local vecLength;
    prunedList = new MutableList from initList; --we are modifying the lists of the greater list so we need a Mutable list
    listLength = #initList;
    vecLength = #initList_0;
    for i from 0 to listLength-1 do (	
	for j from 0 to #pList-1 do (
	    prunedList#i = drop(prunedList#i,{pList_j,pList_j});
	    );
	);
    
    return new List from prunedList
    )


-- minkowskiSum that takes a list of polytopes and recurses 
--appearently can just use fold(minkowskiSum, List) aswell.
myMinkowskiSum = method(TypicalValue => Thing)
myMinkowskiSumHelper = method(TypicalValue => Sequence)
 
myMinkowskiSum (List) := (Thing) => (polyList) -> (
    return (myMinkowskiSumHelper(polyList_0, drop(polyList,1)))_0
)
myMinkowskiSumHelper (Thing, List) := (Thing, List) => (curPoly, polyList) -> (
    if (polyList === {}) then return (curPoly,{});
    if (polyList =!= {}) then return myMinkowskiSumHelper (minkowskiSum (curPoly, polyList_0), drop (polyList,1));
    ) 

--the native mixed volume computation is extremely slow.
--myMixedVolume will take a list of polytope
myMixedVolume  = method(TypicalValue => QQ) --function to call
myMixedVolumeRec = method(TypicalValue => Sequence) --recursive helper function
myMixedVolume (List) := QQ => (initList) -> (
    	tracklist = {};
	for i from 1 to #initList do(
	    tracklist = append(tracklist, i);
	    );
    	return (myMixedVolumeRec (initList, { }, tracklist))#0
    )
myMixedVolumeRec (List ,List, List) := (QQ, List)  => (polyList,accountedList, trackList) -> (--accountedList being used as a lazy redudancy reduction, storage is cheap
    local mvolume; local mixedpoly; local couple; local curpoly; local curtrack;
    	if (any(accountedList, n-> n == trackList)) then return (0, accountedList);
	accountedList = append(accountedList, trackList);
	if (#polyList === 1) then (
	    if (isFullDimensional polyList_0 == false) then return (0, accountedList);
	     return (volume polyList_0, accountedList);		 
		 );
	mvolume = 0;
--	mixedpoly = myMinkowskiSum(polyList);
        mixedpoly = fold(minkowskiSum, polyList);
	if (isFullDimensional mixedpoly == false) then return (0, accountedList);
	 mvolume = volume mixedpoly;
	for i from 0 to #polyList-1 do(
	    curpoly = drop(polyList, {i,i});
	    curtrack = drop(trackList, {i,i});
	    couple = myMixedVolumeRec(curpoly, accountedList, curtrack);
	    mvolume = mvolume - couple#0;
--	   <<couple#0;-- was just here for debugging
--	   <<curtrack;
	    accountedList = couple#1;
	    );
	return (mvolume, accountedList)	
    )

--coneList will take a fan and return the list of matrices defining our cones in each dimension 1 to ambDim(curFan)
coneList = method(TypicalValue => List)
coneList (Thing) := List => (curFan) -> (
    local cList; local adim; local curCone; local curConeList;local curIndex; local raytrix;
    raytrix = rays curFan;
    cList = { };
    adim = ambDim curFan;
    for i from 1 to adim do (	
	curConeList = cones(i, curFan);
      	for j from 0 to #curConeList-1 do ( --honestly a super lazy design that no for loop without a do is implemented....
	    	curIndex = (curConeList_j)_0;
	    	curCone =raytrix_{curIndex};
	    for k from 1 to #curConeList_j-1 do (
		if (#curConeList_j <= 1) then break;
--		<< curCone; lines for debugging
		curIndex = (curConeList_j)_k;
		curCone = curCone | raytrix_{curIndex};
--		<< curIndex;
--		<< curCone;
		);
	    	cList = append(cList, curCone);
	    );
	);
        return cList
     )
 
 --getFaceFunc will take a cone and a polynomial it will return the subpolynomial corresponding to the facet the cone is normal to
--will use exponents to decide which indices are chosed
--then will iterate through a list of terms to create the facial polynomial by successively added the indices chosen
--formate the polynomial so that variables with inverses do not appear with their inverses
--Assumes ring to be of the form R[x_1,...x_n,xi_1,...,xi_n,z_1,...z_m] where xi_j is inverse to i_j and z_i are non invertable vars.
--pass to getFaceFunc(thing, thing) if ring is a laurent polynomial (that is m = 0
--pass to getFaceFunc(thing, thing, ZZ) if ring has non invertable variables, this takes an extra variable m to account for cutting
getFaceFunc = method(TypicalValue => Thing)
getFaceFunc (Thing, Thing) := Thing => (inCone, initf) -> (
    	    return getFaceFunc(inCone, initf, 0)
    );

getFaceFunc (Thing, RingElement, ZZ) := RingElement => (inCone, initf, m) -> (
    local curmin; local minList; local curPoly; local expsf; local termsf; local vartocut; local numinvars;
    numinvars = numRows(inCone) - m;
    vartocut = { };
    for i from 1 to numinvars do (
	    vartocut= append(vartocut, 2*numinvars - i);
	);
  --  << numinvars;
    exps = pruneList(exponents initf, vartocut);
    termsf = terms initf;  
    minList = { };
    curPoly = initf-initf;
    for i from 0 to #exps-1 do (
	if (i == 0) then (
	    checkmin = getIProd(inCone, exps_i);--helper function
	    minList = append(minList, checkmin);
     	    curmin = checkmin;
	    continue;
	    );
	checkmin = getIProd(inCone, exps_i);
	curmin = min {curmin, checkmin};
	minList = append(minList, checkmin);
	);
    for i from 0 to #exps-1 do (
	    if minList_i == curmin then curPoly = curPoly + termsf_i;
	);
    return curPoly
     )
 --getIProd takes a cone/matrix and a vector and computes sum of inner product of the transpose of each column and the vector 
 --helper function to getFaceFunc
 getIProd = method(TypicalValue => QQ)
 getIProd (Matrix, Thing) := QQ => (inCone, expvec) -> (
     local cumsum; local colList;
     cumsum = 0;
     for i from 0 to numColumns(inCone)-1 do(
	 cumsum = cumsum + vecIProd((entries transpose inCone_{i})_0,expvec);
	 );
     return cumsum
     )
 --helper to getIProd, just yields the inner product. takes two lists and computes inner product
 vecIProd = method(TypicalValue => QQ)
 vecIProd (List, List) := QQ => (inCone, expvec) -> (
     local cumsum;
     cumsum = 0;
     if (#inCone =!= #expvec) then (
	 << "Vectors not of equal length nerd";
	 return 0;
	 );
     for i from 0 to #inCone-1 do(
	 cumsum = cumsum + (inCone_i * expvec_i);
	 );
     return cumsum
     )

--will take a coeffiecient and polynomial list and convert it to a polynomial in the given ring
--first input is the polynomial, the second is a list of variables needing pruning, the third is the new ring
--Have not actually tested this yet 
toNewRing = method(TypicalValue => RingElement)
toNewRing (RingElement, List, Ring) := RingElement => (inPoly, toPruneList, desRing) -> (
     local coeffs; local exps; local outPoly; local listFormTemp; local zList;
     listFormTemp = listForm (inPoly);
     coeffs = {};
     exps = {};     
     
     for i from 0 to #listFormTemp-1 do (
	   coeffs = append(coeffs, (listFormTemp_i)_0);
	   exps = append(exps, (listFormTemp_i)_1); 
	 );
     exps = pruneList(exps,toPruneList);
     zList = {};
     for i from 0 to numRows(vars(desRing))-1 do(
	 zList = append(zList, 0);
	 );
     outPoly = desRing_zList;
     for i from 0 to #listFormTemp-1 do(
	 outPoly = outPoly + coeffs_i * desRing_(exps_i);
	 );
     return outPoly
     )


 

--will take a list of two integers m and n 
--m will be the number of actions
--n is the number of nodes
--will return a complete saturated fundamental domain
--by this I mean the graph is a graph that has m complete bipartite subgraphs Kn,n and their intersection gives n nodes that 
--in the graph are a complete graph Kn
--further, besides that intersection the bipartite subgraphs do not intersect
--I think this should fully describe the graphs 
--return the matrix the ring and the quotient ring
 
--takes output of saturated graph system and under a built in specialization yields the output
--next step would probably be to allow for customizable specialization, probably add a parameter at the very least for offset of lambda 



--ended up having to redo graph system within the function since ran into trouble with ringelements not realizing they are in the 
--same ring, rather than working around it with futher abstraction just went with combining the routine
--can probably try to break it up later to get it to work properly.

--also differs from graphsystem in that the matrix produced accounts for subtracting z from each diagonal.


--how it is now it takes m as number of actions/dimensions, n as number of nodes/vertices in fundamental 
--domain, k as offset of w in the projection map

--system bottlenecking at map, will assign edge weights before building the matrix to prevent this.


saturatedPolySystem = method(TypicalValue => (Sequence))   
 --instead of list will take int, int,int
  
saturatedPolySystem (ZZ,ZZ,ZZ) := (Sequence) => (m,n,adjZ) -> (
     local edges; local toteEdges; local invFuncs; local toMatrixList, local toRowList; local zList; local outList; local outMatrix;
     local curIndex;
     local R;
     local W;
     local I;
     invFuncs = {};
     toteEdges = m*n*n + n*(n-1)//2;
     local a;
     a = {0};
     --doing this now since the large number of variables was causing a bottleneck on the map function
     for i from 1 to (m*n*n + n*(n-1)//2) do( --add random aspect maybe
	 a = append(a, random (1, i*107));--no particular reason its 11 just feel thats probably enough
      );

     R = QQ[x_1 .. x_m, y_1 .. y_m,z];
     for i from 1 to m do(
    	invFuncs = append(invFuncs, x_i*y_i -1);	 
	 );     
     I = ideal (toSequence invFuncs);
     W = R/I;
     toMatrixList = new MutableList;
     zList = {};
     for i from 0 to n-1 do(
	 zList = append(zList, 0);
	 );     
     --make a double nested list so that we can convert this to a matrix later, make it nxn
     for i from 1 to n do( -- declare like this so that all rows are unique objects
	 toRowList = new MutableList;
	 for j from 1 to n do (
	 toRowList = append(toRowList, R_zList - R_zList);
	 );
	 toMatrixList = append(toMatrixList, toRowList);
	 );       
     --okay so lets iterate through each bipartite graph
     --we will look at all the edges connected to a node one at a time, going through the nodes not part of the complete subgraph
     -- this will require one loop to go through each bipartite graph and another loop to go through each node and another
     --for the n edges on each node
     --it will then require one more loop after these nested loops to account for the edges in the complete graph
     for i from 1 to m do( --going through bipartite partitions
	 curIndex = (i-1)*n*n + 1;
	 for j from 0 to n-1 do( --going through nodes of a particular partition
	     for k from 0 to n-1 do( --going through the edges
		      if (j == k) then (
	 		  toMatrixList#j#j = toMatrixList#j#j +  (2 - x_i - y_i)*a_curIndex;
			  curIndex =curIndex + 1;
	 		  );
		      if (j != k) then (
     	     	         (toMatrixList#j)#j = (toMatrixList#j)#j + a_curIndex;
		         (toMatrixList#k)#k = (toMatrixList#k)#k + a_curIndex;
		         (toMatrixList#k)#j = (toMatrixList#k)#j - x_i*a_curIndex;
		         (toMatrixList#j)#k = (toMatrixList#j)#k - y_i*a_curIndex;		      
		 	  curIndex =curIndex + 1;
			  );
		 );
	     );
	 );          
	 
	  --need to remember to subtract from each diagonal by w
     for j from 0 to n-1 do( --going through nodes of a particular partition
	     for k from 0 to n-1 do( --going through the edges
		      if (j == k) then (
	 		  toMatrixList#j#j = toMatrixList#j#j - z;
	 		  );
		 );
	     );
     -- also need to finally remember to account for the edges in the complete graph of nodes (accounts for current bug)
      for j from 0 to n-1 do( --
	     for k from 0 to n-1 do(  --only want hit edges once, we can do this by only considering increasing pairs
		      if (j < k) then (
	 		 (toMatrixList#j)#j = (toMatrixList#j)#j + a_curIndex;
		         (toMatrixList#k)#k = (toMatrixList#k)#k + a_curIndex;
		         (toMatrixList#k)#j = (toMatrixList#k)#j - a_curIndex;
		         (toMatrixList#j)#k = (toMatrixList#j)#k - a_curIndex;		      
		 	  curIndex =curIndex + 1;
	 		  );
		 );
	     );   
     
     
	 
	 for j from 0 to n-1 do (
	 toMatrixList#j = new List from toMatrixList#j;
	 );    
    
     toMatrixList = new List from toMatrixList;
     outMatrix = matrix(toMatrixList);
        
     
 --    listvals = {m, n, R, W, outMatrix}
     
     
     --where the old new function starts. recombined the two functions since ring was breaking
     --same RingElement types were being treated as distinct and not multiplying together
     --language type casting is super annoying.
     
     
     local f;     
     f = {};    
     f = append(f, sub(det(outMatrix,Strategy => Cofactor),W));
     for i from 1 to m do (
	 tempI = ideal(x_i^n * f_0);
	 f = append(f, y_i^(n-1)*(diff(x_i, tempI_0) - n*x_i^(n-1) * f_0));
	 );
     local J;
     local Rc;
     local Wc;
     local Ic;
--     J = ideal f;
     Rc = QQ[s_1 .. s_m, t_1 .. t_m, w];
     specList = {}; --building list for variable specialization this part of the code can be modified, will probably build function
     	    	    --to specify this later
     for i from 1 to m do( --no idea what was going on here, why was i appending variables to the inv function list...
       	 specList = append(specList, s_i); --supposed to be speclist and order is wrong
	 );
     for i from 1 to m do(
	 specList = append(specList, t_i);
	 );
     specList = append(specList, w + adjZ); --just weighting edges by numbering 1 to 10 doesnt seem to work well
     --doing this at the start instead to combat bottleneck
--     for i from 1 to (m*n*n + n*(n-1)//2) do( --add random aspect maybe
--	 specList = append(specList, random (i*11));--no particular reason its 11 just feel thats probably enough
--	 );
--     specList = append(specList,1);
     local invFuncs; --making inverse functions to convert certain variables to laurant type
     local adjPoly; --adjusting polynomial that allows us to get rid of the inverse variables for creating convex hull in m+1vars
     invFuncs = {};
     adjPoly = 1;
     for i from 1 to m do(
	 adjPoly = (x_i^n) * adjPoly;
	 );
     for i from 1 to m do(
    	invFuncs = append(invFuncs, s_i*t_i -1);	 
	 );     
     Ic = ideal (toSequence invFuncs);
     Wc = Rc/Ic;
     local SpecMap;
     local fs;
     SpecMap = map(Wc,W,specList);
     fs = apply(f, n -> SpecMap (n*adjPoly)); --list with specialized edges
     local toprunelist; --inverse variables need to be pruned in order to be able to take convex hulls
     toprunelist = {};
     for i from 1 to m do (
	 toprunelist = append(toprunelist, 2*m-i);
	 );     
     local mixedV;
     local fsP; local fsN; local sysf; local clF; local F;
     fsP = apply(fs, n -> pruneList ( exponents(n ) , toprunelist)); --pruned list
     fsN = apply(fsP, n -> convexHull (transpose matrix n)); --list of polytope
--     mixedV = myMixedVolume(fsN);
--     stdio << mixedV;
--     stdio << volume (myMinkowskiSum fsN);
--     fsP_0
     F = normalFan fsN_0; --just need a fan of the main function obtained via the determinate 
     clF = coneList(F); --get a list of all the cones
     
     sysf = apply(fs, n-> apply(clF, m-> getFaceFunc(m,n,1))); --this gives us a list of all the facial system 
     local theideals; --generating output that may be of interest, as in the ideals and relevant info
     local dims;
     local degs;
     local gene; --generators of the grobner basis, necessary since dim of objects with 1 as a generator give -1 dim ideal return
     theideals = {};
     dims = {};
     degs = {};
     gene = {};
     for i from 0 to #clF - 1 do (
	 thelist = {};
	 for j from 0 to m do (
	     thelist = append(thelist, (sysf_j)_i);
	     );
	 theideals = append(theideals,ideal thelist );
	 dims = append(dims, dim theideals_i);
	 degs = append(degs, degree theideals_i);
	 gene = append(gene, gens gb theideals_i);
	 );
     
--     stdio << theideals;
--     stdio << dims;
--     stdio << degs;
     
     return (Wc, sysf,clF, theideals, dims, degs, gene, specList,a)
     )
 
 --this one will take a particular list of edgeweights
 
saturatedPolySystem (ZZ,ZZ,ZZ,List) := (Sequence) => (m,n,adjZ,a) -> (
     local edges; local toteEdges; local invFuncs; local toMatrixList, local toRowList; local zList; local outList; local outMatrix;
     local curIndex;
     local R;
     local W;
     local I;
     invFuncs = {};
     toteEdges = m*n*n + n*(n-1)//2;
     
     if (#a != toteEdges+1) then(
    	return (-1,-1)	 
      );
    

     R = QQ[x_1 .. x_m, y_1 .. y_m,z];
     for i from 1 to m do(
    	invFuncs = append(invFuncs, x_i*y_i -1);	 
	 );     
     I = ideal (toSequence invFuncs);
     W = R/I;
     toMatrixList = new MutableList;
     zList = {};
     for i from 0 to n-1 do(
	 zList = append(zList, 0);
	 );     
     --make a double nested list so that we can convert this to a matrix later, make it nxn
     for i from 1 to n do( -- declare like this so that all rows are unique objects
	 toRowList = new MutableList;
	 for j from 1 to n do (
	 toRowList = append(toRowList, R_zList - R_zList);
	 );
	 toMatrixList = append(toMatrixList, toRowList);
	 );       
     --okay so lets iterate through each bipartite graph
     --we will look at all the edges connected to a node one at a time, going through the nodes not part of the complete subgraph
     -- this will require one loop to go through each bipartite graph and another loop to go through each node and another
     --for the n edges on each node
     --it will then require one more loop after these nested loops to account for the edges in the complete graph
     for i from 1 to m do( --going through bipartite partitions
	 curIndex = (i-1)*n*n + 1;
	 for j from 0 to n-1 do( --going through nodes of a particular partition
	     for k from 0 to n-1 do( --going through the edges
		      if (j == k) then (
	 		  toMatrixList#j#j = toMatrixList#j#j +  (2 - x_i - y_i)*a_curIndex;
			  curIndex =curIndex + 1;
	 		  );
		      if (j != k) then (
     	     	         (toMatrixList#j)#j = (toMatrixList#j)#j + a_curIndex;
		         (toMatrixList#k)#k = (toMatrixList#k)#k + a_curIndex;
		         (toMatrixList#k)#j = (toMatrixList#k)#j - x_i*a_curIndex;
		         (toMatrixList#j)#k = (toMatrixList#j)#k - y_i*a_curIndex;		      
		 	  curIndex =curIndex + 1;
			  );
		 );
	     );
	 );          
	 
	  --need to remember to subtract from each diagonal by w
     for j from 0 to n-1 do( --going through nodes of a particular partition
	     for k from 0 to n-1 do( --going through the edges
		      if (j == k) then (
	 		  toMatrixList#j#j = toMatrixList#j#j - z;
	 		  );
		 );
	     );
     -- also need to finally remember to account for the edges in the complete graph of nodes (accounts for current bug)
      for j from 0 to n-1 do( --
	     for k from 0 to n-1 do(  --only want hit edges once, we can do this by only considering increasing pairs
		      if (j < k) then (
	 		 (toMatrixList#j)#j = (toMatrixList#j)#j + a_curIndex;
		         (toMatrixList#k)#k = (toMatrixList#k)#k + a_curIndex;
		         (toMatrixList#k)#j = (toMatrixList#k)#j - a_curIndex;
		         (toMatrixList#j)#k = (toMatrixList#j)#k - a_curIndex;		      
		 	  curIndex =curIndex + 1;
	 		  );
		 );
	     );   
     
     
	 
	 for j from 0 to n-1 do (
	 toMatrixList#j = new List from toMatrixList#j;
	 );    
    
     toMatrixList = new List from toMatrixList;
     outMatrix = matrix(toMatrixList);
        
     
 --    listvals = {m, n, R, W, outMatrix}
     
     
     --where the old new function starts. recombined the two functions since ring was breaking
     --same RingElement types were being treated as distinct and not multiplying together
     --language type casting is super annoying.
     
     
     local f;     
     f = {};    
     f = append(f, sub(det(outMatrix,Strategy => Cofactor),W));
     for i from 1 to m do (
	 tempI = ideal(x_i^n * f_0);
	 f = append(f, y_i^(n-1)*(diff(x_i, tempI_0) - n*x_i^(n-1) * f_0));
	 );
     local J;
     local Rc;
     local Wc;
     local Ic;
--     J = ideal f;
     Rc = QQ[s_1 .. s_m, t_1 .. t_m, w];
     specList = {}; --building list for variable specialization this part of the code can be modified, will probably build function
     	    	    --to specify this later
     for i from 1 to m do( --no idea what was going on here, why was i appending variables to the inv function list...
       	 specList = append(specList, s_i); --supposed to be speclist and order is wrong
	 );
     for i from 1 to m do(
	 specList = append(specList, t_i);
	 );
     specList = append(specList, w + adjZ); --just weighting edges by numbering 1 to 10 doesnt seem to work well
     --doing this at the start instead to combat bottleneck
--     for i from 1 to (m*n*n + n*(n-1)//2) do( --add random aspect maybe
--	 specList = append(specList, random (i*11));--no particular reason its 11 just feel thats probably enough
--	 );
--     specList = append(specList,1);
     local invFuncs; --making inverse functions to convert certain variables to laurant type
     local adjPoly; --adjusting polynomial that allows us to get rid of the inverse variables for creating convex hull in m+1vars
     invFuncs = {};
     adjPoly = 1;
     for i from 1 to m do(
	 adjPoly = (x_i^n) * adjPoly;
	 );
     for i from 1 to m do(
    	invFuncs = append(invFuncs, s_i*t_i -1);	 
	 );     
     Ic = ideal (toSequence invFuncs);
     Wc = Rc/Ic;
     local SpecMap;
     local fs;
     SpecMap = map(Wc,W,specList);
     fs = apply(f, n -> SpecMap (n*adjPoly)); --list with specialized edges
     local toprunelist; --inverse variables need to be pruned in order to be able to take convex hulls
     toprunelist = {};
     for i from 1 to m do (
	 toprunelist = append(toprunelist, 2*m-i);
	 );     
     local mixedV;
     local fsP; local fsN; local sysf; local clF; local F;
     fsP = apply(fs, n -> pruneList ( exponents(n ) , toprunelist)); --pruned list
     fsN = apply(fsP, n -> convexHull (transpose matrix n)); --list of polytope
--     mixedV = myMixedVolume(fsN);
--     stdio << mixedV;
--     stdio << volume (myMinkowskiSum fsN);
--     fsP_0
     F = normalFan fsN_0; --just need a fan of the main function obtained via the determinate 
     clF = coneList(F); --get a list of all the cones
     
     sysf = apply(fs, n-> apply(clF, m-> getFaceFunc(m,n,1))); --this gives us a list of all the facial system 
     local theideals; --generating output that may be of interest, as in the ideals and relevant info
     local dims;
     local degs;
     local gene; --generators of the grobner basis, necessary since dim of objects with 1 as a generator give -1 dim ideal return
     theideals = {};
     dims = {};
     degs = {};
     gene = {};
     adjPoly = 1;
     for i from 1 to m do(
	 adjPoly = (t_i^n) * adjPoly;
	 );
     for i from 0 to #clF - 1 do (
	 thelist = {};
	 for j from 0 to m do (
	     thelist = append(thelist, adjPoly*(sysf_j)_i);
	     );
	 theideals = append(theideals,ideal thelist );
	 dims = append(dims, dim theideals_i);
	 degs = append(degs, degree theideals_i);
	 gene = append(gene, gens gb theideals_i);
	 );
     
--     stdio << theideals;
--     stdio << dims;
--     stdio << degs;
     
     return (outMatrix, Wc, sysf,clF, theideals, dims, degs, gene, specList,a,f)
     )
 
 

--will take a list of edge weights and specialize to zero in some fashion
--our initial function will delete all edges that are generated by a single action
--0 will remove all single action edges (say nonecrossing), 1 removes all crossing edges
vanishEdges = method(TypicalValue => List)
vanishEdges (ZZ,ZZ,List,ZZ) := List => (m,n,aList,mode) -> (
    local vanishList; local listLength;
    vanishList = new MutableList from aList; --we are modifying the lists of the greater list so we need a Mutable list
    listLength = #aList;
    if (mode == 0) then( --mode 0 keep diagonal edges, graph is  connected since we keep complete graph of inital nodes
    	for j from 1 to m do(
	    for i from 1 to n*n do (	
            	if ( (i - 1) %(n+1) == 0) then(
		    vanishList#((j-1)*n*n+i) = 0; 
		    );
		);
    	    );
    );
    if (mode == 1) then( --mode 0 keep diagonal edges, graph is  connected since we keep complete graph of inital nodes
    	for j from 1 to m do( --mode 1 is necessary to keep the edges needed to saturate the polytope
	    for i from 1 to n*n do (	--so any edge specialization needs these edges.
            	if ( (i - 1) %(n+1) != 0) then( --but these seem to remove too much
		    vanishList#((j-1)*n*n+i) = 0; 
		    );
		);
    	    );
    );
    if (mode == 2) then( --mode works but dunno what graph is
    	for i from 1 to m*n*n do (	
            if ( (i - 1) %(n+1) != 0) then(
	    	vanishList#i = 0; 
		);
	    );
    );
    if (mode == 3) then( --mode works but dunno what graph is
    	for i from 1 to m*n*n do (	
            if ( (i - 1) %(n+1) == 0) then(
	    	vanishList#i = 0; 
		);
	    );
    );
    if (mode == 4) then( --mode 4 keeps edges from mode 1 and will keep other edges too
    	for j from 1 to m do( 
	    for i from 1 to n*n do (
            	if ( (i - 1) %(n+1) != 0 and (i - 1) %(n+1) != 1) then(
		    vanishList#((j-1)*n*n+i) = 0; 
		    );
		);
    	    );
    );
    if (mode == 5) then( --keeps mode one edges as well as a single bipartite graph in one direction. Works exactly as desired
    	for j from 1 to m do( 
	    for i from 1 to n*n do (
            	if ( (i - 1) %(n+1) != 0 and j>1) then(
		    vanishList#((j-1)*n*n+i) = 0; 
		    );
		);
    	    );
    );
    if (mode == 6) then( --keeps mode one edges as well as a single bipartite graph in one direction. Works exactly as desired
    	for j from 1 to m do(  --same as 5 but removes edges in central complete graph, also seems to work
	    for i from 1 to n*n do (
            	if ( (i - 1) %(n+1) != 0 and j>1) then(
		    vanishList#((j-1)*n*n+i) = 0; 
		    );
		);
    	    );
	for j from n*n*m+1 to listLength-1 do(
	    vanishList#j = 0;
	    );
    );
    return new List from vanishList
    ) 








--will take number of actions (m) and number of  vertices (m)
--and return matrix representation for a general laplace beltrami operator
-- on an action adjacent dense graph
AdjacentDensePeriodicMatrix = method(TypicalValue => List)
AdjacentDensePeriodicMatrix(ZZ,ZZ) := (Sequence) => (m,n)-> (
     local edges; local toteEdges; local invFuncs; local toMatrixList, local toRowList; local zList; local outList; local outMatrix;
     local curIndex;
     local R;
     local W;
     local I;
     local adjpoly;
     local yirep;
     invFuncs = {};
     toteEdges = m*n*n + n*(n-1)//2;

     R = QQ[x_1 .. x_m, y_1 .. y_m,z, e_1 .. e_(toteEdges)];
     for i from 1 to m do(
    	invFuncs = append(invFuncs, x_i*y_i -1);	 
	 );     
     I = ideal (toSequence invFuncs);
     toMatrixList = new MutableList;
     zList = {};
     for i from 0 to n-1 do(
	 zList = append(zList, 0);
	 );     
     --make a double nested list so that we can convert this to a matrix later, make it nxn
     for i from 1 to n do( -- declare like this so that all rows are unique objects
	 toRowList = new MutableList;
	 for j from 1 to n do (
	 toRowList = append(toRowList, R_zList - R_zList);
	 );
	 toMatrixList = append(toMatrixList, toRowList);
	 );       
     --okay so lets iterate through each bipartite graph
     --we will look at all the edges connected to a node one at a time, going through the nodes not part of the complete subgraph
     -- this will require one loop to go through each bipartite graph and another loop to go through each node and another
     --for the n edges on each node
     --it will then require one more loop after these nested loops to account for the edges in the complete graph
     
     --yirep to represent y_i
     yirep = new MutableList;
     for i from 1 to m do(
	 yirep = append(yirep,1_R);
	 ); 
     --need yirep in this way so can indentify y_i = adjpoly_i
     for i from 0 to m-1 do(
	 for j from 1 to m do(
	     if ((i +1) != j) then(
	          yirep#i = yirep#i* (x_j)_R;    
		  );
	     );
	 );

--for the rest need to still multiply by all terms, use adjpoly
     adjpoly = 1_R;
     for i from 1 to m do(
	 adjpoly = adjpoly*(x_i)_R;
	 ); 
 
     for i from 1 to m do( --going through bipartite partitions
	 curIndex = (i-1)*n*n + 1;
	 for j from 0 to n-1 do( --going through nodes of a particular partition
	     for k from 0 to n-1 do( --going through the edges
		      if (j == k) then (
	 		  toMatrixList#j#j = toMatrixList#j#j +  ((2 - x_i)*adjpoly - yirep#(i-1))*e_curIndex;
			  curIndex =curIndex + 1;
	 		  );
		      if (j != k) then (
     	     	         (toMatrixList#j)#j = (toMatrixList#j)#j + e_curIndex*adjpoly;
		         (toMatrixList#k)#k = (toMatrixList#k)#k + e_curIndex*adjpoly;
		         (toMatrixList#k)#j = (toMatrixList#k)#j - x_i*e_curIndex*adjpoly;
		         (toMatrixList#j)#k = (toMatrixList#j)#k - yirep#(i-1)*e_curIndex;		      
		 	  curIndex =curIndex + 1;
			  );
		 );
	     );
	 );          
	 
	  --need to remember to subtract from each diagonal by w
     for j from 0 to n-1 do( --going through nodes of a particular partition
	     for k from 0 to n-1 do( --going through the edges
		      if (j == k) then (
	 		  toMatrixList#j#j = toMatrixList#j#j - z*adjpoly;
	 		  );
		 );
	     );
     -- also need to finally remember to account for the edges in the complete graph of nodes (accounts for current bug)
      for j from 0 to n-1 do( --
	     for k from 0 to n-1 do(  --only want hit edges once, we can do this by only considering increasing pairs
		      if (j < k) then (
	 		 (toMatrixList#j)#j = (toMatrixList#j)#j + e_curIndex*adjpoly;
		         (toMatrixList#k)#k = (toMatrixList#k)#k + e_curIndex*adjpoly;
		         (toMatrixList#k)#j = (toMatrixList#k)#j - e_curIndex*adjpoly;
		         (toMatrixList#j)#k = (toMatrixList#j)#k - e_curIndex*adjpoly;		      
		 	  curIndex =curIndex + 1;
	 		  );
		 );
	     );   
     
     
	 
	 for j from 0 to n-1 do (
	 toMatrixList#j = new List from toMatrixList#j;
	 );    
    
     toMatrixList = new List from toMatrixList;
     outMatrix = matrix(toMatrixList);
     
--     stdio << theideals;
--     stdio << dims;
--     stdio << degs;
local DF;
local Ra;
    DF = {};    
    DF = append(DF, sub(det(outMatrix,Strategy => Cofactor),R));
    for i from 1 to m do (
    DF = append(DF,(diff(x_i, DF_0)));
    );
     Ra = QQ[x_1 .. x_m, y_1 .. y_m,z];
     return (outMatrix,R,I,Ra,DF)
     )
 



--will take number of actions (m) and number of  vertices (m)
--and return matrix representation for a general laplace beltrami operator
-- on an action adjacent dense graph
AdjacentDensePeriodicMatrix = method(TypicalValue => List)
AdjacentDensePeriodicMatrix(ZZ,ZZ) := (Sequence) => (m,n)-> (
     local edges; local toteEdges; local invFuncs; local toMatrixList, local toRowList; local zList; local outList; local outMatrix;
     local curIndex;
     local R;
     local W;
     local I;
     local adjpoly;
     local yirep;
     invFuncs = {};
     toteEdges = m*n*n + n*(n-1)//2;

     R = QQ[x_1 .. x_m, y_1 .. y_m,z, e_1 .. e_(toteEdges)];
     for i from 1 to m do(
    	invFuncs = append(invFuncs, x_i*y_i -1);	 
	 );     
     I = ideal (toSequence invFuncs);
     toMatrixList = new MutableList;
     zList = {};
     for i from 0 to n-1 do(
	 zList = append(zList, 0);
	 );     
     --make a double nested list so that we can convert this to a matrix later, make it nxn
     for i from 1 to n do( -- declare like this so that all rows are unique objects
	 toRowList = new MutableList;
	 for j from 1 to n do (
	 toRowList = append(toRowList, R_zList - R_zList);
	 );
	 toMatrixList = append(toMatrixList, toRowList);
	 );       
     --okay so lets iterate through each bipartite graph
     --we will look at all the edges connected to a node one at a time, going through the nodes not part of the complete subgraph
     -- this will require one loop to go through each bipartite graph and another loop to go through each node and another
     --for the n edges on each node
     --it will then require one more loop after these nested loops to account for the edges in the complete graph
     
     --yirep to represent y_i
     yirep = new MutableList;
     for i from 1 to m do(
	 yirep = append(yirep,1_R);
	 ); 
     --need yirep in this way so can indentify y_i = adjpoly_i
     for i from 0 to m-1 do(
	 for j from 1 to m do(
	     if ((i +1) != j) then(
	          yirep#i = yirep#i* (x_j)_R;    
		  );
	     );
	 );

--for the rest need to still multiply by all terms, use adjpoly
     adjpoly = 1_R;
     for i from 1 to m do(
	 adjpoly = adjpoly*(x_i)_R;
	 ); 
 
     for i from 1 to m do( --going through bipartite partitions
	 curIndex = (i-1)*n*n + 1;
	 for j from 0 to n-1 do( --going through nodes of a particular partition
	     for k from 0 to n-1 do( --going through the edges
		      if (j == k) then (
	 		  toMatrixList#j#j = toMatrixList#j#j +  ((2 - x_i)*adjpoly - yirep#(i-1))*e_curIndex;
			  curIndex =curIndex + 1;
	 		  );
		      if (j != k) then (
     	     	         (toMatrixList#j)#j = (toMatrixList#j)#j + e_curIndex*adjpoly;
		         (toMatrixList#k)#k = (toMatrixList#k)#k + e_curIndex*adjpoly;
		         (toMatrixList#k)#j = (toMatrixList#k)#j - x_i*e_curIndex*adjpoly;
		         (toMatrixList#j)#k = (toMatrixList#j)#k - yirep#(i-1)*e_curIndex;		      
		 	  curIndex =curIndex + 1;
			  );
		 );
	     );
	 );          
	 
	  --need to remember to subtract from each diagonal by w
     for j from 0 to n-1 do( --going through nodes of a particular partition
	     for k from 0 to n-1 do( --going through the edges
		      if (j == k) then (
	 		  toMatrixList#j#j = toMatrixList#j#j - z*adjpoly;
	 		  );
		 );
	     );
     -- also need to finally remember to account for the edges in the complete graph of nodes (accounts for current bug)
      for j from 0 to n-1 do( --
	     for k from 0 to n-1 do(  --only want hit edges once, we can do this by only considering increasing pairs
		      if (j < k) then (
	 		 (toMatrixList#j)#j = (toMatrixList#j)#j + e_curIndex*adjpoly;
		         (toMatrixList#k)#k = (toMatrixList#k)#k + e_curIndex*adjpoly;
		         (toMatrixList#k)#j = (toMatrixList#k)#j - e_curIndex*adjpoly;
		         (toMatrixList#j)#k = (toMatrixList#j)#k - e_curIndex*adjpoly;		      
		 	  curIndex =curIndex + 1;
	 		  );
		 );
	     );   
     
     
	 
	 for j from 0 to n-1 do (
	 toMatrixList#j = new List from toMatrixList#j;
	 );    
    
     toMatrixList = new List from toMatrixList;
     outMatrix = matrix(toMatrixList);
     
--     stdio << theideals;
--     stdio << dims;
--     stdio << degs;
local DF;
local Ra;
    DF = {};    
    DF = append(DF, sub(det(outMatrix,Strategy => Cofactor),R));
    for i from 1 to m do (
    DF = append(DF,(diff(x_i, DF_0)));
    );
     Ra = QQ[x_1 .. x_m, y_1 .. y_m,z];
     return (outMatrix,R,I,Ra,DF)
     )
 


--will take number of actions (m) and number of  vertices (m)
--and return matrix representation for a general laplace beltrami operator
-- on an action adjacent dense graph
AdjacentDensePeriodicMatrixInvL = method(TypicalValue => List)
AdjacentDensePeriodicMatrixInvL(ZZ,ZZ) := (Sequence) => (m,n)-> (
     local edges; local toteEdges; local invFuncs; local toMatrixList, local toRowList; local zList; local outList; local outMatrix;
     local curIndex;
     local R;
     local W;
     local I;
     local adjpoly;
     local yirep;
     invFuncs = {};
     toteEdges = m*n*n + n*(n-1)//2;

     R = QQ[x_1 .. x_m,z, y_1 .. y_m,zi, e_1 .. e_(toteEdges)];
     for i from 1 to m do(
    	invFuncs = append(invFuncs, x_i*y_i -1);	 
	 );     
     invFuncs = append(invFuncs, z*zi-1);
     I = ideal (toSequence invFuncs);
     toMatrixList = new MutableList;
     zList = {};
     for i from 0 to n-1 do(
	 zList = append(zList, 0);
	 );     
     --make a double nested list so that we can convert this to a matrix later, make it nxn
     for i from 1 to n do( -- declare like this so that all rows are unique objects
	 toRowList = new MutableList;
	 for j from 1 to n do (
	 toRowList = append(toRowList, R_zList - R_zList);
	 );
	 toMatrixList = append(toMatrixList, toRowList);
	 );       
     --okay so lets iterate through each bipartite graph
     --we will look at all the edges connected to a node one at a time, going through the nodes not part of the complete subgraph
     -- this will require one loop to go through each bipartite graph and another loop to go through each node and another
     --for the n edges on each node
     --it will then require one more loop after these nested loops to account for the edges in the complete graph
     
     --yirep to represent y_i
     yirep = new MutableList;
     for i from 1 to m do(
	 yirep = append(yirep,1_R);
	 ); 
     --need yirep in this way so can indentify y_i = adjpoly_i
     for i from 0 to m-1 do(
	 for j from 1 to m do(
	     if ((i +1) != j) then(
	          yirep#i = yirep#i* (x_j)_R;    
		  );
	     );
	 );

--for the rest need to still multiply by all terms, use adjpoly
     adjpoly = 1_R;
     for i from 1 to m do(
	 adjpoly = adjpoly*(x_i)_R;
	 ); 
 
     for i from 1 to m do( --going through bipartite partitions
	 curIndex = (i-1)*n*n + 1;
	 for j from 0 to n-1 do( --going through nodes of a particular partition
	     for k from 0 to n-1 do( --going through the edges
		      if (j == k) then (
	 		  toMatrixList#j#j = toMatrixList#j#j +  ((2 - x_i)*adjpoly - yirep#(i-1))*e_curIndex;
			  curIndex =curIndex + 1;
	 		  );
		      if (j != k) then (
     	     	         (toMatrixList#j)#j = (toMatrixList#j)#j + e_curIndex*adjpoly;
		         (toMatrixList#k)#k = (toMatrixList#k)#k + e_curIndex*adjpoly;
		         (toMatrixList#k)#j = (toMatrixList#k)#j - x_i*e_curIndex*adjpoly;
		         (toMatrixList#j)#k = (toMatrixList#j)#k - yirep#(i-1)*e_curIndex;		      
		 	  curIndex =curIndex + 1;
			  );
		 );
	     );
	 );          
	 
	  --need to remember to subtract from each diagonal by w
     for j from 0 to n-1 do( --going through nodes of a particular partition
	     for k from 0 to n-1 do( --going through the edges
		      if (j == k) then (
	 		  toMatrixList#j#j = toMatrixList#j#j - z*adjpoly;
	 		  );
		 );
	     );
     -- also need to finally remember to account for the edges in the complete graph of nodes (accounts for current bug)
      for j from 0 to n-1 do( --
	     for k from 0 to n-1 do(  --only want hit edges once, we can do this by only considering increasing pairs
		      if (j < k) then (
	 		 (toMatrixList#j)#j = (toMatrixList#j)#j + e_curIndex*adjpoly;
		         (toMatrixList#k)#k = (toMatrixList#k)#k + e_curIndex*adjpoly;
		         (toMatrixList#k)#j = (toMatrixList#k)#j - e_curIndex*adjpoly;
		         (toMatrixList#j)#k = (toMatrixList#j)#k - e_curIndex*adjpoly;		      
		 	  curIndex =curIndex + 1;
	 		  );
		 );
	     );   
     
     
	 
	 for j from 0 to n-1 do (
	 toMatrixList#j = new List from toMatrixList#j;
	 );    
    
     toMatrixList = new List from toMatrixList;
     outMatrix = matrix(toMatrixList);
     
--     stdio << theideals;
--     stdio << dims;
--     stdio << degs;
local DF;
local Ra;
    DF = {};    
    DF = append(DF, sub(det(outMatrix,Strategy => Cofactor),R));
    for i from 1 to m do (
    DF = append(DF,(diff(x_i, DF_0)));
    );
     Ra = QQ[x_1 .. x_m,z,y_1 .. y_m,zi];
     return (outMatrix,R,I,Ra,DF)
     )
 

--a by b vertex  fundamental domains for the 2x2 matrix
--and return matrix representation for a general laplace beltrami operator
-- on an action adjacent dense graph
BlockPeriodicMatrix = method(TypicalValue => List)
BlockDensePeriodicMatrix(ZZ,ZZ) := (Sequence) => (a,b)-> (
     local edges; local toteEdges; local invFuncs; local toMatrixList, local toRowList; local zList; local outList; local outMatrix;
     local curIndex;
     local R;
     local W;
     local I;
     local adjpoly;
     local yirep;
     invFuncs = {};

     R = QQ[x_1,x_2, y_1,y_2,z, eo_1 .. eo_(a+b),er_1 .. er_((a-1)*b),ed_1 .. ed_((b-1)*a),q_1 .. q_(a*b)];
     for i from 1 to 2 do(
        invFuncs = append(invFuncs, x_i*y_i -1);     
     );     
     I = ideal (toSequence invFuncs);
     toMatrixList = new MutableList;
     zList = {};
     for i from 0 to n-1 do(
     zList = append(zList, 0);
     );     
     --make a double nested list so that we can convert this to a matrix later, make it nxn
     for i from 1 to n do( -- declare like this so that all rows are unique objects
     toRowList = new MutableList;
     for j from 1 to n do (
     toRowList = append(toRowList, R_zList - R_zList);
     );
     toMatrixList = append(toMatrixList, toRowList);
     );       
     --okay so lets iterate through each bipartite graph
     --we will look at all the edges connected to a node one at a time, going through the nodes not part of the complete subgraph
     -- this will require one loop to go through each bipartite graph and another loop to go through each node and another
     --for the n edges on each node
     --it will then require one more loop after these nested loops to account for the edges in the complete graph
     
     --yirep to represent y_i
     yirep = new MutableList;
     for i from 1 to m do(
     yirep = append(yirep,1_R);
     ); 
     --need yirep in this way so can indentify y_i = adjpoly_i
     for i from 0 to m-1 do(
     for j from 1 to m do(
         if ((i +1) != j) then(
              yirep#i = yirep#i* (x_j)_R;    
          );
         );
     );

--for the rest need to still multiply by all terms, use adjpoly
     adjpoly = 1_R;
     for i from 1 to m do(
     adjpoly = adjpoly*(x_i)_R;
     ); 
 
     for j from 0 to a-1 do( --vertical edges leaving FD
  k = a*b-a;
 (toMatrixList#j)#j = (toMatrixList#j)#j + eo_(j+1)*adjpoly;
                 (toMatrixList#k)#k = (toMatrixList#k)#k + eo_(j+1)*adjpoly;
                 (toMatrixList#j)#k = (toMatrixList#j)#k - x_2*eo_(j+1)*adjpoly;
                 (toMatrixList#k)#j = (toMatrixList#k)#j - yirep#1*eo_(j+1)*adjpoly;      
              ); 
     for j from 0 to b-1 do( --horizontal edges leaving FD
  l = 1 + a*j;
  k = a*(j+1);
          curIndex = a+j+1;
         (toMatrixList#l)#l = (toMatrixList#l)#l + eo_curIndex*adjpoly;
                 (toMatrixList#k)#k = (toMatrixList#k)#k + eo_curIndex*adjpoly;
                 (toMatrixList#k)#l = (toMatrixList#k)#l - x_1*eo_curIndex*adjpoly;
                 (toMatrixList#l)#k = (toMatrixList#l)#k - yirep#0*eo_curIndex*adjpoly;
              );
      --need to remember to subtract from each diagonal by w
     for j from 0 to a*b-1 do( --going through nodes of a particular partition
              toMatrixList#j#j = toMatrixList#j#j - z*adjpoly + q_(j+1);
         );
     -- right internal edges
curIndex = 1;
      for j from 0 to b-1 do( --
for k from 0 to a-2 do( --going through the edges
cur = j*a + k;
adj = j*a + k + 1;
(toMatrixList#cur)#cur = (toMatrixList#cur)#cur + er_curIndex*adjpoly;
                 (toMatrixList#adj)#adj = (toMatrixList#adj)#adj + er_curIndex*adjpoly;
(toMatrixList#cur)#adj = (toMatrixList#cur)#adj - er_curIndex*adjpoly;
                 (toMatrixList#adj)#cur = (toMatrixList#adj)#cur - er_curIndex*adjpoly;
curIndex = curIndex+1;
);
         );   

     -- down internal edges
curIndex  = 0;
for j from 0 to a-1 do( --
for k from 0 to b-2 do( --going through the edges
cur = j*b + k;
adj = j*b + k + a;
(toMatrixList#cur)#cur = (toMatrixList#cur)#cur + ed_curIndex*adjpoly;
                (toMatrixList#adj)#adj = (toMatrixList#adj)#adj + ed_curIndex*adjpoly;
(toMatrixList#cur)#adj = (toMatrixList#cur)#adj - ed_curIndex*adjpoly;
                (toMatrixList#adj)#cur = (toMatrixList#adj)#cur - ed_curIndex*adjpoly;
curIndex = curIndex +1;
);
         );        
     
     
     for j from 0 to a*b-1 do (
     toMatrixList#j = new List from toMatrixList#j;
     );    
    
     toMatrixList = new List from toMatrixList;
     outMatrix = matrix(toMatrixList);
     
--     stdio << theideals;
--     stdio << dims;
--     stdio << degs;
local DF;
local Ra;
    DF = {};    
    DF = append(DF, sub(det(outMatrix,Strategy => Cofactor),R));
    for i from 1 to 2 do (
    DF = append(DF,(diff(x_i, DF_0)));
    );
     Ra = QQ[x_1,x_2, y_1,y_2,z];
     return (outMatrix,R,I,Ra,DF)
     )
 



--will take number of actions (m) and number of  vertices (m)
--and return matrix representation for a general laplace beltrami operator
-- on an action adjacent dense graph
AdjacentDensePeriodicMatrix = method(TypicalValue => List)
AdjacentDensePeriodicMatrix(ZZ,ZZ) := (Sequence) => (m,n)-> (
     local edges; local toteEdges; local invFuncs; local toMatrixList, local toRowList; local zList; local outList; local outMatrix;
     local curIndex;
     local R;
     local W;
     local I;
     local adjpoly;
     local yirep;
     invFuncs = {};
     toteEdges = m*n*n + n*(n-1)//2;

     R = QQ[x_1 .. x_m, y_1 .. y_m,z, e_1 .. e_(toteEdges)];
     for i from 1 to m do(
        invFuncs = append(invFuncs, x_i*y_i -1);     
     );     
     I = ideal (toSequence invFuncs);
     toMatrixList = new MutableList;
     zList = {};
     for i from 0 to n-1 do(
     zList = append(zList, 0);
     );     
     --make a double nested list so that we can convert this to a matrix later, make it nxn
     for i from 1 to n do( -- declare like this so that all rows are unique objects
     toRowList = new MutableList;
     for j from 1 to n do (
     toRowList = append(toRowList, R_zList - R_zList);
     );
     toMatrixList = append(toMatrixList, toRowList);
     );       
     --okay so lets iterate through each bipartite graph
     --we will look at all the edges connected to a node one at a time, going through the nodes not part of the complete subgraph
     -- this will require one loop to go through each bipartite graph and another loop to go through each node and another
     --for the n edges on each node
     --it will then require one more loop after these nested loops to account for the edges in the complete graph
     
     --yirep to represent y_i
     yirep = new MutableList;
     for i from 1 to m do(
     yirep = append(yirep,1_R);
     ); 
     --need yirep in this way so can indentify y_i = adjpoly_i
     for i from 0 to m-1 do(
     for j from 1 to m do(
         if ((i +1) != j) then(
              yirep#i = yirep#i* (x_j)_R;    
          );
         );
     );

--for the rest need to still multiply by all terms, use adjpoly
     adjpoly = 1_R;
     for i from 1 to m do(
     adjpoly = adjpoly*(x_i)_R;
     ); 
 
     for i from 1 to m do( --going through bipartite partitions
     curIndex = (i-1)*n*n + 1;
     for j from 0 to n-1 do( --going through nodes of a particular partition
         for k from 0 to n-1 do( --going through the edges
              if (j == k) then (
              toMatrixList#j#j = toMatrixList#j#j +  ((2 - x_i)*adjpoly - yirep#(i-1))*e_curIndex;
              curIndex =curIndex + 1;
              );
              if (j != k) then (
                         (toMatrixList#j)#j = (toMatrixList#j)#j + e_curIndex*adjpoly;
                 (toMatrixList#k)#k = (toMatrixList#k)#k + e_curIndex*adjpoly;
                 (toMatrixList#k)#j = (toMatrixList#k)#j - x_i*e_curIndex*adjpoly;
                 (toMatrixList#j)#k = (toMatrixList#j)#k - yirep#(i-1)*e_curIndex;            
              curIndex =curIndex + 1;
              );
         );
         );
     );          
     
      --need to remember to subtract from each diagonal by w
     for j from 0 to n-1 do( --going through nodes of a particular partition
         for k from 0 to n-1 do( --going through the edges
              if (j == k) then (
              toMatrixList#j#j = toMatrixList#j#j - z*adjpoly;
              );
         );
         );
     -- also need to finally remember to account for the edges in the complete graph of nodes (accounts for current bug)
      for j from 0 to n-1 do( --
         for k from 0 to n-1 do(  --only want hit edges once, we can do this by only considering increasing pairs
              if (j < k) then (
             (toMatrixList#j)#j = (toMatrixList#j)#j + e_curIndex*adjpoly;
                 (toMatrixList#k)#k = (toMatrixList#k)#k + e_curIndex*adjpoly;
                 (toMatrixList#k)#j = (toMatrixList#k)#j - e_curIndex*adjpoly;
                 (toMatrixList#j)#k = (toMatrixList#j)#k - e_curIndex*adjpoly;            
              curIndex =curIndex + 1;
              );
         );
         );   
     
     
     
     for j from 0 to n-1 do (
     toMatrixList#j = new List from toMatrixList#j;
     );    
    
     toMatrixList = new List from toMatrixList;
     outMatrix = matrix(toMatrixList);
     
--     stdio << theideals;
--     stdio << dims;
--     stdio << degs;
local DF;
local Ra;
    DF = {};    
    DF = append(DF, sub(det(outMatrix,Strategy => Cofactor),R));
    for i from 1 to m do (
    DF = append(DF,(diff(x_i, DF_0)));
    );
     Ra = QQ[x_1 .. x_m, y_1 .. y_m,z];
     return (outMatrix,R,I,Ra,DF)
     )
 

