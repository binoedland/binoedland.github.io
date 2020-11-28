newPackage("Parliament",
    Version => "0.1",
    Date => "November 29, 2016",
    Authors => {{Name=> "Bernt Ivar Utstoel Noedland"}},
    PackageExports => {"ToricVectorBundles"},
    Headline => "Computes parliament of polytopes for toric vector bundles and checks positivity properties"
    );

export{"globallyGenerated",
       "Parliament",
       "parliament",
       "getCharacters",
       "getPolytopes",
       "getMatroid"
    }

Parliament = new Type of HashTable



globallyGenerated = method();
globallyGenerated Parliament :=  (P) ->    (
    globalResult := true;
    T := P#"tvb";
    n := rank T;
    polytopeList := P#"Polytopes";
    m := maxCones T;
    A := rays T;
    Usigma := P#"Usigma";   
    for i from 0 to (length m -1) do (
	result := false;
	sigma := m_i;
	usigma := Usigma#sigma;
	for j from 0 to (length usigma - 1) do (
	    current := usigma_j;
	    counter := 0;
	    for k from 0 to (n-1) do (
    	    	v := ((current_k)_0);
		if  contains(polytopeList#((current_k)_1)_0,v) then (
		    counter = counter+1;
		    );
		);
	    if counter == n then (
		result = true;
	    	);
	    );
	if result == false then (
	    globalResult = false;
	    );
	);
    globalResult
    );
		    





coneFindVectors = method();
coneFindVectors (Matrix,ZZ,ZZ,List) := (A,m,counter,L) -> (
    T := L_0;
    sigma := L_1;
    E := L_2;
    n := L_3;
    GM := L_4;
    ray := (rays sigma)_{counter};
    d := numgens source (filtration T)#ray-1;   
    s := nextIndex((filtration T)#ray,m);
    v := modulo(A,(base T)#ray_(s_0)); 
    B := A*v;
    l := numgens source B;
    if l == n then (
	matrixList := changeBasis(B,GM);
	F := {};
	bool := false;
	for i from 0 to (length E -1) do (
	    for k from 0 to (length matrixList -1) do (	
		for j from 0 to (numgens source matrixList_k -1) do (
	    	    if modulo(E_i,(matrixList_k)_{j}) == 0 then (
			F = F | {E_i | (matrixList_k)_{j}};
			bool = true;
			);
		    );
		);
	    );
	if bool == true then (
	    E = F;
	    );
	if A == id_(QQ^(rank T)) then (
	    if counter < (numgens source rays sigma -1) then (
		    newRay := (rays sigma)_{counter + 1};
	    	    E = coneFindVectors(A,(minValue(((filtration T)#newRay))-1),counter+1,{T,sigma,E,n,GM});
		    );
		);
	    )
    else if l > n then (
	if counter < ((numgens source (rays sigma)) - 1 ) then (
	    newRay = (rays sigma)_{counter + 1};
	    E = coneFindVectors(B,(minValue(((filtration T)#newRay))-1),counter+1,L);
	    );
	);
    if (length s_0 - 1) < d then (
	E = coneFindVectors(A,s_1,counter,{T,sigma,E,n,GM});
	);
    E);






findVectors = method();
findVectors (Matrix,ZZ,ZZ,ToricVectorBundleKlyachko) := (A,m,counter,T) -> (
    E := matrix(QQ,{{}});
    for i from 2 to rank T do E = E || matrix(QQ,{{}});
    ray := (rays T)_counter;
    d := numgens source (filtration T)#ray-1;   
    s := nextIndex((filtration T)#ray,m);
    v := modulo(A,(base T)#ray_(s_0)); 
    B := A*v;
    l := numgens source B;
    if l == 1 then (
	E = E | B;
	)
    else if l > 1 then (
	if counter < ((length rays T) - 1 ) then (
	    newRay := (rays T)_(counter + 1);
	    E = E | findVectors(B,(minValue(((filtration T)#newRay))-1),counter+1,T);
	    );
	);
    if (length s_0 -1) < d then (
	E = E | findVectors(A,s_1,counter,T);
	);
    E);


findVectors2 = method();
findVectors2 (Matrix,ZZ,ZZ,List) := (A,m,counter,L) -> (
    T := L_0;
    n := L_1;
    Ematrix := L_2;
    E := matrix(QQ,{{}});
    for i from 2 to rank T do E = E || matrix(QQ,{{}});
    ray := (rays T)_counter;
    d := numgens source (filtration T)#ray-1;   
    s := nextIndex((filtration T)#ray,m);
    v := modulo(A,(base T)#ray_(s_0)); 
    B := A*v;
    l := numgens source B;  
    bool := false;
    if l == n then (
	for i from 0 to (l-1) do (
	    if bool == false then (
	    	c := numgens source Ematrix -1;
	    	colList := 0..c;
		if c+1 >= l-i then (
	    	    subSets :=subsets(colList,l-i); 
	    	    for j from 0 to (length subSets-1) do (
		    	w :=modulo(B,Ematrix_(subSets_j));
		    	if (numgens source w)== l-i then (
			    if complement w != 0 then (
    		    	    	Bnew := B*complement w;
		    	    	E = E | Bnew;
			    );
		    	    bool = true;
		    	    break;
		    	    );
		    	);
		    );
		);
	    );
	if bool == false then (
	    E = E | B;
	    );
	if A == id_(QQ^(rank T)) then (
	    if counter < (length rays T -1) then (
		newRay := (rays T)_(counter + 1);
		E = E | findVectors2(A,(minValue(((filtration T)#newRay))-1),counter+1,L);
		);
	    );
	)
    else if l > n then (
	if counter < ((length rays T) - 1 ) then (
	    newRay = (rays T)_(counter + 1);
	    E = E | findVectors2(B,(minValue(((filtration T)#newRay))-1),counter+1,L);
	    );
	);
    if (length s_0 -1) < d then (
	E = E | findVectors2(A,s_1,counter,L);
	);
    E);



findBases = method();
findBases (ToricVectorBundleKlyachko,Matrix) := (T,GM) -> (
    basisList := {};
    m := maxCones T;
    for i from 0 to (length m - 1) do (
	sigma := m_i;
	ray := (rays sigma)_{0};
	start := minValue(((filtration T)#ray));
	emp := emptyMat(rank T);
	cList := coneFindVectors((id_(QQ^(rank T))),start-1,0,{T,sigma,{emp},1,GM});
	bList := {};
	for i from 0 to (length cList -1) do (
	    bList = bList | {deleteDuplicates(cList_i)};
	    );
	n := 2;
	while (numgens source bList_0 < rank T) do (
	    cList =  coneFindVectors((id_(QQ^(rank T))),start-1,0,{T,sigma,bList,n,GM});	    
	    dList := {};
	    for j from 0 to (length cList -1) do (
		dList = dList | {sort(deleteDuplicates(cList_j))};
		); 
	    dList = deleteDuplicatesList(dList);
	    bList = dList;
	    n = n+1;
	    );
    	basisList = basisList | {(sigma,bList)};
    	);
    H :=  hashTable(basisList);
    H);






findMatroid = method();
findMatroid (ToricVectorBundleKlyachko) := (T) -> (
    basisList := {};
    A := rays T;
    ray := A_0;
    start := minValue(((filtration T)#ray));
    b:=findVectors((id_(QQ^(rank T))),start-1,0,T);
    basisList = basisList | {b};
    H :=  deleteDuplicates(matrix{basisList});
    n := 2;
    while (n < (rank T +1)) do (
	k := H | findVectors2((id_(QQ^(rank T))),start-1,0,{T,n,H});
	H = deleteDuplicates ( H |findVectors2((id_(QQ^(rank T))),start-1,0,{T,n,H}));
	n = n+1;	
	);
    H); 





parliament = method();
parliament (ToricVectorBundleKlyachko) := T -> (
    m :=maxCones T;
    sigmaList := {};
    globalList := {};
    filtList := {};
    a := (rays T)_0;
    G := findMatroid(T);
    H := findBases(T,G);
    for i from 0 to (length m -1) do (
    	sigmaList = {};
    	sigma := m#i;
	A := rays sigma;
    	Elist := H#sigma;
	filtList = {};
	for j from 0 to (length Elist -1 ) do (
	    E := Elist_j;
	    sigmaList = {};
	    for k from 0 to (rank(T)-1) do (
	    	rho := A_{0};
	    	e := E_{k};
	    	rhoE := (base T)#rho;	    
	    	b := matrix{{spanIndex(rhoE,(filtration T)#rho,e)}};
	    	for t from 1 to (dim sigma-1) do (
                    rho =  A_{t};
		    rhoE = (base T)#rho;
		    b = b || matrix{{spanIndex(rhoE,(filtration T)#rho,e)}};
		    );
	    	usigma := solve(transpose A,b);
	    	sigmaList = sigmaList | {{-usigma,e}};
	    	);
	    filtList = filtList  | {sigmaList};
	    );
	globalList = globalList | {(sigma,filtList)};
    	);
    U := hashTable(globalList);
    polytopeList := findPolytopes(T,U,G); 
    pp := hashTable(polytopeList);
    P := new Parliament from {
	"tvb" => T,
	"Usigma" =>U,
	"Gmatroid" =>G,
	"Polytopes" =>pp};
    P);

getMatroid = method();
getMatroid Parliament := (P) -> (
    T:=P#"tvb";
    G := P#"Gmatroid";
    G);

getPolytopes = method();
getPolytopes Parliament := (P) -> (
    T:=P#"tvb";
    polyList := P#"Polytopes";
    G:=P#"Gmatroid";
    result := {};
    for i from 0 to (numgens source G -1) do (
	matPoly := polyList#(G_{i});
	v := (vertices (matPoly_0));
	if (numgens source v == 0 ) then v = {};
	result = result | {{G_{i},v}};
	);
    hashTable(result));


getCharacters = method();
getCharacters Parliament := (P) -> (
    T:=P#"tvb";
    Usigma:=P#"Usigma";
    mCones := maxCones(T);
    F := fan(T);
    result := {};
    for i from 0 to (length mCones-1) do (
	N := emptyMat(dim F);
	sigma := mCones_i;
	uList := (Usigma#sigma)_0;
	for j from 0 to (length uList -1) do (
	    N = N | (uList_j)_0;
	    );
	result = result | {{rays sigma,N}};
	);
    hashTable(result));



findPolytopes = method();
findPolytopes (ToricVectorBundleKlyachko, HashTable,Matrix) := (T,U,G) -> (
    A := matrix{rays T};
    polytopeList := {};
    for i from 0 to (numgens source G -1) do (
	l := {};
	b := G_{i};
	for k from 0 to (numgens source A -1) do (
	    w := A_{k};
	    E := (base T)#w;
	    d := spanIndex(E,(filtration T)#w,b);
	    l = l | {d};
	    );
	v := -transpose matrix{l};
	polytopeList = polytopeList | {{b,{intersection(transpose A,v),(A,v)}}};
	);
    polytopeList);


spanIndex = method();
spanIndex(Matrix,Matrix,Matrix) := (A,F,v) -> (
    start := minValue F -1;
    s := nextIndex(F,start);
    while (length s_0) < (numgens source A) do (
	if modulo(A_(s_0),v) != 0 then (break);
	s = nextIndex(F,s_1);
	);
    s_1);
	


isColumn = method();
isColumn(Matrix,Matrix) := (A,v) -> (
result := false;
for i from 0 to (numgens source A - 1) do (
if A_{i} == v then result = true;
);
result
);

deleteDuplicates = method();
deleteDuplicates(Matrix) := (A) ->  (
    B := matrix{{}};
    if A != 0 then (
    	B = A_{0};
    	for i from 1 to (numgens source A - 1) do (
	    if isColumn(B,A_{i}) == false then (
	    	if isColumn(B,-A_{i}) == false then (
		    B = B | A_{i};
		    );
	    	);
	    );
	)
    else (
	B = A;
	);
    B);


	
columnNumber = method();
columnNumber(Matrix,Matrix) := (ray,M) -> (
    n :=-1;
    for i from 0 to (numgens source M -1) do (
	if M_{i} == ray then (n=i);
	);
    n);
		    

toMatrix = method();
toMatrix (Module) := (A) -> (
    n := rank ambient A;
    a := A_0;
    v := matrix{{a_0}};
    for i from 1 to (n-1) do (
	v = v || matrix{{a_i}};
	);
    v);

    
mult = method();
mult (Matrix,Matrix) := (v,w) -> (
    d := 0;
    for i from 0 to (numgens target v -1 ) do (
	if w_(i,0) != 0 then (
	    d = v_(i,0)/w_(i,0);
	    break;
	    );
	);
    d);

isPositive = method();
isPositive (Matrix) := (v) -> (
    result := true;
    for i from 0 to (numgens target v -1 ) do(
	if v_(i,0) < 0 then (
	    result = false;
	    );
	);
    result);
    


deleteDuplicatesList = method();
deleteDuplicatesList(List) := L -> (
    start := sort(L_0);
    newList := {start};
    for i from 1 to (length L-1) do (
	bool := false;
	for j from 0 to (length newList - 1) do (
	    if sort(L_i) == newList_j then bool = true;
	    );
	if bool == false then (
	    newList = newList | {L_i};
	    );
	);
    newList);

minValue = method();
minValue (Matrix) := (A) -> (
    l := {};
    for i from 0 to (numgens source A -1) do (
	l = l | {A_(0,i)};
	);
    n := min l;
    n);

nextIndex = method();
nextIndex (Matrix,ZZ) := (A,n) -> (
    d := 0;
    l := {};
    m := numgens source A - 1;
    Alist := {};
    for i from 0 to m do (
	Alist  = Alist | {A_(0,i)};
	);
    B := sort(Alist);
    if (n >= max(B)) then (
	l = splice {0..(length B-1)};
	d = n;
	)
    else (
    	for i from 0 to m do (
	    if B_i > n then (
	    	d = B_i;
	    	break;
	    	);
	    );
   	for i from 0 to m do (
       	    if A_(0,i) <= d then (
	   	l = l | {i} ;
	   	);
       	    );
	);
   {l,d}); 


emptyMat = method();
emptyMat (ZZ) := (n) -> (
    E := matrix(QQ,{{}});
    for i from 2 to n do E = E || matrix(QQ,{{}});
    E);

emptyMatZ = method();
emptyMatZ (ZZ) := (n) -> (
    E := matrix{{}};
    for i from 2 to n do E = E || matrix{{}};
    E);

containsM = method();
containsM (List,Matrix) := (ll,B) -> (
    result := false;
    contIndex := -1;
    for i from 0 to (length ll -1) do (
	if (ll_i)_0 == B then (
	    result = true;
	    contIndex = i;
	    break;
	    );
	);
    {result,contIndex});

changeBasis = method();
changeBasis (Matrix,Matrix) := (C,G) -> (
    d:= numgens source G -1;
    l := 0..d;
    e := numgens source C;
    subList := subsets(l,e);
    indexList := {};
    for i from 0 to (length subList -1) do (
	j := subList_i;
	if (modulo(C,G_j) == id_(QQ^e)) then (
    	    indexList = indexList | {j};
	    );
    	);
    matrixList := {};
    for i from 0 to (length indexList -1) do (
	matrixList = matrixList | {G_(indexList_i)};
	);
    matrixList);
   
